#!/usr/bin/env python

from __future__ import print_function
import os
import tempfile
import subprocess
import shutil
import sys
import glob
import timeit
import errno
import atexit
import collections
import yaml
import venv
import logging
from typing import Optional, Mapping

SCRIPT_ROOT = os.path.abspath(os.path.dirname(__file__))
FABM_BASE = os.path.join(SCRIPT_ROOT, "../..")

DEFAULT_FABM_URL = "https://github.com/fabm-model/fabm.git"
DEFAULT_GOTM_URL = "https://github.com/gotm-model/code.git"

def hack_pyyaml():
    # Do not convert on/off to bool
    # [done by pyyaml according to YAML 1.1, dropped from YAML 1.2]
    del yaml.SafeLoader.yaml_implicit_resolvers["o"]
    del yaml.SafeLoader.yaml_implicit_resolvers["O"]

    def none_representer(self, _):
        return self.represent_scalar("tag:yaml.org,2002:null", "")
    yaml.add_representer(type(None), none_representer, Dumper=yaml.SafeDumper)

hack_pyyaml()

def run(phase: str, args, verbose: bool = False, **kwargs):
    indent = "  " * phase.count("/")
    current_phase = phase.rsplit("/", 1)[-1]
    print(f"{indent}{current_phase}... ", end="")
    sys.stdout.flush()
    proc = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        **kwargs,
    )
    stdoutdata, _ = proc.communicate()
    if proc.returncode != 0:
        log_path = f"{phase.replace('/', '_')}.log"
        with open(log_path, "w") as f:
            f.write(stdoutdata)
        logs.append(log_path)
        print(f"FAILED (return code {proc.returncode}, log written to {log_path})")
    else:
        print("SUCCESS")
    if verbose:
        line = 80 * "-"
        print(f"Output:\n{line}\n{stdoutdata}\n{line}")
    return proc.returncode


def git_clone(phase: str, url: str, workdir: str, branch: Optional[str] = None):
    run(f"{phase}/clone", ["git", "clone", url, workdir])
    if branch is not None:
        run(f"{phase}/checkout", ["git", "checkout", branch], cwd=workdir)
    run(
        f"{phase}/submodule",
        ["git", "submodule", "update", "--init", "--recursive"],
        cwd=workdir,
    )


def run_gotm(setup_dir: str, gotm_exe: str):
    start = timeit.default_timer()
    p = subprocess.Popen(
        [gotm_exe],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        cwd=setup_dir,
    )
    stdoutdata, _ = p.communicate()
    duration = timeit.default_timer() - start
    if p.returncode != 0:
        print("FAILED - last output:")
        last_lines = stdoutdata.rsplit("\n", 10)[1:]
        print("\n".join(last_lines))
    else:
        print(f"ok ({duration:.3f} s)")
    return p.returncode == 0, duration


def cmake(
    phase: str,
    build_dir: str,
    source_dir: str,
    cmake_path: str = "cmake",
    target: Optional[str] = None,
    cmake_arguments=[],
):
    # Create and change to build directory
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)
    os.mkdir(build_dir)

    if os.name == "nt" and os.environ.get("CMAKE_GENERATOR", "") != "Ninja":
        x64 = sys.maxsize > 2**32
        cmake_arguments = ["-A", "x64" if x64 else "Win32"] + cmake_arguments

    previous_cache = os.path.join(source_dir, "CMakeCache.txt")
    if os.path.isfile(previous_cache):
        print(
            f'\n\nCannot use "{source_dir}" as source directory for cmake'
            f" because it has previously been used as build directory.\n"
            f'Please delete "{os.path.abspath(previous_cache)}" before continuing.'
        )
        sys.exit(1)

    # Build
    try:
        ret = run(
            f"{phase}/configure",
            [cmake_path, source_dir] + cmake_arguments,
            cwd=build_dir,
        )
    except EnvironmentError as e:
        if e.errno != errno.ENOENT:
            raise
        print(
            f'\n\ncmake executable ("{cmake_path}") not found.'
            f" Specify its location on the command line with --cmake."
        )
        sys.exit(2)

    if ret == 0:
        args = ["--config", "Debug"]
        if target is not None:
            args = args + ["--target", target]
        ret = run(f"{phase}/build", [cmake_path, "--build", "."] + args, cwd=build_dir)

    return ret == 0


def compare_netcdf(path: str, ref_path: str):
    import numpy
    import netCDF4

    perfect = True
    nc = netCDF4.Dataset(path)
    nc_ref = netCDF4.Dataset(ref_path)
    for varname in nc.variables.keys():
        if varname not in nc_ref.variables or varname in (
            "lon",
            "lat",
            "h",
            "z",
            "time",
        ):
            continue
        ncvar = nc.variables[varname]
        ncvar_ref = nc_ref.variables[varname]
        dat = ncvar[...]
        valid = numpy.isfinite(dat)
        if not valid.all():
            nbad = valid.size - valid.sum()
            print(f"    {varname}: {nbad} of {valid.size} values are invalid")
            perfect = False
        else:
            delta = dat - ncvar_ref[...]
            maxdelta = numpy.abs(delta).max()
            perfect = perfect and maxdelta == 0.0
            print(f"    {varname}: max abs difference = {maxdelta}")
    nc.close()
    nc_ref.close()
    return perfect


def test(
    gotm_setup_dir: str,
    work_root: str,
    testcases: Mapping[str, str],
    cmake_path: str = "cmake",
    cmake_arguments=[],
    fabm_url: str = DEFAULT_FABM_URL,
    gotm_url: str = DEFAULT_GOTM_URL,
    fabm_branch: Optional[str] = None,
    gotm_branch: Optional[str] = None,
):
    gotm_base = os.path.join(work_root, "code/gotm")
    build_dir = os.path.join(work_root, "build")

    # Get latest GOTM [public]
    git_clone("test/gotm", gotm_url, gotm_base, gotm_branch)

    cmake(
        "test_gotm",
        build_dir,
        gotm_base,
        cmake_path,
        cmake_arguments=[f"-DFABM_BASE={FABM_BASE}"] + cmake_arguments,
    )
    exe = os.path.join(build_dir, "Debug/gotm.exe" if os.name == "nt" else "gotm")
    with open(os.path.join(gotm_setup_dir, "gotm.yaml"), "r") as f:
        gotm_yaml = yaml.safe_load(f)
    gotm_yaml["fabm"] = {"use": True}
    gotm_yaml["output"] = {
        "result": {
            "time_unit": "day",
            "time_step": 1,
            "time_method": "mean",
            "variables": [{"source": "fabm/*"}, {"source": "h"}],
        }
    }
    custom_gotm_yaml_file = os.path.join(build_dir, "gotm.yaml")
    with open(custom_gotm_yaml_file, "w") as f:
        yaml.dump(gotm_yaml, f, default_flow_style=False)
    for name, path in testcases.items():
        shutil.copyfile(path, os.path.join(gotm_setup_dir, "fabm.yaml"))
        run(f"test/gotm/{name}", [exe, custom_gotm_yaml_file], cwd=gotm_setup_dir)


def compare(
    gotm_setup_dir: str,
    work_root: str,
    testcases: Mapping[str, str],
    cmake_path: str = "cmake",
    cmake_arguments=[],
    fabm_url: str = DEFAULT_FABM_URL,
    gotm_url: str = DEFAULT_GOTM_URL,
    fabm_branch: Optional[str] = None,
    gotm_branch: Optional[str] = None,
    fabm_ref_branch: Optional[str] = None,
    gotm_ref_branch: Optional[str] = None,
):
    assert fabm_branch != fabm_ref_branch or gotm_branch != gotm_ref_branch
    fabm_base = os.path.join(work_root, "code/fabm")
    gotm_base = os.path.join(work_root, "code/gotm")
    git_clone(fabm_url, fabm_base, fabm_branch)
    git_clone(gotm_url, gotm_base, gotm_branch)
    cmake(
        "test_gotm",
        os.path.join(work_root, "build"),
        gotm_base,
        cmake_path,
        f"-DFABM_BASE={fabm_base}",
        cmake_arguments=cmake_arguments,
    )

    ref_fabm_base = os.path.join(work_root, "ref/code/fabm")
    ref_gotm_base = os.path.join(work_root, "ref/code/gotm")
    git_clone(fabm_url, ref_fabm_base, fabm_ref_branch)
    git_clone(gotm_url, ref_gotm_base, gotm_ref_branch)
    cmake(
        "ref_gotm",
        os.path.join(work_root, "ref/build"),
        ref_gotm_base,
        cmake_path,
        f"-DFABM_BASE={ref_fabm_base}",
        cmake_arguments=cmake_arguments,
    )

    faster, slower = [], []
    failed, success, crashed = [], [], []
    for name in enumerate_testcases(
        testcase_dir, os.path.join(work_root, "code/fabm/testcases/*.yaml")
    ):
        print(f"TESTING {name}...")
        print("  reference...", end="")
        valid_ref, duration_ref = run_gotm(
            testcase_dir, os.path.join(work_root, "ref/build/gotm")
        )
        os.rename(
            os.path.join(testcase_dir, "result.nc"),
            os.path.join(testcase_dir, "result_ref.nc"),
        )
        print("  target...", end="")
        valid, duration = run_gotm(testcase_dir, os.path.join(work_root, "build/gotm"))
        if valid and valid_ref:
            if compare_netcdf(
                os.path.join(testcase_dir, "result_ref.nc"),
                os.path.join(testcase_dir, "result.nc"),
            ):
                success.append(name)
            else:
                failed.append(name)
        else:
            crashed.append(name)
        if duration < duration_ref:
            faster.append(name)
        else:
            slower.append(name)

    print(f"{len(success)} perfect matches: {', '.join(success)}")
    print(f"{len(failed)} mismatches: {', '.join(failed)}")
    print(f"{len(crashed)} failed to run: {', '.join(crashed)}")
    print(
        f"Faster than reference? {len(faster)} out of {len(faster) + len(slower)} times."
    )


def test_gotm(args, testcases: Mapping[str, str]):
    assert (
        args.gotm_setup is not None
    ), "You must specify --gotm_setup when testing GOTM"
    if args.fabm_ref is not None:
        print("Running in comparison mode.")
        compare(
            args.gotm_setup,
            args.work_root,
            testcases,
            args.cmake,
            cmake_arguments=args.cmake_arguments,
            fabm_ref_branch=args.fabm_ref,
        )
    else:
        test(
            args.gotm_setup,
            args.work_root,
            testcases,
            args.cmake,
            cmake_arguments=args.cmake_arguments,
            gotm_branch=args.gotm_branch
        )


def test_pyfabm(args, testcases: Mapping[str, str]):
    if not args.inplace:
        env_root = os.path.join(args.work_root, "python")
        print(f"Setting up virtual environment in {env_root}...")
        builder = venv.EnvBuilder(with_pip=True)
        builder.create(env_root)
        context = builder.ensure_directories(env_root)
        sys.stdout.flush()
        subprocess.check_call([context.env_exe, "-m", "pip", "install", "pyyaml"])
        return subprocess.call(
            [context.env_exe, os.path.abspath(sys.argv[0])]
            + sys.argv[1:]
            + ["--inplace"],
            cwd=args.work_root,
        )
    with open(os.path.join(FABM_BASE, "setup.cfg"), "w") as cfg:
        cfg.write("[build_ext]\n")
        cfg.write("debug=1\n")
        cfg.write("force=1\n")
        if len(args.cmake_arguments) > 0:
            cfg.write(f"cmake_opts={' '.join(args.cmake_arguments)}\n")
    if (
        run(
            "test/pyfabm/install",
            [sys.executable, "-m", "pip", "install", ".", "-v"],
            cwd=FABM_BASE,
        )
        != 0
    ):
        return
    with open(os.path.join(SCRIPT_ROOT, "environment.yaml")) as f:
        environment = yaml.safe_load(f)
    import pyfabm

    pyfabm.logger = logging.getLogger()
    dependency_names = set()
    print("Running FABM testcases with pyfabm:")
    for case, path in testcases.items():
        print(f"  {case}... ", end="")
        sys.stdout.flush()
        m0d = pyfabm.Model(path)
        counts = collections.Counter(v.name for v in m0d.variables).items()
        dup = [v for v, c in counts if c > 1]
        assert not dup, f"Duplicate variable names in 0D: {dup}"
        m1d = pyfabm.Model(path, shape=(5,))
        counts = collections.Counter(v.name for v in m1d.variables).items()
        dup = [v for v, c in counts if c > 1]
        assert not dup, f"Duplicate variable names in 1D: {dup}"
        for m in (m0d, m1d):
            m.cell_thickness = environment["cell_thickness"]
            for d in m.dependencies:
                dependency_names.add(d.name)
                if d.required:
                    d.value = environment[d.name]
            m.start()
        r0d = m0d.getRates(t=0.0)
        r1d = m1d.getRates(t=0.0)
        if (r1d != r1d[:, :1]).any():
            ran = r1d.max(axis=1) - r1d.min(axis=1)
            bad = {}
            for var, val in zip(m1d.state_variables, ran):
                if val != 0.0:
                    bad[var.name] = val
            assert False, f"Variability among 1D results: {r1d} (range: {bad})"
        if (r1d[:, 0] != r0d).any():
            diff = r1d[:, 0] - r0d
            bad = {}
            for var, val in zip(m1d.state_variables, diff):
                if val != 0.0:
                    bad[var.name] = val
            assert False, f"Mismatch between 0D and 1D results: {r0d} vs {r1d[:, 0]}. Difference: {diff}. {bad}"
        print("SUCCESS")
    pyfabm_libs = ", ".join([f"{n}={l._name}" for n, l in pyfabm.name2lib.items()])
    print(f"pyfabm {pyfabm.__version__} loaded from {pyfabm.__file__} ({pyfabm_libs})")
    try:
        pyfabm.unload()
    except Exception as e:
        print(f"Failed to unload pyfabm: {e}")
    if args.verbose:
        dependencies = "\n".join(sorted(dependency_names))
        print(f"Combined dependency list:\n{dependencies}")


def test_0d(args, testcases: Mapping[str, str], gotm_url=DEFAULT_GOTM_URL):
    build_dir = os.path.join(args.work_root, "build")
    gotm_dir = os.path.join(args.work_root, "code/gotm")
    run_dir = os.path.join(args.work_root, "run")
    shutil.copytree(os.path.join(FABM_BASE, "testcases/0d"), run_dir)
    with open(os.path.join(SCRIPT_ROOT, "environment.yaml")) as f:
        var2data = yaml.safe_load(f)
    with open(os.path.join(run_dir, "input.yaml"), "w") as f:
        yaml.dump(
            dict([(n, {"constant_value": v}) for n, v in var2data.items()]),
            f,
            default_flow_style=False,
        )
    git_clone("test/0d", gotm_url, gotm_dir)
    cmake(
        "test/0d",
        build_dir,
        os.path.join(FABM_BASE, "src/drivers/0d"),
        args.cmake,
        cmake_arguments=[f"-DGOTM_BASE={gotm_dir}"] + args.cmake_arguments,
    )
    exe = os.path.join(build_dir, "Debug/fabm0d.exe" if os.name == "nt" else "fabm0d")
    assert os.path.isfile(exe), f"{exe} not found"
    print("Running FABM testcases with 0d driver:")
    for case, path in testcases.items():
        shutil.copy(path, os.path.join(run_dir, "fabm.yaml"))
        run(f"test/0d/{case}", [exe], cwd=run_dir)


def test_harness(args, testcases: Mapping[str, str]):
    print("Running FABM testcases with testing harness:")
    for host in sorted(os.listdir(os.path.join(FABM_BASE, "src/drivers"))):
        print(f"  host {host}")
        build_dir = os.path.join(args.work_root, f"build_{host}")
        success = cmake(
            f"test_harness/{host}",
            build_dir,
            FABM_BASE,
            args.cmake,
            cmake_arguments=[f"-DFABM_HOST={host}", "-DCMAKE_BUILD_TYPE=debug"]
            + args.cmake_arguments,
            target="test_host",
        )
        if not success:
            continue
        exe = "test_host"
        if os.name == "nt":
            exe += ".exe"
            if not os.path.isfile(os.path.join(build_dir, exe)):
                exe = os.path.join("Debug", exe)
        for case, path in testcases.items():
            run(
                f"test_harness/{host}/{case}",
                [os.path.join(build_dir, exe), "--simulate", "-n", "10", path],
                cwd=SCRIPT_ROOT,
            )


def clean(workdir: str):
    print(f"Clean-up: deleting {workdir}")
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    host2function = {
        "gotm": test_gotm,
        "pyfabm": test_pyfabm,
        "0d": test_0d,
        "harness": test_harness,
    }

    import argparse

    parser = argparse.ArgumentParser(
        description="This script runs all FABM testcases in either pyfabm or GOTM."
    )
    parser.add_argument(
        "host", choices=host2function.keys(), help="Host to use for testing."
    )
    parser.add_argument(
        "--work_root", help="Path to use for code, testcases, results.", default=None
    )
    parser.add_argument(
        "--fabm_ref",
        help="Name of FABM branch/commit to compare results against.",
        default=None,
    )
    parser.add_argument(
        "--gotm_branch",
        help="Name of GOTM branch",
        default=None,
    )
    parser.add_argument(
        "--gotm_setup",
        help="Path to directory with GOTM setup (gotm.yaml, etc.)",
        default=None,
    )
    parser.add_argument("--cmake", help="path to cmake executable", default="cmake")
    parser.add_argument("--compiler", help="Fortran compiler executable")
    parser.add_argument(
        "--show_logs",
        action="store_true",
        help="Show contents of log files for failures at end of testing",
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Use current python environment instead of new virtual environment for pyfabm testing",
    )
    parser.add_argument(
        "--ext",
        nargs=2,
        action="append",
        help="Additional institute (name + dir) to include",
        default=[],
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable more detailed output"
    )
    args, cmake_arguments = parser.parse_known_args()
    if args.compiler is not None:
        cmake_arguments.append(f"-DCMAKE_Fortran_COMPILER={args.compiler}")

    # Add built-in test cases
    testcases = collections.OrderedDict()
    for path in sorted(glob.glob(os.path.join(FABM_BASE, "testcases/*.yaml"))):
        testcases[os.path.basename(path)[:-5]] = path

    # Add test cases from submodules (any file named fabm*.yaml)
    extern_dir = os.path.join(FABM_BASE, "extern")
    if os.path.isdir(extern_dir):
        for source in sorted(os.listdir(extern_dir)):
            extra_cases = glob.glob(os.path.join(extern_dir, source, "**/fabm*.yaml"))
            if extra_cases:
                print(f"Found {len(extra_cases)} additional test cases in extern/{source}")
            for path in sorted(extra_cases):
                testcases[f"{source}_{os.path.basename(path)[:-5]}"] = path

    # Add additional institutes specified on command line (--ext),
    # along with any contained test cases (any file named fabm*.yaml)
    if len(args.ext) > 0:
        cmake_arguments.append(
            f"-DFABM_EXTRA_INSTITUTES={';'.join([e[0] for e in args.ext])}"
        )
    for name, basedir in args.ext:
        cases = sorted(glob.glob(os.path.join(basedir, "**/fabm*.yaml")))
        print(
            f"Adding external sources of {name} from {basedir}"
            f" ({len(cases)} test cases)"
        )
        for path in cases:
            testcases[f"{name}/{os.path.basename(path)[:-5]}"] = path
        cmake_arguments.append(f"-DFABM_{name.upper()}_BASE={basedir}")

    args.cmake_arguments = cmake_arguments

    tmp = args.work_root is None
    if tmp:
        args.work_root = tempfile.mkdtemp()
        atexit.register(clean, args.work_root)
    args.work_root = os.path.abspath(args.work_root)
    print(f"Root of test directory: {args.work_root}")
    if len(args.cmake_arguments) > 0:
        print(f"Additional cmake arguments: {' '.join(args.cmake_arguments)}")

    logs = []
    retcode = host2function[args.host](args, testcases)
    if logs:
        log_list = "\n".join(logs)
        print(f"{len(logs)} ERRORS! Check the logs:\n{log_list}")
        if args.show_logs:
            print()
            for path in logs:
                print("=" * 80)
                print(path)
                print("=" * 80)
                with open(path) as f:
                    print(f.read())
                print()
        retcode = 1
    if not retcode:
        print("ALL TESTS SUCCEEDED")
    sys.exit(retcode)
