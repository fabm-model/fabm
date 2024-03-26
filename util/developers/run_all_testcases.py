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

script_root = os.path.abspath(os.path.dirname(__file__))
fabm_base = os.path.join(script_root, "../..")

default_fabm_url = "https://github.com/fabm-model/fabm.git"
default_gotm_url = "https://github.com/gotm-model/code.git"


def run(phase, args, verbose=False, **kwargs):
    print("%s%s... " % ("  " * phase.count("/"), phase.rsplit("/", 1)[-1]), end="")
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
        log_path = "%s.log" % phase.replace("/", "_")
        with open(log_path, "w") as f:
            f.write(stdoutdata)
        logs.append(log_path)
        print(
            "FAILED (return code %i, log written to %s)" % (proc.returncode, log_path)
        )
    else:
        print("SUCCESS")
    if verbose:
        print("Output:\n%s\n%s\n%s" % (80 * "-", stdoutdata, 80 * "-"))
    return proc.returncode


def git_clone(phase, url, workdir, branch=None):
    run("%s/clone" % phase, ["git", "clone", url, workdir])
    olddir = os.getcwd()
    if branch is not None:
        run("%s/checkout" % phase, ["git", "checkout", branch], cwd=workdir)
    run(
        "%s/submodule" % phase,
        ["git", "submodule", "update", "--init", "--recursive"],
        cwd=workdir,
    )


def run_gotm(setup_dir, gotm_exe):
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
        print("ok (%.3f s)" % duration)
    return p.returncode == 0, duration


def cmake(
    phase, build_dir, source_dir, cmake_path="cmake", target=None, cmake_arguments=[]
):
    # Create and change to build directory
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)
    os.mkdir(build_dir)

    if os.name == "nt":
        x64 = sys.maxsize > 2**32
        cmake_arguments = ["-A", "x64" if x64 else "Win32"] + cmake_arguments

    if os.path.isfile(os.path.join(source_dir, "CMakeCache.txt")):
        print(
            '\n\nCannot use "%s" as source directory for cmake because it has previously been used as build directory.\nPlease delete "%s" before continuing.'
            % (source_dir, os.path.abspath(os.path.join(source_dir, "CMakeCache.txt")))
        )
        sys.exit(1)

    # Build
    try:
        ret = run(
            "%s/configure" % phase,
            [cmake_path, source_dir] + cmake_arguments,
            cwd=build_dir,
        )
    except EnvironmentError as e:
        if e.errno != errno.ENOENT:
            raise
        print(
            '\n\ncmake executable ("%s") not found. Specify its location on the command line with --cmake.'
            % cmake_path
        )
        sys.exit(2)

    if ret == 0:
        args = ["--config", "Debug"]
        if target is not None:
            args = args + ["--target", target]
        ret = run(
            "%s/build" % phase, [cmake_path, "--build", "."] + args, cwd=build_dir
        )

    return ret == 0


def compare_netcdf(path, ref_path):
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
            print(
                "    %s: %i of %i values are invalid"
                % (varname, valid.size - valid.sum(), valid.size)
            )
            perfect = False
        else:
            delta = dat - ncvar_ref[...]
            maxdelta = numpy.abs(delta).max()
            perfect = perfect and maxdelta == 0.0
            print("    %s: max abs difference = %s" % (varname, maxdelta))
    nc.close()
    nc_ref.close()
    return perfect


def test(
    gotm_setup_dir,
    work_root,
    testcases,
    cmake_path="cmake",
    cmake_arguments=[],
    fabm_url=default_fabm_url,
    gotm_url=default_gotm_url,
    fabm_branch=None,
    gotm_branch=None,
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
        cmake_arguments=["-DFABM_BASE=%s" % fabm_base] + cmake_arguments,
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
    with open(os.path.join(gotm_setup_dir, "gotm.yaml"), "w") as f:
        yaml.dump(gotm_yaml, f, default_flow_style=False)
    for name, path in testcases.items():
        shutil.copyfile(path, os.path.join(gotm_setup_dir, "fabm.yaml"))
        run("test/gotm/%s" % name, [exe], cwd=gotm_setup_dir)


def compare(
    gotm_setup_dir,
    work_root,
    testcases,
    cmake_path="cmake",
    cmake_arguments=[],
    fabm_url=default_fabm_url,
    gotm_url=default_gotm_url,
    fabm_branch=None,
    gotm_branch=None,
    fabm_ref_branch=None,
    gotm_ref_branch=None,
):
    assert fabm_branch != fabm_ref_branch or gotm_branch != gotm_ref_branch
    git_clone(fabm_url, os.path.join(work_root, "code/fabm"), fabm_branch)
    git_clone(gotm_url, os.path.join(work_root, "code/gotm"), gotm_branch)
    cmake(
        "test_gotm",
        os.path.join(work_root, "build"),
        os.path.join(work_root, "code/gotm"),
        cmake_path,
        "-DFABM_BASE=%s" % os.path.join(work_root, "code/fabm"),
        cmake_arguments=cmake_arguments,
    )

    git_clone(fabm_url, os.path.join(work_root, "ref/code/fabm"), fabm_ref_branch)
    git_clone(gotm_url, os.path.join(work_root, "ref/code/gotm"), gotm_ref_branch)
    cmake(
        "ref_gotm",
        os.path.join(work_root, "ref/build"),
        os.path.join(work_root, "ref/code/gotm"),
        cmake_path,
        "-DFABM_BASE=%s" % os.path.join(work_root, "ref/code/fabm"),
        cmake_arguments=cmake_arguments,
    )

    faster, slower = [], []
    failed, success, crashed = [], [], []
    for name in enumerate_testcases(
        testcase_dir, os.path.join(work_root, "code/fabm/testcases/*.yaml")
    ):
        print("TESTING %s..." % name)
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

    print("%i perfect matches: %s" % (len(success), ", ".join(success)))
    print("%i mismatches: %s" % (len(failed), ", ".join(failed)))
    print("%i failed to run: %s" % (len(crashed), ", ".join(crashed)))
    print(
        "Faster than reference? %i out of %i times."
        % (len(faster), len(faster) + len(slower))
    )


def test_gotm(args, testcases):
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
        )


def test_pyfabm(args, testcases):
    if not args.inplace:
        env_root = os.path.join(args.work_root, "python")
        print("Setting up virtual environment in %s..." % env_root)
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
    with open(os.path.join(fabm_base, "setup.cfg"), "w") as cfg:
        cfg.write("[build_ext]\n")
        cfg.write("debug=1\n")
        cfg.write("force=1\n")
        if len(args.cmake_arguments) > 0:
            cfg.write("cmake_opts=%s\n" % " ".join(args.cmake_arguments))
    if (
        run(
            "test/pyfabm/install",
            [sys.executable, "-m", "pip", "install", "."],
            cwd=fabm_base,
        )
        != 0
    ):
        return
    with open(os.path.join(script_root, "environment.yaml")) as f:
        environment = yaml.safe_load(f)
    import pyfabm

    pyfabm.logger = logging.getLogger()
    dependency_names = set()
    print("Running FABM testcases with pyfabm:")
    for case, path in testcases.items():
        print("  %s... " % case, end="")
        sys.stdout.flush()
        m0d = pyfabm.Model(path)
        m1d = pyfabm.Model(path, shape=(5,))
        for m in (m0d, m1d):
            m.cell_thickness = environment["cell_thickness"]
            for d in m.dependencies:
                dependency_names.add(d.name)
                if d.required:
                    d.value = environment[d.name]
            m.start()
        r0d = m0d.getRates()
        r1d = m1d.getRates()
        if (r1d != r1d[:, :1]).any():
            ran = r1d.max(axis=1) - r1d.min(axis=1)
            bad = {}
            for var, val in zip(m1d.state_variables, ran):
                if val != 0.0:
                    bad[var.name] = val
            assert False, "Variability among 1D results: %s (range: %s)" % (r1d, bad)
        assert (
            r1d[:, 0] == r0d
        ).all(), "Mismatch between 0D and 1D results: %s vs %s. Difference: %s" % (
            r0d,
            r1d[:, 0],
            r1d[:, 0] - r0d,
        )
        print("SUCCESS")
    print(
        "pyfabm %s loaded from %s (%s)"
        % (
            pyfabm.get_version(),
            pyfabm.__file__,
            ", ".join(["%s=%s" % (n, l._name) for n, l in pyfabm.name2lib.items()]),
        )
    )
    try:
        pyfabm.unload()
    except Exception as e:
        print("Failed to unload pyfabm: %s" % e)
    if args.verbose:
        print("Combined dependency list:\n%s" % "\n".join(sorted(dependency_names)))


def test_0d(args, testcases, gotm_url=default_gotm_url):
    build_dir = os.path.join(args.work_root, "build")
    gotm_dir = os.path.join(args.work_root, "code/gotm")
    run_dir = os.path.join(args.work_root, "run")
    shutil.copytree(os.path.join(fabm_base, "testcases/0d"), run_dir)
    with open(os.path.join(script_root, "environment.yaml")) as f:
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
        os.path.join(fabm_base, "src/drivers/0d"),
        args.cmake,
        cmake_arguments=["-DGOTM_BASE=%s" % gotm_dir] + args.cmake_arguments,
    )
    exe = os.path.join(build_dir, "Debug/fabm0d.exe" if os.name == "nt" else "fabm0d")
    assert os.path.isfile(exe), "%s not found" % exe
    print("Running FABM testcases with 0d driver:")
    for case, path in testcases.items():
        shutil.copy(path, os.path.join(run_dir, "fabm.yaml"))
        run("test/0d/%s" % case, [exe], cwd=run_dir)


def test_harness(args, testcases):
    run_dir = os.path.join(args.work_root, "run")
    os.mkdir(run_dir)
    shutil.copy(os.path.join(script_root, "environment.yaml"), run_dir)
    print("Running FABM testcases with testing harness:")
    for host in sorted(os.listdir(os.path.join(fabm_base, "src/drivers"))):
        print("  host %s" % host)
        build_dir = os.path.join(args.work_root, "build_%s" % host)
        success = cmake(
            "test_harness/%s" % host,
            build_dir,
            fabm_base,
            args.cmake,
            cmake_arguments=["-DFABM_HOST=%s" % host, "-DCMAKE_BUILD_TYPE=debug"]
            + args.cmake_arguments,
            target="test_host",
        )
        if not success:
            continue
        exe = os.path.join(
            build_dir, "Debug/test_host.exe" if os.name == "nt" else "test_host"
        )
        for case, path in testcases.items():
            shutil.copy(path, os.path.join(run_dir, "fabm.yaml"))
            run(
                "test_harness/%s/%s" % (host, case),
                [exe, "--simulate", "-n", "10"],
                cwd=run_dir,
            )


def clean(workdir):
    print("Clean-up: deleting %s" % workdir)
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
        cmake_arguments.append("-DCMAKE_Fortran_COMPILER=%s" % args.compiler)

    testcases = collections.OrderedDict()
    for path in sorted(glob.glob(os.path.join(fabm_base, "testcases/*.yaml"))):
        testcases[os.path.basename(path)[:-5]] = path

    for path in sorted(glob.glob(os.path.join(fabm_base, "extern/*/testcases/*.yaml"))):
        relpath = os.path.relpath(path, start=os.path.join(fabm_base, "extern"))
        source = os.path.dirname(os.path.dirname(relpath))
        testcases[source + "_" + os.path.basename(path)] = path

    for name, basedir in args.ext:
        cases = sorted(glob.glob(os.path.join(basedir, "testcases/*.yaml")))
        print(
            "Adding external sources of %s from %s (%i test cases)"
            % (name, basedir, len(cases))
        )
        for path in cases:
            testcases["%s/%s" % (name, os.path.basename(path)[:-5])] = path
        cmake_arguments.append("-DFABM_%s_BASE=%s" % (name.upper(), basedir))

    if len(args.ext) > 0:
        cmake_arguments.append(
            "-DFABM_EXTRA_INSTITUTES=%s" % ";".join([e[0] for e in args.ext])
        )

    args.cmake_arguments = cmake_arguments

    tmp = args.work_root is None
    if tmp:
        args.work_root = tempfile.mkdtemp()
        atexit.register(clean, args.work_root)
    args.work_root = os.path.abspath(args.work_root)
    print("Root of test directory: %s" % args.work_root)
    if len(args.cmake_arguments) > 0:
        print("Additional cmake arguments: %s" % " ".join(args.cmake_arguments))

    logs = []
    retcode = host2function[args.host](args, testcases)
    if logs:
        print("%i ERRORS! Check the logs:\n%s" % (len(logs), "\n".join(logs)))
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
