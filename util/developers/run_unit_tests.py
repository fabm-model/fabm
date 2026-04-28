#!/usr/bin/env python

import sys
import tempfile
import subprocess
import shutil
import argparse
import timeit
import errno
from pathlib import Path
from typing import Dict, List

script_root = Path(__file__).parent
root = script_root.parent.parent
allowed_hosts = sorted(p.name for p in (root / "src/drivers").iterdir())

parser = argparse.ArgumentParser()
parser.add_argument(
    "--host",
    action="append",
    dest="hosts",
    choices=allowed_hosts,
    help="host to test (may appear multiple times)",
)
parser.add_argument("--cmake", default="cmake", help="path to cmake executable")
parser.add_argument("--compiler", help="Fortran compiler executable")
parser.add_argument(
    "--performance",
    action="store_true",
    help="test performance with a specific model and environment (see --config and --env)",
)
parser.add_argument(
    "--config",
    default="fabm.yaml",
    help="model configuration for performance testing, default: fabm.yaml",
    type=Path,
)
parser.add_argument(
    "--env",
    default="environment.yaml",
    help="model environment for performance testing (YAML file containing a dictionary with variable: value combinations), default: environment.yaml",
    type=Path,
)
parser.add_argument(
    "--report",
    default=None,
    help="file to write performance report to (only used with --performance), default: performance_<BRANCH>_<COMMIT>.log",
    type=Path,
)
parser.add_argument(
    "--repeat",
    type=int,
    default=5,
    help="number of times to run each performance test. Increase this to reduce the noise in timings",
)
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help="show test results even if completed successfully",
)
args, cmake_arguments = parser.parse_known_args()
if args.performance:
    if not args.config.is_file():
        print(
            f"Model configuration {args.config} does not exist. Specify (or change) --config."
        )
        sys.exit(2)
    if not args.env.is_file():
        print(
            f"Model environment {args.env} does not exist. Specify (or change) --env."
        )
        sys.exit(2)
    if args.report is None:
        git_branch = (
            subprocess.check_output(["git", "name-rev", "--name-only", "HEAD"])
            .decode("ascii")
            .strip()
        )
        git_commit = (
            subprocess.check_output(["git", "describe", "--always", "--dirty"])
            .decode("ascii")
            .strip()
        )
        args.report = Path(f"performance_{git_branch}_{git_commit}.log")
    print(f"Performance report will be written to {args.report}")

if args.compiler is not None:
    cmake_arguments.append(f"-DCMAKE_Fortran_COMPILER={args.compiler}")

generates: Dict[str, int] = {}
builds: Dict[str, int] = {}
tests: Dict[str, int] = {}

if not args.hosts:
    args.hosts = allowed_hosts
print(f"Selected hosts: {', '.join(args.hosts)}")

logs: List[Path] = []


def run(phase: str, args: List[str], verbose: bool = False, **kwargs) -> int:
    proc = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        **kwargs,
    )
    stdoutdata, _ = proc.communicate()
    if proc.returncode != 0:
        log_path = Path(f"{phase}.log")
        with log_path.open("w") as f:
            f.write(stdoutdata)
        logs.append(log_path)
        print(f"FAILED (return code {proc.returncode}, log written to {log_path})")
    else:
        print("SUCCESS")
    if verbose:
        print(f"Output:\n{80 * '-'}\n{stdoutdata}\n{80 * '-'}")
    return proc.returncode


build_root = Path(tempfile.mkdtemp())
try:
    vsconfig = "Release" if args.performance else "Debug"
    host2exe: Dict[str, Path] = {}
    for host in args.hosts:
        print(host)
        build_dir = build_root / host
        build_dir.mkdir()
        print("  generating...", end="")
        sys.stdout.flush()
        try:
            generates[host] = run(
                f"{host}_generate",
                [args.cmake, str(root), f"-DFABM_HOST={host}"] + cmake_arguments,
                cwd=build_dir,
            )
        except EnvironmentError as e:
            if e.errno != errno.ENOENT:
                raise
            print(
                "\n\ncmake executable not found. Specify its location on the command line with --cmake."
            )
            sys.exit(2)
        if generates[host] != 0:
            continue
        print("  building...", end="")
        sys.stdout.flush()
        builds[host] = run(
            f"{host}_build",
            [
                args.cmake,
                "--build",
                str(build_dir),
                "--target",
                "test_host",
                "--config",
                vsconfig,
            ],
        )
        if builds[host] != 0:
            continue
        print("  testing...", end="")
        sys.stdout.flush()
        for exepath in (
            build_dir / vsconfig / "test_host.exe",
            build_dir / "test_host",
        ):
            if exepath.is_file():
                break
        else:
            raise Exception(
                f"Could not find test executable for host {host} in expected locations"
            )
        host2exe[host] = exepath
        tests[host] = run(
            f"{host}_test",
            [str(host2exe[host]), "-n", f"{1 if args.performance else args.repeat}"],
            verbose=args.verbose,
        )

    if args.performance:
        print("Measuring runtime")
        timings: Dict[str, List[float]] = {}
        shutil.copy(args.config, build_root / "fabm.yaml")
        shutil.copy(args.env, build_root / "environment.yaml")
        for i in range(args.repeat):
            print(f"  replicate {i}")
            for host in args.hosts:
                if tests.get(host, 1) != 0:
                    continue
                start = timeit.default_timer()
                print(f"    {host}...", end="")
                run(
                    f"{host}_perfrun_{i}",
                    [str(host2exe[host]), "--simulate"],
                    cwd=build_root,
                )
                timings.setdefault(host, []).append(timeit.default_timer() - start)

finally:
    shutil.rmtree(build_root)

if logs:
    print(f"All tests complete - {len(logs)} FAILED")
    print(f"See the following log files:\n{'\n'.join(map(str, logs))}")
else:
    print("All tests complete - no failures")

if args.performance:
    print("Timings:")
    for host in args.hosts:
        ts = timings.get(host, ())
        timing = "NA" if not ts else f"{sum(ts) / len(ts):.3f} s"
        print(f"  {host}: {timing}")
    with args.report.open("w") as f:
        header = (
            ["host"] + [f"run {i} (s)" for i in range(args.repeat)] + ["average (s)"]
        )
        f.write("\t".join(header) + "\n")
        for host in args.hosts:
            if host in timings:
                ts = timings[host]
                items = [host] + [f"{t:.3f}" for t in ts] + [f"{sum(ts) / len(ts):.3f}"]
                f.write("\t".join(items) + "\n")
