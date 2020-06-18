#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import tempfile
import subprocess
import shutil
import argparse
import timeit
import io

script_root = os.path.abspath(os.path.dirname(__file__))
root = os.path.join(script_root, '../..')
allowed_hosts = sorted(os.listdir(os.path.join(root, 'src/drivers')))

parser = argparse.ArgumentParser()
parser.add_argument('--host', action='append', dest='hosts', choices=allowed_hosts, help='host to test (may appear multiple times)')
parser.add_argument('--cmake', default='cmake', help='path to cmake executable')
parser.add_argument('--compiler', help='Fortran compiler executable')
parser.add_argument('--performance', action='store_true', help='test performance with a specific model and environment (see --config and --env)')
parser.add_argument('--config', default='fabm.yaml', help='model configuration for performance testing, default: fabm.yaml')
parser.add_argument('--env', default='environment.yaml', help='model environment for performance testing (YAML file containing a dictionary with variable: value combinations), default: environment.yaml')
parser.add_argument('--report', default=None, help='file to write performance report to (only used with --performance), default: performance_<BRANCH>_<COMMIT>.log')
parser.add_argument('--repeat', type=int, default=5, help='number of times to run each performance test. Increase this to reduce the noise in timings')
parser.add_argument('-v', '--verbose', action='store_true', help='show test results even if completed successfully')
args, cmake_arguments = parser.parse_known_args()
if args.performance:
    if not os.path.isfile(args.config):
        print('Model configuration %s does not exist. Specify (or change) --config.' % args.config)
        sys.exit(2)
    if not os.path.isfile(args.env):
        print('Model environment %s does not exist. Specify (or change) --env.' % args.env)
        sys.exit(2)
    if args.report is None:
        git_branch = subprocess.check_output(['git', 'name-rev', '--name-only', 'HEAD']).decode('ascii').strip()
        git_commit = subprocess.check_output(['git', 'describe', '--always', '--dirty']).decode('ascii').strip()
        args.report = 'performance_%s_%s.log' % (git_branch, git_commit)
    print('Performance report will be written to %s' % args.report)

if args.compiler is not None:
    cmake_arguments.append('-DCMAKE_Fortran_COMPILER=%s' % args.compiler)

generates = {}
builds = {}
tests = {}

if not args.hosts:
    args.hosts = allowed_hosts
print('Selected hosts: %s' % ', '.join(args.hosts))

logs = []
def run(phase, args, verbose=False, **kwargs):
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, **kwargs)
    stdoutdata, _ = proc.communicate()
    if proc.returncode != 0:
        log_path = '%s.log' % phase
        with io.open(log_path, 'w') as f:
            f.write(stdoutdata)
        logs.append(log_path)
        print('FAILED (return code %i, log written to %s)' % (proc.returncode, log_path))
    else:
        print('SUCCESS')
    if verbose:
        print('Output:\n%s\n%s\n%s' % (80 * '-', stdoutdata, 80 * '-'))
    return proc.returncode

build_root = tempfile.mkdtemp()
try:
    vsconfig = 'Release' if args.performance else 'Debug'
    host2exe = {}
    for host in args.hosts:
        print(host)
        build_dir = os.path.join(build_root, host)
        os.mkdir(build_dir)
        print('  generating...', end='', flush=True)
        try:
            generates[host] = run('%s_generate' % host, [args.cmake, os.path.join(root, 'src'), '-DFABM_HOST=%s' % host] + cmake_arguments, cwd=build_dir)
        except FileNotFoundError:
            print('\n\ncmake executable not found. Specify its location on the command line with --cmake.')
            sys.exit(2)
        if generates[host] != 0:
            continue
        print('  building...', end='', flush=True)
        builds[host] = run('%s_build' % host, [args.cmake, '--build', build_dir, '--target', 'test_host', '--config', vsconfig])
        if builds[host] != 0:
            continue
        print('  testing...', end='', flush=True)
        for exename in ('%s/test_host.exe' % vsconfig, 'test_host'):
            exepath = os.path.join(build_dir, exename)
            if os.path.isfile(exepath):
                host2exe[host] = exepath
        tests[host] = run('%s_test' % host, [host2exe[host], '-n', '%i' % (1 if args.performance else args.repeat)], verbose=args.verbose)

    if args.performance:
        print('Measuring runtime')
        timings = {}
        shutil.copy(args.config, os.path.join(build_root, 'fabm.yaml'))
        shutil.copy(args.env, os.path.join(build_root, 'environment.yaml'))
        for i in range(args.repeat):
            print('  replicate %i' % i)
            for host in args.hosts:
                if tests.get(host, 1) != 0:
                    continue
                start = timeit.default_timer()
                print('    %s...' % (host,), end='')
                run('%s_perfrun_%i' % (host, i), [host2exe[host], '--simulate'], cwd=build_root)
                timings.setdefault(host, []).append(timeit.default_timer() - start)

finally:
    shutil.rmtree(build_root)

if logs:
    print('All tests complete - %i FAILED' % len(logs))
    print('See the following log files:\n%s' % '\n'.join(logs))
else:
    print('All tests complete - no failures')

if args.performance:
    print('Timings:')
    for host in args.hosts:
        ts = timings.get(host, ())
        timing = 'NA' if not ts else '%.3f s' % (sum(ts) / len(ts))
        print('  %s: %s' % (host, timing))
    with open(args.report, 'w') as f:
        f.write('host\t%s\taverage (s)\n' % '\t'.join(['run %i (s)' % i for i in range(args.repeat)]))
        for host in args.hosts:
            if host in timings:
                ts = timings[host]
                f.write('%s\t%s\t%.3f\n' % (host, '\t'.join(['%.3f' % t for t in ts]), sum(ts) / len(ts)))
