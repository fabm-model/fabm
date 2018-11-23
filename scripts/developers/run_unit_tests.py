#!/usr/bin/env python

from __future__ import print_function
import os.path
import tempfile
import subprocess
import shutil
import argparse

script_root = os.path.abspath(os.path.dirname(__file__))
root = os.path.join(script_root, '../..')
allowed_hosts = os.listdir(os.path.join(root, 'src/drivers'))

parser = argparse.ArgumentParser()
parser.add_argument('--host', action='append', dest='hosts', choices=allowed_hosts)
args = parser.parse_args()

cmake_arguments = []

generates = {}
builds = {}
tests = {}

if not args.hosts:
    args.hosts = allowed_hosts
print('Selected hosts: %s' % ', '.join(args.hosts))

logs = []
def run(phase, args, **kwargs):
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)
    stdoutdata, stderrdata = proc.communicate()
    if proc.returncode != 0:
        log_path = '%s.log' % phase
        with open(log_path, 'w') as f:
            f.write(stdoutdata)
        logs.append(log_path)
        print('FAILED (log written to %s)' % log_path)
    else:
        print('SUCCESS')
    return proc.returncode

build_root = tempfile.mkdtemp()
try:
    for host in args.hosts:
        print(host)
        build_dir = os.path.join(build_root, host)
        os.mkdir(build_dir)
        print('  generating...', end='')
        generates[host] = run('%s_generate' % host, ['cmake', os.path.join(root, 'src'), '-DFABM_HOST=%s' % host] + cmake_arguments, cwd=build_dir)
        if generates[host] != 0:
            continue
        print('  building...', end='')
        builds[host] = run('%s_build' % host, ['cmake', '--build', build_dir, '--target', 'test_host'])
        if builds[host] != 0:
            continue
        print('  testing...', end='')
        if os.path.isfile(os.path.join(build_dir, 'Debug/test_host.exe')):
            tests[host] = run('%s_test' % host, [os.path.join(build_dir, 'Debug/test_host.exe')])
        else:
            tests[host] = run('%s_test' % host, [os.path.join(build_dir, 'test_host')])
finally:
    shutil.rmtree(build_root)

if logs:
    print('All tests complete - %i FAILED' % len(logs))
    print('See the following log files:\n%s' % '\n'.join(logs))
else:
    print('All tests complete - no failures')
