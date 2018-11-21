#!/usr/bin/env python

from __future__ import print_function
import os.path
import tempfile
import subprocess
import shutil

script_root = os.path.abspath(os.path.dirname(__file__))
cmake_arguments = []

generates = {}
builds = {}
tests = {}

root = os.path.join(script_root, '../..')
hosts = os.listdir(os.path.join(root, 'src/drivers'))

build_root = tempfile.mkdtemp()
try:
    for host in hosts:
        builddir = os.path.join(build_root, host)
        os.mkdir(builddir)
        generates[host] = subprocess.call(['cmake', os.path.join(root, 'src'), '-DFABM_HOST=%s' % host] + cmake_arguments, cwd=builddir)
        if generates[host] != 0:
            continue
        builds[host] = subprocess.call(['cmake', '--build', builddir, '--target', 'test_host'])
        if builds[host] != 0:
            continue
        if os.path.isfile(os.path.join(builddir, 'Debug/test_host.exe')):
            tests[host] = subprocess.call([os.path.join(builddir, 'Debug/test_host.exe')])
finally:
    shutil.rmtree(build_root)