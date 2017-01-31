#!/usr/bin/env python

from __future__ import print_function
import os.path
import tempfile
import subprocess
import shutil
import sys
import glob
import timeit

script_root = os.path.abspath(os.path.dirname(__file__))

default_fabm_url = 'https://github.com/fabm-model/fabm.git'
default_gotm_url = 'https://github.com/gotm-model/code.git'

def run(*args):
    returncode = subprocess.call(args)
    if returncode != 0:
        print('Command failed: %s' % (args,))
        sys.exit(1)

def git_clone(url, workdir, branch=None):
    run('git', 'clone', url, workdir)
    if branch is not None:
        olddir = os.getcwd()
        os.chdir(workdir)
        run('git', 'checkout', branch)
        os.chdir(olddir)

def run_gotm(testcase_dir, gotm_exe):
    start = timeit.default_timer()
    p = subprocess.Popen([gotm_exe], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=testcase_dir)
    stdoutdata, stderrdata = p.communicate()
    duration = timeit.default_timer() - start
    if p.returncode != 0:
        print('FAILED - last output:')
        last_lines = stdoutdata.rsplit('\n', 10)[1:]
        print('\n'.join(last_lines))
    else:
        print('ok (%.3f s)' % duration)
    return p.returncode == 0, duration

def build(build_dir, fabm_base, gotm_base, cmake_arguments=()):
    # Save current working directory
    olddir = os.getcwd()

    # Create and change to build directory
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)
    os.mkdir(build_dir)
    os.chdir(build_dir)

    # Build GOTM-FABM
    run('cmake', os.path.join(gotm_base, 'src'), '-DFABM_BASE=%s' % fabm_base, *cmake_arguments)
    run('make', '-j8')

    # Restore original working directory
    os.chdir(olddir)

def enumerate_testcases(testcase_dir, fabm_testcases):
    with open(os.path.join(testcase_dir, 'output.yaml'), 'w') as f:
        f.write("""
result:
  time_unit: day
  time_step: 1
  time_method: mean
  variables:
    - source: fabm/*
    - source: h
""")

    for current_fabm_yaml in sorted(glob.glob(fabm_testcases)):
        # Copy fabm.yaml
        shutil.copyfile(current_fabm_yaml, os.path.join(testcase_dir, 'fabm.yaml'))

        yield os.path.basename(current_fabm_yaml)

def compare_netcdf(path, ref_path):
    import numpy
    import netCDF4
    perfect = True
    nc = netCDF4.Dataset(path)
    nc_ref = netCDF4.Dataset(ref_path)
    for varname in nc.variables.keys():
        if varname not in nc_ref.variables or varname in ('lon', 'lat', 'h', 'z', 'time'):
            continue
        ncvar = nc.variables[varname]
        ncvar_ref = nc_ref.variables[varname]
        dat = ncvar[...]
        valid = numpy.isfinite(dat)
        if not valid.all():
            print('    %s: %i of %i values are invalid' % (varname, valid.size - valid.sum(), valid.size))
            perfect = False
        else:
            delta = dat - ncvar_ref[...]
            maxdelta = numpy.abs(delta).max()
            perfect = perfect and maxdelta == 0.0
            print('    %s: max abs difference = %s' % (varname, maxdelta))
    nc.close()
    nc_ref.close()
    return perfect

def test(testcase_dir, work_root, cmake_arguments=(), fabm_url=default_fabm_url, gotm_url=default_gotm_url, fabm_branch=None, gotm_branch=None):
    fabm_base = os.path.join(work_root, 'code/fabm')
    gotm_base = os.path.join(work_root, 'code/gotm')
    build_dir = os.path.join(work_root, 'build')
    gotm_exe = os.path.join(build_dir, 'gotm')

    # Get latest FABM [public]
    git_clone(fabm_url, fabm_base, fabm_branch)

    # Get latest GOTM [public]
    git_clone(gotm_url, gotm_base, gotm_branch)

    build(build_dir, fabm_base, gotm_base, cmake_arguments)

    for name in enumerate_testcases(testcase_dir, os.path.join(fabm_base, 'testcases/*.yaml')):
        print('TESTING %s...' % name, end='')
        run_gotm(testcase_dir, gotm_exe)

def compare(testcase_dir, work_root=None, cmake_arguments=(), fabm_url=default_fabm_url, gotm_url=default_gotm_url, fabm_branch=None, gotm_branch=None, fabm_ref_branch=None, gotm_ref_branch=None):
    assert fabm_branch != fabm_ref_branch or gotm_branch != gotm_ref_branch
    git_clone(fabm_url, os.path.join(work_root, 'code/fabm'), fabm_branch)
    git_clone(gotm_url, os.path.join(work_root, 'code/gotm'), gotm_branch)
    build(os.path.join(work_root, 'build'), os.path.join(work_root, 'code/fabm'), os.path.join(work_root, 'code/gotm'), cmake_arguments)

    git_clone(fabm_url, os.path.join(work_root, 'ref/code/fabm'), fabm_ref_branch)
    git_clone(gotm_url, os.path.join(work_root, 'ref/code/gotm'), gotm_ref_branch)
    build(os.path.join(work_root, 'ref/build'), os.path.join(work_root, 'ref/code/fabm'), os.path.join(work_root, 'ref/code/gotm'), cmake_arguments)

    faster, slower = [], []
    failed, success, crashed = [], [], []
    for name in enumerate_testcases(testcase_dir, os.path.join(work_root, 'code/fabm/testcases/*.yaml')):
        print('TESTING %s...' % name)
        print('  reference...', end='')
        valid_ref, duration_ref = run_gotm(testcase_dir, os.path.join(work_root, 'ref/build/gotm'))
        os.rename(os.path.join(testcase_dir, 'result.nc'), os.path.join(testcase_dir, 'result_ref.nc'))
        print('  target...', end='')
        valid, duration = run_gotm(testcase_dir, os.path.join(work_root, 'build/gotm'))
        if valid and valid_ref:
            if compare_netcdf(os.path.join(testcase_dir, 'result_ref.nc'), os.path.join(testcase_dir, 'result.nc')):
                success.append(name)
            else:
                failed.append(name)
        else:
            crashed.append(name)
        if duration < duration_ref:
            faster.append(name)
        else:
            slower.append(name)

    print('%i perfect matches: %s' % (len(success), ', '.join(success)))
    print('%i mismatches: %s' % (len(failed), ', '.join(failed)))
    print('%i failed to run: %s' % (len(crashed), ', '.join(crashed)))
    print('Faster than reference? %i out of %i times.' % (len(faster), len(faster) + len(slower)))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script runs GOTM-FABM-ERSEM tests.')
    parser.add_argument('scenario_dir', help='Path to directory with GOTM scenario (gotmrun.nml, etc.)')
    parser.add_argument('--work_root', help='Path to use for code, testcases, results.', default=None)
    parser.add_argument('--fabm_branch', help='Name of FABM branch to test.', default=None)
    parser.add_argument('--fabm_ref_branch', help='Name of FABM branch to compare against.', default=None)
    parser.add_argument('cmake_arguments', nargs=argparse.REMAINDER, help='Additional arguments to pass to cmake.', default=[])
    args = parser.parse_args()

    if args.work_root is None:
        args.work_root = tempfile.mkdtemp()
    args.work_root = os.path.abspath(args.work_root)
    print('Root of test directory: %s' % args.work_root)

    if args.fabm_ref_branch is not None:
        print('Running in comparison mode.')
        compare(args.scenario_dir, args.work_root, args.cmake_arguments, fabm_branch=args.fabm_branch, fabm_ref_branch=args.fabm_ref_branch)
    else:
        test(args.scenario_dir, args.work_root, args.cmake_arguments, fabm_branch=args.fabm_branch)

