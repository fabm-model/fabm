#!/usr/bin/env python

import os.path
import tempfile
import subprocess
import shutil
import sys
import glob

script_root = os.path.abspath(os.path.dirname(__file__))

def run(*args):
    returncode = subprocess.call(args)
    if returncode != 0:
        print('Command failed: %s' % (args,))
        sys.exit(1)

def git_clone(url,workdir,branch=None):
    run('git','clone',url,workdir)
    if branch is not None:
        olddir = os.getcwd()
        os.chdir(workdir)
        run('git','checkout',branch)
        os.chdir(olddir)

def build(work_root,fabm_url,gotm_url,fabm_branch=None,gotm_branch=None,cmake_arguments=()):
    # Save current working directory
    olddir = os.getcwd()

    # Change to test directory
    os.chdir(work_root)
    
    # Get latest FABM [public]
    git_clone(fabm_url,'code/fabm',fabm_branch)
    
    # Get latest GOTM [public]
    git_clone(gotm_url,'code/gotm',gotm_branch)
    
    # Build GOTM-FABM
    os.mkdir('build')
    os.chdir('build')
    run('cmake','../code/gotm/src','-DFABM_BASE=../code/fabm',*cmake_arguments)
    run('make')

    # Restore original working directory
    os.chdir(olddir)

def test(testcase_dir,work_root=None,cmake_arguments=(),fabm_url='git://git.code.sf.net/p/fabm/code',gotm_url='https://github.com/gotm-model/code.git',fabm_branch=None,gotm_branch=None):
    if work_root is None: work_root = tempfile.mkdtemp()
    work_root = os.path.abspath(work_root)
    print 'Root of test directory: %s' % work_root

    build(work_root,fabm_url,gotm_url,fabm_branch,gotm_branch,cmake_arguments)
    
    # Run L4 test case
    os.chdir(testcase_dir)
    with open('output.yaml','w') as f:
       f.write("""
result:
  time_unit: day
  time_step: 1
  time_method: mean
  variables:
    - source: fabm/*
    - source: h
""")

    for current_fabm_yaml in sorted(glob.glob(os.path.join(work_root,'code/fabm/testcases/*.yaml'))):
        print 'TESTING %s...' % os.path.basename(current_fabm_yaml),

        # Copy fabm.yaml
        shutil.copyfile(current_fabm_yaml, 'fabm.yaml') 

        # Run GOTM
        p = subprocess.Popen([os.path.join(work_root,'build/gotm')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        stdoutdata, stderrdata = p.communicate()
        if p.returncode!=0:
            print 'FAILED - last output:'
            last_lines = stdoutdata.rsplit('\n',10)[1:]
            print '\n'.join(last_lines)
        else:
            print 'ok'

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This script runs GOTM-FABM-ERSEM tests.')
    parser.add_argument('testcase_dir',help='Path to directory with GOTM test case')
    parser.add_argument('--work_root',help='Path to use for code, testcases, results.',default=None)
    parser.add_argument('--cmake_arg',action='append',dest='cmake_arguments',help='Additional arguments to pass to cmake.',default=[])
    args = parser.parse_args()
    test(args.testcase_dir,args.work_root,args.cmake_arguments)

