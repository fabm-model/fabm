import os.path
import glob
import argparse

import pyfabm.complete_yaml

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true', help='Show process id and wait before accessign FABM so a debugger can be attached.')
args = parser.parse_args()

if args.debug:
    print('Process id: %i. Press return to continue...' % os.getpid())
    input()

for path in glob.glob(os.path.join(os.path.dirname(__file__), '../../testcases/*.yaml')):
    print('Processing %s...' % path)
    pyfabm.complete_yaml.processFile(path, path, add_missing=False)