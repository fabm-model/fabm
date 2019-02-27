#!/usr/bin/env python

"""
This script takes a model configuration file and rewrites it with standardized indentation and comments describing each of the model\'s variables and parameters.
"""

from __future__ import print_function

import sys

try:
   import pyfabm.complete_yaml
except ImportError:
   print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
   sys.exit(1)

def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',help='Path to a YAML file with the model configuration that needs to be completed.')
    parser.add_argument('output_path',help='Path to save the completed YAML file to. If not provided, this defaults to the file from which the model configuration is read.',nargs='?')
    parser.add_argument('--add_missing',action='store_true',default=False,help='Whether to add missing parameter that have a default set in the code.')
    args = parser.parse_args()

    if args.output_path is None: args.output_path = args.path

    pyfabm.complete_yaml.processFile(args.path,args.output_path,add_missing=args.add_missing)

if __name__ == "__main__":
    # execute only if run as a script
    main()

