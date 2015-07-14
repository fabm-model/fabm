#!/usr/bin/env python

import sys

try:
   import pyfabm.complete_yaml
except ImportError:
   print 'Unable to load pyfabm. Please build and install python-fabm by running cmake+make for source directory $FABMDIR/src/drivers/python.'
   sys.exit(1)

if len(sys.argv)!=3:
   print 'This script takes two argument: the path to the input YAML file and the path to the output YAML file.'
   sys.exit(2)
infile = sys.argv[1]
outfile = sys.argv[2]
pyfabm.complete_yaml.processFile(infile,outfile)
