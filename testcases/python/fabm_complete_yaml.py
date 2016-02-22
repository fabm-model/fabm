#!/usr/bin/env python

import sys

try:
   import pyfabm.complete_yaml
except ImportError:
   print 'Unable to load pyfabm. Please build and install python-fabm by running cmake+make for source directory $FABMDIR/src/drivers/python.'
   sys.exit(1)

if len(sys.argv)<2:
   print 'This script takes at least one argument: the path to the YAML file that needs to be completed. Optionally, the path to the output file can be given as second argument. By default, the input file is modified in-place.'
   sys.exit(2)
infile = sys.argv[1]
outfile = infile
if len(sys.argv)>2: outfile = sys.argv[2]
pyfabm.complete_yaml.processFile(infile,outfile)
