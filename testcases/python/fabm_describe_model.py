#!/usr/bin/env python

import sys

yamlfile = '../fabm-npzd-carbonate.yaml'
if len(sys.argv)>1: yamlfile = sys.argv[1]

try:
   import pyfabm
except ImportError:
   print 'Unable to load pyfabm. Please build and install FABM with FABMHOST=python.'
   sys.exit(1)

# Create model object from YAML file.
model = pyfabm.Model(yamlfile)

# List pelagic state variables
for variable in model.bulk_state_variables:
   print variable.name,variable.units,variable.long_name
