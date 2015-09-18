#!/usr/bin/env python

import sys

if len(sys.argv)!=2:
   print 'This script takes one argument: the path to a YAML file with FABM settings (typically fabm.yaml).'
   sys.exit(2)
yamlfile = sys.argv[1]

try:
   import pyfabm
except ImportError:
   print 'Unable to load pyfabm. Please build and install FABM with FABMHOST=python.'
   sys.exit(1)

# Create model object from YAML file.
model = pyfabm.Model(yamlfile)

print 'Interior state variables:'
for variable in model.bulk_state_variables:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)

print 'Surface-attached state variables:'
for variable in model.surface_state_variables:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)

print 'Bottom-attached state variables:'
for variable in model.bottom_state_variables:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)

print 'Interior diagnostic variables:'
for variable in model.bulk_diagnostic_variables:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)

print 'Horizontal diagnostic variables:'
for variable in model.horizontal_diagnostic_variables:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)

print 'Conserved quantities:'
for variable in model.conserved_quantities:
   print '  %s (%s)' % (variable.name,variable.units)

print 'Dependencies:'
for variable in model.dependencies:
   print '  %s = %s (%s)' % (variable.name,variable.long_name,variable.units)
