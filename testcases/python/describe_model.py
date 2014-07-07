#!/usr/bin/env python

yamlfile = '../fabm-npzd-carbonate.yaml'

import sys, os

# Try to locate FABM library.
# Retrieve FABM installation directory from enviornemtn variable FABM_PREFIX if set.
# Otherwise, try platform-specific default installation locations.
if 'FABM_PREFIX' in os.environ:
    prefix = os.environ['FABM_PREFIX']
elif os.name=='nt':
    prefix = '%s%s' % (os.environ['HOMEDRIVE'],os.environ['HOMEPATH'])
else:
    prefix = os.path.join(os.environ['HOME'],'local')
sys.path.append(os.path.join(prefix,'fabm/python/python'))
import pyfabm

# Create model object from YAML file.
model = pyfabm.Model(os.path.join(os.path.dirname(__file__),yamlfile))

# List pelagic state variables
for variable in model.bulk_state_variables:
   print variable.name,variable.units,variable.long_name
