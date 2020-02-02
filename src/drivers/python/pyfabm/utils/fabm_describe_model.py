#!/usr/bin/env python

"""
This script lists all state variables, diagnostic variables, conserved quantities and environmental dependencies of a biogeochemical model.
"""

from __future__ import print_function
import sys

try:
   import pyfabm
except ImportError:
   print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
   sys.exit(1)

def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', help='Path to a YAML file with the model configuration (typically fabm.yaml)', nargs='?', default='fabm.yaml')
    parser.add_argument('--all', '-a', action='store_true', help='Show diagnostics that are by default excluded from output')
    args = parser.parse_args()

    # Create model object from YAML file.
    model = pyfabm.Model(args.path)

    print('Interior state variables:')
    for variable in model.interior_state_variables:
       print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

    print('Surface-attached state variables:')
    for variable in model.surface_state_variables:
       print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

    print('Bottom-attached state variables:')
    for variable in model.bottom_state_variables:
       print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

    print('Interior diagnostic variables:')
    for variable in model.interior_diagnostic_variables:
       if variable.output or args.all:
          print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

    print('Horizontal diagnostic variables:')
    for variable in model.horizontal_diagnostic_variables:
       if variable.output or args.all:
          print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

    print('Conserved quantities:')
    for variable in model.conserved_quantities:
       print('  %s (%s)' % (variable.name, variable.units))

    print('Dependencies:')
    for variable in model.dependencies:
       print('  %s = %s (%s)' % (variable.name, variable.long_name, variable.units))

if __name__ == "__main__":
    # execute only if run as a script
    main()
