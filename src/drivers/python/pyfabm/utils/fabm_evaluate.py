#!/usr/bin/env python
"""This script tests a FABM configuration (fabm.yaml) by reading the value of
the model state and environmental dependencies from one or more files (e.g.,
outputs of a 3D model), and then evaluating the model's rate of change for that
state and environment. The absolute and relative rates of change per state variable
are then shown, allowing one to assess model behaviour and time step constraints.

The values of state variables and environmental dependencies can be read from NetCDF
or YAML files (the YAML file should be a simple dictionary with variable_name: value pairs).
You can set additional variables on the command line with -v/--values.
"""

import sys
import os

try:
    input = raw_input
except NameError:
    pass

import numpy
import netCDF4
import yaml

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
    sys.exit(1)

def evaluate(yaml_path, sources=(), location={}, assignments={}, verbose=True, ignore_missing=False, surface=True, bottom=True):
    # Create model object from YAML file.
    model = pyfabm.Model(yaml_path)

    allvariables = list(model.state_variables) + list(model.dependencies)
    name2variable = {}
    for variable in allvariables:
        name2variable[variable.name] = variable
        if hasattr(variable, 'output_name'):
            name2variable[variable.output_name] = variable
    lcname2variable = dict([(name.lower(), variable) for (name, variable) in name2variable.items()])

    def set_state(**dim2index):
        missing = set(allvariables)
        variable2source = {}

        def set_variable(variable, value, source):
            missing.discard(variable)
            if variable in variable2source:
                print('WARNING: %s = %s set by %s is overwritten with %s set by %s' % (variable.name, variable.value, variable2source[variable], value, source))
            variable2source[variable] = source
            variable.value = value

        for path in sources:
            if path.endswith('yaml'):
                with open(path) as f:
                    data = yaml.load(f)
                for name, value in data.items():
                    variable = name2variable.get(name)
                    if variable is None:
                        variable = lcname2variable.get(name.lower())
                    if variable is None:
                        print('ERROR: variable "%s" specified in %s not found in model' % (name, path))
                        sys.exit(1)
                    set_variable(variable, float(value), path)
            else:
                with netCDF4.Dataset(path) as nc:
                    for variable in allvariables:
                        if variable.output_name not in nc.variables:
                            continue
                        ncvar = nc.variables[variable.output_name]
                        indices = []
                        for dim, length in zip(ncvar.dimensions, ncvar.shape):
                            index = 0
                            if length > 1:
                                if dim not in dim2index:
                                    print('ERROR: Dimension %s of %s has length > 1; an index must be specified with %s=INDEX' % (dim, variable.output_name, dim))
                                    sys.exit(1)
                                index = dim2index[dim]
                            indices.append(index)
                        set_variable(variable, float(ncvar[tuple(indices)]), path)

        for name, value in assignments.items():
            if name not in name2variable:
                print('Explicitly specified variable "%s" not found in model.' % name)
                sys.exit(2)
            variable = name2variable[name]
            missing.discard(variable)
            variable2source[variable] = 'command line'
            variable.value = float(value)

        if verbose:
            print()
            print('State:')
            for variable in sorted(model.state_variables, key=lambda x: x.name.lower()):
                print('  %s: %s [%s]' % (variable.name, variable.value, variable2source.get(variable)))
            print('Environment:')
            for variable in sorted(model.dependencies, key=lambda x: x.name.lower()):
                print('  %s: %s [%s]' % (variable.name, variable.value, variable2source.get(variable)))

        if missing:
            print('The following variables are still missing:')
            for variable in sorted(missing, key=lambda x: x.name.lower()):
                print('- %s' % variable.name,)
                if variable.name != variable.output_name:
                    print('(NetCDF: %s)' % variable.output_name,)
                print()

        return missing

    missing = set_state(**location)
    if missing and not ignore_missing:
        sys.exit(1)

    print('State variables with largest value:')
    for variable in sorted(model.state_variables, key=lambda x: abs(x.value), reverse=True)[:3]:
        print('  %s: %s %s' % (variable.name, variable.value, variable.units))

    # Get model rates
    rates = model.getRates(surface=surface, bottom=bottom)
    assert len(rates) == len(model.state_variables), 'Length of array with rates does not match number of state variables'

    if verbose:
        print('Diagnostics:')
        for variable in sorted(model.diagnostic_variables, key=lambda x: x.name.lower()):
            if variable.output:
                print('  %s: %s %s' % (variable.name, variable.value, variable.units))

    # Check whether rates of change are valid numbers
    valids = numpy.isfinite(rates)
    if not valids.all():
        print('The following state variables have an invalid rate of change:')
        for variable, rate, valid in zip(model.state_variables, rates, valids):
            if not valid:
                print('  %s: %s' % (variable.name, rate))

    eps = 1e-30
    relative_rates = numpy.array([rate/(variable.value+eps) for variable, rate in zip(model.state_variables, rates)])

    if verbose:
        # Show all rates of change, odered by their value relative to the state variable's value.
        print('Relative rates of change (low to high):')
        for variable, rate, relative_rate in sorted(zip(model.state_variables, rates, relative_rates), key=lambda x: x[2]):
            print('  %s: %s d-1' % (variable.name, 86400*relative_rate))

    print('Largest relative rates of change:')
    for variable, rate, relative_rate in sorted(zip(model.state_variables, rates, relative_rates), key=lambda x: abs(x[2]), reverse=True)[:3]:
        print('  %s: %s d-1' % (variable.name, 86400*relative_rate))

    i = relative_rates.argmin()
    print('Minimum time step = %.3f s due to decrease in %s' % (-1./relative_rates[i], model.state_variables[i].name))

def main():
    import argparse
    parser = argparse.ArgumentParser(description='This script evaluates a biogeochemical model for a state and environment specified in one or more NetCDF files, yaml files, and command line arguments.')
    parser.add_argument('model_path', help='Path to a YAML file with the model configuration (typically fabm.yaml)')
    parser.add_argument('sources', nargs='+', help='Path to NetCDF or yaml file with the model state and environment')
    parser.add_argument('-l', '--location', nargs='+', help='NetCDF dimension to fix at particular index (specify: DIMENSION_NAME=INDEX)', default=[])
    parser.add_argument('-v', '--values', nargs='+', help='Additional state variable/environmental dependency values (specify: VARIABLE_NAME=VALUE)', default=[])
    parser.add_argument('--ignore_missing', action='store_true', help='Whether to ignore missing values for state variables and dependencies (the model will be evaluated with a default value of 0 for such missing variables)', default=False)
    parser.add_argument('--no_surface', dest='surface', action='store_false', help='Whether to omit surface processes (do_surface calls)', default=True)
    parser.add_argument('--no_bottom', dest='bottom', action='store_false', help='Whether to omit surface processes (do_bottom calls)', default=True)
    parser.add_argument('--pause', action='store_true', help='Whether to pause before model evaluation to manually attach a debugger.', default=False)
    args = parser.parse_args()

    if args.pause:
        input('Attach the debugger (process id = %i) and then press Enter.' % os.getpid())

    evaluate(args.model_path, args.sources, location=dict([dimension2index.split('=') for dimension2index in args.location]), assignments=dict([name2value.split('=') for name2value in args.values]), ignore_missing=args.ignore_missing, surface=args.surface, bottom=args.bottom)

if __name__ == '__main__':
    main()
