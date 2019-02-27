#!/usr/bin/env python

"""
This script verifies that a biogeochemical model returns valid derivatives under a wide variety of inputs (state variables, environmental dependencies). Different tests can be run, using either random values or extremes for inputs, and running randomized or exhaustive tests.
"""

from __future__ import print_function
import sys
import yaml

import numpy

try:
   import pyfabm
except ImportError:
   print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
   sys.exit(1)

def testRangePresence(ranges, variables, vary, found_variables):
    for variable in variables:
        if variable.name not in ranges:
            print('No range specified for variable %s.' % (variable.name,))
            sys.exit(1)
        variable_range = ranges[variable.name]
        if isinstance(variable_range, (int, float)):
            variable.value = float(variable_range)
        else:
            if len(variable_range) != 2:
                print('Invalid range %s specified for %s. It must be either a constant value, or a [minimum,maximum] array.' % (variable_range, variable.name))
                sys.exit(1)
            try:
                variable_range[0] = float(variable_range[0])
            except ValueError:
                print('Invalid minimum value "%s" specified for %s.' % (variable_range[0], variable.name))
                sys.exit(1)
            try:
                variable_range[1] = float(variable_range[1])
            except ValueError:
                print('Invalid maximum value "%s" specified for %s.' % (variable_range[1], variable.name))
                sys.exit(1)
            value = float(variable.value)
            if value < variable_range[0]:
                print('Raising default value of %s to prescribed minimum %s (previous default was %s).' % (variable.name, variable_range[0], value))
                variable.value = variable_range[0]
            elif value > variable_range[1]:
                print('Lowering default value of %s to prescribed maximum %s (previous default was %s).' % (variable.name, variable_range[1], value))
                variable.value = variable_range[1]
            if variable_range[0] == variable_range[1]:
               variable.value = float(variable_range[0])
            else:
               vary.append((variable, variable_range))
        found_variables.add(variable.name)

ndone = 0

def check(model):
    rates = model.getRates()
    assert len(rates) == len(model.state_variables)
    valid = numpy.isfinite(rates)
    global ndone
    ndone += 1
    if not valid.all():
       print('Test %i FAILED!' % ndone)
       for variable, value in zip(model.state_variables, rates):
          if not numpy.isfinite(value):
             print('Change in %s has invalid value %s' % (variable.name, value))
       values = {}
       print('MODEL STATE:')
       for variable in model.state_variables:
           print('- %s = %s' % (variable.name, variable.value))
           values[variable.name] = variable.value
       print('ENVIRONMENT:')
       for variable in model.dependencies:
           print('- %s = %s' % (variable.name, variable.value))
           values[variable.name] = variable.value
       with open('last_error.yaml', 'w') as f:
           yaml.dump(values, f, default_flow_style=False)
       print('This model state and environment has been saved in last_error.yaml.')
       print('To retest with these exact inputs (e.g., after introducing model fixes),')
       print('specify this file as the "ranges" argument.')
       sys.exit(1)

def testRandomized(model, vary):
    # Perpetual random test:
    # for each model input, pick a value from its valid range [minimum,maximum]
    while 1:
       random_values = numpy.random.rand(len(vary))
       for (variable, (minimum, maximum)), random_value in zip(vary, random_values):
          variable.value = minimum + (maximum - minimum) * random_value
       check(model)
       if ndone % 1000 == 0:
           print('Test %i completed.' % ndone)

def testRandomizedExtremes(model, vary):
    # Perpetual random test:
    # for each model input, pick either its minimum or its maximum value.
    while 1:
       pick_maxs = numpy.random.rand(len(vary)) > 0.5
       for (variable, (minimum, maximum)), pick_max in zip(vary, pick_maxs):
          variable.value = maximum if pick_max else minimum
       check(model)
       if ndone % 1000 == 0:
           print('Test %i completed.' % ndone)

def testExtremes(model, vary):
    # Finite deterministic test:
    # for each model input, test minimum and maximum, leaving all other inputs at their default value.
    for variable, (minimum, maximum) in vary:
       print('Testing %s = %s...' % (variable.name, minimum))
       oldvalue = float(variable.value)
       variable.value = minimum
       check(model)
       print('Testing %s = %s...' % (variable.name, maximum))
       variable.value = maximum
       check(model)
       variable.value = oldvalue

def testExtremesRecursive(model, vary, ntot=None):
    # Finite deterministic test:
    # test all possible combinations of minimum and maximm for each model input.
    if ntot is None:
        ntot = 2**len(vary)
    if len(vary) == 0:
       check(model)
       if ndone % 1000 == 0:
           print('Completed %i of %i tests' % (ndone, ntot))
       return
    variable, (minimum, maximum) = vary[0]
    oldvalue = float(variable.value)
    variable.value = minimum
    testExtremesRecursive(model, vary[1:], ntot)
    variable.value = maximum
    testExtremesRecursive(model, vary[1:], ntot)
    variable.value = oldvalue

def main():
    import argparse

    tests = {'randomized': testRandomized, 'extremes_per_variable': testExtremes, 'extremes_randomized': testRandomizedExtremes, 'extremes_all': testExtremesRecursive}

    parser = argparse.ArgumentParser(description='This script verifies that a biogeochemical model returns valid derivatives under a wide variety of inputs (state variables, environmental dependencies). Different tests can be run, using either random values or extremes for inputs, and running randomized or exhaustive tests.')
    parser.add_argument('model_path', help='Path to a YAML file with the model configuration (typically fabm.yaml)')
    parser.add_argument('ranges_path', help='Path to a YAML file with ranges for all model inputs (state variables, environmental dependencies)')
    parser.add_argument('--test', choices=tests.keys(), action='append', help='Path to a YAML file with ranges for all model inputs (state variables, environmental dependencies)')
    parser.add_argument('--write-ranges', action='store_true', help='Write ranges file with default variable ranges')
    args = parser.parse_args()

    # Create model object from YAML file.
    model = pyfabm.Model(args.model_path)

    if args.write_ranges:
        with open(args.ranges_path, 'w') as f:
            def writeRanges(variables):
                for variable in variables:
                    f.write('%s: [0,%s]\n' % (variable.name, 10*variable.value))
            writeRanges(model.state_variables)
            writeRanges(model.dependencies)
        print('Default ranges have been written to %s.' % args.ranges_path)
        sys.exit(0)

    with open(args.ranges_path, 'rU') as f:
        ranges = yaml.load(f)
    if not isinstance(ranges, dict):
        print('Range file %s should contain a dictionary mapping each variable to its range (or constant value).' % args.ranges_path)
        sys.exit(1)

    vary = []
    found_variables = set()
    testRangePresence(ranges, model.state_variables, vary, found_variables)
    testRangePresence(ranges, model.dependencies, vary, found_variables)

    for variable_name in ranges.keys():
        if variable_name not in found_variables:
            print('WARNING: range specification for unknown variable %s in %s will be ignored.' % (variable_name, args.ranges_path))

    if args.test is None:
        check(model)
    else:
        for test in args.test:
            tests[test](model, vary)

if __name__ == '__main__':
    main()

