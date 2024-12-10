#!/usr/bin/env python

"""
This script verifies that a biogeochemical model returns valid derivatives
under a wide variety of inputs (state variables, environmental dependencies).
Different tests can be run, using either random values or extremes for inputs,
and running randomized or exhaustive tests.
"""

from typing import (
    MutableSequence,
    MutableSet,
    Tuple,
    Callable,
    Mapping,
    Sequence,
    Iterable,
    Union,
    Dict,
    cast,
    Optional,
)
import sys
import yaml

import numpy as np

try:
    import pyfabm
except ImportError:
    print("Unable to load pyfabm. See https://fabm.net/python.")
    sys.exit(1)


VAR_TYPE = Union[pyfabm.StateVariable, pyfabm.Dependency]
VARY_TYPE = Sequence[Tuple[VAR_TYPE, Sequence[float]]]
MUTABLE_VARY_TYPE = MutableSequence[Tuple[VAR_TYPE, MutableSequence[float]]]


def testRangePresence(
    ranges: Mapping[str, Union[float, int, MutableSequence[Union[int, float]]]],
    variables: Iterable[VAR_TYPE],
    vary: MUTABLE_VARY_TYPE,
    found_variables: MutableSet[str],
):
    for variable in variables:
        if variable.name not in ranges:
            print(f"No range specified for variable {variable.name}.")
            sys.exit(1)
        variable_range = ranges[variable.name]
        if isinstance(variable_range, (int, float)):
            variable.value = cast(np.ndarray, float(variable_range))
        else:
            if len(variable_range) != 2:
                print(
                    f"Invalid range {variable_range!r} specified for"
                    f" {variable.name}. It must be either a constant value,"
                    f" or a [minimum, maximum] array."
                )
                sys.exit(1)
            try:
                variable_range[0] = float(variable_range[0])
            except ValueError:
                print(
                    f"Invalid minimum value {variable_range[0]!r} specified for"
                    f" {variable.name}."
                )
                sys.exit(1)
            try:
                variable_range[1] = float(variable_range[1])
            except ValueError:
                print(
                    f"Invalid maximum value {variable_range[1]!r} specified for"
                    f" {variable.name}."
                )
                sys.exit(1)
            if variable.value is not None:
                value = float(variable.value)
                if value < variable_range[0]:
                    print(
                        f"Raising default value of {variable.name} to prescribed"
                        f" minimum {variable_range[0]} (previous default was {value})."
                    )
                    variable.value = cast(np.ndarray, variable_range[0])
                elif value > variable_range[1]:
                    print(
                        f"Lowering default value of {variable.name} to prescribed"
                        f" maximum {variable_range[1]} (previous default was {value})."
                    )
                    variable.value = cast(np.ndarray, variable_range[1])
            else:
                center = 0.5 * (variable_range[0] + variable_range[1])
                variable.value = cast(np.ndarray, center)
                print(
                    f"Setting default value of {variable.name} to center {center}"
                    f" of prescribed range {variable_range[0]} - {variable_range[1]}"
                    f" (no previous default was set)."
                )
            if variable_range[0] == variable_range[1]:
                variable.value = cast(np.ndarray, variable_range[0])
            else:
                vary.append((variable, variable_range))
        found_variables.add(variable.name)


ndone = 0


def check(model: pyfabm.Model):
    rates = model.getRates()
    assert len(rates) == len(model.state_variables)
    valid = np.isfinite(rates)
    global ndone
    ndone += 1
    if not valid.all():
        print(f"Test {ndone} FAILED!")
        variable: VAR_TYPE
        for variable, value in zip(model.state_variables, rates):
            if not np.isfinite(value):
                print(f"Change in {variable.name} has invalid value {value}")
        values: Dict[str, np.ndarray] = {}
        print("MODEL STATE:")
        for variable in model.state_variables:
            print(f"- {variable.name} = {variable.value}")
            values[variable.name] = np.array(variable.value)
        print("ENVIRONMENT:")
        for variable in model.dependencies:
            print(f"- {variable.name} = {variable.value}")
            values[variable.name] = np.array(variable.value)
        with open("last_error.yaml", "w") as f:
            yaml.safe_dump(values, f, default_flow_style=False)
        print("This model state and environment has been saved in last_error.yaml.")
        print(
            "To retest with these exact inputs (e.g., after introducing model fixes),"
            " specify this file in the ranges_path argument."
        )
        sys.exit(1)


def testRandomized(model: pyfabm.Model, vary: VARY_TYPE):
    # Perpetual random test:
    # for each model input, pick a value from its valid range [minimum,maximum]
    while 1:
        random_values = np.random.rand(len(vary))
        for (variable, (minimum, maximum)), random_value in zip(vary, random_values):
            value = minimum + (maximum - minimum) * random_value
            variable.value = cast(np.ndarray, value)
        check(model)
        if ndone % 1000 == 0:
            print(f"Test {ndone} completed.")


def testRandomizedExtremes(model: pyfabm.Model, vary: VARY_TYPE):
    # Perpetual random test:
    # for each model input, pick either its minimum or its maximum value.
    while 1:
        pick_maxs = np.random.rand(len(vary)) > 0.5
        for (variable, (minimum, maximum)), pick_max in zip(vary, pick_maxs):
            value = maximum if pick_max else minimum
            variable.value = cast(np.ndarray, value)
        check(model)
        if ndone % 1000 == 0:
            print(f"Test {ndone} completed.")


def testExtremes(model: pyfabm.Model, vary: VARY_TYPE):
    # Finite deterministic test:
    # for each model input, test minimum and maximum, leaving all other inputs
    # at their default value.
    for variable, (minimum, maximum) in vary:
        print(f"Testing {variable.name} = {minimum}...")
        oldvalue = np.array(variable.value)
        variable.value = cast(np.ndarray, minimum)
        check(model)
        print(f"Testing {variable.name} = {maximum}...")
        variable.value = cast(np.ndarray, maximum)
        check(model)
        variable.value = oldvalue


def testExtremesRecursive(
    model: pyfabm.Model,
    vary: VARY_TYPE,
    ntot: Optional[int] = None,
):
    # Finite deterministic test:
    # test all possible combinations of minimum and maximm for each model input.
    if ntot is None:
        ntot = 2 ** len(vary)
    if len(vary) == 0:
        check(model)
        if ndone % 1000 == 0:
            print(f"Completed {ndone} of {ntot} tests")
        return
    variable, (minimum, maximum) = vary[0]
    oldvalue = np.array(variable.value)
    variable.value = cast(np.ndarray, minimum)
    testExtremesRecursive(model, vary[1:], ntot)
    variable.value = cast(np.ndarray, maximum)
    testExtremesRecursive(model, vary[1:], ntot)
    variable.value = oldvalue


def main() -> None:
    import argparse

    tests: Mapping[str, Callable[[pyfabm.Model, VARY_TYPE], None]] = {
        "randomized": testRandomized,
        "extremes_per_variable": testExtremes,
        "extremes_randomized": testRandomizedExtremes,
        "extremes_all": testExtremesRecursive,
    }

    parser = argparse.ArgumentParser(
        description=(
            "This script verifies that a biogeochemical model returns valid"
            " derivatives under a wide variety of inputs (state variables,"
            " environmental dependencies). Different tests can be run, using"
            " either random values or extremes for inputs, and running"
            " randomized or exhaustive tests."
        )
    )
    parser.add_argument(
        "model_path",
        help="Path to a YAML file with the model configuration (typically fabm.yaml)",
    )
    parser.add_argument(
        "ranges_path",
        help=(
            "Path to a YAML file with ranges for all model inputs (state"
            " variables, environmental dependencies)"
        ),
    )
    parser.add_argument(
        "--test",
        choices=tests.keys(),
        action="append",
        help=(
            "Path to a YAML file with ranges for all model inputs (state variables,"
            " environmental dependencies)"
        ),
    )
    parser.add_argument(
        "--write-ranges",
        action="store_true",
        help="Write ranges file with default variable ranges",
    )
    args = parser.parse_args()

    # Create model object from YAML file.
    model = pyfabm.Model(args.model_path)

    if args.write_ranges:
        with open(args.ranges_path, "w") as f:

            def writeRanges(variables: Iterable[VAR_TYPE]):
                for variable in variables:
                    strmax = "?" if variable.value is None else 10 * variable.value
                    f.write(f"{variable.name}: [0,{strmax}]\n")

            writeRanges(model.state_variables)
            writeRanges(model.dependencies)
        print(f"Default ranges have been written to {args.ranges_path}.")
        sys.exit(0)

    with open(args.ranges_path, "r") as f:
        ranges = yaml.safe_load(f)
    if not isinstance(ranges, dict):
        print(
            f"Range file {args.ranges_path} should contain a dictionary mapping"
            f" each variable to its range (or constant value)."
        )
        sys.exit(1)

    vary: MUTABLE_VARY_TYPE = []
    found_variables: MutableSet[str] = set()
    testRangePresence(ranges, model.state_variables, vary, found_variables)
    testRangePresence(ranges, model.dependencies, vary, found_variables)

    for variable_name in ranges.keys():
        if variable_name not in found_variables:
            print(
                f"WARNING: range specification for unknown variable {variable_name}"
                f" in {args.ranges_path} will be ignored."
            )

    model.cell_thickness = 1
    model.start()

    if args.test is None:
        check(model)
    else:
        for test in args.test:
            tests[test](model, vary)


if __name__ == "__main__":
    main()
