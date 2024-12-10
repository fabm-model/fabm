#!/usr/bin/env python
"""This script tests a FABM configuration (fabm.yaml) by reading the value of
the model state and environmental dependencies from one or more files (e.g.,
outputs of a 3D model), and then evaluating the model's rate of change for that
state and environment. The absolute and relative rates of change per state variable
are then shown, allowing one to assess model behaviour and time step constraints.

The values of state variables and environmental dependencies can be read from NetCDF
or YAML files (the YAML file should be a simple dictionary with variable_name: value
pairs). You can set additional variables on the command line with -v/--values.
"""

import sys
import os
from typing import Union, MutableMapping, Mapping, Iterable, cast

import numpy as np
import netCDF4
import yaml

try:
    import pyfabm
except ImportError:
    print("Unable to load pyfabm. See https://fabm.net/python.")
    sys.exit(1)


def evaluate(
    yaml_path: str,
    sources: Iterable[str] = (),
    location: Mapping[str, int] = {},
    assignments: Mapping[str, float] = {},
    verbose: bool = True,
    ignore_missing: bool = False,
    surface: bool = True,
    bottom: bool = True,
):
    # Create model object from YAML file.
    model = pyfabm.Model(yaml_path)

    allvariables: pyfabm.NamedObjectList[
        Union[pyfabm.StateVariable, pyfabm.Dependency]
    ] = pyfabm.NamedObjectList(model.state_variables, model.dependencies)
    name2variable: MutableMapping[
        str, Union[pyfabm.StateVariable, pyfabm.Dependency]
    ] = {}
    for variable in allvariables:
        name2variable[variable.name] = variable
        if hasattr(variable, "output_name"):
            name2variable[variable.output_name] = variable
    lcname2variable = dict(
        [(name.lower(), variable) for (name, variable) in name2variable.items()]
    )

    def set_state(**dim2index: int):
        missing = set(allvariables)
        variable2source: MutableMapping[
            Union[pyfabm.StateVariable, pyfabm.Dependency], str
        ] = {}

        def set_variable(
            variable: Union[pyfabm.StateVariable, pyfabm.Dependency],
            value: float,
            source: str,
        ):
            missing.discard(variable)
            if variable in variable2source:
                print(
                    f"WARNING: {variable.name} = {variable.value} set by"
                    f" {variable2source[variable]} is overwritten with {value}"
                    f" set by {source}"
                )
            variable2source[variable] = source
            variable.value = cast(np.ndarray, value)

        for path in sources:
            if path.endswith("yaml"):
                with open(path) as f:
                    data = yaml.safe_load(f)
                for name, value in data.items():
                    variable = name2variable.get(name)
                    if variable is None:
                        variable = lcname2variable.get(name.lower())
                    if variable is None:
                        print(
                            f"ERROR: variable {name!r} specified in {path}"
                            " not found in model"
                        )
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
                                    print(
                                        f"ERROR: Dimension {dim} of"
                                        f" {variable.output_name} has length > 1;"
                                        f" an index must be specified with {dim}=INDEX"
                                    )
                                    sys.exit(1)
                                index = dim2index[dim]
                            indices.append(index)
                        set_variable(variable, float(ncvar[tuple(indices)]), path)

        for name, value in assignments.items():
            if name not in name2variable:
                print(f"Explicitly specified variable {name!r} not found in model.")
                sys.exit(2)
            variable = name2variable[name]
            missing.discard(variable)
            variable2source[variable] = "command line"
            variable.value = cast(np.ndarray, float(value))

        if verbose:
            print()
            print("State:")
            for sv in sorted(model.state_variables, key=lambda x: x.name.lower()):
                print(f"  {sv.name}: {sv.value}" f" [{variable2source.get(sv)}]")
            print("Environment:")
            for d in sorted(model.dependencies, key=lambda x: x.name.lower()):
                print(f"  {d.name}: {d.value} [{variable2source.get(d)}]")

        if missing:
            print("The following variables are still missing:")
            for variable in sorted(missing, key=lambda x: x.name.lower()):
                print(f"- {variable.name}", end="")
                if variable.name != variable.output_name:
                    print(f" (NetCDF: {variable.output_name})")
                print()

        return missing

    missing = set_state(**location)
    if missing and not ignore_missing:
        sys.exit(1)

    print("State variables with largest value:")
    for sv in sorted(
        model.state_variables, key=lambda x: abs(float(x.value)), reverse=True
    )[:3]:
        print(f"  {sv.name}: {sv.value} {sv.units}")

    # Get model rates
    rates = model.getRates(surface=surface, bottom=bottom)
    assert len(rates) == len(
        model.state_variables
    ), "Length of array with rates does not match number of state variables"

    if verbose:
        print("Diagnostics:")
        for dv in sorted(model.diagnostic_variables, key=lambda x: x.name.lower()):
            if dv.output:
                print(f"  {dv.name}: {dv.value} {dv.units}")

    # Check whether rates of change are valid numbers
    valids = np.isfinite(rates)
    if not valids.all():
        print("The following state variables have an invalid rate of change:")
        for sv, rate, valid in zip(model.state_variables, rates, valids):
            if not valid:
                print(f"  {sv.name}: {rate}")

    eps = 1e-30
    relative_rates = np.array(
        [
            rate / (variable.value + eps)
            for variable, rate in zip(model.state_variables, rates)
        ]
    )

    if verbose:
        # Show all rates of change, odered by their value relative to the state
        # variable's value.
        print("Relative rates of change (low to high):")
        for variable, rate, relative_rate in sorted(
            zip(model.state_variables, rates, relative_rates), key=lambda x: x[2]
        ):
            print(f"  {variable.name}: {86400 * relative_rate} d-1")

    print("Largest relative rates of change:")
    for variable, rate, relative_rate in sorted(
        zip(model.state_variables, rates, relative_rates),
        key=lambda x: abs(x[2]),
        reverse=True,
    )[:3]:
        print(f"  {variable.name}: {86400 * relative_rate} d-1")

    i = int(relative_rates.argmin())
    print(
        f"Minimum time step = {-1.0 / relative_rates[i]:%.3f} s due to decrease"
        f" in {model.state_variables[i].name}"
    )


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "This script evaluates a biogeochemical model for a state and"
            " environment specified in one or more NetCDF files, yaml files,"
            " and command line arguments."
        )
    )
    parser.add_argument(
        "model_path",
        help="Path to a YAML file with the model configuration (typically fabm.yaml)",
    )
    parser.add_argument(
        "sources",
        nargs="+",
        help="Path to NetCDF or yaml file with the model state and environment",
    )
    parser.add_argument(
        "-l",
        "--location",
        nargs="+",
        help=(
            "NetCDF dimension to fix at particular index"
            " (specify: DIMENSION_NAME=INDEX)"
        ),
        default=[],
    )
    parser.add_argument(
        "-v",
        "--values",
        nargs="+",
        help=(
            "Additional state variable/environmental dependency values"
            " (specify: VARIABLE_NAME=VALUE)"
        ),
        default=[],
    )
    parser.add_argument(
        "--ignore_missing",
        action="store_true",
        help=(
            "Whether to ignore missing values for state variables and"
            " dependencies (the model will be evaluated with a default value"
            " of 0 for such missing variables)"
        ),
        default=False,
    )
    parser.add_argument(
        "--no_surface",
        dest="surface",
        action="store_false",
        help="Whether to omit surface processes (do_surface calls)",
        default=True,
    )
    parser.add_argument(
        "--no_bottom",
        dest="bottom",
        action="store_false",
        help="Whether to omit surface processes (do_bottom calls)",
        default=True,
    )
    parser.add_argument(
        "--pause",
        action="store_true",
        help="Whether to pause before model evaluation to manually attach a debugger.",
        default=False,
    )
    args = parser.parse_args()

    if args.pause:
        input(f"Attach the debugger (process id = {os.getpid()}) and then press Enter.")

    evaluate(
        args.model_path,
        args.sources,
        location=dict(
            [dimension2index.split("=") for dimension2index in args.location]
        ),
        assignments=dict([name2value.split("=") for name2value in args.values]),
        ignore_missing=args.ignore_missing,
        surface=args.surface,
        bottom=args.bottom,
    )


if __name__ == "__main__":
    main()
