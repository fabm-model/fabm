# Fortran-YAML

This is a lightweight [YAML](http://yaml.org) parser written in object-oriented Fortran.

The latest source code can be found at [GitHub](https://github.com/BoldingBruggeman/fortran-yaml).

## Features

This parser handles a subset of YAML only, currently subject to the following limitations:

* Only block style is supported for mappings and sequences; flow style is not.
* All scalars are left as strings; interpreting them as native data types is left to the caller. Thus, keys in the (key : value) pairs of a mapping are strings by definition.
* Use of quotes around strings is not supported.
* Multi-line strings are not supported.

Comments (starting with #) are allowed. As per the YAML specification, indentation must consist of spaces only (no tabs!)

For instance:

```yaml
instances:
  P1:
    model: pml/ersem/vphyt
    parameters:
      mu_max: 2.2
      K: 1
    coupling:
      R: pom/R6
```

While this is not an attempt to write a complete YAML parser in Fortran, it _is_ meant
to only accept documents that are valid YAML. If you find that this parser
permits constructs that the YAML specification disallows (but not vice versa),
please contact the author.

## Copyright and license

This software is copyright 2013-2016 Bolding & Bruggeman ApS.

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by [the Free Software
Foundation](https://www.gnu.org/licenses/gpl.html).

It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. A copy of the license is provided in the COPYING file.