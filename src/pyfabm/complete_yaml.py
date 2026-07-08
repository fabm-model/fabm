from typing import Union
import os

import pyfabm


def processFile(
    infile: Union[str, os.PathLike],
    outfile: Union[str, os.PathLike],
    subtract_background: bool = False,
    add_missing: bool = False,
):
    # Create model object from YAML file.
    model = pyfabm.Model(infile)
    model.save_settings(
        outfile, pyfabm.DISPLAY_NORMAL if add_missing else pyfabm.DISPLAY_MINIMUM
    )
