import pyfabm


def processFile(
    infile: str,
    outfile: str,
    subtract_background: bool = False,
    add_missing: bool = False,
):
    # Create model object from YAML file.
    model = pyfabm.Model(infile)
    model.save_settings(
        outfile, pyfabm.DISPLAY_NORMAL if add_missing else pyfabm.DISPLAY_MINIMUM
    )
