import pyfabm


def processFile(infile, outfile, subtract_background=False, add_missing=False):
    # Create model object from YAML file.
    model = pyfabm.Model(infile)
    model.save_settings(
        outfile, pyfabm.DISPLAY_NORMAL if add_missing else pyfabm.DISPLAY_MINIMUM
    )
