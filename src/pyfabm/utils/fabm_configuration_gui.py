#!/usr/bin/env python

"""
This script open an interactive tree view with all settings of a biogeochemical model.
"""

import sys

try:
    import pyfabm
except ImportError:
    print("Unable to load pyfabm. See https://fabm.net/python.")
    sys.exit(1)
import pyfabm.gui_qt

QtCore = pyfabm.gui_qt.QtCore
QtWidgets = pyfabm.gui_qt.QtWidgets


def main() -> None:
    import argparse

    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "path",
        help="Path to a YAML file with the model configuration",
        nargs="?",
        default="fabm.yaml",
    )
    args = parser.parse_args()

    # Create model object
    model = pyfabm.Model(args.path)

    # Create Qt application object
    app = QtWidgets.QApplication([" "])

    # Create dialog box with model configuration tree.
    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle("Configure model")
    layout = QtWidgets.QHBoxLayout()
    tree = pyfabm.gui_qt.TreeView(model, dialog)
    layout.addWidget(tree)
    dialog.setLayout(layout)
    dialog.show()

    # Show dialog
    app.exec()


if __name__ == "__main__":
    # execute only if run as a script
    main()
