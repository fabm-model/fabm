#!/usr/bin/env python

"""
This script open an interactive tree view with all settings of a biogeochemical model.
"""

import sys

try:
    import pyfabm
except ImportError:
    print 'Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.'
    sys.exit(1)
import pyfabm.gui_qt

from PySide import QtCore,QtGui

def main():
    import argparse
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path',help='Path to a YAML file with the model configuration (typically fabm.yaml)',nargs='?',default='fabm.yaml')
    args = parser.parse_args()

    # Create model object
    model = pyfabm.Model(args.path)

    # Create Qt application object
    app = QtGui.QApplication([' '])

    # Create dialog box with model configuration tree.
    dialog = QtGui.QDialog()
    dialog.setWindowTitle('Configure model')
    layout = QtGui.QHBoxLayout()
    tree = pyfabm.gui_qt.TreeView(model,dialog)
    layout.addWidget(tree)
    dialog.setLayout(layout)
    dialog.show()

    # Show dialog
    ret = app.exec_()

if __name__ == "__main__":
    # execute only if run as a script
    main()
