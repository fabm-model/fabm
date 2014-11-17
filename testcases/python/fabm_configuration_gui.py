#!/usr/bin/env python

import sys,optparse

try:
    import pyfabm
except ImportError:
    print 'Unable to load pyfabm. Please build and install FABM with FABM_HOST=python.'
    sys.exit(1)
import pyfabm.gui_qt

from PySide import QtCore,QtGui

# Parse command line arguments.
parser = optparse.OptionParser()
options, args = parser.parse_args()
if len(sys.argv)!=2:
   print 'This script takes one argument: the path to a YAML file with FABM settings (typically fabm.yaml).'
   sys.exit(2)
yamlfile = sys.argv[1]

# Create model object
model = pyfabm.Model(yamlfile)

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
