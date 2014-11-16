#!/usr/bin/env python

import sys,optparse

parser = optparse.OptionParser()
options, args = parser.parse_args()

if len(sys.argv)!=2:
   print 'This script takes one argument: the path to a YAML file with FABM settings (typically fabm.yaml).'
   sys.exit(2)
yamlfile = sys.argv[1]

try:
    import pyfabm
except ImportError:
    print 'Unable to load pyfabm. Please build and install FABM with FABM_HOST=python and make sure your Python installation includes numpy.'
    sys.exit(1)

try:
    import pyfabm.gui_qt
except ImportError:
    print 'Unable to load pyfabm.gui_qt. Do you have PySide installed?'
    sys.exit(1)


# Create model object
model = pyfabm.Model(yamlfile)

from PySide import QtCore,QtGui

app = QtGui.QApplication([' '])
dialog = QtGui.QDialog()
layout = QtGui.QHBoxLayout()

tree = pyfabm.gui_qt.TreeView(model,dialog)

layout.addWidget(tree)
dialog.setLayout(layout)
dialog.show()

ret = app.exec_()
