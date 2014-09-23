#!/usr/bin/env python

import sys,optparse

parser = optparse.OptionParser()
options, args = parser.parse_args()

yamlfile = '../fabm-npzd-carbonate.yaml'
if args: yamlfile = args[0]

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

tree = QtGui.QTreeView()

itemmodel = pyfabm.gui_qt.ItemModel(model,dialog)
tree.setItemDelegate(pyfabm.gui_qt.Delegate(dialog))
tree.setModel(itemmodel)
tree.setUniformRowHeights(True)
tree.expandAll()
def onTreeViewContextMenu(pos):
    index = tree.indexAt(pos)
    if index.isValid() and index.column()==1:
        data = index.internalPointer().object
        if isinstance(data,pyfabm.Parameter) and data.value!=data.default:
            def reset():
                data.reset()
                itemmodel.rebuild()
            contextMenu = QtGui.QMenu(tree)
            contextMenu.addAction('Reset to default (%s)' % data.default,reset)
            contextMenu.exec_(tree.mapToGlobal(pos))
tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
tree.customContextMenuRequested.connect(onTreeViewContextMenu)
for i in range(3): tree.resizeColumnToContents(i)

layout.addWidget(tree)
dialog.setLayout(layout)
dialog.show()

ret = app.exec_()
