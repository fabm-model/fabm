#!/usr/bin/env python

import sys, os

yamlfile = '../fabm-npzd-carbonate.yaml'
if len(sys.argv)>1: yamlfile = sys.argv[1]

# Try to locate FABM library.
# Retrieve FABM installation directory from enviornemtn variable FABM_PREFIX if set.
# Otherwise, try platform-specific default installation locations.
if 'FABM_PREFIX' in os.environ:
    prefix = os.environ['FABM_PREFIX']
elif os.name=='nt':
    prefix = '%s%s' % (os.environ['HOMEDRIVE'],os.environ['HOMEPATH'])
else:
    prefix = os.path.join(os.environ['HOME'],'local')
sys.path.append(os.path.join(prefix,'fabm/python'))
import pyfabm
import pyfabm.gui_qt

# Create model object
model = pyfabm.Model(yamlfile)

from PySide import QtCore,QtGui

app = QtGui.QApplication([' '])
dialog = QtGui.QDialog()
layout = QtGui.QHBoxLayout()

tree = QtGui.QTreeView()

itemmodel = pyfabm.gui_qt.ItemModel(model,dialog)

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
            contextMenu.addAction('Reset to %s' % data.default,reset)
            contextMenu.exec_(tree.mapToGlobal(pos))
tree.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
tree.customContextMenuRequested.connect(onTreeViewContextMenu)
for i in range(3): tree.resizeColumnToContents(i)

layout.addWidget(tree)
dialog.setLayout(layout)
dialog.show()

ret = app.exec_()
