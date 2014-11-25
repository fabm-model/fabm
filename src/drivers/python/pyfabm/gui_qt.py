import sys
import pyfabm

try:
    from PySide import QtCore,QtGui
except ImportError:
    print 'Unable to load PySide. Is it installed?'
    sys.exit(1)

class Delegate(QtGui.QStyledItemDelegate):
    def __init__(self,parent=None):
        QtGui.QStyledItemDelegate.__init__(self,parent)
    def createEditor(self,parent,option,index):
        assert index.isValid()
        data = index.internalPointer().object
        if not isinstance(data,basestring):
            options = data.getOptions()
            if options is not None:
                widget = QtGui.QComboBox(parent)
                widget.addItems(options)
                return widget
            elif isinstance(data.value,float):
                widget = ScientificDoubleEditor(parent)
                if data.units: widget.setSuffix(u' %s' % data.units_unicode)
                return widget
        return QtGui.QStyledItemDelegate.createEditor(self,parent,option,index)
    def setEditorData(self,editor,index):
        if isinstance(editor,QtGui.QComboBox):
            data = index.internalPointer().object
            if not isinstance(data,basestring):
                options = data.getOptions()
                if options is not None:
                    editor.setCurrentIndex(list(options).index(data.value))
                    return
        elif isinstance(editor,ScientificDoubleEditor):
            data = index.internalPointer().object
            editor.setValue(data.value)
            return
        return QtGui.QStyledItemDelegate.setEditorData(self,editor,index)
    def setModelData(self,editor,model,index):
        if isinstance(editor,QtGui.QComboBox):
            data = index.internalPointer().object
            if not isinstance(data,basestring):
                options = data.getOptions()
                if options is not None:
                    i = editor.currentIndex()
                    model.setData(index,options[i],QtCore.Qt.EditRole)
                    return
        elif isinstance(editor,ScientificDoubleEditor):
            model.setData(index,editor.value(),QtCore.Qt.EditRole)
            return
        return QtGui.QStyledItemDelegate.setModelData(self,editor,model,index)

class Entry(object):
    def __init__(self,object=None,name=''):
        if name=='' and object is not None: name = object
        self.object = object
        self.name = name
        self.parent = None
        self.children = []
        assert isinstance(self.name,basestring)

    def addChild(self,child):
        child.parent = self
        self.children.append(child)

    def insertChild(self,index,child):
        child.parent = self
        self.children.insert(index,child)

    def removeChild(self,index):
        self.children.pop(index).parent = None

    def findChild(self,name):
        for child in self.children:
            if isinstance(child.object,basestring) and child.object==name:
                return child
        child = Entry(name)
        self.addChild(child)
        return child

    def addTree(self,arr,category=None):
        for variable in arr:
            pathcomps = variable.path.split('/')
            if len(pathcomps)<2:
                #print 'Skipping root level entity %s' % pathcomps[0]
                continue
            if category is not None: pathcomps = [pathcomps[0],category] + list(pathcomps[1:])
            parent = self
            for component in pathcomps[:-1]:
                parent = parent.findChild(component)
            entry = Entry(variable,pathcomps[-1])
            parent.addChild(entry)

class Submodel:
    def __init__(self,long_name):
        self.units = None
        self.units_unicode = None
        self.value = None
        self.long_name = long_name

class ItemModel(QtCore.QAbstractItemModel):
    def __init__(self,model,parent):
        QtCore.QAbstractItemModel.__init__(self,parent)
        self.root = None
        self.model = model
        self.rebuild()

    def rebuild(self):
        root = Entry()
        env = Entry('environment')
        for d in self.model.dependencies: env.addChild(Entry(d,d.name))
        root.addTree(self.model.parameters,'parameters')
        root.addTree(self.model.state_variables,'initialization')
        root.addTree(self.model.couplings,'coupling')
        root.addChild(env)

        # For all models, create an object that returns appropriate metadata.
        def processNode(n,path=()):
            for i in range(len(n.children)-1,-1,-1):
                child = n.children[i]
                childpath = path + (child.name,)
                if isinstance(child.object,basestring):
                    if childpath[-1] in ('parameters','initialization','environment','coupling'):
                        childpath = childpath[:-1]
                    else:
                        submodel = self.model.getSubModel('/'.join(childpath))
                        if submodel.user_created:
                            child.object = Submodel(submodel.long_name)
                        else:
                            n.removeChild(i)
                            child = None
                if child is not None: processNode(child,childpath)
        processNode(root)

        if self.root is not None:
            # We already have an old tree - compare and amend model.
            def processChange(newnode,oldnode,parent):
                oldnode.object = newnode.object
                if not newnode.children: return
                ioldstart = 0
                for node in newnode.children:
                    iold = -1
                    for i in range(ioldstart,len(oldnode.children)):
                        if node.name==oldnode.children[i].name:
                            iold = i
                            break
                    if iold!=-1:
                        # New node was found among old nodes; remove any unused old nodes that precede it.
                        if iold>ioldstart:
                            self.beginRemoveRows(parent,ioldstart,iold-1)
                            for i in range(iold-1,ioldstart-1,-1): oldnode.removeChild(i)
                            self.endRemoveRows()
                        # Process changes to children of node.
                        processChange(node,oldnode.children[ioldstart],self.createIndex(ioldstart,0,oldnode.children[ioldstart]))
                    else:
                        # New node not found; insert it.
                        self.beginInsertRows(parent,ioldstart,ioldstart)
                        oldnode.insertChild(ioldstart,node)
                        self.endInsertRows()
                    ioldstart = ioldstart + 1
                if ioldstart<len(oldnode.children):
                    # Remove any trailing unused old nodes.
                    self.beginRemoveRows(parent,ioldstart,len(oldnode.children)-1)
                    for i in range(len(oldnode.children)-1,ioldstart-1,-1): oldnode.removeChild(i)
                    self.endRemoveRows()

            processChange(root,self.root,QtCore.QModelIndex())
        else:
            # First time a tree was created - store it and move on.
            self.root = root
        #self.reset()

    def rowCount(self,index):
        if not index.isValid():
            return len(self.root.children)
        elif index.column()==0:
            return len(index.internalPointer().children)
        return 0

    def columnCount(self,index):
        return 4

    def index(self,row,column,parent):
        if not parent.isValid():
            # top-level
            children = self.root.children
        else:
            # not top-level
            children = parent.internalPointer().children
        if row<0 or row>=len(children) or column<0 or column>=4: return QtCore.QModelIndex()
        return self.createIndex(row,column,children[row])

    def parent(self,index):
        if not index.isValid(): return QtCore.QModelIndex()
        parent = index.internalPointer().parent
        if parent is self.root: return QtCore.QModelIndex()   
        irow = parent.parent.children.index(parent)
        return self.createIndex(irow,0,parent)

    def data(self,index,role):
        if not index.isValid(): return
        entry = index.internalPointer()
        data = entry.object
        if role==QtCore.Qt.DisplayRole:
            if index.column()==0:
                return entry.name if isinstance(data,basestring) else data.long_name
            if not isinstance(data,(basestring,Submodel)):
                if index.column()==1:
                    value = data.value
                    if not isinstance(value,bool):
                        if data.units:
                            return u'%s %s' % (value,data.units_unicode)
                        else:
                            return unicode(value)
                elif index.column()==2 and data.units:
                    return data.units_unicode
                elif index.column()==3:
                    return entry.name
        elif role==QtCore.Qt.ToolTipRole and index.parent().isValid():
           if not isinstance(data,basestring): return data.long_path
        elif role==QtCore.Qt.EditRole:
           if not isinstance(data,(basestring,Submodel)):
              #print data.getOptions()
              return data.getValue()
        elif role==QtCore.Qt.FontRole and index.column()==1:
            if isinstance(data,pyfabm.Parameter) and data.value!=data.default:
                font = QtGui.QFont()
                font.setBold(True)
                return font
        elif role==QtCore.Qt.CheckStateRole and index.column()==1 and not isinstance(data,basestring) and isinstance(data.value,bool):
            return QtCore.Qt.Checked if data.value else QtCore.Qt.Unchecked
        return None

    def setData(self,index,value,role):
        if role==QtCore.Qt.CheckStateRole: value = value==QtCore.Qt.Checked
        if role in (QtCore.Qt.EditRole,QtCore.Qt.CheckStateRole):
            data = index.internalPointer().object
            data.setValue(value)
            if isinstance(data,pyfabm.Parameter): self.rebuild()
            return True
        return False

    def flags(self,index):
        flags = QtCore.Qt.NoItemFlags
        if not index.isValid(): return flags
        if index.column()==1:
            entry = index.internalPointer().object
            if not isinstance(entry,(basestring,Submodel)):
                if isinstance(entry.value,bool):
                    flags |= QtCore.Qt.ItemIsUserCheckable
                else:
                    flags |= QtCore.Qt.ItemIsEditable
        return flags | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self,section,orientation,role):
        if orientation==QtCore.Qt.Horizontal and role==QtCore.Qt.DisplayRole and section>=0 and section<4:
            return ('name','value','units','symbol')[section]

class TreeView(QtGui.QTreeView):
    def __init__(self,model,parent):
        QtGui.QTreeView.__init__(self,parent)
        itemmodel = pyfabm.gui_qt.ItemModel(model,parent)
        self.setItemDelegate(Delegate(parent))
        self.setModel(itemmodel)
        self.setUniformRowHeights(True)
        self.expandAll()
        def onTreeViewContextMenu(pos):
            index = self.indexAt(pos)
            if index.isValid() and index.column()==1:
                data = index.internalPointer().object
                if isinstance(data,pyfabm.Parameter) and data.value!=data.default and data.default is not None:
                    def reset():
                        data.reset()
                        itemmodel.rebuild()
                    contextMenu = QtGui.QMenu(self)
                    default = data.default
                    if data.units: default = u'%s %s' % (data.default,data.units_unicode)
                    contextMenu.addAction(u'Reset to default: %s' % default,reset)
                    contextMenu.exec_(self.mapToGlobal(pos))
        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(onTreeViewContextMenu)
        for i in range(3): self.resizeColumnToContents(i)

# From xmlstore (GOTM-GUI)

class ScientificDoubleValidator(QtGui.QValidator):
    """Qt validator for floating point values
    Less strict than the standard QDoubleValidator, in the sense that is
    also accepts values in scientific format (e.g. 1.2e6)
    Also has properties 'minimum' and 'maximum', used for validation and
    fix-up.
    """
    def __init__(self,parent=None):
        QtGui.QValidator.__init__(self,parent)
        self.minimum = None
        self.maximum = None
        self.suffix = ''

    def validate(self,input,pos):
        assert isinstance(input,basestring),'input argument is not a string (old PyQt4 API?)'

        # Check for suffix (if ok, cut it off for further value checking)
        if not input.endswith(self.suffix): return (QtGui.QValidator.Invalid,input,pos)
        vallength = len(input)-len(self.suffix)

        # Check for invalid characters
        rx = QtCore.QRegExp('[^\d\-+eE,.]')
        if rx.indexIn(input[:vallength])!=-1: return (QtGui.QValidator.Invalid,input,pos)
        
        # Check if we can convert it into a floating point value
        try:
            v = float(input[:vallength])
        except ValueError:
            return (QtGui.QValidator.Intermediate,input,pos)

        # Check for minimum and maximum.
        if self.minimum is not None and v<self.minimum: return (QtGui.QValidator.Intermediate,input,pos)
        if self.maximum is not None and v>self.maximum: return (QtGui.QValidator.Intermediate,input,pos)
        
        return (QtGui.QValidator.Acceptable,input,pos)

    def fixup(self,input):
        assert isinstance(input,basestring),'input argument is not a string (old PyQt4 API?)'
        if not input.endswith(self.suffix): return input

        try:
            v = float(input[:len(input)-len(self.suffix)])
        except ValueError:
            return input

        if self.minimum is not None and v<self.minimum: input = u'%s%s' % (self.minimum,self.suffix)
        if self.maximum is not None and v>self.maximum: input = u'%s%s' % (self.maximum,self.suffix)
        print u'"%s"' % input
        return input

    def setSuffix(self,suffix):
        self.suffix = suffix

class ScientificDoubleEditor(QtGui.QLineEdit):
    """Editor for a floating point value.
    """
    def __init__(self,parent):
        QtGui.QLineEdit.__init__(self,parent)

        self.curvalidator = ScientificDoubleValidator(self)
        self.setValidator(self.curvalidator)
        self.suffix = ''
        #self.editingFinished.connect(self.onPropertyEditingFinished)

    def setSuffix(self,suffix):
        value = self.value()
        self.suffix = suffix
        self.curvalidator.setSuffix(suffix)
        self.setValue(value)

    def value(self):
        text = self.text()
        text = text[:len(text)-len(self.suffix)]
        if text=='': return 0
        return float(text)

    def setValue(self,value,format=None):
        if value is None:
            strvalue = ''
        else:  
            if format is None:
                strvalue = unicode(value)
            else:
                strvalue = format % value
        self.setText(u'%s%s' % (strvalue,self.suffix))

    def focusInEvent(self,e):
        QtGui.QLineEdit.focusInEvent(self,e)
        self.selectAll()

    def selectAll(self):
        QtGui.QLineEdit.setSelection(self,0,len(self.text())-len(self.suffix))

    def setMinimum(self,minimum):
        self.curvalidator.minimum = minimum

    def setMaximum(self,maximum):
        self.curvalidator.maximum = maximum

    def interpretText(self):
        if not self.hasAcceptableInput():
            text = self.text()
            textnew = self.curvalidator.fixup(text)
            if textnew is None: textnew = text
            self.setText(textnew)
