import pyfabm
from PySide import QtCore,QtGui

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
            if category is not None: pathcomps = [pathcomps[0],category] + list(pathcomps[1:])
            parent = self
            for component in pathcomps[:-1]:
                parent = parent.findChild(component)
            entry = Entry(variable,pathcomps[-1])
            parent.addChild(entry)

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
        root.addChild(env)
        class Submodel:
            def __init__(self,long_name):
                self.units = None
                self.units_unicode = None
                self.value = None
                self.long_name = long_name
        def processNode(n,prefix=()):
            if isinstance(n.object,basestring) and prefix and prefix[-1] not in ('parameters','initialization','environment'):
                n.object = Submodel(self.model.getModelLongName('/'.join(prefix)))
            for child in n.children: processNode(child,prefix+(child.name,))
        processNode(root)
        if self.root is not None:
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
            # We already have an old tree - compare and amend model.
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
            if not isinstance(data,basestring):
                if index.column()==1:
                    value = data.value
                    if not isinstance(value,bool): return value
                elif index.column()==2 and data.units:
                    return data.units_unicode
                elif index.column()==3:
                    return entry.name
        elif role==QtCore.Qt.ToolTipRole and index.parent().isValid():
           if not isinstance(data,basestring): return '%s (%s)' % (data.long_name,entry.name)
        elif role==QtCore.Qt.EditRole:
           if not isinstance(data,basestring): return data.getValue()
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
            if not isinstance(entry,basestring):
                if isinstance(entry.value,bool):
                    flags |= QtCore.Qt.ItemIsUserCheckable
                else:
                    flags |= QtCore.Qt.ItemIsEditable
        return flags | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def headerData(self,section,orientation,role):
        if orientation==QtCore.Qt.Horizontal and role==QtCore.Qt.DisplayRole and section>=0 and section<4:
            return ('name','value','units','symbol')[section]
