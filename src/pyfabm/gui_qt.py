import sys
from typing import Iterable, Union, List, Optional, Tuple, cast
import re

import numpy as np
import pyfabm

try:
    from PySide6 import QtCore, QtGui, QtWidgets
except ImportError as e1:
    try:
        from PyQt5 import QtCore, QtGui, QtWidgets  # type: ignore[no-redef]
    except ImportError as e2:
        print(e1)
        print(e2)
        print("Unable to load PyQt5 or PySide6. Is either installed?")
        sys.exit(1)


class Delegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, parent: Optional[QtCore.QObject] = None):
        QtWidgets.QStyledItemDelegate.__init__(self, parent)

    def createEditor(
        self,
        parent: QtWidgets.QWidget,
        option: QtWidgets.QStyleOptionViewItem,
        index: QtCore.QModelIndex,  # type: ignore[override]
    ):
        assert index.isValid()
        entry = cast(Entry, index.internalPointer())
        data = entry.object
        if isinstance(data, pyfabm.Variable):
            if data.options is not None:
                combobox = QtWidgets.QComboBox(parent)
                combobox.addItems(data.options)
                return combobox
            elif data.value is None or isinstance(data.value, (float, np.ndarray)):
                editor = ScientificDoubleEditor(parent)
                if data.units:
                    editor.setSuffix(f" {data.units_unicode}")
                return editor
        return QtWidgets.QStyledItemDelegate.createEditor(self, parent, option, index)

    def setEditorData(
        self,
        editor: QtWidgets.QWidget,
        index: QtCore.QModelIndex,  # type: ignore[override]
    ):
        entry = cast(Entry, index.internalPointer())
        data = entry.object
        if isinstance(editor, QtWidgets.QComboBox):
            assert isinstance(data, pyfabm.Variable) and data.options is not None
            if data.value is not None:
                editor.setCurrentIndex(data.options.index(data.value))
                return
        elif isinstance(editor, ScientificDoubleEditor):
            assert isinstance(data, pyfabm.Variable)
            editor.setValue(data.value if data.value is None else float(data.value))
            return
        return QtWidgets.QStyledItemDelegate.setEditorData(self, editor, index)

    def setModelData(
        self,
        editor: QtWidgets.QWidget,
        model: QtCore.QAbstractItemModel,
        index: QtCore.QModelIndex,  # type: ignore[override]
    ):
        if isinstance(editor, QtWidgets.QComboBox):
            entry = cast(Entry, index.internalPointer())
            data = entry.object
            assert isinstance(data, pyfabm.Variable) and data.options is not None
            i = editor.currentIndex()
            model.setData(index, data.options[i], QtCore.Qt.ItemDataRole.EditRole)
            return
        elif isinstance(editor, ScientificDoubleEditor):
            model.setData(index, editor.value(), QtCore.Qt.ItemDataRole.EditRole)
            return
        return QtWidgets.QStyledItemDelegate.setModelData(self, editor, model, index)


class Entry:
    def __init__(
        self,
        name: str = "",
        object: Union[
            None,
            "Entry",
            pyfabm.Parameter,
            pyfabm.StateVariable,
            pyfabm.Dependency,
            pyfabm.Coupling,
            "Submodel",
        ] = None,
    ):
        self.name = name
        self.object = object
        self.parent: Optional["Entry"] = None
        self.children: List["Entry"] = []
        assert isinstance(self.name, str)

    def addChild(self, child: "Entry"):
        child.parent = self
        self.children.append(child)

    def insertChild(self, index: int, child: "Entry"):
        child.parent = self
        self.children.insert(index, child)

    def removeChild(self, index: int):
        self.children.pop(index).parent = None

    def findChild(self, name: str):
        for child in self.children:
            if child.object is None and child.name == name:
                return child
        child = Entry(name)
        self.addChild(child)
        return child

    def addTree(
        self,
        arr: Iterable[
            Union[
                pyfabm.Parameter,
                pyfabm.StateVariable,
                pyfabm.Dependency,
                pyfabm.Coupling,
            ]
        ],
        category: Optional[str] = None,
    ):
        for variable in arr:
            pathcomps = variable.path.split("/")
            if len(pathcomps) < 2:
                # print 'Skipping root level entity %s' % pathcomps[0]
                continue
            if category is not None:
                pathcomps = [pathcomps[0], category] + list(pathcomps[1:])
            parent = self
            for component in pathcomps[:-1]:
                parent = parent.findChild(component)
            entry = Entry(pathcomps[-1], variable)
            parent.addChild(entry)


class Submodel:
    def __init__(self, long_name: str):
        self.units = None
        self.units_unicode = None
        self.value = None
        self.long_name = long_name


class ItemModel(QtCore.QAbstractItemModel):
    def __init__(self, model: pyfabm.Model, parent: Optional[QtCore.QObject]):
        QtCore.QAbstractItemModel.__init__(self, parent)
        self.model = model
        self.root = self._build_tree()

    def _build_tree(self) -> Entry:
        root = Entry()
        env = Entry("environment")
        for d in self.model.dependencies:
            env.addChild(Entry(d.name, d))
        root.addTree(self.model.parameters, "parameters")
        root.addTree(self.model.state_variables, "initialization")
        root.addTree(self.model.couplings, "coupling")
        root.addChild(env)

        # For all models, create an object that returns appropriate metadata.
        def processNode(n: Entry, path: Tuple[str, ...] = ()):
            for i in range(len(n.children) - 1, -1, -1):
                child = n.children[i]
                childpath = path + (child.name,)
                if child.object is None:
                    # Child without underlying object: remove redundant subheaders
                    if childpath[-1] in (
                        "parameters",
                        "initialization",
                        "environment",
                        "coupling",
                    ):
                        childpath = childpath[:-1]
                    else:
                        submodel = self.model.getSubModel("/".join(childpath))
                        if submodel.user_created:
                            child.object = Submodel(submodel.long_name)
                        else:
                            n.removeChild(i)
                            continue
                processNode(child, childpath)

        processNode(root)
        return root

    def rebuild(self) -> None:
        # We already have an old tree - compare and amend model.
        def processChange(newnode: Entry, oldnode: Entry, parent: QtCore.QModelIndex):
            oldnode.object = newnode.object
            if not newnode.children:
                return
            ioldstart = 0
            for node in newnode.children:
                iold = -1
                for i in range(ioldstart, len(oldnode.children)):
                    if node.name == oldnode.children[i].name:
                        iold = i
                        break
                if iold != -1:
                    # New node was found among old nodes;
                    # remove any unused old nodes that precede it.
                    if iold > ioldstart:
                        self.beginRemoveRows(parent, ioldstart, iold - 1)
                        for i in range(iold - 1, ioldstart - 1, -1):
                            oldnode.removeChild(i)
                        self.endRemoveRows()
                    # Process changes to children of node.
                    processChange(
                        node,
                        oldnode.children[ioldstart],
                        self.createIndex(ioldstart, 0, oldnode.children[ioldstart]),
                    )
                else:
                    # New node not found; insert it.
                    self.beginInsertRows(parent, ioldstart, ioldstart)
                    oldnode.insertChild(ioldstart, node)
                    self.endInsertRows()
                ioldstart = ioldstart + 1
            if ioldstart < len(oldnode.children):
                # Remove any trailing unused old nodes.
                self.beginRemoveRows(parent, ioldstart, len(oldnode.children) - 1)
                for i in range(len(oldnode.children) - 1, ioldstart - 1, -1):
                    oldnode.removeChild(i)
                self.endRemoveRows()

        processChange(self._build_tree(), self.root, QtCore.QModelIndex())

    def rowCount(
        self, index: QtCore.QModelIndex = QtCore.QModelIndex()  # type: ignore[override]
    ) -> int:
        if not index.isValid():
            return len(self.root.children)
        elif index.column() == 0:
            entry = cast(Entry, index.internalPointer())
            return len(entry.children)
        return 0

    def columnCount(
        self, index: QtCore.QModelIndex = QtCore.QModelIndex()  # type: ignore[override]
    ) -> int:
        return 4

    def index(
        self,
        row: int,
        column: int,
        parent: QtCore.QModelIndex = QtCore.QModelIndex(),  # type: ignore[override]
    ) -> QtCore.QModelIndex:
        if not parent.isValid():
            # top-level
            children = self.root.children
        else:
            # not top-level
            entry = cast(Entry, parent.internalPointer())
            children = entry.children
        if row < 0 or row >= len(children) or column < 0 or column >= 4:
            return QtCore.QModelIndex()
        return self.createIndex(row, column, children[row])

    def parent(self, index: QtCore.QModelIndex = QtCore.QModelIndex()) -> QtCore.QModelIndex:  # type: ignore[override]
        if not index.isValid():
            return QtCore.QModelIndex()
        entry = cast(Entry, index.internalPointer())
        parent = entry.parent
        assert parent is not None
        if parent.parent is None:
            return QtCore.QModelIndex()
        irow = parent.parent.children.index(parent)
        return self.createIndex(irow, 0, parent)

    def data(
        self,
        index: QtCore.QModelIndex,  # type: ignore[override]
        role: int = QtCore.Qt.ItemDataRole.DisplayRole,
    ):
        if not index.isValid():
            return
        entry = cast(Entry, index.internalPointer())
        data = entry.object
        assert not isinstance(data, Entry)
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            if index.column() == 0:
                return entry.name if data is None else data.long_name
            if isinstance(data, pyfabm.Variable):
                if index.column() == 1:
                    value = data.value
                    if value is None:
                        return "<not set>"
                    if not isinstance(value, bool):
                        if data.units and not isinstance(value, str):
                            return f"{value} {data.units_unicode}"
                        else:
                            return f"{value}"
                elif index.column() == 2 and data.units:
                    return data.units_unicode
                elif index.column() == 3:
                    return entry.name
        elif role == QtCore.Qt.ItemDataRole.ToolTipRole and index.parent().isValid():
            if isinstance(data, pyfabm.Variable):
                return data.long_path
        elif role == QtCore.Qt.ItemDataRole.EditRole:
            assert isinstance(data, pyfabm.Variable)
            return data.value
        elif role == QtCore.Qt.ItemDataRole.FontRole and index.column() == 1:
            if isinstance(data, pyfabm.Parameter) and data.value != data.default:
                font = QtGui.QFont()
                font.setBold(True)
                return font
        elif (
            role == QtCore.Qt.ItemDataRole.CheckStateRole
            and index.column() == 1
            and isinstance(data, pyfabm.Variable)
            and isinstance(data.value, bool)
        ):
            return (
                QtCore.Qt.CheckState.Checked
                if data.value
                else QtCore.Qt.CheckState.Unchecked
            )
        return None

    def setData(
        self,
        index: QtCore.QModelIndex,  # type: ignore[override]
        value: object,
        role: int = QtCore.Qt.ItemDataRole.EditRole,
    ) -> bool:
        if role == QtCore.Qt.ItemDataRole.CheckStateRole:
            assert isinstance(value, int)
            value = QtCore.Qt.CheckState(value) == QtCore.Qt.CheckState.Checked
        if role in (
            QtCore.Qt.ItemDataRole.EditRole,
            QtCore.Qt.ItemDataRole.CheckStateRole,
        ):
            assert isinstance(value, (float, int, bool, str))
            entry: Entry = index.internalPointer()
            data = entry.object
            assert isinstance(data, pyfabm.Variable)
            data.value = value
            if isinstance(data, pyfabm.Parameter):
                self.rebuild()
            return True
        return False

    def flags(self, index: QtCore.QModelIndex):  # type: ignore[override]
        flags = QtCore.Qt.ItemFlags()  # type: ignore[attr-defined]
        if not index.isValid():
            return flags
        if index.column() == 1:
            entry = cast(Entry, index.internalPointer())
            data = entry.object
            if isinstance(data, pyfabm.Variable):
                if isinstance(data.value, bool):
                    flags |= QtCore.Qt.ItemFlag.ItemIsUserCheckable
                else:
                    flags |= QtCore.Qt.ItemFlag.ItemIsEditable
        return (
            flags
            | QtCore.Qt.ItemFlag.ItemIsEnabled
            | QtCore.Qt.ItemFlag.ItemIsSelectable
        )

    def headerData(
        self,
        section: int,
        orientation: QtCore.Qt.Orientation,
        role: int = QtCore.Qt.ItemDataRole.EditRole,
    ):
        if (
            orientation == QtCore.Qt.Orientation.Horizontal
            and role == QtCore.Qt.ItemDataRole.DisplayRole
            and section >= 0
            and section < 4
        ):
            return ("name", "value", "units", "symbol")[section]


class TreeView(QtWidgets.QTreeView):
    def __init__(self, model: pyfabm.Model, parent: Optional[QtWidgets.QWidget]):
        QtWidgets.QTreeView.__init__(self, parent)
        itemmodel = pyfabm.gui_qt.ItemModel(model, parent)
        self.setItemDelegate(Delegate(parent))
        self.setModel(itemmodel)
        self.setUniformRowHeights(True)
        self.expandAll()

        def onTreeViewContextMenu(pos: QtCore.QPoint):
            index = self.indexAt(pos)
            if index.isValid() and index.column() == 1:
                entry = cast(Entry, index.internalPointer())
                data = entry.object
                if (
                    isinstance(data, pyfabm.Parameter)
                    and data.value != data.default
                    and data.default is not None
                ):

                    def reset() -> None:
                        data.reset()
                        itemmodel.rebuild()

                    contextMenu = QtWidgets.QMenu(self)
                    default = data.default
                    if data.units:
                        default = f"{data.default} {data.units_unicode}"
                    contextMenu.addAction(f"Reset to default: {default}", reset)
                    contextMenu.exec(self.mapToGlobal(pos))

        self.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(onTreeViewContextMenu)
        for i in range(3):
            self.resizeColumnToContents(i)


# From xmlstore (GOTM-GUI)


class ScientificDoubleValidator(QtGui.QValidator):
    """Qt validator for floating point values
    Less strict than the standard QDoubleValidator, in the sense that is
    also accepts values in scientific format (e.g. 1.2e6)
    Also has properties 'minimum' and 'maximum', used for validation and
    fix-up.
    """

    def __init__(self, parent: Optional[QtCore.QObject] = None):
        QtGui.QValidator.__init__(self, parent)
        self.minimum: Optional[float] = None
        self.maximum: Optional[float] = None
        self.suffix = ""

    def validate(self, input: str, pos: int) -> Tuple[QtGui.QValidator.State, str, int]:
        assert isinstance(input, str), "input argument is not a string (old PyQt4 API?)"

        # Check for suffix (if ok, cut it off for further value checking)
        if not input.endswith(self.suffix):
            return (QtGui.QValidator.State.Invalid, input, pos)
        vallength = len(input) - len(self.suffix)

        # Check for invalid characters
        if re.match(r"[^\d\-+eE,.]", input[:vallength]):
            return (QtGui.QValidator.State.Invalid, input, pos)

        # Check if we can convert it into a floating point value
        try:
            v = float(input[:vallength])
        except ValueError:
            return (QtGui.QValidator.State.Intermediate, input, pos)

        # Check for minimum and maximum.
        if self.minimum is not None and v < self.minimum:
            return (QtGui.QValidator.State.Intermediate, input, pos)
        if self.maximum is not None and v > self.maximum:
            return (QtGui.QValidator.State.Intermediate, input, pos)

        return (QtGui.QValidator.State.Acceptable, input, pos)

    def fixup(self, input: str):
        assert isinstance(input, str), "input argument is not a string (old PyQt4 API?)"
        if not input.endswith(self.suffix):
            return input

        try:
            v = float(input[: len(input) - len(self.suffix)])
        except ValueError:
            return input

        if self.minimum is not None and v < self.minimum:
            input = f"{self.minimum}{self.suffix}"
        if self.maximum is not None and v > self.maximum:
            input = f"{self.maximum}{self.suffix}"
        print(repr(input))
        return input

    def setSuffix(self, suffix: str):
        self.suffix = suffix


class ScientificDoubleEditor(QtWidgets.QLineEdit):
    """Editor for a floating point value."""

    def __init__(self, parent: Optional[QtWidgets.QWidget]):
        QtWidgets.QLineEdit.__init__(self, parent)

        self.curvalidator = ScientificDoubleValidator(self)
        self.setValidator(self.curvalidator)
        self.suffix = ""
        # self.editingFinished.connect(self.onPropertyEditingFinished)

    def setSuffix(self, suffix: str):
        value = self.value()
        self.suffix = suffix
        self.curvalidator.setSuffix(suffix)
        self.setValue(value)

    def value(self) -> float:
        text = self.text()
        text = text[: len(text) - len(self.suffix)]
        if text == "":
            return 0
        return float(text)

    def setValue(self, value: Optional[float], format: Optional[str] = None):
        if value is None:
            strvalue = ""
        else:
            if format is None:
                strvalue = str(value)
            else:
                strvalue = format % value
        self.setText(f"{strvalue}{self.suffix}")

    def focusInEvent(self, e: QtGui.QFocusEvent):
        QtWidgets.QLineEdit.focusInEvent(self, e)
        self.selectAll()

    def selectAll(self) -> None:
        QtWidgets.QLineEdit.setSelection(self, 0, len(self.text()) - len(self.suffix))

    def setMinimum(self, minimum: Optional[float]):
        self.curvalidator.minimum = minimum

    def setMaximum(self, maximum: Optional[float]):
        self.curvalidator.maximum = maximum

    def interpretText(self) -> None:
        if not self.hasAcceptableInput():
            text = self.text()
            textnew = self.curvalidator.fixup(text)
            if textnew is None:
                textnew = text
            self.setText(textnew)
