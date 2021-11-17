#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
import sys
import os
import ctypes
import re

try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

try:
   import numpy
except ImportError:
   print('Unable to import NumPy. Please ensure it is installed.')
   sys.exit(1)

def find_library(basedir, names):
    for name in names:
        path = os.path.join(basedir, name)
        if os.path.isfile(path):
            return path

name2lib = {}
def get_lib(name):
    if name in name2lib:
       return name2lib[name]

    # Determine potential names of dynamic library.
    if os.name == 'nt':
       names = ('%s.dll' % name, 'lib%s.dll' % name)
    elif os.name == 'posix' and sys.platform == 'darwin':
       names = ('lib%s.dylib' % name,)
    else:
       names = ('lib%s.so' % name,)

    # Find FABM dynamic library.
    # Look first in pyfabm directory, then in Python path.
    path = find_library(os.path.dirname(os.path.abspath(__file__)), names)
    if not path:
        for basedir in sys.path:
            path = find_library(basedir, names)
            if path:
                break
        else:
            raise Exception('Unable to locate dynamic library %s (tried %s).' % (name, ', '.join(names),))

    # Load FABM library.
    lib = ctypes.CDLL(str(path))
    lib.dtype = ctypes.c_double
    lib.numpy_dtype = numpy.dtype(lib.dtype).newbyteorder('=')

    # Driver settings (number of spatial dimensions, depth index)
    lib.get_driver_settings.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.get_driver_settings.restype = ctypes.c_void_p

    ndim_c = ctypes.c_int()
    idepthdim_c = ctypes.c_int()
    ihas_mask = ctypes.c_int()
    lib.get_driver_settings(ctypes.byref(ndim_c), ctypes.byref(idepthdim_c), ctypes.byref(ihas_mask))
    ndim_int = ndim_c.value
    ndim_hz = ndim_int if idepthdim_c.value == -1 else ndim_int - 1
    lib.ndim_int = ndim_int
    lib.ndim_hz = ndim_hz
    lib.idepthdim = idepthdim_c.value
    lib.has_mask = ihas_mask.value != 0

    CONTIGUOUS = str('CONTIGUOUS')
    arrtype0D = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=0, flags=CONTIGUOUS)
    arrtype1D = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=1, flags=CONTIGUOUS)
    arrtypeInterior = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_int, flags=CONTIGUOUS)
    arrtypeHorizontal = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_hz, flags=CONTIGUOUS)
    arrtypeInteriorExt = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_int + 1, flags=CONTIGUOUS)
    arrtypeHorizontalExt = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_hz + 1, flags=CONTIGUOUS)
    arrtypeInteriorExt2 = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_int + 2, flags=CONTIGUOUS)
    arrtypeHorizontalExt2 = numpy.ctypeslib.ndpointer(dtype=lib.dtype, ndim=ndim_hz + 2, flags=CONTIGUOUS)

    # Initialization
    lib.create_model.argtypes = [ctypes.c_char_p] + [ctypes.c_int] * ndim_int
    lib.create_model.restype = ctypes.c_void_p

    # Access to model objects (variables, parameters, dependencies, couplings, model instances)
    lib.get_counts.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.get_counts.restype = None
    lib.get_variable_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.get_variable_metadata.restype = None
    lib.set_variable_save.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
    lib.set_variable_save.restype = None
    lib.get_variable.argtypes = [ctypes.c_void_p, ctypes.c_int,ctypes.c_int]
    lib.get_variable.restype = ctypes.c_void_p
    lib.get_parameter_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
    lib.get_parameter_metadata.restype = None
    lib.get_model_metadata.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
    lib.get_model_metadata.restype = None
    lib.get_coupling.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_void_p)]
    lib.get_coupling.restype = None
    lib.get_error_state.argtypes = []
    lib.get_error_state.restype = ctypes.c_int
    lib.get_error.argtypes = [ctypes.c_int, ctypes.c_char_p]
    lib.get_error.restype = None
    lib.reset_error_state.argtypes = []
    lib.reset_error_state.restype = None
    lib.configure.argtypes = [ctypes.c_int]
    lib.configure.restype = None
    if lib.has_mask:
        lib.set_mask.restype = None
        lib.set_mask.argtypes = [ctypes.c_void_p, numpy.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=ndim_hz, flags=CONTIGUOUS)]

    # Read access to variable attributes
    lib.variable_get_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.variable_get_metadata.restype = None
    lib.variable_get_background_value.argtypes = [ctypes.c_void_p]
    lib.variable_get_background_value.restype = lib.dtype
    lib.variable_get_missing_value.argtypes = [ctypes.c_void_p]
    lib.variable_get_missing_value.restype = lib.dtype
    lib.variable_get_long_path.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
    lib.variable_get_long_path.restype = None
    lib.variable_get_output_name.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
    lib.variable_get_output_name.restype = None
    lib.variable_get_suitable_masters.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.variable_get_suitable_masters.restype = ctypes.c_void_p
    lib.variable_get_output.argtypes = [ctypes.c_void_p]
    lib.variable_get_output.restype = ctypes.c_int
    lib.variable_is_required.argtypes = [ctypes.c_void_p]
    lib.variable_is_required.restype = ctypes.c_int
    lib.variable_get_property_type.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.variable_get_property_type.restype = ctypes.c_int
    lib.variable_get_real_property.argtypes = [ctypes.c_void_p, ctypes.c_char_p, lib.dtype]
    lib.variable_get_real_property.restype = lib.dtype
    lib.variable_get_integer_property.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.variable_get_integer_property.restype = ctypes.c_int
    lib.variable_get_logical_property.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.variable_get_logical_property.restype = ctypes.c_int

    # Read/write/reset access to parameters.
    lib.get_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_real_parameter.restype = lib.dtype
    lib.get_integer_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_integer_parameter.restype = ctypes.c_int
    lib.get_logical_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
    lib.get_logical_parameter.restype = ctypes.c_int
    lib.get_string_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_char_p]
    lib.get_string_parameter.restype = None
    lib.reset_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.reset_parameter.restype = None
    lib.set_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, lib.dtype]
    lib.set_real_parameter.restype = None
    lib.set_integer_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.set_integer_parameter.restype = None
    lib.set_logical_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
    lib.set_logical_parameter.restype = None
    lib.set_string_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
    lib.set_string_parameter.restype = None

    # Read access to lists of variables (e.g., suitable coupling targets).
    lib.link_list_count.argtypes = [ctypes.c_void_p]
    lib.link_list_count.restype = ctypes.c_int
    lib.link_list_index.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.link_list_index.restype = ctypes.c_void_p
    lib.link_list_finalize.argtypes = [ctypes.c_void_p]
    lib.link_list_finalize.restype = None

    # Routines for sending pointers to state and dependency data.
    lib.link_interior_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, arrtypeInterior]
    lib.link_interior_state_data.restype = None
    lib.link_surface_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, arrtypeHorizontal]
    lib.link_surface_state_data.restype = None
    lib.link_bottom_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, arrtypeHorizontal]
    lib.link_bottom_state_data.restype = None
    lib.link_interior_data.argtypes = [ctypes.c_void_p, ctypes.c_void_p, arrtypeInterior]
    lib.link_interior_data.restype = None
    lib.link_horizontal_data.argtypes = [ctypes.c_void_p, ctypes.c_void_p, arrtypeHorizontal]
    lib.link_horizontal_data.restype = None
    lib.link_scalar.argtypes = [ctypes.c_void_p, ctypes.c_void_p, arrtype0D]
    lib.link_scalar.restype = None

    # Read access to diagnostic data.
    lib.get_interior_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.get_interior_diagnostic_data.restype = ctypes.POINTER(lib.dtype)
    lib.get_horizontal_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.get_horizontal_diagnostic_data.restype = ctypes.POINTER(lib.dtype)

    lib.start.argtypes = [ctypes.c_void_p]
    lib.start.restype = None

    # Routine for retrieving source-sink terms for the interior domain.
    lib.get_sources.argtypes = [ctypes.c_void_p, lib.dtype, arrtypeInteriorExt, arrtypeHorizontalExt, arrtypeHorizontalExt, ctypes.c_int, ctypes.c_int, arrtypeInterior]
    lib.get_sources.restype = None
    lib.check_state.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib.check_state.restype = ctypes.c_int

    # Routine for getting git repository version information.
    lib.get_version.argtypes = (ctypes.c_int, ctypes.c_char_p)
    lib.get_version.restype = None

    if ndim_int == 0:
        lib.integrate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, arrtype1D, arrtypeInteriorExt, arrtypeInteriorExt2, lib.dtype, ctypes.c_int, ctypes.c_int, arrtypeInterior]
        lib.integrate.restype = None

    name2lib[name] = lib
    return lib

INTERIOR_STATE_VARIABLE        = 1
SURFACE_STATE_VARIABLE         = 2
BOTTOM_STATE_VARIABLE          = 3
INTERIOR_DIAGNOSTIC_VARIABLE   = 4
HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
CONSERVED_QUANTITY             = 6
INTERIOR_DEPENDENCY            = 7
HORIZONTAL_DEPENDENCY          = 8
SCALAR_DEPENDENCY              = 9
ATTRIBUTE_LENGTH               = 256

unicodesuperscript = {'1':'\u00B9','2':'\u00B2','3':'\u00B3',
                      '4':'\u2074','5':'\u2075','6':'\u2076',
                      '7':'\u2077','8':'\u2078','9':'\u2079',
                      '0':'\u2070','-':'\u207B'}
unicodesubscript = {'1':'\u2081','2':'\u2082','3':'\u2083',
                    '4':'\u2084','5':'\u2085','6':'\u2086',
                    '7':'\u2087','8':'\u2088','9':'\u2089',
                    '0':'\u2080'}
supnumber = re.compile(r'(?<=\w)(-?\d+)(?=[ \*+\-/]|$)')
supenumber = re.compile(r'(?<=\d)e(-?\d+)(?=[ \*+\-/]|$)')
oldsupminus = re.compile(r'/(\w+)(?:\*\*|\^)(\d+)(?=[ \*+\-/]|$)')
oldsup = re.compile(r'(?<=\w)(?:\*\*|\^)(-?\d+)(?=[ \*+\-/]|$)')
oldsub = re.compile(r'(?<=\w)_(-?\d+)(?=[ \*+\-/]|$)')
def createPrettyUnit(unit):
    def replace_superscript(m):
        return ''.join([unicodesuperscript[n] for n in m.group(1)])
    def replace_subscript(m):
        return ''.join([unicodesubscript[n] for n in m.group(1)])
    def reple(m):
        return '\u00D710%s' % ''.join([unicodesuperscript[n] for n in m.group(1)])
    def reploldminus(m):
        return ' %s\u207B%s' % (m.group(1),''.join([unicodesuperscript[n] for n in m.group(2)]))
    #def replold(m):
    #    return u'%s%s' % (m.group(1),u''.join([unicodesuperscript[n] for n in m.group(2)]))
    unit = oldsup.sub(replace_superscript,unit)
    unit = oldsub.sub(replace_subscript,unit)
    unit = supenumber.sub(reple,unit)
    unit = supnumber.sub(replace_superscript,unit)
    #unit = oldsupminus.sub(reploldminus,unit)
    return unit

class FABMException(Exception):
    pass

def hasError():
   for lib in name2lib.values():
       if lib.get_error_state() != 0:
           return True
   return False

def getError():
   for lib in name2lib.values():
       if lib.get_error_state() != 0:
           strmessage = ctypes.create_string_buffer(1024)
           lib.get_error(1024, strmessage)
           return strmessage.value.decode('ascii')

def printTree(root, stringmapper, indent=''):
    """Print an indented tree of objects, encoded by dictionaries linking the names of children to
    their subtree, or to their object. Objects are finally printed as string obtained by
    calling the provided stringmapper method."""
    for name, item in root.items():
        if isinstance(item, dict):
            print('%s%s' % (indent, name))
            printTree(item, stringmapper, indent + '   ')
        else:
            print('%s%s = %s' % (indent, name, stringmapper(item)))

class VariableProperties:
    def __init__(self, model, variable_pointer):
        self.model = model
        self.variable_pointer = variable_pointer

    def __getitem__(self, key):
        typecode = self.model.fabm.variable_get_property_type(self.variable_pointer, key.encode('ascii'))
        if typecode == 1:
            return self.model.fabm.variable_get_real_property(self.variable_pointer, key.encode('ascii'), -1.)
        elif typecode == 2:
            return self.model.fabm.variable_get_integer_property(self.variable_pointer, key.encode('ascii'), 0)
        elif typecode == 3:
            return self.model.fabm.variable_get_logical_property(self.variable_pointer, key.encode('ascii'), 0) != 0
        raise KeyError

class Variable(object):
    def __init__(self, model, name=None, units=None, long_name=None, path=None, variable_pointer=None):
        self.model = model
        self.variable_pointer = variable_pointer
        if variable_pointer:
           strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           self.model.fabm.variable_get_metadata(variable_pointer, ATTRIBUTE_LENGTH, strname, strunits, strlong_name)
           name = strname.value.decode('ascii')
           units = strunits.value.decode('ascii')
           long_name = strlong_name.value.decode('ascii')

        self.name = name
        self.units = units
        self.units_unicode = None if units is None else createPrettyUnit(units)
        self.long_name = long_name or name
        self.path = path or name
        self.properties = VariableProperties(self.model, self.variable_pointer)

    @property
    def long_path(self):
        if self.variable_pointer is None:
            return self.long_name
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_long_path(self.variable_pointer, ATTRIBUTE_LENGTH, strlong_name)
        return strlong_name.value.decode('ascii')

    @property
    def output_name(self):
        if self.variable_pointer is None:
            return self.name
        stroutput_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_output_name(self.variable_pointer, ATTRIBUTE_LENGTH, stroutput_name)
        return stroutput_name.value.decode('ascii')

    @property
    def missing_value(self):
        if self.variable_pointer is not None:
            return self.model.fabm.variable_get_missing_value(self.variable_pointer)

    def getOptions(self):
        pass

    def getRealProperty(self, name, default=-1.0):
        return self.model.fabm.variable_get_real_property(self.variable_pointer, name.encode('ascii'), default)

    def __repr__(self):
        return '<%s=%s>' % (self.name, self.value)

class Dependency(Variable):
    def __init__(self, model, variable_pointer, shape, link_function):
        Variable.__init__(self, model, variable_pointer=variable_pointer)
        self.is_set = False
        self.link_function = link_function
        self.shape = shape

    def getValue(self):
        return None if not self.is_set else self.data

    def setValue(self, value):
        if not self.is_set:
            self.link(numpy.empty(self.shape, dtype=float))
        self.data[...] = value

    def link(self, data):
        assert data.shape == self.shape, '%s: shape of provided array %s does not match the shape required %s' % (self.name, data.shape, self.shape)
        self.data = data
        self.link_function(self.model.pmodel, self.variable_pointer, self.data)
        self.is_set = True

    value = property(getValue, setValue)

    @property
    def required(self):
        return self.model.fabm.variable_is_required(self.variable_pointer) != 0

class StateVariable(Variable):
    def __init__(self, model, variable_pointer, data):
        Variable.__init__(self, model, variable_pointer=variable_pointer)
        self.data = data

    def getValue(self):
        return self.data

    def setValue(self, value):
        self.data[...] = value

    value = property(getValue, setValue)

    @property
    def background_value(self):
        return self.model.fabm.variable_get_background_value(self.variable_pointer)

    @property
    def output(self):
        return self.model.fabm.variable_get_output(self.variable_pointer) != 0

class DiagnosticVariable(Variable):
    def __init__(self, model, variable_pointer, index, horizontal):
        Variable.__init__(self, model, variable_pointer=variable_pointer)
        self.data = None
        self.horizontal = horizontal
        self.index = index

    def getValue(self):
        return self.data
    value = property(getValue)

    @property
    def output(self):
        return self.model.fabm.variable_get_output(self.variable_pointer) != 0

    def setSave(self, value):
        self.model.fabm.set_variable_save(self.model.pmodel, HORIZONTAL_DIAGNOSTIC_VARIABLE if self.horizontal else INTERIOR_DIAGNOSTIC_VARIABLE, self.index, 1 if value else 0)
    save = property(fset=setSave)

class Parameter(Variable):
    def __init__(self, model, name, index, units=None, long_name=None, type=None, has_default=False):
        Variable.__init__(self, model, name, units, long_name)
        self.type = type
        self.index = index + 1
        self.has_default = has_default

    def getValue(self, default=False):
        default = 1 if default else 0
        if self.type == 1:
            return self.model.fabm.get_real_parameter(self.model.pmodel, self.index, default)
        elif self.type == 2:
            return self.model.fabm.get_integer_parameter(self.model.pmodel, self.index, default)
        elif self.type == 3:
            return self.model.fabm.get_logical_parameter(self.model.pmodel, self.index, default) != 0
        elif self.type == 4:
            result = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
            self.model.fabm.get_string_parameter(self.model.pmodel, self.index, default, ATTRIBUTE_LENGTH, result)
            return result.value.decode('ascii')

    def setValue(self, value):
        settings = self.model.saveSettings()

        if self.type == 1:
            self.model.fabm.set_real_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type == 2:
            self.model.fabm.set_integer_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type == 3:
            self.model.fabm.set_logical_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type == 4:
            self.model.fabm.set_string_parameter(self.model.pmodel, self.name.encode('ascii'), value)

        # Update the model configuration (arrays with variables and parameters have changed)
        self.model.updateConfiguration(settings)

    def getDefault(self):
        if not self.has_default:
            return None
        return self.getValue(True)

    def reset(self):
        settings = self.model.saveSettings()
        self.model.fabm.reset_parameter(self.model.pmodel, self.index)
        self.model.updateConfiguration(settings)

    value = property(getValue, setValue)
    default = property(getDefault)

class NamedObjectList(Sequence):
    def __init__(self, *data):
        self._data = []
        for d in data:
            self._data.extend(d)
        self._lookup = None
        self._lookup_ci = None

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.find(key)
        return self._data[key]

    def __repr__(self):
        return repr(self._data)

    def __add__(self, other):
        return NamedObjectList(self._data, other._data)

    def find(self, name, case_insensitive=False):
        if self._lookup is None:
            self._lookup_ci = dict([(obj.name.lower(), obj) for obj in self._data])
            self._lookup = dict([(obj.name, obj) for obj in self._data])
        if case_insensitive:
            return self._lookup_ci[name.lower()]
        return self._lookup[name]

class Coupling(Variable):
    def __init__(self, model, index):
        self.model = model
        self.master = ctypes.c_void_p()
        self.slave = ctypes.c_void_p()
        self.model.fabm.get_coupling(self.model.pmodel, index, ctypes.byref(self.slave), ctypes.byref(self.master))
        Variable.__init__(self, model, variable_pointer=self.slave)

    def getValue(self):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_long_path(self.master, ATTRIBUTE_LENGTH, strlong_name)
        return strlong_name.value.decode('ascii')

    def setValue(self,value):
        print('New coupling specified: %s' % value)
        pass

    def getOptions(self):
        options = []
        list = self.model.fabm.variable_get_suitable_masters(self.model.pmodel, self.slave)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        for i in range(self.model.fabm.link_list_count(list)):
            variable = self.model.fabm.link_list_index(list, i + 1)
            self.model.fabm.variable_get_long_path(variable, ATTRIBUTE_LENGTH, strlong_name)
            options.append(strlong_name.value.decode('ascii'))
        self.model.fabm.link_list_finalize(list)
        return options

    @property
    def long_path(self):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        self.model.fabm.variable_get_long_path(self.slave, ATTRIBUTE_LENGTH, strlong_name)
        return strlong_name.value.decode('ascii')

    value = property(getValue, setValue)

class SubModel(object):
    def __init__(self, model, name):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        iuser = ctypes.c_int()
        model.fabm.get_model_metadata(model.pmodel, name.encode('ascii'), ATTRIBUTE_LENGTH, strlong_name, iuser)
        self.long_name = strlong_name.value.decode('ascii')
        self.user_created = iuser.value != 0

class Model(object):
    def __init__(self, path='fabm.yaml', shape=(), libname=None):
        delete = False
        if isinstance(path, dict):
            import tempfile
            import yaml
            import io
            with tempfile.NamedTemporaryFile(suffix='.yaml', prefix='fabm', delete=False) as f:
                yaml.safe_dump(path, io.TextIOWrapper(f, encoding='ascii'))
                path = f.name
            delete = True

        if libname is None:
            if len(shape) > 1:
                raise FABMException('Invalid domain shape %s. Domain must have 0 or 1 dimensions.' % (shape,))
            libname = {0: 'fabm_0d', 1: 'fabm_1d'}[len(shape)]
        self.fabm = get_lib(libname)
        self.fabm.reset_error_state()
        self._cell_thickness = None
        self.pmodel = self.fabm.create_model(path.encode('ascii'), *shape[::-1])
        if hasError():
            raise FABMException('An error occurred while parsing %s:\n%s' % (path, getError()))
        self.interior_domain_shape = tuple(shape)
        self.horizontal_domain_shape = tuple([l for i, l in enumerate(self.interior_domain_shape) if i != self.fabm.idepthdim])
        if delete:
            os.remove(path)
        self.updateConfiguration()
        self._mask = None

    def link_mask(self, data):
        assert data.shape == self.horizontal_domain_shape and data.dtype == numpy.intc and data.flags['C_CONTIGUOUS']
        self._mask = data
        self.fabm.set_mask(self.pmodel, self._mask)

    def _get_mask(self):
        return self._mask
    def _set_mask(self, value):
        if self._mask is None:
            self.link_mask(numpy.ones(self.horizontal_domain_shape, dtype=numpy.intc))
        if value is not self._mask:
            self._mask[...] = value
    mask = property(_get_mask, _set_mask)

    def _get_state(self):
        return self._state
    def _set_state(self, value):
        if value is not self._state:
            self._state[...] = value
    state = property(_get_state, _set_state)

    def _get_interior_state(self):
        return self._interior_state
    def _set_interior_state(self, value):
        if value is not self._interior_state:
            self._interior_state[...] = value
    interior_state = property(_get_interior_state, _set_interior_state)

    def _get_surface_state(self):
        return self._surface_state
    def _set_surface_state(self, value):
        if value is not self._surface_state:
            self._surface_state[...] = value
    surface_state = property(_get_surface_state, _set_surface_state)

    def _get_bottom_state(self):
        return self._bottom_state
    def _set_bottom_state(self, value):
        if value is not self._bottom_state:
            self._bottom_state[...] = value
    bottom_state = property(_get_bottom_state, _set_bottom_state)

    def link_cell_thickness(self, data):
        assert data.shape == self.interior_domain_shape and data.dtype == self.fabm.numpy_dtype and data.flags['C_CONTIGUOUS']
        self._cell_thickness = data

    def setCellThickness(self, value):
        if self._cell_thickness is None:
            self.link_cell_thickness(numpy.empty(self.interior_domain_shape))
        self._cell_thickness[...] = value

    cell_thickness = property(fset=setCellThickness)

    def getSubModel(self,name):
        return SubModel(self, name)

    def saveSettings(self):
        environment = dict([(dependency.name, dependency.value) for dependency in self.dependencies])
        state = dict([(variable.name,variable.value) for variable in self.state_variables])
        return environment,state

    def restoreSettings(self, data):
        environment,state = data
        for dependency in self.dependencies:
            if dependency.name in environment:
                dependency.value = environment[dependency.name]
        for variable in self.state_variables:
            if variable.name in state:
                variable.value = state[variable.name]

    def updateConfiguration(self, settings=None):
        # Get number of model variables per category
        nstate_interior = ctypes.c_int()
        nstate_surface = ctypes.c_int()
        nstate_bottom = ctypes.c_int()
        ndiag_interior = ctypes.c_int()
        ndiag_horizontal = ctypes.c_int()
        ndependencies_interior = ctypes.c_int()
        ndependencies_horizontal = ctypes.c_int()
        ndependencies_scalar = ctypes.c_int()
        nconserved = ctypes.c_int()
        nparameters = ctypes.c_int()
        ncouplings = ctypes.c_int()
        self.fabm.get_counts(self.pmodel,
            ctypes.byref(nstate_interior), ctypes.byref(nstate_surface), ctypes.byref(nstate_bottom),
            ctypes.byref(ndiag_interior), ctypes.byref(ndiag_horizontal),
            ctypes.byref(ndependencies_interior), ctypes.byref(ndependencies_horizontal), ctypes.byref(ndependencies_scalar),
            ctypes.byref(nconserved), ctypes.byref(nparameters), ctypes.byref(ncouplings)
        )

        # Allocate memory for state variable values, and send ctypes.pointer to this memory to FABM.
        if self.fabm.idepthdim == -1:
            self._state = numpy.empty((nstate_interior.value + nstate_surface.value + nstate_bottom.value,) + self.interior_domain_shape, dtype=self.fabm.numpy_dtype)
            self._interior_state = self._state[:nstate_interior.value, ...]
            self._surface_state = self._state[nstate_interior.value:nstate_interior.value + nstate_surface.value, ...]
            self._bottom_state = self._state[nstate_interior.value + nstate_surface.value:, ...]
        else:
            self._interior_state = numpy.empty((nstate_interior.value,) + self.interior_domain_shape, dtype=self.fabm.numpy_dtype)
            self._surface_state = numpy.empty((nstate_surface.value,) + self.horizontal_domain_shape, dtype=self.fabm.numpy_dtype)
            self._bottom_state = numpy.empty((nstate_bottom.value,) + self.horizontal_domain_shape, dtype=self.fabm.numpy_dtype)
        for i in range(nstate_interior.value):
            self.fabm.link_interior_state_data(self.pmodel, i + 1, self._interior_state[i, ...])
        for i in range(nstate_surface.value):
            self.fabm.link_surface_state_data(self.pmodel, i + 1, self._surface_state[i, ...])
        for i in range(nstate_bottom.value):
            self.fabm.link_bottom_state_data(self.pmodel, i + 1, self._bottom_state[i, ...])

        # Retrieve variable metadata
        strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strpath = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        typecode = ctypes.c_int()
        has_default = ctypes.c_int()
        required = ctypes.c_int()
        self.interior_state_variables = NamedObjectList()
        self.surface_state_variables = NamedObjectList()
        self.bottom_state_variables = NamedObjectList()
        self.interior_diagnostic_variables = NamedObjectList()
        self.horizontal_diagnostic_variables = NamedObjectList()
        self.conserved_quantities = NamedObjectList()
        self.parameters = NamedObjectList()
        self.interior_dependencies = NamedObjectList()
        self.horizontal_dependencies = NamedObjectList()
        self.scalar_dependencies = NamedObjectList()
        for i in range(nstate_interior.value):
            ptr = self.fabm.get_variable(self.pmodel, INTERIOR_STATE_VARIABLE, i + 1)
            self.interior_state_variables._data.append(StateVariable(self, ptr, self._interior_state[i, ...]))
        for i in range(nstate_surface.value):
            ptr = self.fabm.get_variable(self.pmodel, SURFACE_STATE_VARIABLE, i + 1)
            self.surface_state_variables._data.append(StateVariable(self, ptr, self._surface_state[i, ...]))
        for i in range(nstate_bottom.value):
            ptr = self.fabm.get_variable(self.pmodel, BOTTOM_STATE_VARIABLE, i + 1)
            self.bottom_state_variables._data.append(StateVariable(self, ptr, self._bottom_state[i, ...]))
        for i in range(ndiag_interior.value):
            ptr = self.fabm.get_variable(self.pmodel, INTERIOR_DIAGNOSTIC_VARIABLE, i + 1)
            self.interior_diagnostic_variables._data.append(DiagnosticVariable(self, ptr, i + 1, False))
        for i in range(ndiag_horizontal.value):
            ptr = self.fabm.get_variable(self.pmodel, HORIZONTAL_DIAGNOSTIC_VARIABLE, i + 1)
            self.horizontal_diagnostic_variables._data.append(DiagnosticVariable(self, ptr, i + 1, True))
        for i in range(ndependencies_interior.value):
            ptr = self.fabm.get_variable(self.pmodel, INTERIOR_DEPENDENCY, i + 1)
            self.interior_dependencies._data.append(Dependency(self, ptr, self.interior_domain_shape, self.fabm.link_interior_data))
        for i in range(ndependencies_horizontal.value):
            ptr = self.fabm.get_variable(self.pmodel, HORIZONTAL_DEPENDENCY, i + 1)
            self.horizontal_dependencies._data.append(Dependency(self, ptr, self.horizontal_domain_shape, self.fabm.link_horizontal_data))
        for i in range(ndependencies_scalar.value):
            ptr = self.fabm.get_variable(self.pmodel, SCALAR_DEPENDENCY, i + 1)
            self.scalar_dependencies._data.append(Dependency(self, ptr, (), self.fabm.link_scalar))
        for i in range(nconserved.value):
            self.fabm.get_variable_metadata(self.pmodel, CONSERVED_QUANTITY, i + 1, ATTRIBUTE_LENGTH, strname, strunits, strlong_name, strpath)
            self.conserved_quantities._data.append(Variable(self, strname.value.decode('ascii'), strunits.value.decode('ascii'), strlong_name.value.decode('ascii'), strpath.value.decode('ascii')))
        for i in range(nparameters.value):
            self.fabm.get_parameter_metadata(self.pmodel, i + 1, ATTRIBUTE_LENGTH, strname, strunits, strlong_name, ctypes.byref(typecode), ctypes.byref(has_default))
            self.parameters._data.append(Parameter(self, strname.value.decode('ascii'), i, type=typecode.value, units=strunits.value.decode('ascii'), long_name=strlong_name.value.decode('ascii'), has_default=has_default.value != 0))

        self.couplings = NamedObjectList([Coupling(self, i + 1) for i in range(ncouplings.value)])

        # Arrays that combine variables from pelagic and boundary domains.
        self.state_variables = self.interior_state_variables + self.surface_state_variables + self.bottom_state_variables
        self.diagnostic_variables = self.interior_diagnostic_variables + self.horizontal_diagnostic_variables
        self.dependencies = self.interior_dependencies + self.horizontal_dependencies + self.scalar_dependencies

        if settings is not None:
            self.restoreSettings(settings)

        # For backward compatibility
        self.bulk_state_variables = self.interior_state_variables
        self.bulk_diagnostic_variables = self.interior_diagnostic_variables

        self.itime = -1.

    def getRates(self, t=None, surface=True, bottom=True):
        """Returns the local rate of change in state variables,
        given the current state and environment.
        """
        assert self.fabm.idepthdim == -1
        if t is None:
            t = self.itime
        sources = numpy.empty_like(self._state)
        sources_interior = sources[:len(self.interior_state_variables), ...]
        sources_surface = sources[len(self.interior_state_variables):len(self.interior_state_variables) + len(self.surface_state_variables), ...]
        sources_bottom = sources[len(self.interior_state_variables) + len(self.surface_state_variables):, ...]
        assert not ((surface or bottom) and self._cell_thickness is None), 'You must assign model.cell_thickness to use getRates'
        self.fabm.get_sources(self.pmodel, t, sources_interior, sources_surface, sources_bottom, surface, bottom, self._cell_thickness)
        if hasError():
            raise FABMException(getError())
        return sources

    def get_sources(self, t=None, out=None):
        if t is None:
            t = self.itime
        if out is None:
            sources_interior = numpy.empty_like(self._interior_state)
            sources_surface = numpy.empty_like(self._surface_state)
            sources_bottom = numpy.empty_like(self._bottom_state)
        else:
            sources_interior, sources_surface, sources_bottom = out
        assert self._cell_thickness is not None, 'You must assign model.cell_thickness to use get_sources'
        self.fabm.get_sources(self.pmodel, t, sources_interior, sources_surface, sources_bottom, True, True, self._cell_thickness)
        if hasError():
            raise FABMException(getError())
        return sources_interior, sources_surface, sources_bottom

    def checkState(self, repair=False):
        valid = self.fabm.check_state(self.pmodel, repair) != 0
        if hasError():
            raise FABMException(getError())
        return valid

    def getJacobian(self,pert=None):
        # Define perturbation per state variable.
        y_pert = numpy.empty_like(self.state)
        if pert is None: pert = 1e-6
        y_pert[:] = pert

        # Compute dy for original state (used as reference for finite differences later on)
        dy_ori = self.getRates()

        # Create memory for Jacobian
        Jac = numpy.empty((len(self.state), len(self.state)), dtype=self.state.dtype)

        for i in range(len(self.state)):
            # Save original state variable value, create perturbed one.
            y_ori = self.state[i]
            self.state[i] += y_pert[i]

            # Compute dy for perturbed state, compute Jacobian elements using finite difference.
            dy_pert = self.getRates()
            Jac[:,i] = (dy_pert - dy_ori) / y_pert[i]

            # Restore original state variable value.
            self.state[i] = y_ori

        return Jac

    def findParameter(self,name,case_insensitive=False):
        return self.parameters.find(name, case_insensitive)

    def findDependency(self,name,case_insensitive=False):
        return self.dependencies.find(name, case_insensitive)

    def findStateVariable(self,name,case_insensitive=False):
        return self.state_variables.find(name, case_insensitive)

    def findDiagnosticVariable(self,name,case_insensitive=False):
        return self.diagnostic_variables.find(name, case_insensitive)

    def findCoupling(self,name,case_insensitive=False):
        return self.couplings.find(name, case_insensitive)

    def getParameterTree(self):
        root = {}
        for parameter in self.parameters:
            pathcomps = parameter.name.split('/')
            parent = root
            for component in pathcomps[:-1]:
                parent = root.setdefault(component,{})
            parent[pathcomps[-1]] = parameter
        return root

    def start(self, verbose=True, stop=False):
        ready = True
        if self.fabm.has_mask:
            if self._mask is None:
                print('Mask not yet assigned')
                ready = False

        def process_dependencies(dependencies, link_function):
            ready = True
            for dependency in dependencies:
                if dependency.required and not dependency.is_set:
                    print('Value for dependency %s is not set.' % dependency.name)
                    ready = False
            return ready

        ready = process_dependencies(self.interior_dependencies, self.fabm.link_interior_data) and ready
        ready = process_dependencies(self.horizontal_dependencies, self.fabm.link_horizontal_data) and ready
        ready = process_dependencies(self.scalar_dependencies, self.fabm.link_scalar) and ready
        assert ready or not stop, 'Not all dependencies have been fulfilled.'

        self.fabm.start(self.pmodel)
        if hasError():
            return False
        for i, variable in enumerate(self.interior_diagnostic_variables):
            pdata = self.fabm.get_interior_diagnostic_data(self.pmodel, i + 1)
            variable.data = None if not pdata else numpy.ctypeslib.as_array(pdata, self.interior_domain_shape).newbyteorder('=')
        for i, variable in enumerate(self.horizontal_diagnostic_variables):
            pdata = self.fabm.get_horizontal_diagnostic_data(self.pmodel, i + 1)
            variable.data = None if not pdata else numpy.ctypeslib.as_array(pdata, self.horizontal_domain_shape).newbyteorder('=')
        return ready
    checkReady = start

    def updateTime(self, nsec):
       self.itime = nsec

    def printInformation(self):
        """Show information about the model."""
        def printArray(classname, array):
            if not array:
                return
            print(' %i %s:' % (len(array), classname))
            for variable in array:
                print('    %s = %s %s' % (variable.name, variable.value, variable.units))

        print('FABM model contains the following:')
        printArray('interior state variables', self.interior_state_variables)
        printArray('bottom state variables', self.bottom_state_variables)
        printArray('surface state variables', self.surface_state_variables)
        printArray('interior diagnostic variables', self.interior_diagnostic_variables)
        printArray('horizontal diagnostic variables', self.horizontal_diagnostic_variables)
        printArray('external variables', self.dependencies)
        print(' %i parameters:' % len(self.parameters))
        printTree(self.getParameterTree(), lambda x:'%s %s' % (x.value, x.units), '    ')

class Simulator(object):
    def __init__(self, model):
        assert model._cell_thickness is not None, 'You must assign model.cell_thickness to use Simulator'
        self.model = model

    def integrate(self, y0, t, dt, surface=True, bottom=True):
        y = numpy.empty((t.size, self.model.state.size))
        self.model.fabm.integrate(self.model.pmodel, t.size, self.model.state.size, t, y0, y, dt, surface, bottom, ctypes.byref(self.model._cell_thickness))
        if hasError():
            raise FABMException(getError())
        return y

def unload():
    global name2lib, ctypes

    for lib in name2lib.values():
        handle = lib._handle
        if os.name == 'nt':
            import ctypes.wintypes
            ctypes.windll.kernel32.FreeLibrary.argtypes = [ctypes.wintypes.HMODULE]
            ctypes.windll.kernel32.FreeLibrary(handle)
        else:
            dlclose = ctypes.CDLL(None).dlclose
            dlclose.argtypes = [ctypes.c_void_p]
            dlclose.restype = ctypes.c_int
            dlclose(handle)
    name2lib = {}

def get_version():
    for lib in name2lib.values():
        version_length = 256
        strversion = ctypes.create_string_buffer(version_length)
        lib.get_version(version_length, strversion)
        return strversion.value.decode('ascii')

