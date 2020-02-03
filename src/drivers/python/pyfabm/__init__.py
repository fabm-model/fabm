#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
import sys
import os
import ctypes
import re

try:
   import numpy
except ImportError:
   print('Unable to import NumPy. Please ensure it is installed.')
   sys.exit(1)

# Determine potential names of FABM dynamic library.
if os.name == 'nt':
   dllpaths = ('python_fabm.dll', 'libpython_fabm.dll')
elif os.name == "posix" and sys.platform == "darwin":
   dllpaths = ('libpython_fabm.dylib',)
else:
   dllpaths = ('libpython_fabm.so',)

def find_library(basedir):
    for dllpath in dllpaths:
        dllpath = os.path.join(basedir, dllpath)
        if os.path.isfile(dllpath):
            return dllpath

# Find FABM dynamic library.
# Look first in pyfabm directory, then in Python path.
dllpath = find_library(os.path.dirname(os.path.abspath(__file__)))
if not dllpath:
    for basedir in sys.path:
        dllpath = find_library(basedir)
        if dllpath:
            break

if not dllpath:
   print('Unable to locate FABM dynamic library %s.' % (' or '.join(dllpaths),))
   sys.exit(1)

# Load FABM library.
fabm = ctypes.CDLL(str(dllpath))

# Initialization
fabm.create_model.argtypes = [ctypes.c_char_p]
fabm.create_model.restype = ctypes.c_void_p

CONTIGUOUS = str('CONTIGUOUS')

# Access to model objects (variables, parameters, dependencies, couplings, model instances)
fabm.get_counts.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
fabm.get_counts.restype = None
fabm.get_variable_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
fabm.get_variable_metadata.restype = None
fabm.get_variable.argtypes = [ctypes.c_void_p, ctypes.c_int,ctypes.c_int]
fabm.get_variable.restype = ctypes.c_void_p
fabm.get_parameter_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int,ctypes.c_int,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
fabm.get_parameter_metadata.restype = None
fabm.get_dependency_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int,ctypes.c_int,ctypes.c_char_p,ctypes.c_char_p]
fabm.get_dependency_metadata.restype = None
fabm.get_model_metadata.argtypes = [ctypes.c_void_p, ctypes.c_char_p,ctypes.c_int,ctypes.c_char_p,ctypes.POINTER(ctypes.c_int)]
fabm.get_model_metadata.restype = None
fabm.get_coupling.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_void_p)]
fabm.get_coupling.restype = None
fabm.get_error_state.argtypes = []
fabm.get_error_state.restype = ctypes.c_int
fabm.get_error.argtypes = [ctypes.c_int, ctypes.c_char_p]
fabm.get_error.restype = None

# Read access to variable attributes
fabm.variable_get_metadata.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
fabm.variable_get_metadata.restype = None
fabm.variable_get_background_value.argtypes = [ctypes.c_void_p]
fabm.variable_get_background_value.restype = ctypes.c_double
fabm.variable_get_long_path.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
fabm.variable_get_long_path.restype = None
fabm.variable_get_output_name.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p]
fabm.variable_get_output_name.restype = None
fabm.variable_get_suitable_masters.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
fabm.variable_get_suitable_masters.restype = ctypes.c_void_p
fabm.variable_get_output.argtypes = [ctypes.c_void_p]
fabm.variable_get_output.restype = ctypes.c_int
fabm.variable_get_real_property.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double]
fabm.variable_get_real_property.restype = ctypes.c_double

# Read/write/reset access to parameters.
fabm.get_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
fabm.get_real_parameter.restype = ctypes.c_double
fabm.get_integer_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
fabm.get_integer_parameter.restype = ctypes.c_int
fabm.get_logical_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int]
fabm.get_logical_parameter.restype = ctypes.c_int
fabm.get_string_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_char_p]
fabm.get_string_parameter.restype = None
fabm.reset_parameter.argtypes = [ctypes.c_void_p, ctypes.c_int]
fabm.reset_parameter.restype = None
fabm.set_real_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double]
fabm.set_real_parameter.restype = None
fabm.set_integer_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
fabm.set_integer_parameter.restype = None
fabm.set_logical_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
fabm.set_logical_parameter.restype = None
fabm.set_string_parameter.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
fabm.set_string_parameter.restype = None

# Read access to lists of variables (e.g., suitable coupling targets).
fabm.link_list_count.argtypes = [ctypes.c_void_p]
fabm.link_list_count.restype = ctypes.c_int
fabm.link_list_index.argtypes = [ctypes.c_void_p,ctypes.c_int]
fabm.link_list_index.restype = ctypes.c_void_p
fabm.link_list_finalize.argtypes = [ctypes.c_void_p]
fabm.link_list_finalize.restype = None

# Routines for sending pointers to state and dependency data.
fabm.link_interior_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
fabm.link_interior_state_data.restype = None
fabm.link_surface_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
fabm.link_surface_state_data.restype = None
fabm.link_bottom_state_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
fabm.link_bottom_state_data.restype = None
fabm.link_dependency_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
fabm.link_dependency_data.restype = None

# Read access to diagnostic data.
fabm.get_interior_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
fabm.get_interior_diagnostic_data.restype = None
fabm.get_horizontal_diagnostic_data.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
fabm.get_horizontal_diagnostic_data.restype = None

fabm.start.argtypes = [ctypes.c_void_p]
fabm.start.restype = None

# Routine for retrieving source-sink terms for the interior domain.
fabm.get_rates.argtypes = [ctypes.c_void_p, ctypes.c_double, numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=1, flags=CONTIGUOUS), ctypes.c_int, ctypes.c_int]
fabm.get_rates.restype = None
fabm.check_state.argtypes = [ctypes.c_void_p, ctypes.c_int]
fabm.check_state.restype = ctypes.c_int
fabm.integrate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=1, flags=CONTIGUOUS), numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=1, flags=CONTIGUOUS), numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, flags=CONTIGUOUS), ctypes.c_double, ctypes.c_int, ctypes.c_int]
fabm.integrate.restype = None

# Routine for getting git repository version information.
fabm.get_version.argtypes = (ctypes.c_int,ctypes.c_char_p)
fabm.get_version.restype = None

INTERIOR_STATE_VARIABLE        = 1
SURFACE_STATE_VARIABLE         = 2
BOTTOM_STATE_VARIABLE          = 3
INTERIOR_DIAGNOSTIC_VARIABLE   = 4
HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
CONSERVED_QUANTITY             = 6
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

def hasError():
   return fabm.get_error_state() != 0

def getError():
    if hasError():
        strmessage = ctypes.create_string_buffer(1024)
        fabm.get_error(1024, strmessage)
        return strmessage.value.decode('ascii')

def printTree(root,stringmapper,indent=''):
    """Print an indented tree of objects, encoded by dictionaries linking the names of children to
    their subtree, or to their object. Objects are finally printed as string obtained by
    calling the provided stringmapper method."""
    for name, item in root.items():
        if isinstance(item, dict):
            print('%s%s' % (indent, name))
            printTree(item, stringmapper, indent + '   ')
        else:
            print('%s%s = %s' % (indent, name, stringmapper(item)))

class Variable(object):
    def __init__(self, name=None, units=None, long_name=None, path=None, variable_pointer=None):
        self.variable_pointer = variable_pointer
        if variable_pointer:
           strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
           fabm.variable_get_metadata(variable_pointer, ATTRIBUTE_LENGTH, strname, strunits, strlong_name)
           name = strname.value.decode('ascii')
           units = strunits.value.decode('ascii')
           long_name = strlong_name.value.decode('ascii')

        self.name = name
        self.units = units
        if units is None:
           self.units_unicode = None
        else:
           self.units_unicode = createPrettyUnit(units)
        if long_name is None:
            long_name = name
        if path is None:
            path = name
        self.long_name = long_name
        self.path = path

    @property
    def long_path(self):
        if self.variable_pointer is None: return self.long_name
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        fabm.variable_get_long_path(self.variable_pointer,ATTRIBUTE_LENGTH,strlong_name)
        return strlong_name.value.decode('ascii')

    @property
    def output_name(self):
        if self.variable_pointer is None: return self.name
        stroutput_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        fabm.variable_get_output_name(self.variable_pointer,ATTRIBUTE_LENGTH,stroutput_name)
        return stroutput_name.value.decode('ascii')

    def getOptions(self):
        pass

    def getRealProperty(self, name, default=-1.0):
        return fabm.variable_get_real_property(self.variable_pointer, name.encode('ascii'), default)

class Dependency(Variable):
    def __init__(self, name, units=None, long_name=None):
        if long_name is None:
            long_name = name
        Variable.__init__(self, name, units, long_name.replace('_',' '))
        self.data = ctypes.c_double(0.)
        self.is_set = False

    def getValue(self):
        return self.data.value

    def setValue(self, value):
        self.is_set = True
        self.data.value = value

    value = property(getValue, setValue)

class StateVariable(Variable):
    def __init__(self,variable_pointer,statearray,index):
        Variable.__init__(self,variable_pointer=variable_pointer)
        self.index = index
        self.statearray = statearray

    def getValue(self):
        return float(self.statearray[self.index])

    def setValue(self,value):
        self.statearray[self.index] = value

    value = property(getValue, setValue)

    @property
    def background_value(self):
        return fabm.variable_get_background_value(self.variable_pointer)

    @property
    def output(self):
        return fabm.variable_get_output(self.variable_pointer) != 0

class DiagnosticVariable(Variable):
    def __init__(self, variable_pointer, index, horizontal):
        Variable.__init__(self, variable_pointer=variable_pointer)
        self.data = None

    def getValue(self):
        return None if self.data is None else self.data.value

    @property
    def output(self):
        return fabm.variable_get_output(self.variable_pointer) != 0

    value = property(getValue)

class Parameter(Variable):
    def __init__(self,name,index,units=None,long_name=None,type=None,model=None,has_default=False):
        Variable.__init__(self,name,units,long_name)
        self.type = type
        self.index = index+1
        self.model = model
        self.has_default = has_default

    def getValue(self,default=False):
        default = 1 if default else 0
        if self.type==1:
            return fabm.get_real_parameter(self.model.pmodel, self.index, default)
        elif self.type==2:
            return fabm.get_integer_parameter(self.model.pmodel, self.index, default)
        elif self.type==3:
            return fabm.get_logical_parameter(self.model.pmodel, self.index, default) != 0
        elif self.type==4:
            result = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
            fabm.get_string_parameter(self.model.pmodel, self.index, default, ATTRIBUTE_LENGTH, result)
            return result.value.decode('ascii')

    def setValue(self,value):
        settings = self.model.saveSettings()

        if self.type==1:
            fabm.set_real_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type==2:
            fabm.set_integer_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type==3:
            fabm.set_logical_parameter(self.model.pmodel, self.name.encode('ascii'), value)
        elif self.type==4:
            fabm.set_string_parameter(self.model.pmodel, self.name.encode('ascii'), value)

        # Update the model configuration (arrays with variables and parameters have changed)
        self.model.updateConfiguration(settings)

    def getDefault(self):
       if not self.has_default: return None
       return self.getValue(True)

    def reset(self):
        settings = self.model.saveSettings()
        fabm.reset_parameter(self.model.pmodel, self.index)
        self.model.updateConfiguration(settings)

    value = property(getValue, setValue)
    default = property(getDefault)

class Coupling(Variable):
    def __init__(self, pmodel, index):
        self.pmodel = pmodel
        self.master = ctypes.c_void_p()
        self.slave = ctypes.c_void_p()
        fabm.get_coupling(self.pmodel, index, ctypes.byref(self.slave), ctypes.byref(self.master))
        Variable.__init__(self, variable_pointer=self.slave)

    def getValue(self):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        fabm.variable_get_long_path(self.master,ATTRIBUTE_LENGTH,strlong_name)
        return strlong_name.value.decode('ascii')

    def setValue(self,value):
        print('New coupling specified: %s' % value)
        pass

    def getOptions(self):
        options = []
        list = fabm.variable_get_suitable_masters(self.pmodel, self.slave)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        for i in range(fabm.link_list_count(list)):
           variable = fabm.link_list_index(list,i+1)
           fabm.variable_get_long_path(variable,ATTRIBUTE_LENGTH,strlong_name)
           options.append(strlong_name.value.decode('ascii'))
        fabm.link_list_finalize(list)
        return options

    @property
    def long_path(self):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        fabm.variable_get_long_path(self.slave,ATTRIBUTE_LENGTH,strlong_name)
        return strlong_name.value.decode('ascii')

    value = property(getValue, setValue)

class SubModel(object):
    def __init__(self, pmodel, name):
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        iuser = ctypes.c_int()
        fabm.get_model_metadata(pmodel, name.encode('ascii'), ATTRIBUTE_LENGTH, strlong_name, iuser)
        self.long_name = strlong_name.value.decode('ascii')
        self.user_created = iuser.value!=0

class Model(object):
    def __init__(self,path='fabm.yaml'):
        self.lookup_tables = {}
        self.pmodel = fabm.create_model(path.encode('ascii'))
        assert not hasError(), 'An error occurred while parsing %s:\n%s' % (path, getError())
        self.updateConfiguration()

    def getSubModel(self,name):
        return SubModel(self.pmodel, name)

    def saveSettings(self):
        environment = dict([(dependency.name, dependency.value) for dependency in self.dependencies])
        state = dict([(variable.name,variable.value) for variable in self.state_variables])
        return environment,state

    def restoreSettings(self,data):
        environment,state = data
        for dependency in self.dependencies:
            if dependency.name in environment:
                dependency.value = environment[dependency.name]
        for variable in self.state_variables:
            if variable.name in state:
                variable.value = state[variable.name]

    def updateConfiguration(self,settings=None):
        # Get number of model variables per category
        nstate_interior = ctypes.c_int()
        nstate_surface = ctypes.c_int()
        nstate_bottom = ctypes.c_int()
        ndiag_interior = ctypes.c_int()
        ndiag_horizontal = ctypes.c_int()
        nconserved = ctypes.c_int()
        ndependencies = ctypes.c_int()
        nparameters = ctypes.c_int()
        ncouplings = ctypes.c_int()
        fabm.get_counts(self.pmodel, ctypes.byref(nstate_interior), ctypes.byref(nstate_surface), ctypes.byref(nstate_bottom),
                                 ctypes.byref(ndiag_interior), ctypes.byref(ndiag_horizontal),
                                 ctypes.byref(nconserved), ctypes.byref(ndependencies), ctypes.byref(nparameters), ctypes.byref(ncouplings))

        # Allocate memory for state variable values, and send ctypes.pointer to this memory to FABM.
        self.state = numpy.empty((nstate_interior.value + nstate_surface.value + nstate_bottom.value,), dtype=float)
        for i in range(nstate_interior.value):
            fabm.link_interior_state_data(self.pmodel, i + 1, self.state[i:].ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        for i in range(nstate_surface.value):
            fabm.link_surface_state_data(self.pmodel, i + 1, self.state[i+nstate_interior.value:].ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        for i in range(nstate_bottom.value):
            fabm.link_bottom_state_data(self.pmodel, i + 1, self.state[i+nstate_interior.value+nstate_surface.value:].ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

        # Retrieve variable metadata
        strname = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strunits = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strlong_name = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        strpath = ctypes.create_string_buffer(ATTRIBUTE_LENGTH)
        typecode = ctypes.c_int()
        has_default = ctypes.c_int()
        self.interior_state_variables = []
        self.surface_state_variables = []
        self.bottom_state_variables = []
        self.interior_diagnostic_variables = []
        self.horizontal_diagnostic_variables = []
        self.conserved_quantities = []
        self.parameters = []
        self.dependencies = []
        for i in range(nstate_interior.value):
            ptr = fabm.get_variable(self.pmodel, INTERIOR_STATE_VARIABLE, i + 1)
            self.interior_state_variables.append(StateVariable(ptr, self.state, i))
        for i in range(nstate_surface.value):
            ptr = fabm.get_variable(self.pmodel, SURFACE_STATE_VARIABLE, i + 1)
            self.surface_state_variables.append(StateVariable(ptr, self.state, nstate_interior.value + i))
        for i in range(nstate_bottom.value):
            ptr = fabm.get_variable(self.pmodel, BOTTOM_STATE_VARIABLE, i + 1)
            self.bottom_state_variables.append(StateVariable(ptr, self.state, nstate_interior.value + nstate_surface.value + i))
        for i in range(ndiag_interior.value):
            ptr = fabm.get_variable(self.pmodel, INTERIOR_DIAGNOSTIC_VARIABLE, i + 1)
            self.interior_diagnostic_variables.append(DiagnosticVariable(ptr, i, False))
        for i in range(ndiag_horizontal.value):
            ptr = fabm.get_variable(self.pmodel, HORIZONTAL_DIAGNOSTIC_VARIABLE, i + 1)
            self.horizontal_diagnostic_variables.append(DiagnosticVariable(ptr, i, True))
        for i in range(nconserved.value):
            fabm.get_variable_metadata(self.pmodel, CONSERVED_QUANTITY, i + 1, ATTRIBUTE_LENGTH, strname, strunits, strlong_name, strpath)
            self.conserved_quantities.append(Variable(strname.value.decode('ascii'), strunits.value.decode('ascii'), strlong_name.value.decode('ascii'), strpath.value.decode('ascii')))
        for i in range(nparameters.value):
            fabm.get_parameter_metadata(self.pmodel, i + 1, ATTRIBUTE_LENGTH, strname, strunits, strlong_name, ctypes.byref(typecode), ctypes.byref(has_default))
            self.parameters.append(Parameter(strname.value.decode('ascii'), i, type=typecode.value, units=strunits.value.decode('ascii'), long_name=strlong_name.value.decode('ascii'), model=self, has_default=has_default.value!=0))
        for i in range(ndependencies.value):
            fabm.get_dependency_metadata(self.pmodel, i + 1, ATTRIBUTE_LENGTH, strname, strunits)
            dependency = Dependency(strname.value.decode('ascii'), units=strunits.value.decode('ascii'))
            fabm.link_dependency_data(self.pmodel, i + 1, ctypes.byref(dependency.data))
            self.dependencies.append(dependency)

        self.couplings = [Coupling(self.pmodel, i + 1) for i in range(ncouplings.value)]

        # Arrays that combine variables from pelagic and boundary domains.
        self.state_variables = self.interior_state_variables + self.surface_state_variables + self.bottom_state_variables
        self.diagnostic_variables = self.interior_diagnostic_variables + self.horizontal_diagnostic_variables

        if settings is not None: self.restoreSettings(settings)

        # For backward compatibility
        self.bulk_state_variables = self.interior_state_variables
        self.bulk_diagnostic_variables = self.interior_diagnostic_variables

        self.itime = -1.

    def getRates(self, t=None, surface=True, bottom=True):
        """Returns the local rate of change in state variables,
        given the current state and environment.
        """
        if t is None:
            t = self.itime
        localrates = numpy.empty_like(self.state)
        fabm.get_rates(self.pmodel, t, localrates, surface, bottom)
        return localrates

    def checkState(self, repair=False):
        return fabm.check_state(self.pmodel, repair) != 0

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

    def findObject(self, name, objecttype, case_insensitive=False):
        tablename = str(objecttype)
        if case_insensitive:
            tablename += '_ci'
        table = self.lookup_tables.get(tablename, None)
        if table is None:
           data = getattr(self, objecttype)
           if case_insensitive:
              table = dict([(obj.name.lower(), obj) for obj in data])
           else:
              table = dict([(obj.name, obj) for obj in data])
           self.lookup_tables[tablename] = table
        if case_insensitive:
            name = name.lower()
        return table[name]

    def findParameter(self,name,case_insensitive=False):
        return self.findObject(name, 'parameters', case_insensitive)

    def findDependency(self,name,case_insensitive=False):
        return self.findObject(name, 'dependencies', case_insensitive)

    def findStateVariable(self,name,case_insensitive=False):
        return self.findObject(name, 'state_variables', case_insensitive)

    def findDiagnosticVariable(self,name,case_insensitive=False):
        return self.findObject(name, 'diagnostic_variables', case_insensitive)

    def findCoupling(self,name,case_insensitive=False):
        return self.findObject(name, 'couplings', case_insensitive)

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
       for dependency in self.dependencies:
          if not dependency.is_set:
             print('Value for dependency %s is not set.' % dependency.name)
             ready = False
       assert ready or not stop, 'Not all dependencies have been fulfilled.'
       fabm.start(self.pmodel)
       if hasError():
           return False
       for i, variable in enumerate(self.interior_diagnostic_variables):
            pdata = ctypes.POINTER(ctypes.c_double)()
            fabm.get_interior_diagnostic_data(self.pmodel, i + 1, ctypes.byref(pdata))
            variable.data = None if not pdata else pdata.contents
       for i, variable in enumerate(self.horizontal_diagnostic_variables):
            pdata = ctypes.POINTER(ctypes.c_double)()
            fabm.get_horizontal_diagnostic_data(self.pmodel, i + 1, ctypes.byref(pdata))
            variable.data = None if not pdata else pdata.contents
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
        self.model = model

    def integrate(self, y0, t, dt, surface=True, bottom=True):
        y = numpy.empty((t.size, self.model.state.size))
        fabm.integrate(self.model.pmodel, t.size, self.model.state.size, t, y0, y, dt, surface, bottom)
        return y

def get_version():
    version_length = 256
    strversion = ctypes.create_string_buffer(version_length)
    fabm.get_version(version_length,strversion)
    return strversion.value.decode('ascii')

