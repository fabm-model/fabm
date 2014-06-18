#!/usr/bin/env python

import ctypes,ctypes.util,sys
from ctypes import CDLL, POINTER, c_int, c_double, c_char_p, byref, create_string_buffer, pointer
from numpy.ctypeslib import ndpointer

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'../../../compilers/vs2010/x64-Debug/fabm-python'))
sys.path.append(os.path.dirname(__file__))
import numpy

if sys.platform=='win32':
   dllpath = '../../../compilers/vs2010/x64-Debug/fabm-python/fabm-python.dll'
else:
   dllpath = 'fabm-python.so'
fabm = CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)),dllpath))
fabm.initialize.argtypes = [c_char_p]
fabm.get_variable_counts.argtypes = [POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int),POINTER(c_int)]
fabm.get_variable_metadata.argtypes = [c_int,c_int,c_int,c_char_p,c_char_p,c_char_p]
fabm.get_parameter_metadata.argtypes = [c_int,c_int,c_char_p,c_char_p,c_char_p,POINTER(c_int)]
fabm.get_dependency_metadata.argtypes = [c_int,c_int,c_char_p,c_char_p]
fabm.get_real_parameter.argtypes = [c_char_p,c_double]
fabm.get_real_parameter.restype = c_double
fabm.get_integer_parameter.argtypes = [c_char_p,c_int]
fabm.get_integer_parameter.restype = c_int
fabm.get_logical_parameter.argtypes = [c_char_p,c_int]
fabm.get_logical_parameter.restype = c_int
#fabm.get_string_parameter.argtypes = [c_char_p,c_char_p]
#fabm.get_string_parameter.restype = c_char_p
fabm.set_real_parameter.argtypes = [c_char_p,c_double]
fabm.set_integer_parameter.argtypes = [c_char_p,c_int]
fabm.set_logical_parameter.argtypes = [c_char_p,c_int]
fabm.link_bulk_state_data.argtypes = [c_int,POINTER(c_double)]
fabm.link_surface_state_data.argtypes = [c_int,POINTER(c_double)]
fabm.link_bottom_state_data.argtypes = [c_int,POINTER(c_double)]
fabm.link_dependency_data.argtypes = [c_int,POINTER(c_double)]
fabm.get_bulk_diagnostic_data.argtypes = [c_int,POINTER(POINTER(c_double))]
fabm.get_horizontal_diagnostic_data.argtypes = [c_int,POINTER(POINTER(c_double))]
fabm.get_rates.argtypes = [ndpointer(dtype=c_double, ndim=1, flags='CONTIGUOUS')]

BULK_STATE_VARIABLE            = 1
SURFACE_STATE_VARIABLE         = 2
BOTTOM_STATE_VARIABLE          = 3
BULK_DIAGNOSTIC_VARIABLE       = 4
HORIZONTAL_DIAGNOSTIC_VARIABLE = 5
CONSERVED_QUANTITY             = 6
ATTRIBUTE_LENGTH               = 256

def printTree(root,stringmapper,indent=''):
    """Print an indented tree of objects, encoded by dictionaries linking the names of children to
    their subtree, or to their object. Objects are finally printed as string obtained by
    calling the provided stringmapper method."""
    for name,item in root.iteritems():
        if isinstance(item,dict):
            print '%s%s' % (indent,name)
            printTree(item,stringmapper,indent+'   ')
        else:
            print '%s%s = %s' % (indent,name,stringmapper(item))

class Variable(object):
    def __init__(self,name,units=None,long_name=None):
        self.name = name
        self.units = units
        if long_name is None: long_name = name
        self.long_name = long_name

class Dependency(Variable):
    def __init__(self,name,index,units=None,long_name=None):
        Variable.__init__(self,name,units,long_name)
        self.data = c_double(0.)
        fabm.link_dependency_data(index+1,byref(self.data))

    def getValue(self):
        return self.data.value

    def setValue(self,value):
        self.data.value = value

    value = property(getValue, setValue)

class StateVariable(Variable):
    def __init__(self,statearray,name,index,units=None,long_name=None):
        Variable.__init__(self,name,units,long_name)
        self.index = index
        self.statearray = statearray

    def getValue(self):
        return float(self.statearray[self.index])

    def setValue(self,value):
        self.statearray[self.index] = value

    value = property(getValue, setValue)

class DiagnosticVariable(Variable):
    def __init__(self,name,index,horizontal,units=None,long_name=None):
        Variable.__init__(self,name,units,long_name)
        pdata = POINTER(c_double)()
        if horizontal:
            fabm.get_horizontal_diagnostic_data(index+1,byref(pdata))
        else:
            fabm.get_bulk_diagnostic_data(index+1,byref(pdata))
        self.data = pdata.contents

    def getValue(self):
        return self.data.value

    value = property(getValue)

class Parameter(object):
    def __init__(self,name,units=None,long_name=None,type=None,model=None):
        self.name = name
        self.units = units
        self.type = type
        self.model = model
        if long_name is None: long_name = name
        self.long_name = long_name

    def getValue(self):
        if self.type==1:
            return fabm.get_real_parameter(self.name,0.)
        elif self.type==2:
            return fabm.get_integer_parameter(self.name,0)
        elif self.type==3:
            return fabm.get_logical_parameter(self.name,0)!=0
        elif self.type==4:
            return fabm.get_string_parameter(self.name,'').rstrip()

    def setValue(self,value):
        settings = self.model.saveSettings()

        if self.type==1:
            fabm.set_real_parameter(self.name,value)
        elif self.type==2:
            fabm.set_integer_parameter(self.name,value)
        elif self.type==3:
            fabm.set_logical_parameter(self.name,value)
        elif self.type==4:
            fabm.set_string_parameter(self.name,value)

        # Update the model configuration (arrays with variables and parameters have changed)
        self.model.updateConfiguration(settings)

    def reset(self):
        fabm.reset_parameter(self.name)

    value = property(getValue, setValue)

class Model(object):
    def __init__(self,path='fabm.yaml'):
        fabm.initialize(path)
        self.updateConfiguration()

    def saveSettings(self):
        environment = dict([(dependency.name,dependency.value) for dependency in self.dependencies])
        state = dict([(variable.name,variable.value) for variable in self.state_variables])
        return environment,state

    def restoreSettings(self,data):
        environment,state = data
        for dependency in self.dependencies:
            if dependency.name in environment: dependency.value = environment[dependency.name]
        for variable in self.state_variables:
            if variable.name in state: variable.value = state[variable.name]

    def updateConfiguration(self,settings=None):
        # Get number of model variables per category
        nstate_bulk = c_int()
        nstate_surface = c_int()
        nstate_bottom = c_int()
        ndiag_bulk = c_int()
        ndiag_horizontal = c_int()
        nconserved = c_int()
        ndependencies = c_int()
        nparameters = c_int()
        fabm.get_variable_counts(byref(nstate_bulk),byref(nstate_surface),byref(nstate_bottom),
                                 byref(ndiag_bulk),byref(ndiag_horizontal),
                                 byref(nconserved),byref(ndependencies),byref(nparameters))

        # Allocate memory for state variable values, and send pointer to this memory to FABM.
        self.state = numpy.empty((nstate_bulk.value+nstate_surface.value+nstate_bottom.value,),dtype=float)
        for i in range(nstate_bulk.value):
            fabm.link_bulk_state_data(i+1,self.state[i:].ctypes.data_as(POINTER(c_double)))
        for i in range(nstate_surface.value):
            fabm.link_surface_state_data(i+1,self.state[i+nstate_bulk.value:].ctypes.data_as(POINTER(c_double)))
        for i in range(nstate_bottom.value):
            fabm.link_bottom_state_data(i+1,self.state[i+nstate_bulk.value+nstate_surface.value:].ctypes.data_as(POINTER(c_double)))

        # Retrieve variable metadata
        strname = create_string_buffer(ATTRIBUTE_LENGTH)
        strunits = create_string_buffer(ATTRIBUTE_LENGTH)
        strlong_name = create_string_buffer(ATTRIBUTE_LENGTH)
        typecode = c_int()
        self.bulk_state_variables = []
        self.surface_state_variables = []
        self.bottom_state_variables = []
        self.bulk_diagnostic_variables = []
        self.horizontal_diagnostic_variables = []
        self.conserved_quantities = []
        self.parameters = []
        self.dependencies = []
        for i in range(nstate_bulk.value):
            fabm.get_variable_metadata(BULK_STATE_VARIABLE,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.bulk_state_variables.append(StateVariable(self.state,strname.value,i,strunits.value,strlong_name.value))
        for i in range(nstate_surface.value):
            fabm.get_variable_metadata(SURFACE_STATE_VARIABLE,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.surface_state_variables.append(StateVariable(self.state,strname.value,nstate_bulk.value+i,strunits.value,strlong_name.value))
        for i in range(nstate_bottom.value):
            fabm.get_variable_metadata(BOTTOM_STATE_VARIABLE,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.bottom_state_variables.append(StateVariable(self.state,strname.value,nstate_bulk.value+nstate_surface.value+i,strunits.value,strlong_name.value))
        for i in range(ndiag_bulk.value):
            fabm.get_variable_metadata(BULK_DIAGNOSTIC_VARIABLE,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.bulk_diagnostic_variables.append(DiagnosticVariable(strname.value,i,False,strunits.value))
        for i in range(ndiag_horizontal.value):
            fabm.get_variable_metadata(HORIZONTAL_DIAGNOSTIC_VARIABLE,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.horizontal_diagnostic_variables.append(DiagnosticVariable(strname.value,i,True,strunits.value))
        for i in range(nconserved.value):
            fabm.get_variable_metadata(CONSERVED_QUANTITY,i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name)
            self.conserved_quantities.append(Variable(strname.value,strunits.value))
        for i in range(nparameters.value):
            fabm.get_parameter_metadata(i+1,ATTRIBUTE_LENGTH,strname,strunits,strlong_name,byref(typecode))
            self.parameters.append(Parameter(strname.value,type=typecode.value,units=strunits.value,long_name=strlong_name.value,model=self))
        for i in range(ndependencies.value):
            fabm.get_dependency_metadata(i+1,ATTRIBUTE_LENGTH,strname,strunits)
            self.dependencies.append(Dependency(strname.value,i,units=strunits.value))

        # Arrays that combine variables from pelagic and boudnary domains.
        self.state_variables = self.bulk_state_variables + self.surface_state_variables + self.bottom_state_variables
        self.diagnostic_variables = self.bulk_diagnostic_variables + self.horizontal_diagnostic_variables

        if settings is not None: self.restoreSettings(settings)

    def getRates(self):
        """Returns the local rate of change in state variables,
        given the current state and environment.
        """
        localrates = numpy.empty_like(self.state)
        consrates = numpy.empty((len(self.conserved_quantities),))
        fabm.get_rates(localrates)
        return localrates,consrates

    def findParameter(self,name):
        for parameter in self.parameters:
            if parameter.name==name: return parameter
        raise Exception('Parameter "%s" was not found.' % name)

    def findDependency(self,name):
        for dependency in self.dependencies:
            if dependency.name==name: return dependency
        raise Exception('Dependency "%s" was not found.' % name)

    def findStateVariable(self,name):
        for variable in self.state_variables:
            if variable.name==name: return variable
        raise Exception('State variable "%s" was not found.' % name)

    def getParameterTree(self):
        root = {}
        for parameter in self.parameters:
            pathcomps = parameter.name.split('/')
            parent = root
            for component in pathcomps[:-1]:
                parent = root.setdefault(component,{})
            parent[pathcomps[-1]] = parameter
        return root

    def printInformation(self):
        """Show information about the model."""
        print 'FABM model contains the following:'
        print ' %i state variables:' % len(self.state_variables)
        for variable in self.state_variables: print '    %s = %s %s' % (variable.name,variable.value,variable.units)
        print ' %i diagnostic variables:' % len(self.diagnostic_variables)
        for variable in self.diagnostic_variables: print '    %s = %s %s' % (variable.name,variable.value,variable.units)
        print ' %i external variables:' % len(self.dependencies)
        for variable in self.dependencies: print '    %s = %s %s' % (variable.name,variable.value,variable.units)
        print ' %i parameters:' % len(self.parameters)
        printTree(self.getParameterTree(),lambda x:'%s %s' % (x.value,x.units),'    ')
