#!/usr/bin/env python

import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'../../../compilers/vs2010/x64-Debug/fabm-python-c++'))
sys.path.append(os.path.join(os.path.dirname(__file__),'fabm'))
import numpy
import fabm

# Shortcut to FABM module living on the Fortran side
fabm = fabm.fabm_python

def getStringArray(d):
    """Maps a F2Py-exposed Fortran string array to a Python list of strings."""
    shape = d.shape
    return [''.join(d.T.reshape(shape)[i,:]).rstrip() for i in range(d.shape[0])]

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
    def __init__(self,name,units=None):
        self.name = name
        self.units = units

class Dependency(Variable):
    def __init__(self,name,index,units=None):
        Variable.__init__(self,name,units)
        self.index = index

    def getValue(self):
        return fabm.environment[self.index]

    def setValue(self,value):
        fabm.environment[self.index] = value

    value = property(getValue, setValue)

class StateVariable(Variable):
    def __init__(self,name,index,units=None):
        Variable.__init__(self,name,units)
        self.index = index

    def getValue(self):
        return fabm.state[self.index]

    def setValue(self,value):
        fabm.state[self.index] = value

    value = property(getValue, setValue)

class Parameter(object):
    def __init__(self,name,units=None,type=None,model=None):
        self.name = name
        self.units = units
        self.type = type
        self.model = model

    def getValue(self):
        if self.type==1:
            return fabm.get_real_parameter(self.name,0.)
        elif self.type==2:
            return fabm.get_integer_parameter(self.name,0)
        elif self.type==3:
            return bool(fabm.get_logical_parameter(self.name,False))
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
    def __init__(self):
        fabm.initialize()
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
            if variable.name in state: variable.value = state[dependency.name]

    def updateConfiguration(self,settings=None):
        # Map string arrays in Fortran to variable objects with members decribing metadata.
        self.state_variables = tuple([StateVariable(name,i,units) for i,(name,units) in enumerate(zip(getStringArray(fabm.state_names),getStringArray(fabm.state_units)))])
        self.dependencies = tuple([Dependency(name,i,units) for i,(name,units) in enumerate(zip(getStringArray(fabm.environment_names),getStringArray(fabm.environment_units)))])
        self.parameters = tuple([Parameter(name,type=type,model=self,units=units) for name,units,type in zip(getStringArray(fabm.parameter_names),getStringArray(fabm.parameter_units),fabm.parameter_types)])
        self.conserved_quantities = tuple([Variable(name,units) for name,units in zip(getStringArray(fabm.conserved_quantity_names),getStringArray(fabm.conserved_quantity_units))])

        # Make module-level arrays in Fortran accessible as member variables of the model.
        self.dependency_values = fabm.environment
        self.state = fabm.state

        if settings is not None: self.restoreSettings(settings)

    def getRates(self):
        """Returns the local rate of change in state variables,
        given the current state and environment.
        """
        localrates = numpy.empty_like(self.state)
        consrates = numpy.empty((len(self.conserved_quantities),))
        fabm.get_rates(localrates,consrates)
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
        print ' %i external variables:' % len(self.dependencies)
        for variable in self.dependencies: print '    %s = %s %s' % (variable.name,variable.value,variable.units)
        print ' %i parameters:' % len(self.parameters)
        printTree(self.getParameterTree(),lambda x:'%s %s' % (x.value,x.units),'    ')
