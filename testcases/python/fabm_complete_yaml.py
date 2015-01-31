#!/usr/bin/env python

import sys

if len(sys.argv)!=3:
   print 'This script takes two argument: the path to the input YAML file and the path to the output YAML file.'
   sys.exit(2)
infile = sys.argv[1]
outfile = sys.argv[2]

try:
   import pyfabm
except ImportError:
   print 'Unable to load pyfabm. Please build and install FABM with FABMHOST=python.'
   sys.exit(1)

import yaml # http://pyyaml.org/wiki/PyYAML

# ------------------------------------------
# Hook into PyYAML to make it preserve the order of dictionary elements.
try:
    import collections
    def dict_representer(dumper, data):                                                            
        return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())                                                                                         
    def dict_constructor(loader, node):                                                            
        return collections.OrderedDict(loader.construct_pairs(node))                               
    yaml.add_representer(collections.OrderedDict, dict_representer)                                
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
except ImportError:
    pass
# ------------------------------------------

# Create model object from YAML file.
model = pyfabm.Model(infile)

# Load the old configuration
with open(infile,'rU') as f:
   config = yaml.load(f)

def findMaximumDepth(d):
   n = 0
   for key,value in d.iteritems():
      if isinstance(value,dict):
         l = 2+findMaximumDepth(value)
      else:
         l = len(key)+2+len(str(value)) # key, folowed by ": ", followed by value
      n = max(n,l)
   return n

def reorderParameters(modelname,parameters):
   newparameters = collections.OrderedDict()
   parameters_lower = dict([(key.lower(),value) for key,value in parameters.iteritems()])
   modelname = modelname.lower()
   for parameter in model.parameters:
      if parameter.name.lower().startswith(modelname+'/'):
         name = parameter.name[len(modelname)+1:].lower()
         if name in parameters_lower: newparameters[name] = parameters_lower[name]
   assert len(newparameters)>=len(parameters)
   return newparameters

def processDict(f,d,path=[]):
   # If processing parameter list, reorder according to their registration by the model.
   if len(path)==3 and path[-1]=='parameters':
      d = reorderParameters(path[1],d)

   # If processing a model dictionary, reorder according to prescribed order.
   if len(path)==2 and path[0]=='instances':
      newd = collections.OrderedDict()
      order = ('long_name','model','parameters','initialization','coupling')
      for key in order:
         if key in d: newd[key] = d.pop(key)
      assert not d,'Model "%s" contains unknown keys %s' % (path[1],', '.join(d.keys()))
      d = newd

   nspace = len(path)*2
   for key,value in d.iteritems():
      f.write(' '*nspace)
      if isinstance(value,dict):
         f.write('%s:\n' % key)
         processDict(f,value,path=path+[key])
      else:
         metadata = None
         if len(path)==3:
            if path[-1]=='parameters':
               metadata = model.findParameter(path[1]+'/'+key,case_insensitive=True)
            elif path[-1]=='initialization':
               metadata = model.findStateVariable(path[1]+'_'+key)
         if metadata is not None:
            f.write('%s: %s' % (metadata.name[len(path[1])+1:],value))
            l = nspace+len(key)+2+len(str(value))
            f.write('%s# %s' % (' '*(icommentstart-l),metadata.long_name,))
            if metadata.units: f.write(' (%s)' % (metadata.units,))
         else:
            f.write('%s: %s' % (key,value))
         f.write('\n')

icommentstart = findMaximumDepth(config)+3
with open(outfile,'w') as f:   
   processDict(f,config)
      
