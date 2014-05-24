#!/usr/bin/env python

import sys,os.path
import urllib
import xml.etree.ElementTree

import yaml # http://pyyaml.org

variablespath = 'variables.yaml'
localpath = 'cf-standard-name-table.xml'
url = 'https://raw.githubusercontent.com/cf-convention/cf-documents/master/cf-standard-names/cf-standard-name-table-25.xml'
output_F90 = '../../include/standard_variables.h'
output_F90_assignments = '../../include/standard_variable_assignments.h'
output_wiki = 'standard_variables.wiki'

domain2type = {'bulk':'type_bulk_standard_variable',
               'horizontal':'type_horizontal_standard_variable',
               'global':'type_global_standard_variable',
               'conserved':'type_bulk_standard_variable'}

print 'Parsing %s...' % variablespath
stream = file(variablespath, 'rU')
selection = yaml.load(stream)
stream.close()
name2data = {}
for cls,items in selection.iteritems():
    for item in items: name2data[item['name']] = item

def strip_start(text, prefix):
    if not text.startswith(prefix): return text
    return text[len(prefix):]

def strip_end(text, suffix):
    if not text.endswith(suffix): return text
    return text[:len(text)-len(suffix)]

def CF2FABM(id):
    id = strip_end(id,'_in_sea_water')
    id = strip_start(id,'sea_water_')
    if id.endswith('_at_sea_floor'):
        id = 'bottom_'+strip_end(id,'_at_sea_floor')
    if id.startswith('sea_floor_'):
        id = 'bottom_'+strip_start(id,'sea_floor_')
    return id

if not os.path.isfile(localpath):
    print 'Downloading %s...' % url,
    xmlsrc = urllib.urlretrieve(url,localpath)
    print 'Done.'

print 'Parsing CF standard names XML...'
tree = xml.etree.ElementTree.parse(localpath)

print 'Linking CF variables to FABM variables...'
cfid2fabmid = {}
for element in tree.findall('entry'):
    cfid = element.get('id')
    cfunits = element.find('canonical_units').text
    fabmid = CF2FABM(cfid)
    if fabmid in name2data:
        if cfunits!=name2data[fabmid]['units']:
            print '   %s:\n      WARNING: mismatch between FABM units (%s) and CF units (%s)' % (fabmid,name2data[fabmid]['units'],cfunits)
        name2data[fabmid].setdefault('cf_names',[]).append(cfid)
        cfid2fabmid[cfid] = fabmid
for element in tree.findall('alias'):
    alias = element.get('id')
    fabmid = CF2FABM(cfid)
    basename = element.find('entry_id').text
    if fabmid in name2data:
        print 'ERROR: id %s is a synonym of %s' % (fabmid,basename)
        sys.exit(1)
    elif basename in cfid2fabmid:
        name2data[cfid2fabmid[basename]]['cf_names'].append(alias)

print 'The following FABM variables were not found in the CF convention:'
for name in sorted(name2data.keys()):
    if 'cf_names' not in name2data[name]: print '   %s' % name

print yaml.dump(selection,default_flow_style=False)

fout = open(output_F90,'w')
fout_assignments = open(output_F90_assignments,'w')
fwiki = open(output_wiki,'w')
for domain,items in selection.iteritems():
    fwiki.write('== %s variables ==\n\n{|\n|-\n! Variable\n! Units\n! Corresponding name in [http://cf-pcmdi.llnl.gov/documents/cf-standard-names/ CF convention]\n' % (domain[0].upper()+domain[1:]))
    fout.write('! %s variables\n' % (domain[0].upper()+domain[1:]))
    for i,item in enumerate(sorted(items,cmp=lambda x,y:cmp(x['name'],y['name']))):
        # Collect variable attribute]
        data = [('name',"'%s'" % item['name']),('units',"'%s'" % item['units'])]
        if 'cf_names' in item: data.append(('cf_names',"'%s'" % ','.join(item['cf_names'])))
        if item.get('aggregate_variable',domain=='conserved'): data.append(('aggregate_variable','.true.'))

        # Declare standard variable in Fortran.
        fout.write('type (%s) :: %s\n' % (domain2type[domain],item['name']))

        # Assign variable attributes in Fortran.
        for k,v in data: fout_assignments.write('standard_variables%%%s%%%s = %s\n' % (item['name'],k,v))
        fout_assignments.write('\n')

        # Create wiki entry for this variable.
        fwiki.write('|-\n| %s\n| %s\n' % (item['name'],item['units']))
        if 'cf_names' in item:
            fwiki.write('| %s\n' % item['cf_names'][0])
        else:
            fwiki.write('| \n')

    # Close entries for this variable category
    fwiki.write('|}\n\n')
    fout.write('\n')
