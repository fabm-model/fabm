#!/usr/bin/env python

# This script operates on outputs of run_unit_tests.py --performance

from __future__ import print_function

import os.path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('reference')
parser.add_argument('other', nargs='+')
args = parser.parse_args()

def loadResult(path):
    result = {}
    with open(path, 'rU') as f:
        f.readline()
        for l in f:
            items = l.rstrip('\n').split('\t')
            if items[-1] == 'NA':
                continue
            host, duration = items[0], float(items[-1])
            result[host] = duration
    return result

reference = loadResult(args.reference)
others = [loadResult(path) for path in args.other]
print('host\tchange in runtime (%)')
for host, reflength in sorted(reference.items(), key=lambda x: x[0].lower()):
    print(host, end='')
    for other in others:
        print('\t%.0f %%' % (100 * (other[host] / reflength - 1)), end='')
    print()