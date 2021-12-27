#!/usr/bin/env python
from argparse import ArgumentParser
import os

from setupTest import LMTestFromNBs

this_dir = os.path.dirname(os.path.realpath(__file__))

def getInitKwargs(kw):
    initKeys = (
        'lm_sbml_import',
        'lm_setp',
        'lmes',
        'r',
        'gr',
        'c',
        'useavx',
    )

    return {key:kw[key] for key in initKeys}

def getRunKwargs(kw):
    runKeys = (
        'dryrun',
        'verbose',
    )

    return {key:kw[key] for key in runKeys}

def main():
    parser = ArgumentParser()

    parser.add_argument('root_pth', default=this_dir, nargs='?')
    parser.add_argument('-out', '--outroot_pth')

    parser.add_argument('-ls', action='store_true')
    parser.add_argument('-tid', '--testIDs', nargs='*')

    parser.add_argument('-d', '--dryrun', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')

    parser.add_argument('-lm_sbml_import')
    parser.add_argument('-lm_setp')
    parser.add_argument('-lmes')

    parser.add_argument('-r')
    parser.add_argument('-gr')
    parser.add_argument('-c')
    parser.add_argument('-useavx', action='store_true')

    kwargs = vars(parser.parse_args())

    initKwargs = getInitKwargs(kwargs)
    runKwargs = getRunKwargs(kwargs)

    root_pth = kwargs['root_pth']
    outroot_pth = kwargs['outroot_pth']

    lmTests = LMTestFromNBs(root_pth=root_pth, outroot_pth=outroot_pth, **initKwargs)

    if kwargs['testIDs']:
        for testID in kwargs['testIDs']:
            lmTests.runTestByID(testID, **runKwargs)
    else:
        if kwargs['ls']:
            lmTests.ls()
        else:
            lmTests.runTests(**runKwargs)

if __name__=='__main__':
    main()
