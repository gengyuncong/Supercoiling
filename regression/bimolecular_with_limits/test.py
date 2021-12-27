#!/usr/bin/env python
import os
import sys; sys.path.append("..")
from setupTest import addOrderParameter, runTest, LMTestFromNB

this_dir = os.path.dirname(os.path.realpath(__file__))

def _addOP(file_pth):
    addOrderParameter(file_pth, 0, 2, (0,3), (1/1000,1/100))

def executeWithTimeLimit(kwBase):
    kwTest = dict([
        ('dir_pth', 'data/1'),
        ('parameterDict', dict([
            ('writeInterval', 1e-2),
            ('maxTime', 1e-1),
        ])),
        ('r', 4),
    ])

    runTest(kwBase, kwTest)

def executeWithSpeciesUpperLimit(kwBase):
    kwTest = dict([
        ('dir_pth', 'data/2'),
        ('parameterDict', dict([
            ('writeInterval', 1e3),
            ('maxTime', 1e4),
            ('speciesUpperLimitList', '2:550,5:65')
        ])),
        ('r', 1000),
    ])

    runTest(kwBase, kwTest)

def executeWithSpeciesLowerLimit(kwBase):
    kwTest = dict([
        ('dir_pth', 'data/3'),
        ('parameterDict', dict([
            ('writeInterval', 1e3),
            ('maxTime', 1e4),
            ('speciesLowerLimitList', '2:450,5:35')
        ])),
        ('r', 1000),
    ])

    runTest(kwBase, kwTest)

def executeWithOrderParameterUpperLimit(kwBase):
    kwTest = dict([
        ('dir_pth', 'data/4'),
        ('parameterDict', dict([
            ('writeInterval', 1e3),
            ('maxTime', 1e4),
            ('orderParameterUpperLimitList', '0:1.2')
        ])),
        ('r', 48),
    ])

    runTest(kwBase, kwTest)

def executeWithOrderParameterLowerLimit(kwBase):
    kwTest = dict([
        ('dir_pth', 'data/5'),
        ('parameterDict', dict([
            ('writeInterval', 1e3),
            ('maxTime', 1e4),
            ('orderParameterLowerLimitList', '0:0.8')
        ])),
        ('r', 48),
    ])

    runTest(kwBase, kwTest)

def main():
    kwBase = dict([
        ('file_name', 'bimolecular_with_limits.lm'),
        ('sbml_pth', 'bimolecular_with_limits.sbml'),
        ('addOP', _addOP),
        ('gr', 0),
        ('c', 4),
        ('useavx', False),
    ])

    executeWithTimeLimit(kwBase)

    executeWithSpeciesUpperLimit(kwBase)
    executeWithSpeciesLowerLimit(kwBase)

    executeWithOrderParameterUpperLimit(kwBase)
    executeWithOrderParameterLowerLimit(kwBase)

def mainFromNB():
    lmtest = LMTestFromNB(os.path.join(this_dir, 'analysis.ipynb'))
    lmtest.runTests()

if __name__=='__main__':
    # main()
    mainFromNB()
