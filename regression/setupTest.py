from collections import OrderedDict
import contextlib
import h5py
import multiprocessing
import nbformat
import numpy as np
import os
import re
from six import print_, StringIO, string_types
from string import Template
from subprocess import PIPE, Popen
import sys

__all__ = ['addOrderParameter', 'executeBash', 'LMTest', 'runTest', 'LMTestFromNB', 'LMTestFromNBs']

this_dir = os.path.dirname(os.path.realpath(__file__))
build_dir = os.path.realpath(os.path.join(this_dir, '..', 'build'))

def addOrderParameter(filename, id, type, speciesIDs, speciesCoefficients):
    print_("Adding order parameters to %s" % filename)
    fp = h5py.File(filename, "r+")
    if "OrderParameters" not in fp.keys():
        opg = fp.create_group("/OrderParameters/%07d" % id)
        opg.attrs.create("ID", id, dtype=np.uint32)
        opg.attrs.create("Type", type, dtype=np.uint32)
        d1=opg.create_dataset("SpeciesIDs", (len(speciesIDs),), dtype=np.uint32)
        d1[:]=np.array(speciesIDs)
        d2=opg.create_dataset("SpeciesCoefficients", (len(speciesIDs),), dtype=np.double)
        d2[:]=np.array(speciesCoefficients)
    else:
        print_("ERROR: file already has order parameters")
    fp.close()

def banner(s, padlen=4):
    _pad = '#'*padlen

    lshort = '%s ' % _pad
    rshort = ' %s' % _pad
    long = '#'*(len(s) + len(lshort) + len(rshort))

    return '\n'.join([
        long,
        '%s%s%s' % (lshort, s, rshort),
        long,
    ])

@contextlib.contextmanager
def remember_cwd():
    curdir= os.getcwd()
    try: yield
    finally: os.chdir(curdir)

class _BashTemplate(Template):
    delimiter = '$$'

def _setupScript(script, dryrun, verbose, kind='script', templTipe=None, **varDict):
    if not isinstance(script, string_types):
        script = '\n'.join(script)

    if templTipe:
        scriptTmpl = templTipe(script)
        script = scriptTmpl.substitute(varDict)

    if dryrun or verbose:
        print_('####%s:' % kind)
        print_(script)
        if dryrun: return None
        print_()

    return script

def executeBash(script, consoleOutput=False, dryrun=False, verbose=False, **varDict):
    script = _setupScript(script, dryrun=dryrun, verbose=verbose, kind='Bash', templTipe=_BashTemplate, **varDict)
    if script is None: return

    if consoleOutput:
        p = Popen(script, executable='/bin/bash', shell=True)
        p.wait()
    else:
        if sys.version_info > (3,):
            with Popen(script, executable='/bin/bash', shell=True, stdout=PIPE, stderr=PIPE, bufsize=1) as p:
                for line in p.stdout:
                    print_(line.decode('latin-1'), end='', flush=True)
                print_(p.communicate()[0].decode('latin-1'), end='', flush=True)
        else:
            p = Popen(script, executable='/bin/bash', shell=True, stdout=PIPE, stderr=PIPE, bufsize=1)

            for line in iter(p.stdout.readline, ''):
                print_(line, end='', flush=True)
            p.stdout.close()
            p.wait()

def executePython(script, dryrun=False, verbose=False):
    script = _setupScript(script, dryrun=dryrun, verbose=verbose, kind='Python')
    if script is None: return

    exec(script)

class LMTest(object):
    def __init__(self, dir_pth, file_name, sbml_pth, lm_sbml_import=None, lm_setp=None, lmes=None, addOP=None, parameterDict=None, r=10, gr=0, c=4, useavx=False):
        self.dir_pth = dir_pth
        self.file_name = file_name
        self.file_pth = os.path.join(dir_pth, file_name)

        self.lm_sbml_import = os.path.join(build_dir, 'utils', 'c', 'lm_sbml_import') if lm_sbml_import is None else lm_sbml_import
        self.lm_setp = os.path.join(build_dir, 'utils', 'c', 'lm_setp') if lm_setp is None else lm_setp
        self.lmes = os.path.join(build_dir, 'lmes') if lmes is None else lmes

        self.sbml_pth = sbml_pth
        self.addOP = addOP
        self.parameterDict = {} if parameterDict is None else parameterDict

        self.r = r
        self.gr = gr
        self.c = c
        self.useavx = useavx

        self.bashDict = None
        self.initBashDict()

    def initBashDict(self):
        self.bashDict = OrderedDict([
            ('dir_pth',         self.dir_pth),
            ('file_pth',        self.file_pth),
            ('sbml_pth',        self.sbml_pth),
            ('lm_sbml_import', self.lm_sbml_import),
            ('lm_setp',        self.lm_setp),
            ('lmes',          self.lmes)
        ])

    def setupTest(self, dryrun=False, verbose=False):
        parameters = ' '.join(['%s=%s' % (key,val) for key,val in self.parameterDict.items()])

        bashScript = [
            "mkdir -p $${dir_pth}",
            "rm -f $${file_pth}",
            "$${lm_sbml_import} $${file_pth} $${sbml_pth}",
            "$${lm_setp} $${file_pth} $${parameters}",
        ]

        bashDict = OrderedDict([
            ('parameters', parameters),
        ])
        bashDict.update(self.bashDict)

        executeBash(bashScript, dryrun=dryrun, verbose=verbose, **bashDict)

        if self.addOP:
            print_()
            self.addOP(self.file_pth)

        print_()

    def executeTest(self, dryrun=False, verbose=False):
        solver = 'lm::avx::GillespieDSolverAVX' if self.useavx else 'lm::cme::GillespieDSolver'

        lmArgDict = OrderedDict([
            ('r', '1-%s' % self.r),
            ('gr', self.gr),
            ('c',  self.c),
            ('sl', solver)
        ])
        lmArgs = ' '.join('-%s %s' % (key,val) for key,val in lmArgDict.items())

        bashScript = "$${lmes} -f $${file_pth} $${lmArgs}"

        bashDict = OrderedDict([
            ('lmArgs', lmArgs),
        ])
        bashDict.update(self.bashDict)

        executeBash(bashScript, dryrun=dryrun, verbose=verbose, **bashDict)

        print_()

    def runTest(self, dryrun=False, verbose=False):
        print_("Running tests using the following executables:")
        for pthAttr in ('lm_sbml_import', 'lm_setp', 'lmes'):
            print_('%s => %s' % (pthAttr, getattr(self, pthAttr)))

        self.setupTest(dryrun=dryrun, verbose=verbose)
        self.executeTest(dryrun=dryrun, verbose=verbose)

def runTest(*kws, **kwTestRun):
    kwTest = {}
    for kw in kws:
        kwTest.update(kw)

    lmtest = LMTest(**kwTest)
    lmtest.runTest(**kwTestRun)

bashIgnoreRe = re.compile(r'^\s*(?:(?:$)|(?:#)|(?:%%bash))')
redirectRe = re.compile(r'(>.*)')

dirnameRe = re.compile(r'(dirname\s*=\s*)([\'\"]?)(\S+)')

lmSbmlImportRe = re.compile(r'(\S*lm_sbml_import)')
lmSetpRe = re.compile(r'(\S*lm_setp)')
lmesRe = re.compile(r'(\S*lmes)')

rRe = re.compile(r'(-r\s*)(\S+)')
grRe = re.compile(r'(-gr\s*)(\S+)')
cRe = re.compile(r'(-c\s*)(\S+)')
slRe = re.compile(r'(-sl\s*)(\S+)')

class LMTestFromNB(object):
    scriptRunners = dict([
        ('bash', executeBash),
        ('python', executePython),
    ])

    @property
    def dir_pth(self):
        return os.path.dirname(self.nb_pth)

    @property
    def useavx(self):
        return 'AVX' in self.sl

    def __init__(self, nb_pth, redirect=None, lm_sbml_import=None, lm_setp=None, lmes=None, dirname=None, r=None, gr=None, c=None, useavx=None):
        redirect = False if redirect is None else redirect
        useavx = False if useavx is None else useavx

        self.nb_pth = None
        self.testName = None
        self.initPth(nb_pth)

        self.redirect = redirect

        self.lm_sbml_import = os.path.join(build_dir, 'utils', 'c', 'lm_sbml_import') if lm_sbml_import is None else lm_sbml_import
        self.lm_setp = os.path.join(build_dir, 'utils', 'c', 'lm_setp') if lm_setp is None else lm_setp
        self.lmes = os.path.join(build_dir, 'lmes') if lmes is None else lmes

        self.dirname = dirname
        self.r = None if r is None else ('1-%s' % r)
        self.gr = gr
        self.c = multiprocessing.cpu_count() if c is None else c
        self.sl = 'lm::avx::GillespieDSolverAVX' if useavx else 'lm::cme::GillespieDSolver'

        self.testDict = None
        self.initTests()

    def initPth(self, nb_pth):
        if ':' in nb_pth:
            self.nb_pth, self.testName = nb_pth.split(':')
        else:
            self.nb_pth = nb_pth

    def initTests(self):
        self.testDict = OrderedDict()

        workednb = nbformat.read(str(self.nb_pth), as_version=nbformat.NO_CONVERT)
        for cell in workednb.cells:
            if 'test_name' in cell.metadata:
                testName = cell.metadata['test_name']
                self.testDict[testName] = testCells = self.testDict.get(testName, [])

                if '%%bash' in cell.source:
                    testCells.append((
                        'bash',
                        ''.join(self.cleanBashLine(line) for line in StringIO(cell.source)),
                    ))
                else:
                    testCells.append((
                        'python',
                        ''.join(self.cleanPythonLine(line) for line in StringIO(cell.source)),
                    ))

    def cleanBashLine(self, line):
        if bashIgnoreRe.search(line):
            return ''

        if self.dirname is not None:
            search = dirnameRe.search(line)
            if search:
                dirname = os.path.join(str(self.dirname), search.group(3).strip("\"\'"))
                line = dirnameRe.sub(r'\g<1>\g<2>' + dirname + r'\g<2>', line, count=1)

        for subRe,subRepl in zip((lmSbmlImportRe, lmSetpRe, lmesRe),
                                 (self.lm_sbml_import, self.lm_setp, self.lmes)):
            line = subRe.sub(str(subRepl), line, count=1)

        for subRe,subRepl in zip((rRe, grRe, cRe, slRe),
                                 (self.r, self.gr, self.c, self.sl)):
            if subRepl is not None:
                line = subRe.sub(r'\g<1>' + str(subRepl), line, count=1)

        if not self.redirect:
            line = redirectRe.sub('', line)
            
        return line

    def cleanPythonLine(self, line):
        if self.dirname is not None:
            search = dirnameRe.search(line)
            if search:
                dirname = os.path.join(str(self.dirname), search.group(3).strip("\"\'"))
                line = dirnameRe.sub(r'\g<1>\g<2>' + dirname + r'\g<2>', line, count=1)

        return line

    def getTestNameVisible(self, testName):
        if testName is None:
            return None
        else:
            return testName + '_AVX' if self.useavx else testName

    def runTest(self, testName, dryrun=False, verbose=False):
        if not dryrun and self.dirname is not None:
            try:
                os.makedirs(self.dirname)
            except OSError:
                pass

        # restore the current working directory at the end of the test no matter what
        with remember_cwd():
            os.chdir(self.dir_pth)

            print_(banner('%s %s' % (self.dir_pth, self.getTestNameVisible(testName))))
            print_()

            print_("Using the following executables:")
            for pthAttr in ('lm_sbml_import', 'lm_setp', 'lmes'):
                print_('%s => %s' % (pthAttr, getattr(self, pthAttr)))
            print_()

            for testCellType,testCell in self.testDict[testName]:
                self.scriptRunners[testCellType](testCell, dryrun=dryrun, verbose=verbose)
                print_()

    def runTests(self, dryrun=False, verbose=False):
        if self.testName is not None:
            self.runTest(testName=self.testName, dryrun=dryrun, verbose=verbose)
        else:
            for testName in self.testDict.keys():
                self.runTest(testName=testName, dryrun=dryrun, verbose=verbose)

nbNameRe = re.compile('(?:analysis)|(?:test)\.ipynb')

class LMTestFromNBs(object):
    def __init__(self, root_pth, outroot_pth=None, redirect=None, lm_sbml_import=None, lm_setp=None, lmes=None, r=None, gr=None, c=None, useavx=None):
        self.root_pth = root_pth
        self.outroot_pth = outroot_pth
        self.redirect = redirect

        self.lm_sbml_import = lm_sbml_import
        self.lm_setp = lm_setp
        self.lmes = lmes

        self.r = r
        self.gr = gr
        self.c = c
        self.useavx = useavx

        self.testDict = None
        self.getLMTests()

    def getLMTests(self):
        self.testDict = OrderedDict()

        for dir_pth, subdirs, fnames in os.walk(self.root_pth):
            if '.ipynb_checkpoints' in dir_pth:
                # skip the auto-saved backup notebooks
                continue

            for fname in (f for f in fnames if nbNameRe.search(f)):
                nb_pth = os.path.join(dir_pth, fname)

                dirname = os.path.join(self.outroot_pth, os.path.relpath(dir_pth, self.root_pth)) if self.outroot_pth is not None else None

                # lmTest will automatically search for tests in the notebook at nb_pth upon initialization
                lmTest = LMTestFromNB(
                    nb_pth=nb_pth,
                    dirname=dirname,
                    redirect=self.redirect,
                    lm_sbml_import=self.lm_sbml_import,
                    lm_setp=self.lm_setp,
                    lmes=self.lmes,
                    r=self.r,
                    gr=self.gr,
                    c=self.c,
                    useavx=self.useavx,
                )

                if lmTest.testDict:
                    # if lmTest discovered any tests in the notebook, keep it
                    self.testDict[nb_pth] = lmTest

    def _ls(self):
        for lmTest in self.testDict.values():
            nb_pth = lmTest.nb_pth
            for testName in lmTest.testDict:
                yield nb_pth,testName

    def ls(self):
        out = '\n'.join([('%s:%s' % (nb_pth, testName)) for nb_pth,testName in self._ls()])
        print_(out, end='')

    def runTest(self, nb_pth, testName, dryrun=False, verbose=False):
        self.testDict[nb_pth].runTest(testName=testName, dryrun=dryrun, verbose=verbose)

    def runTestByID(self, testID, dryrun=False, verbose=False):
        nb_pth,testName = testID.split(':')
        self.runTest(nb_pth=nb_pth, testName=testName, dryrun=dryrun, verbose=verbose)

    def runTests(self, dryrun=False, verbose=False):
        for nb_pth,testName in self._ls():
            self.runTest(nb_pth=nb_pth, testName=testName, dryrun=dryrun, verbose=verbose)
