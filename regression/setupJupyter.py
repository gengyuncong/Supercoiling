from argparse import ArgumentParser
from IPython.core.magic import register_cell_magic
import matplotlib
import matplotlib.pyplot as plt
import shlex

from setupTest import addOrderParameter, executeBash

__all__ = ['addOrderParameter', 'setupJupyter']

try:
    # shadows/replaces the IPython built-in bash cell magic
    @register_cell_magic
    def bash(line, cell):
        """An improved version of the %%bash IPython cell magic.
        Allows for the use of python variables in bash scripts, using the `$${var}` syntax.
        Also support real time display of stdout.
        """
        # if we have received any arguments on the magic line, parse them
        parser = ArgumentParser()

        parser.add_argument('varDict', nargs='?', help='(optional) dictionary or container of dictionaries, to be used to in place of the default globals()/locals() combo when substituting variables')
        parser.add_argument('--echo', action='store_true', help='do a dry-run, ie print commands instead of running them')
        parser.add_argument('--nosub', action='store_true', help='skip variable substitution on the bash script')
        parser.add_argument('--console_output', dest='consoleOutput', action='store_true', help='pipe stdout/err to the ipython console output instead of the notebook cell output')

        lineTokens = shlex.split(line)
        lineDict = vars(parser.parse_args(lineTokens))

        # figure out the dictionary we're going to use when substituting variables in the bash script
        if lineDict['nosub']:
            # skip variable substitution on the bash script
            varDict = {}
        elif lineDict['varDict'] is not None:
            # assume we were passed a dict or a container of dicts
            varDict = eval(lineDict['varDict'])
        else:
            # no varDict was passed to our magic, so search globals() and locals() for variables to substitute
            varDict = {}
            varDict.update(globals())
            varDict.update(locals())

        # execute the bash script
        optDict = {key:lineDict[key] for key in ['echo', 'consoleOutput']}
        if isinstance(varDict, dict):
            # we were passed a single dict on the magic line, so execute once
            kwargs = dict(optDict)
            kwargs.update(varDict)
            executeBash(cell, **kwargs)
        else:
            # we were passed a container of dicts, so execute once for each dict
            for subVarDict in varDict:
                kwargs = dict(optDict)
                kwargs.update(subVarDict)
                executeBash(cell, **kwargs)
except NameError:
    # NameError: Decorator can only run in context where `get_ipython` exists
    pass

def setupJupyter(figBkg='light', pylab='inline', rerun=False):
    """Initializes the Jupyter environment
    """
    if rerun or not ('_guard_setupJupyter' in globals() and globals()['_guard_setupJupyter']):
        # run some magics to get inline plotting and automatic reloading of edited module code
        mStrs = [
            'pylab %s' % pylab,
            "config InlineBackend.figure_format = 'retina'",
            'load_ext autoreload',
            'autoreload 2',
            ]
        for mStr in mStrs:
            get_ipython().magic(mStr)

        def _setupJupyter():
            light_background = {key:plt.rcParamsDefault[key] for key in plt.style.library['dark_background'].keys()}

            # set the plot dpi to increase the size of the inline plots
            matplotlib.rcParams['figure.dpi'] = 100

            if figBkg == 'light':
                plt.style.use(light_background)
            elif figBkg == 'dark':
                plt.style.use('dark_background')

            # remove the hook used to run this setup func
            try:
                get_ipython().events.unregister('post_run_cell', _setupJupyter)
            except ValueError:
                pass

            # we reached the end of setup without an exception, set the guard to prevent double setup
            globals()['_guard_setupJupyter'] = True

        # some magics run at the end of a cell. Run the rest of our setup after those magics are finished.
        get_ipython().events.register('post_run_cell', _setupJupyter)
