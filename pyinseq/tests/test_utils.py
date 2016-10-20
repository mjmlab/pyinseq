#!/usr/bin/env python3
from distutils import dir_util
import os
from pkg_resources import Requirement, resource_filename, ResolutionError
import shutil
import sys
import traceback

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pytest
'''
These script running functions were taken from the khmer project:
https://github.com/dib-lab/khmer/blob/master/tests/khmer_tst_utils.py
'''

@pytest.fixture(scope='session')
def datadir(tmpdir_factory, request):
    '''
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    '''
    tmpdir = tmpdir_factory.mktemp('data')
    # this is the data/ dir in pyinseq/tests
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    dir_util.copy_tree(data_dir, str(tmpdir))

    def getter(filename='', as_str=True):
        filepath = tmpdir.join(filename)
        if as_str:
            return str(filepath)
        return filepath

    return getter


def scriptpath(scriptname='pyinseq'):
    "Return the path to the scripts, in both dev and install situations."

    path = os.path.join(os.path.dirname(__file__), "../../scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../../EGG-INFO/scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ['PATH'].split(':'):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname):
    """
    Find & run a script with exec (i.e. not via os.system or subprocess).
    """

    import pkg_resources
    ns = {"__name__": "__main__"}
    ns['sys'] = globals()['sys']

    try:
        pkg_resources.get_distribution("pyinseq").run_script(scriptname, ns)
        return 0
    except pkg_resources.ResolutionError as err:
        path = scriptpath()

        scriptfile = os.path.join(path, scriptname)
        if os.path.isfile(scriptfile):
            if os.path.isfile(scriptfile):
                exec(compile(open(scriptfile).read(), scriptfile, 'exec'), ns)
                return 0

    return -1


def runscript(scriptname, args, directory=None,
              fail_ok=False, sandbox=False):
    """Run a Python script using exec().
    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.
    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    sysargs = [scriptname]
    sysargs.extend(args)
    cwd = os.getcwd()

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        if directory:
            os.chdir(directory)
        else:
            directory = cwd

        try:
            print('running:', scriptname, 'in:', directory, file=oldout)
            print('arguments', sysargs, file=oldout)

            status = _runscript(scriptname)
        except SystemExit as e:
            status = e.code
        except:
            traceback.print_exc(file=sys.stderr)
            status = -1
    finally:
        sys.argv = oldargs
        out, err = sys.stdout.getvalue(), sys.stderr.getvalue()
        sys.stdout, sys.stderr = oldout, olderr

        os.chdir(cwd)

    if status != 0 and not fail_ok:
        print('Script Failed:', scriptname,
              'Status:', status,
              'Output:', out,
              'Error:', err,
              sep='\n')
        assert False

    return status, out, err
