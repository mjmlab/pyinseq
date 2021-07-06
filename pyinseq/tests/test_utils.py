#!/usr/bin/env python3
import os
import sys
import filecmp
import traceback
from distutils import dir_util
from pathlib import Path

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pytest


class cd:
    """Context manager to change to the specified directory then back."""

    def __init__(self, new_path):
        self.new_path = os.path.expanduser(new_path)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def run_pytest():
    """Runs pytest directly from source code. CWD should be in pyinseq source dir"""
    rootdir = Path(__file__).parents[2]
    with cd(rootdir):
        pytest.main([])


def get_dump():
    """Helper for creating and returning a dump Path"""
    dump = Path("pyinseq/tests/dump")
    if not dump.exists():
        dump.mkdir()
    return dump


@pytest.fixture()
def load_settings():
    def fetch(command):
        config_file = "../data/input/test-pyinseq-config.yml"
        dump = get_dump()
        with cd(dump):
            from pyinseq.settings import SETTINGS_CONSTRUCTORS

            settings_constructor = SETTINGS_CONSTRUCTORS[command]
            return settings_constructor(config_file)

    return fetch


"""
These script running functions were taken from the khmer project:
https://github.com/dib-lab/khmer/blob/master/tests/khmer_tst_utils.py
"""


@pytest.fixture(scope="session")
def datadir(tmpdir_factory, request):
    """
    Fixture responsible for locating the test data directory and copying it
    into a temporary directory.
    """
    tmpdir = tmpdir_factory.mktemp("data")
    # this is the data/ dir in pyinseq/tests
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    dir_util.copy_tree(data_dir, str(tmpdir))

    # Create dump in temporary directory
    tmpdir.join("dump").mkdir()

    def getter(filename="", as_str=True):
        filepath = tmpdir.join(filename)
        if as_str:
            return str(filepath)
        return filepath

    return getter


def scriptpath(scriptname="pyinseq"):
    "Return the path to the scripts, in both dev and install situations."

    path = os.path.join(os.path.dirname(__file__), "../../scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    path = os.path.join(os.path.dirname(__file__), "../../../EGG-INFO/scripts")
    if os.path.exists(os.path.join(path, scriptname)):
        return path

    for path in os.environ["PATH"].split(":"):
        if os.path.exists(os.path.join(path, scriptname)):
            return path


def _runscript(scriptname):
    """
    Find & run a script with exec (i.e. not via os.system or subprocess).
    """

    import pkg_resources

    ns = {"__name__": "__main__"}
    ns["sys"] = globals()["sys"]

    try:
        pkg_resources.get_distribution("pyinseq").run_script(scriptname, ns)
        return 0
    except pkg_resources.ResolutionError as err:
        path = scriptpath()

        scriptfile = os.path.join(path, scriptname)
        if os.path.isfile(scriptfile):
            if os.path.isfile(scriptfile):
                exec(compile(open(scriptfile).read(), scriptfile, "exec"), ns)
                return 0

    return -1


def runscript(
    scriptname: str, args: list, directory=None, fail_ok=False, sandbox=False
) -> [str, str, str]:
    """Run a Python script using exec().
    Run the given Python script, with the given args, in the given directory,
    using 'exec'.  Mimic proper shell functionality with argv, and capture
    stdout and stderr.
    When using :attr:`fail_ok`=False in tests, specify the expected error.
    """
    sysargs = [scriptname]
    sysargs.extend(args)
    cwd = os.getcwd()
    print(os.getcwd())

    try:
        status = -1
        oldargs = sys.argv
        sys.argv = sysargs

        oldout, olderr = sys.stdout, sys.stderr
        sys.stdout = StringIO()
        sys.stdout.name = "StringIO"
        sys.stderr = StringIO()

        # Use with to change directories
        if not directory:
            directory = cwd

        with cd(directory):
            try:
                print("running:", scriptname, "in:", directory, file=oldout)
                print("arguments", sysargs, file=oldout)
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
        print(out)

    if status != 0 and not fail_ok:
        print(
            "Script Failed:",
            scriptname,
            "Status:",
            status,
            "Output:",
            out,
            "Error:",
            err,
            sep="\n",
        )
        assert False
    return status, out, err


def compare_directories(expected_output, output_dir, ignore=[]):
    """Compare output from pyinseq runs to expected output"""
    dcmp = filecmp.dircmp(expected_output, output_dir, ignore=ignore)
    # check that log file is created from pyinseq
    assert "log.txt" in os.listdir(output_dir)

    # checks that files are same in both directories
    assert not dcmp.left_only and not dcmp.right_only

    # check files to see if content differs
    assert not dcmp.diff_files
    assert not dcmp.funny_files  # Check for files that cannot be compared

    # because subdirs is a dict keyed by subdir name
    # with dircmp objects as values
    for subdcmp in dcmp.subdirs.values():
        assert not subdcmp.diff_files
        assert not subdcmp.funny_files
        assert not subdcmp.left_only and not subdcmp.right_only
