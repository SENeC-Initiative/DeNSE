# -*- coding: utf-8 -*-
#
# conf_helpers.py
#
# This file is part of DeNSE.
#
# Copyright (C) 2019 SeNEC Initiative
#
# DeNSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# DeNSE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DeNSE. If not, see <http://www.gnu.org/licenses/>.

import importlib
import inspect
import os
import re


def gen_autosum(source, target, module, autotype, dtype="all", ignore=None):
    '''
    Automatically write a sphinx-parsable file, adding a list of functions or
    classes to the autosummary method of sphinx, in place of the @autosum@
    keyword.

    Parameters
    ----------
    source : str
        Name of the input file, usually of the form "source.rst.in".
    target : str
        Name of the output file, usually "source.rst".
    module : str
        Name of the module on which autosummary should be performed.
    autotype : str
        Type of summary (normal for all, 'autofunction' or 'autoclass').
    dtype : str, optional (default: all)
        Type of object that should be kept ('func' or 'class'), depending on
        `autotype`.
    ignore : list, optional (default: None)
        Names of the objects that should not be included in the summary.
    '''
    # load module and get content
    mod = importlib.import_module(module)
    mod_dir = dir(mod)
    # set ignored classes
    ignore = [] if ignore is None else ignore
    # list classes and functions
    str_autosum = ''
    for member in mod_dir:
        if not member.startswith('_') and not member in ignore:
            m = getattr(mod, member)
            keep = 1
            if dtype == "func":
                keep *= inspect.isfunction(m)
            elif dtype == "class":
                keep *= inspect.isclass(m)
            else:
                keep *= inspect.isfunction(m) + inspect.isclass(m)
            if keep:
                if autotype == "summary":
                    str_autosum += '    ' + module + '.' + member + '\n'
                else:
                    str_autosum += '\n.. ' + autotype + ':: ' + member + '\n'
    # write to file
    with open(source, "r") as rst_input:
        with open(target, "w") as main_rst:
            for line in rst_input:
                if line.find("%autosum%") != -1:
                    main_rst.write(str_autosum)
                else:
                    main_rst.write(line)


def pygrowth_mocker():
    '''
    Generate a pseudo _pygrowth.py file to mock the cython import.
    '''
    root_path = os.path.abspath("..")
    doc_path = os.path.abspath(".")

    ignore = {"init", "finalize"}


    with open(root_path + "/src/pymodule/dense/_pygrowth.pyx", "r") as pgr:
        with open(doc_path + "/pg_mock.py", "w") as pgw:
            raw_data = pgr.read()

            regexp = re.compile(r"^def\s+(?P<name>\w+)\s*\((?P<arg>[^)]*)\):[\s\n]+'''(?P<doc>.*?)'''",
                                re.MULTILINE|re.DOTALL)

            functions = re.finditer(regexp, raw_data)
            fnames    = []

            for f in functions:
                name = f.group("name")

                if name not in ignore and not name.startswith("_"):
                    pgw.write(f.group())
                    pgw.write("\n    pass\n\n")
                    # store function name
                    fnames.append(name)

            # write __all__
            pgw.write("__all__ = [\n")
            for name in fnames:
                pgw.write("    '" + name + "',\n")
            pgw.write("]\n\n")

            # mock init and finalize functions
            pgw.write("def init(argv):\n    pass\n\n")
            pgw.write("def finalize():\n    pass\n")


def configure_doxygen():
    ''' Write Doxyfile from .in '''
    root_path = os.path.abspath("..")
    doc_path = os.path.abspath(".")

    with open(doc_path + "/Doxyfile.in", "r") as dfin:
        data = dfin.read()
        data = data.replace("@PY_NAME@", "DeNSE")
        data = data.replace("@PROJECT_SOURCE_DIR@", root_path + "/src")

        with open(doc_path + "/Doxyfile", "w") as dfout:
            dfout.write(data)