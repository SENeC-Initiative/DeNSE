# -*- coding: utf-8 -*-
#
# headers.py
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

"""
Check if all files have proper copyright header, add it otherwise.

This script checks the copyright headers of all C/C++/Cython/Python
files in the source code against the corresponding templates defined
in "copyright_header.*". It uses the relative path to determine the
path to DeNSE source directory.
"""

import os
import re
import sys


source_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# Use encoding-aware Py3 open also in Py2
if sys.version_info[0] < 3:
    from io import open


exclude_dirs = [
    '.git',
    'CMakeFiles',
    'build',
    'docs',
    'src/pymodule/dense/environment'
]

# match all file names against these regular expressions. if a match
# is found the file is excluded from the check
exclude_file_patterns = ['\.#.*', '#.*', '.*~', '.*.bak', '.*.txt']
exclude_file_regex = [re.compile(pattern) for pattern in exclude_file_patterns]

exclude_files = [
    'extra/copyright_header.cpp',
    'extra/copyright_header.py',
    'src/libgrowth/getindex.hpp',
    'src/libgrowth/cttrie.hpp',
    'src/libgrowth/cstr.hpp',
    'src/libgrowth/stringview.hpp',
]

templates = {
    ('cpp', 'hpp'): 'cpp',
    ('py', 'pyx', 'pxd'): 'py',
}

template_contents = {}

for extensions, template_ext in templates.items():
    template_name = "{0}/extra/copyright_header.{1}".format(source_dir,
                                                            template_ext)
    with open(template_name) as template_file:
        template = template_file.readlines()
        for ext in extensions:
            template_contents[ext] = template

to_be_corrected = []

for dirpath, _, fnames in os.walk(source_dir):

    if any([exclude_dir in dirpath for exclude_dir in exclude_dirs]):
        continue

    for fname in fnames:
        if any([regex.search(fname) for regex in exclude_file_regex]):
            continue

        extension = os.path.splitext(fname)[1][1:]
        if not (extension and extension in template_contents.keys()):
            continue

        tested_file = os.path.join(dirpath, fname)
        filename = tested_file[tested_file.rfind("/") + 1:]

        if any([exclude_file in tested_file
                for exclude_file in exclude_files]):
            continue

        with open(tested_file, encoding='utf-8') as source_file:
            for template_line in template_contents[extension]:
                try:
                    line_src = source_file.readline()
                except Exception as e:
                    print("[ERROR]: invalid character in {}".format(
                          tested_file))
                    print(line_src)
                    raise e
                line_exp = template_line.replace('{{file_name}}', filename)

                if line_src != line_exp:
                    fname = os.path.relpath(tested_file)
                    print("[ERROR]: {0}: expected '{1}', found '{2}'.".format(
                          filename, line_exp.rstrip('\n'),
                          line_src.rstrip('\n')))
                    print("... {}".format(filename))
                    to_be_corrected.append(tested_file)
                    break

# correct incorrect headers

for incorrect_file in to_be_corrected:
    print("[correcting]: {}".format(incorrect_file))
    extension    = os.path.splitext(incorrect_file)[1][1:]
    last_line_nb = 0

    # save previous content
    with open(incorrect_file, "r", encoding='utf-8') as source_file:
        content = source_file.readlines()

    # count first lines starting with comment
    with open(incorrect_file, "r", encoding='utf-8') as source_file:
        line   = source_file.readline()
        cstart = line.startswith("/*")

        if extension == "py":
            while line.startswith("#"):
                last_line_nb += 1
                line = source_file.readline()
        elif cstart:
            while not line.strip().endswith("*/"):
                last_line_nb += 1
                line = source_file.readline()
            last_line_nb += 1

    # remove them
    content = "".join(content[last_line_nb:])

    # finally, replace file content
    with open(incorrect_file, "w", encoding='utf-8') as target_file:
        data = "".join(template_contents[extension])
        filename = incorrect_file[incorrect_file.rfind("/") + 1:]
        data = data.replace('{{file_name}}', filename)
        target_file.write(data)
        target_file.write(content)