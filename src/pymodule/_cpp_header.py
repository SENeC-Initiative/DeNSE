# -*- coding: utf-8 -*-
#
# _cpp_header.py
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


''' Remove cython header from _pygrowth.cpp '''


def clean_cpp(filename):
    with open(filename, "r+") as f:
        data = [line for line in f]

        # get the header
        i, idx_end = 0, 0
        searching = True
        while not idx_end and i < len(data):
            if "END: Cython Metadata */" in data[i]:
                idx_end = i + 1
            i += 1

        # remove header and any commented line overwrite
        include = True
        output = ""
        for line in data[idx_end:]:
            if line.strip().startswith("/*"):
                if "*/" not in line:
                    include = False
            if include:
                output += line
            if "*/" in line:
                include = True
        f.seek(0)
        f.write(output)
        f.truncate()

