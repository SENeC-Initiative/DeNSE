#!/usr/bin/env python
#-*- coding:utf-8 -*-

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

