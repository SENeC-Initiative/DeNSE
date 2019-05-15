#!/usr/bin/env python
#-*- coding:utf-8 -*-

import hashlib, json


""" Utilities """


def HashID(*args):
    experiment_dict = {}
    for num, dict_ in enumerate(args):
        experiment_dict[num] = dict_
    return _hash_dict(experiment_dict)

def _hash_dict(_dict):
    sha=hashlib.sha1()
    sha.update(str(json.dumps(_dict, sort_keys =True)).encode('utf-8'))
    hash_name = sha.hexdigest()
    return hash_name[:16]
