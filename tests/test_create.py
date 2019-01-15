#-*- coding:utf-8 -*-

""" Testing delete functions """

import dense as ds
from dense.units import *


def test_create():
    '''
    Create neurons and neurites
    '''
    ds.reset_kernel()

    # create one neuron
    pos    = (0., 0.)*um
    neuron = ds.create_neurons(params={"position": pos}, num_neurites=0)

    # create a new neurite
    neuron.create_neurites(names="new_dendrite")

    assert len(neuron.neurites) == 2
    assert "new_dendrite" in neuron.neurites


if __name__ == "__main__":
    test_create()
