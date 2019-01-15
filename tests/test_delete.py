#-*- coding:utf-8 -*-

""" Testing delete functions """

import numpy as np

import dense as ds
from dense.units import *


def test_delete_neurons():
    '''
    Neurons deletion
    '''
    ds.reset_kernel()

    num_neurons = 50
    simtime     = 2.*minute

    ds.set_kernel_status("environment_required", False)

    # create and delete
    neurons = ds.create_neurons(num_neurons, params={"position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um, "growth_cone_model": "res_po_rt"},
                                num_neurites=2)

    ds.delete_neurons(3)
    ds.delete_neurons(neurons[8])

    # check neurons have been deleted
    neurons_got = ds.get_neurons()
    n_to_ints   = [int(n) for n in neurons_got]
    n_to_ints  += ds.get_neurons(True)

    assert len(n_to_ints) == 2*(num_neurons - 2)
    assert 3 not in n_to_ints
    assert 8 not in n_to_ints

    # simulate then delete
    ds.simulate(simtime)
    ds.delete_neurons([5, 7])
    ds.delete_neurons(neurons[40:45])

    # recreate neurons and resimulate
    _ = ds.create_neurons(num_neurons, params={"position": np.random.uniform(-1000, 1000, (num_neurons, 2))*um, "growth_cone_model": "res_po_rt"}, 
                                 num_neurites=2)
    ds.simulate(simtime)

    neurons_got = ds.get_neurons()
    n_to_ints   = [int(n) for n in neurons_got]
    n_to_ints  += ds.get_neurons(True)

    assert len(n_to_ints) == 2*(2*num_neurons - 4 - len(neurons[40:45]))
    assert 5 not in n_to_ints
    assert 7 not in n_to_ints
    for n in neurons[40:45]:
        assert int(n) not in n_to_ints

    # delete all neurons
    ds.delete_neurons()

    # check no neurons are left
    assert not ds.get_neurons()


def test_delete_neurites():
    '''
    Neurons deletion
    '''
    ds.reset_kernel()

    num_neurons = 50
    simtime     = 2.*minute

    ds.set_kernel_status("environment_required", False)

    # create and delete
    pos     = np.random.uniform(-1000, 1000, (num_neurons, 2))*um
    neurons = ds.create_neurons(num_neurons, params={"position": pos},
                                num_neurites=2)

    ds.delete_neurites("axon")

    for n in neurons:
        assert not n.has_axon
        assert "axon" not in n.neurites
        assert len(n.neurites) == 1

    ds.create_neurites(neurons[0], neurite_types="axon")

    assert neurons[0].has_axon
    assert len(neurons[0].neurites) == 2

    ds.delete_neurites("dendrite_1", neurons[1])

    for n in neurons:
        if n == neurons[1]:
            assert "dendrite_1" not in n.neurites
            assert not n.neurites
        elif n == neurons[0]:
            assert len(n.neurites) == 2
        else:
            assert len(n.neurites) == 1

    ds.delete_neurites()

    for n in neurons:
        assert not n.neurites


if __name__ == "__main__":
    test_delete_neurons()
    test_delete_neurites()
