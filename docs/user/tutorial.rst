.. _tuto:

========
Tutorial
========

A quick tour of what can be done with |name|.

.. contents::
    :local:


A single neuron
===============


Creating one neuron: 

    neuron_prop = {
        "position": (5., 12.)*um,
        "description": "my_special_neuron"
    }

    axon_prop = {
        "speed_growth_cone": 1.*um/hour,
        "persistence_length": 300.*um
    }

    neuron = ds.create_neurons(
        params=neuron_prop, num_neurites=3, # axon + 2 dendrites
        axon_params=axon_prop)


Embedding neurons in space
==========================


Complex structures
==================


Generating neuronal networks
============================
