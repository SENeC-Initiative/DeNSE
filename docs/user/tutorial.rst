.. _tuto:

========
Tutorial
========

A quick tour of what can be done with |name|.

.. contents::
    :local:


A single neuron
===============

First, one should import all the modules and variables that are necessary
for the simulation: ::

    import dense as ds
    from dense.units import *

The first line makes the DeNSE simulator available and let us call it through the variable ``ds``. The second line imports all the
units (e.g. time units like ``ms``, ``minutes``...) which will be used to set the properties of the simulation and of the neurons.

Once this is done, we can set the various parameters for the simulation and the neuronal properties:


Two interacting neurons
=======================

Again, import all necessary modules and variables:

.. literalinclude:: ../../examples/tutorials/1_first_steps.py
    :linenos:
    :language: python
    :lines: 33-34


The first line makes the \texttt{DeNSE} simulator available and let us call it through the variable \texttt{ds}. The second line imports all the
units (e.g. time units like \texttt{ms}, \texttt{minute}\ldots) which will be used to set the properties of the simulation and of the neurons.

Once this is done, we can set the various parameters for the simulation and the neuronal properties:

\begin{lstlisting}
num_omp       = 2
num_neurons   = 2

simu_params   = {
    "resolution": 1.*minute,
    "num_local_threads": num_omp,
    "seeds": [0, 1],
    "environment_required": False,
}

neuron_params = {
    "axon_diameter": 4.*um,
    "dendrite_diameter": 3.*um,
    "growth_cone_model": "run-and-tumble",
    "position": [(0., 0.), (100., 100.)]*um,
    "persistence_length": 200.*um,
    "speed_growth_cone": 0.03*um/minute,
    "taper_rate": 1./400.,
    "use_uniform_branching": True,
    "uniform_branching_rate": 0.009*cph,
}
\end{lstlisting}

The first line here declares the number of OpenMP processes that will be used, i.e. how many parallel threads will be used to perform the simulation.
The second line will be used to set the number of neurons that will be simulated.

The third part contains the parameters related to the simulation: the timestep used to run the simulation, the number of threads used in the parallel simulation, the random seeds that will be used to  generate random numbers (one per OpenMP thread), and whether the neurons are embedded in spatial boundaries.

Eventually, the \texttt{neuron\_params} dictionary on line 11 contains the information that will be used to describe the growth of the two neurons, all expressed with their proper units (distances in microns, speed in micron per minutes, and branching events in ``counts per hour'').  As no environment was specified here, the positions of the neurons is also specified; there orientation (the direction of the axon) will be set randomly. Though we used the same parameters for both neurons here, this is not necessary and different parameters can be passed for each neurons through a list, as shown for the positions on line 15.

Once all these parameters are declared, we can configure \texttt{DeNSE} and create the neurons so that everything is ready for the simulation.

\begin{lstlisting}
# configure DeNSE
ds.set_kernel_status(simu_params)

# create neurons
n = ds.create_neurons(n=num_neurons,
                      params=neuron_params,
                      num_neurites=2)
\end{lstlisting}

As can be seen above, one uses the \textit{\texttt{ds}} variable to access
the simulator main function. 
The \texttt{set\_kernel\_status} function is used here to transfer the parameters to the \texttt{DeNSE} kernel (the main simulator units).
Once this is done, the \texttt{create\_neurons} is called in order to obtain two neurons with the specific set of parameters declared in \texttt{neuron\_params}, and two neurites (by default the first is an axon, and the second is a dendrite) which initially possess a single branch protruding from the soma.

Once the neurons are created one can visualize them using \texttt{ds.plot.plot\_neurons()}, which outputs the picture shown on Figure \ref{fig:example}.\textbf{a}. The initial condition of the neurons can thus be visualized before starting the simulation.

Following neuron creation, the simulation can be started and its result can be visualized again and is shown on Figure \ref{fig:example}.\textbf{b}.

\begin{lstlisting}
ds.simulate(7*day)
ds.plot.plot_neurons()
\end{lstlisting}

After this first 7 day simulation, the parameters of the neurons can be changed to account for changes in developmental mechanisms, so that these new parameters can be used to simulate the next part of these cells' growth.

\begin{lstlisting}
# new dendritic and axonal parameters
axon_params = {
    "speed_growth_cone": 0.02*um/minute,
}

dend_params = {
    "use_uniform_branching": False,
    "speed_growth_cone": 0.01*um/minute,
}

# update the properties of the neurons
ds.set_object_properties(n, dendrites_params=dend_params,
                         axon_params=axon_params)

# simulate and plot again
ds.simulate(7*day)
ds.plot.plot_neurons()
\end{lstlisting}

\begin{figure}
\centering\includegraphics[scale=0.4]{example}

\caption{Visualization of the different steps in the example. \textbf{a.} Initial neurons with slight protrusions representing the positions of the neurites. \textbf{b.} Neurons after 7 days of growth. \textbf{c.} Neurons after 14 days of growth. Axons are in red, dendrites in blue, and somas in black. Scale bars are 50 $\mu$m.}

\label{fig:example}
\end{figure}

Here we changed separately the dendritic and axonal parameters using the \texttt{set\_object\_properties} function on the two neurons which are stored in the \texttt{n} variable.
The different changes in the properties of the dendrites and axons, notably through the suppression of lateral (or interstitial) branching in the dendrites, can be seen on Figure \ref{fig:example}.\textbf{c}.

Once the simulation is over, the shapes of the neurons that were obtained can be saved in standard morphology formats such as \texttt{SWC} or \texttt{MorphoML} (\texttt{NeuroML}).

\begin{lstlisting}
ds.io.save_to_swc(n, "neurons.swc")
ds.io.save_to_neuroml(n, "neurons.nml")
\end{lstlisting}

Eventually, once the neuronal growth has been simulated, one can ask \texttt{DeNSE} to generate connections between the neurons in order to obtain a set of synapses.
This information can then be used to simulate the activities of these two neurons and their interactions by inputting it into activity simulators such as \texttt{NEURON} or \texttt{NEST}.

\begin{lstlisting}
synapses = ds.morphology.get_synapses()
network  = ds.morphology.generate_network()
\end{lstlisting}

Embedding neurons in space
==========================


Complex structures
==================


Generating neuronal networks
============================
