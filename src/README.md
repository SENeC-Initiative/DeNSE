# Neural network growth project

## Organization

* Main objects are written in c++
    + Headers: ``Branch.h``, ``Neurite.h``, ``Neuron.h``
    + **Neuron:** main container, this is where you define the average properties (axonal and dendritic) that will be used during growth.
    + **Neurite:** this container holds a link to a random number generator and the specific parameters for its growing process; it is where each of new segments characteritics (length, angle) are generated.
    + **Branch:** container (created at each ``branching`` event of a ``Neurite``) that holds a vector of all points of a ``Neurite``'s line.

* Reasons for the objects' existence:
    + ``Neuron`` is obvious,
    + ``Neurite`` is an element composed of parts that will all grow accordong to the same parameters,
    + ``Branch`` is necessary because it is the only time-sequential container.

* Environment
    + Environment can be read from a `.svg` or a `.dfx` file
    + Geometric operations are implemented using Boost::geometry


## Current limitations

* ``Branch`` currently stores a ``std::vector<Point>`` where ``Point`` is itself a ``std::vector<double>``; this might be too heavy.
    + **Idea:** define at the beginning whether we are in 2 or 3D, then define ``typedef std::tuple<double, double> Point`` in 2D and with an additional ``double`` in 3D.
    + **Or** even lighter: use a ``std::tuple<std::vector<double>, std::vector<double>>`` which is the equivalent of a 2D array (make it 3D if necessary).
* I don't know how to make sure the threads do not accidentaly modify variables in other threads in ``cython``.
* The random generators/seeds are not well implemented (too many of them).
* Branching process is simplistic.
* The soma is a point.


## TODO

### Code

* [x] Create a ``Configuration`` object which will contain the parameters of the simulation (e.g. the number of mpi/omp processes, the timestep...)
* [x] Make one random number generator per ``omp`` process inside ``Area`` and pass a pointer to the neurites or use one RNG per neuron.
* Implement the obstacles and environment with Boost::geometry
* set branching proba to be timestep independant (square root of timestep)

### Model

* Rules for the end of growth
* Interactions
* Synapse formation

We need to have a look at simple existing homeostatic models to be able to build from them.
NEST apparently includes examples:

* [structural_plasticity](http://www.nest-simulator.org/py_sample/structural-plasticity-example/) + [.py](https://github.com/nest/nest-simulator/blob/master/pynest/examples/structural_plasticity.py)
* Linked to following articles: [PlosOne 2013](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003259) and [Front. Synaptic Neuro. 2014](http://journal.frontiersin.org/article/10.3389/fnsyn.2014.00007/full)


## Problems, caveats and tests

### Doxygen

Bug combining XML and ``friend class`` in 1.8.13, so documentation cannot be
generated with this version of Doxygen.

### Python 3

There is a difference between strings (``str``) and bytes (``bytes``); all characters passed to C++ ``std::string`` apparently need to be ``bytes``.

### C++

Careful when passing arguments to function:
* ``func(std::vector<double> vec)`` means that a (temporary) local copy of ``vec`` will be used inside ``func``,
* ``func(std::vector<double>& vec)`` means that the real object will be used (and potentially modified) inside ``func``,
* ``func(const std::vector<double>& vec)`` means that the real object will be used but that it cannot be modified,
* ``func(std::vector<double>* vec)`` means that a pointer to the object will be passed; using the pointer, the real object can also be modified.

### Tests
