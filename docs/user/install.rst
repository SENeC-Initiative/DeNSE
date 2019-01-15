.. _install:

======================
Installing the library
======================

The configuration and installation of |name| is provided by a ``cmake_``
script which provides a convenient multiplatform experience.

**NB** if you run against some bugs after the installation, please first have
a look at the `Known bugs`_ section, in case it is not already listed there.


Requirements
============

Main library:

* cmake_ 3.6+, python_ 2.7+, GEOS_ 3.5+ (install both binary and header/dev files via apt on debian/ubuntu or rpm/pacman on other linux distributions)
* numpy_, scipy_, shapely_, pint_ (install via ``pip install --user numpy scipy shapely pint``)

Optional:

* cython_ (``pip install --user cython``)
* doxygen_ (other than 1.8.13), sphinx_, and breathe_ for local documentation,
* matplotlib_ for graphical representation (``pip install --user matplotlib``).
* seaborn_  for graphical representation (``pip2 install --user seaborn``)
* neurom_ for graphical representation (``pip2 install --user neurom``)
* h5py_ for neurom ( ``apt install libhdf5-dev`` then ``pip2 install --user h5py``)


On a debian-based linux install
-------------------------------

Non-python part:
* Need g++ + python-dev to compile matplotlib (maybe also automake)
* say that doxygen is via apt
* say that we need to install python and python-dev
* sudo apt install libgeos-X.Y.Z et libgeos-dev libgeos++-dev

Python part:

* Install pip
* Install setuptools, then python libraries
* say that sphinx and breathe are via pip
* use sudo apt install python-tk python-matplotlib

On a debian-based linux install
-------------------------------

Non-python part:
* Need g++ + python-dev to compile matplotlib (maybe also automake)
* say that doxygen is via apt
* say that we need to install python and python-dev
* sudo apt install libgeos-X.Y.Z et libgeos-dev libgeos++-dev

Python part:

* Install pip
* Install setuptools, then python libraries
* say that sphinx and breathe are via pip
* use sudo apt install python-tk python-matplotlib

Install
=======

Commands
--------

I will consider that you downloaded the sources into a folder called
``DeNSE`` and that you created a ``dense_install`` folder somewhere.::

    cd DeNSE
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/absolute/path/to/dense_install
    make && make install

.. warning::
    ``-DCMAKE_INSTALL_PREFIX:PATH`` is **required** and **must** be given as
    an absolute path!


Options
-------

You can provide the following options to cmake using ``-Doption_name``, where
``option_name`` can be:

* ``with-python``, which should be either ``2`` or ``3`` depending on the
  version of Python you want to use,
* ``with-docs`` -- ``ON`` (default) or ``OFF`` if you do not want to install
  the documentation on your computer,
* ``with-openmp`` -- ``ON`` (default) or ``OFF`` if you do not want to use
  multithreading.
* ``with-debug`` -- ``ON`` to compile with debug symbols or ``OFF`` (default).


Known bugs
==========

* Any random error message
    - you forgot to empty the ``build`` folder before reexecuting cmake/make?
* Script runs into
  ``Assertion failed: (!static_cast<bool>("should never be reached"))``.
    - This is due to a bug in Shapely (#553_) and can be fixed by installing it
      without its binaries: ``pip install --no-binary shapely, shapely``
      (note that you first have to uninstall it before reinstalling it).
* Segmentation fault at the end of the script, after all the code has run:
    - This often seems to be a matplotlib bug for which the fix is unknown, but
      it does not worry us too much because, as mentioned, all the code is
      executed properly...
    - Alternatively, for people importing NNGT and with NEST installed, it can
      be caused by an error in the finalize function of NEST, so again, nothing
      to worry about.


.. References

.. _cmake: https://cmake.org/
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _GEOS: https://trac.osgeo.org/geos/
.. _shapely: http://toblerity.org/shapely/manual.html
.. _pint: https://pint.readthedocs.io/en/latest/
.. _cython: http://cython.org/
.. _doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _sphinx: http://www.sphinx-doc.org/
.. _breathe: http://breathe.readthedocs.io/en/latest/
.. _matplotlib: http://matplotlib.org/
.. _python: https://www.python.org/
.. _seaborn:  https://seaborn.pydata.org/
.. _neurom: https://github.com/BlueBrain/NeuroM
.. _#553: https://github.com/Toblerity/Shapely/issues/553
