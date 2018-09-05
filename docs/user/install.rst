======================
Installing the library
======================

The configuration and installation of |name| is provided by a ``cmake_``
script which provides a convenient multiplatform experience.


Requirements
============

Main library:

* cmake_
* python_
* numpy_, scipy_

Optional:

* cython_
* doxygen_ (other than 1.8.13), sphinx_, and breathe_ for local documentation,
* matplotlib_ for graphical representation.


Install
=======

Commands
--------

I will consider that you downloaded the sources into a folder called
``DeNSE`` and that you created a ``dense_install`` folder somewhere.::

    cd DeNSE
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/abolute/path/to/dense_install
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


.. References

.. _cmake: https://cmake.org/
.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _cython: http://cython.org/
.. _doxygen: http://www.stack.nl/~dimitri/doxygen/
.. _sphinx: http://www.sphinx-doc.org/
.. _breathe: http://breathe.readthedocs.io/en/latest/
.. _matplotlib: http://matplotlib.org/
.. _python: https://www.python.org/
