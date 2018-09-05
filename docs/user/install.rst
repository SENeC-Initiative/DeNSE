======================
Installing the library
======================

The configuration and installation of |name| is provided by a ``cmake_``
script which provides a convenient multiplatform experience.


Requirements
============

Main library:

* cmake_, python_, GEOS_ (install both binary and header/dev files via apt on debian/ubuntu or rpm/pacman on other linux distributions)
* numpy_, scipy_, shapely_, pint_ (install via ``pip install --user numpy scipy shapely pint``)

Optional:

* cython_ (``pip install --user cython``)
* doxygen_ (other than 1.8.13), sphinx_, and breathe_ for local documentation,
* matplotlib_ for graphical representation (``pip install --user matplotlib``).


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
