.. _python:

The Python Examples
===================

These example works with a python wrapper that calls the various library 
functions.  These require `Cython <http://www.cython.org/>`_ to be installed.  
The best thing to do is to install the `yt analysis toolkit 
<http://yt-project.org>`_, which includes Cython.

Installing the Python Wrappers
------------------------------

After building the grackle library, some additional environment variables must 
be set.  Move into the **src/python** directory and run the **set_libs** script 
and follow the instructions.  Then, run *python setup.py install* to build and 
install the python wrappers.

.. highlight:: none

::

    ~/grackle $ cd src/python
    ~/grackle/src/python $ ./set_libs
    Issue the following commands:
    export PYTHONPATH=$PYTHONPATH:/grackle/src/python
    export LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/grackle/src/python/../clib
    You can also set your LD_LIBRARY_PATH to point to where you installed libgrackle.
    ~/grackle/src/python $ python setup.py install
    running install
    running build
    running config_cc
    [clipped]
    running install_clib
    customize UnixCCompiler

Running the Example Scripts
---------------------------

The python examples are located in the **src/python/examples** directory.  Before 
running them, make sure to copy the data files from the **input** directory into this 
directory.


Cooling Rate Figure Example
+++++++++++++++++++++++++++

This sets up a one-dimensional grid at a constant density with 
logarithmically spaced temperatures from 10 K to 10\ :sup:`9` K.  Radiative cooling 
is disabled and the chemistry solver is iterated until the species fractions 
have converged.  The cooling time is then calculated and used to compute the cooling 
rate.  This script also provides a good example for setting up cosmological unit 
systems.  

.. highlight:: none

::

    python cooling_rate.py

.. image:: _images/cooling_rate.png
   :width: 500

Cooling Cell Example
++++++++++++++++++++

This sets up a single grid cell with an initial density and temperature and solves 
the chemistry and cooling for a given amount of time.  The script will output a table 
of temperature vs. time.  This can be used to test implementations of the Grackle in 
simulation codes.

.. highlight:: none

::

    python cooling_rate.py

.. image:: _images/cooling_cell.png
   :width: 500

Free-Fall Collapse Example
++++++++++++++++++++++++++

This sets up a single grid cell with an initial number density of 1 cm\ :sup:`-3`.  
The density increases with time following a free-fall collapse model.  As the density 
increases, thermal energy is added to model adiabatic compression heating.  This can be 
useful for testing chemistry networks over a large range in density.

.. highlight:: none

::

    python freefall.py

.. image:: _images/freefall.png
   :width: 500

Simulation Dataset Example
++++++++++++++++++++++++++

This provides an example of using the grackle library for calculating chemistry and 
cooling quantities for a pre-existing simulation dataset.  To run this example, you 
must have `yt <http://yt-project.org>`_ installed and must also download the 
*IsolatedGalaxy* dataset from the `yt sample data page <http://yt-project.org/data/>`_.

.. highlight:: none

::

    python run_from_yt.py
