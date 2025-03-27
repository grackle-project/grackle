.. grackle documentation master file, created by
   sphinx-quickstart on Sun Mar 24 09:23:44 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grackle's documentation!
===================================

Grackle is a chemistry and radiative cooling library for astrophysical 
simulations and models.  Grackle has interfaces for C, C++, Fortran, and
Python codes and provides:

- two options for primordial chemistry and cooling:

   1. non-equilibrium primordial chemistry network for atomic H, D, and He
      as well as H\ :sub:`2`\  and HD, including H\ :sub:`2`\  formation on
      dust grains.

   2. tabulated H and He cooling rates calculated with the photo-ionization
      code, `Cloudy <http://nublado.org>`_.

- tabulated metal cooling rates calculated with `Cloudy <http://nublado.org>`_.

- photo-heating and photo-ionization from two UV backgrounds:

   1. `Faucher-Giguere et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`_.

   2. `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_.

- support for user-provided arrays of volumetric and specific heating rates.

The Grackle provides functions to update chemistry species; solve radiative 
cooling and update internal energy; and calculate cooling time, temperature, 
pressure, and ratio of specific heats (gamma).

.. toctree::
   :maxdepth: 1
   :hidden:

   Introduction <self>

.. toctree::
   :caption: Getting Started
   :maxdepth: 2
   :hidden:

   Installation.rst
   Python.rst
   Testing.rst
   Getting Help <Help.rst>

.. toctree::
   :caption: Using Grackle in Your Code
   :hidden:
   :maxdepth: 2

   Integration.rst
   Interaction.rst
   Parameters.rst
   RateFunctions.rst
   Reference.rst

.. toctree::
   :caption: About the Project
   :hidden:

   Versioning.rst
   Conduct.rst
   Contributing.rst

.. toctree::
   :caption: Useful links
   :hidden:

   GitHub <https://github.com/grackle-project/grackle>
   Tracker <https://github.com/grackle-project/grackle/issues>

.. include:: Help.rst

Contributing
------------

Development of Grackle happens in the open on GitHub `here
<https://github.com/grackle-project/grackle>`__.  We welcome new
contributors. Please, see the :ref:`conduct`.  For a guide to developing
Grackle, see :ref:`contributing-code`.

Citing grackle
--------------

The Grackle method paper can be found
`here <http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`__.

The Grackle library was born out of the chemistry and cooling routines of the 
`Enzo <http://enzo-project.org/>`__ simulation code.  As such, all of those who 
have contributed to Enzo development, and especially to the chemistry and 
cooling, have contributed to the Grackle.

If you used the Grackle library in your work, please cite it as "the Grackle 
chemistry and cooling library (`Smith et al. 2017
<http://adsabs.harvard.edu/abs/2017MNRAS.466.2217S>`_)."  Also, please add a
footnote to `https://grackle.readthedocs.io/ <https://grackle.readthedocs.io/>`_.


Search
------

* :ref:`search`

