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

Contents:

.. toctree::
   :maxdepth: 2

   Installation.rst
   Testing.rst
   Integration.rst
   Parameters.rst
   Reference.rst
   Python.rst
   Conduct.rst
   Contributing.rst           

Help
----

If you have any questions, please join the `Grackle Users Google Group
<https://groups.google.com/forum/#!forum/grackle-cooling-users>`_.  Feel 
free to post any questions or ideas for development.

Contributing
------------

Development of Grackle happens in the open on Bitbucket `here
<https://bitbucket.org/grackle/grackle>`_.  We welcome new contributors.
Please, see the :ref:`conduct`.  For a guide to developing Grackle, see
:ref:`contributing-code`.

Citing grackle
--------------

The Grackle method paper can be found
`here <https://arxiv.org/abs/1610.09591>`__.

The Grackle library was born out of the chemistry and cooling routines of the 
`Enzo <http://enzo-project.org/>`__ simulation code.  As such, all of those who 
have contributed to Enzo development, and especially to the chemistry and 
cooling, have contributed to the Grackle.

If you used the Grackle library in your work, please cite it as "the Grackle 
chemistry and cooling library (`Smith et al. 2016
<https://arxiv.org/abs/1610.09591>`_)."  Also, please add a footnote to
`https://grackle.readthedocs.org/ <https://grackle.readthedocs.org/>`_.



Search
------

* :ref:`search`

