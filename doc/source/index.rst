.. grackle documentation master file, created by
   sphinx-quickstart on Sun Mar 24 09:23:44 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to grackle's documentation!
===================================

Grackle is a chemistry and radiative cooling library for astrophysical 
simulations.  It is a generalized and trimmed down version of the chemistry 
network of the `Enzo <http://enzo-project.org>`_ simulation code.  It 
provides a non-equilibrium primordial chemistry network for atomic H and He 
as well as H\ :sub:`2`\  and HD; tabulated metal cooling rates calculated 
with the photo-ionization code, `Cloudy <http://nublado.org>`_; and 
photo-heating and photo-ionization from the UV backgrounds of 
`Faucher-Giguere et al. (2009) 
<http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`_ or 
`Haardt & Madau (2012) 
<http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`_.

Contents:

.. toctree::
   :maxdepth: 2

   Installation.rst
   Integration.rst
   Parameters.rst
   Python.rst

Help
----

If you have any questions, please join the `Grackle Users Google Group
<https://groups.google.com/forum/#!forum/grackle-cooling-users>`_.  Feel 
free to post any questions or ideas for development.

Citing grackle
--------------

The Grackle library was born out of the chemistry and cooling routines of the 
`Enzo <http://enzo-project.org/>`_ simulation code.  As such, all of those who 
have contributed to Enzo development, and especially to the chemistry and 
cooling, have contributed to the Grackle.  There is currently no paper that 
specifically presents the Grackle library on its own, but the functionality 
was fully described in the `Enzo method paper 
<http://adsabs.harvard.edu/abs/2013arXiv1307.2265T>`_.  The Grackle was 
originally designed for the `AGORA Project 
<https://sites.google.com/site/santacruzcomparisonproject/>`_ and first referred 
to by name in the `AGORA method paper 
<http://adsabs.harvard.edu/abs/2013arXiv1308.2669K>`_.

If you used the Grackle library in your work, please cite it as "the Grackle 
chemistry and cooling library (`The Enzo Collaboration et al. 2013 
<http://adsabs.harvard.edu/abs/2013arXiv1307.2265T>`_; `Kim, J. et al. 2013 
<http://adsabs.harvard.edu/abs/2013arXiv1308.2669K>`_)."

  The Enzo Collaboration, Bryan, G. L., Norman, M. L., et al. 2013, arXiv:1307.2265

  Kim, J.-h., Abel, T., Agertz, O., et al. 2013, arXiv:1308.2669 







Search
------

* :ref:`search`

