
.. _capability:

What does Grackle do?
=====================

Grackle is conceptually simple.
At its core, the library provides machinery to specify the :ref:`"physical model" <what-is-a-physical-model>` for an astrophysical fluid (i.e. gas or plasma), (ii) a handful of functions using the model to compute the fluid's physical properties (e.g. temperature), and (iii) a function useing the model to solve the fluid's time evolution (from chemical/radiative/thermodynamic processes).


The fine-grained configurability of the physical model is Grackle's "killer feature."
Users can choose from a rich selection of battle-tested options, which allows Grackle to:

- support state of the art calculations capturing lots of detailed physics
- allows users to :ref:`tailor the model <why-tailor>` to suit the needs of a given calculation

.. _what-is-a-physical-model:

What is a "physical model"?
---------------------------

In the current discussion we define the "physical model" as the collection of assumptions and physics that affect the properties (e.g. thermal properties/composition) of the fluid, itself.
Grackle primarily cares about details relevant to a fluid element's **local micro-scale evolution**; it doesn't need to know any details about interactions between fluid elements (gravity, fluid dynamics, viscosity, conduction, etc.).

Fluid calculations **almost always** have a physical model.
There's usually an assumption about the equation of state.
An assumption is also made about radiative effects (or the lack thereof).
Grackle becomes useful to calculations as soon as you want to start considering how radiative cooling (and heating) affect your calculation (**and it can do much more**).

.. _why-tailor:

Why you tailor a physical model?
--------------------------------

Part of Grackle's value-proposition is that it lets users adopt a physical model that is "sophisticated enough" to suit their needs. [#f1]_

It may be instructive to consider a few scenarios that Grackle could support.
We order these in terms of increasing complexity:

- There are lots of cases with idealized simulations, where you are trying to tease out the impacts of other physics and you just want a realistic-enough cooling curve that is easy to reason about.
  You might want to artificially stop cooling below some temperature in certain cases.
  More concretely, this might include cloud-wind interactions or galactic wind simulations.

- You might be running a cosmological simulation, where you want to solve a basic chemical network of atomic primordial species, but you aren't running at high enough resolution to benefit from modelling molecular species.

- You might be running high resolution simulations focused on the formation of Population III stars, where a more complete chemical network of primordial species is useful.

- Then, of course there are the simulations where you want as much sophistication as possible (e.g. you're trying to model self-consistent star formation and feedback).
  You might even couple a radiative transfer solver to Grackle.
  You could also couple external sources of heating.

.. COMMENT-BLOCK

   If we add more to the following section, maybe we should spin it off as a separate part of the user-guide?

.. _capabilities:

Grackle's Capabilities
----------------------

We provide a brief overview of Grackle's capabilities.
You should read through the :ref:`available parameters (and data tables) <parameters>` for more details.

Primordial Chemistry
^^^^^^^^^^^^^^^^^^^^

Grackle provides several options for modelling the compositionof the gas.
This is controlled by the integer value of the :c:var:`primordial_chemistry` parameter.
We summarize the different values down below:

.. table::
   :widths: 1 8 12 15

   +-----+-------------------+----------------------------------------------+-----------------------------------------------+
   | val | Description       | Chemical species                             | Required fields                               |
   +=====+===================+==============================================+===============================================+
   | 0   | Tabulated (Equilibirum) mode                                     | :c:var:`density`, :c:var:`internal_energy`    |
   |     |                                                                  |                                               |
   |     | Assumes gas is composed of atomic e\ :sup:`-`, H, H\ :sup:`+`,   |                                               |
   |     | He, He\ :sup:`+`, He\ :sup:`++` that are always in equilibrium.  |                                               |
   |     | Cooling rates & mean molecular weights come from tables          |                                               |
   |     | calculated with the photo-ionization code, `Cloudy`_.            |                                               |
   |     |                                                                  |                                               |
   +-----+-------------------+----------------------------------------------+-----------------------------------------------+
   | 1   | Non-Equilibrium   | e\ :sup:`-`,                                 | above fields + :c:var:`e_density`,            |
   |     | Chemical Network  | H, H\ :sup:`+`,                              | :c:var:`HI_density`, :c:var:`HII_density`,    |
   |     | Solver            | He, He\ :sup:`+`, He\ :sup:`++`              | :c:var:`HeI_density`, :c:var:`HeII_density`   |
   |     |                   |                                              | :c:var:`HeIII_density`                        |
   +-----+                   +----------------------------------------------+-----------------------------------------------+
   | 2   |                   | above species + H\ :sup:`-`, H\ :sub:`2`,    | above fields + :c:var:`HM_density`,           |
   |     |                   | H\ :sub:`2`\ :sup:`+`                        | :c:var:`H2I_density`, :c:var:`H2II_density`   |
   +-----+                   +----------------------------------------------+-----------------------------------------------+
   | 3   |                   | above species + D,                           | above fields + :c:var:`DI_density`,           |
   |     |                   | D\ :sup:`+`, HD                              | :c:var:`DII_density`, :c:var:`HDI_density`    |
   +-----+-------------------+----------------------------------------------+-----------------------------------------------+

.. note:: In Development: A more extended primordial network

Dust Chemistry
^^^^^^^^^^^^^^
Grackle has basic dust chemistry to model the influence of a non-evolving dust population on the chemothermal evolution of the gas. This includes H\ :sub:`2` formation on dust grains, photo-electric heating, and dust recombination cooling.
This is controlled by the :c:data:`dust_chemistry` parameter.

Grackle optionally gives users fine-grained control over the interstellar radiation field used to compute the dust temperature and photoelectric heating.

.. note:: In Development: Dust Grain Chemistry Network

Metal Cooling
^^^^^^^^^^^^^
Grackle supports metallicity-dependant cooling from heavy elements (using tabulated rates calculate with `Cloudy`_).
This is enabled by the :c:data:`metal_cooling` parameter.

.. note:: In Development: Metal Chemistry Network

UV Background
^^^^^^^^^^^^^

Grackle supports photo-heating and photo-ionization from two models of the time-evolving, spatially-independent UV metal-galactic background:

   1. `Faucher-Giguere et al. (2009) <http://adsabs.harvard.edu/abs/2009ApJ...703.1416F>`__.

   2. `Haardt & Madau (2012) <http://adsabs.harvard.edu/abs/2012ApJ...746..125H>`__.

Configuration of this functionality is tied to the :ref:`data files <data-files>` that Grackle provides.
Additionally, the cooling and heating rates can be modified to include the effects of atomic self-shielding in a few different ways.
User-Extendability
^^^^^^^^^^^^^^^^^^

Grackle also provides support for extension beyond the above capabilities.
For examples users can provide arrays of volumetric and specific heating rates
Users can also couple Grackle to external radiative transfer solvers.


.. _Cloudy: http://nublado.org

.. rubric:: Footnotes

.. [#f1] If your inclination is to always use the most detailed model that captures the most physics, that's a valid choice, but it certainly isn't the only option.
         While adding more detailed physics theoretically improves accuracy, it also increases complexity.
         There are obvious performance considerations (more complex models are often more computationally expensive).
         It's also important to understand that the benefits you get from introducing physics with higher order effects may be marginal if some other part of your calculation is less accurate or uncertain.
         Furthermore, there are certain context where a simpler model can make it easier to tease out the underlying physics in a numerical simulation or generalize your results.
