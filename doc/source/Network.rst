.. _Network:

Chemical Network
==================
The table below shows the reactions in Grackle`s chemical network. All reactions are of the form left-hand side (LHS) -> right-hand side (RHS). The rate column corresponds to the interval table responsable for the reaction rate. The citation colomns denotates the publication in which reaction rate appears in the form used within Grackle. Note, some reaction rate have multiple options depending on the setting of various parameters.

.. list-table::
   :widths: auto
   :header-rows: 1

   * - LHS
     - RHS
     - Rate
     - Citation
     - Enabled By
   * - H + e\ :sup:`-`
     - H\ :sup:`+` + e\ :sup:`-` + e\ :sup:`-`
     - k1
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >0     
   * - H\ :sup:`+` +  e\ :sup:`-`
     - H +  :math:`{\gamma}`
     - k2
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__  & `Hui & Gnedin (1997) <https://ui.adsabs.harvard.edu/abs/1997MNRAS.292...27H/abstract>`__
     - :c:data:`primordial_chemistry` >0   
   * - He +  e\ :sup:`-`
     - He\ :sup:`+` +  e\ :sup:`-` +  e\ :sup:`-`
     - k3
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >0  
   * - He\ :sup:`+` + e\ :sup:`-`
     - He\ :sup:`+` + :math:`{\gamma}`
     - k4
     - `Hui & Gnedin (1997) <https://ui.adsabs.harvard.edu/abs/1997MNRAS.292...27H/abstract>`_  & `Aldrovandi & Pequignot (1973) <https://ui.adsabs.harvard.edu/abs/1973A%26A....25..137A/abstract>`_ & `Black (1981) <https://ui.adsabs.harvard.edu/abs/1981MNRAS.197..553B/abstract>`_
     - :c:data:`primordial_chemistry` >0   
   * - He\ :sup:`+` +  e\ :sup:`-`
     - He\ :sup:`+`\ :sup:`+` + e\ :sup:`-` + e\ :sup:`-`
     - k5
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >0     
   * - He\ :sup:`+`\ :sup:`+` + e\ :sup:`-`
     - He\ :sup:`+` +  :math:`{\gamma}`
     - k6
     - `Hui & Gnedin (1997) <https://ui.adsabs.harvard.edu/abs/1997MNRAS.292...27H/abstract>`__  & `Cen (1992) <https://ui.adsabs.harvard.edu/abs/1992ApJS...78..341C/abstract>`__
     - :c:data:`primordial_chemistry` >0   
   * - H + H
     - H\ :sup:`+` + e\ :sup:`-` + H
     - k57
     - `Lenzuni, Chernoff & Salpeter (1991) <https://ui.adsabs.harvard.edu/abs/1991ApJS...76..759L/abstract>`__
     - :c:data:`primordial_chemistry` >0   
   * - H + He
     - H\ :sup:`+` + e\ :sup:`-` + He
     - k58
     - `Lenzuni, Chernoff & Salpeter (1991) <https://ui.adsabs.harvard.edu/abs/1991ApJS...76..759L/abstract>`__
     - :c:data:`primordial_chemistry` >0  
   * - H + :math:`{\gamma}`
     - H\ :sup:`+` + e\ :sup:`-`
     - \-
     - none
     - :c:data:`primordial_chemistry` >0      
   * - He +  :math:`{\gamma}`
     - He\ :sup:`+` + e\ :sup:`-`
     - \-
     - none
     - :c:data:`primordial_chemistry` >0              
   * - He\ :sup:`+` + :math:`{\gamma}`
     - He\ :sup:`+` :sup:`+` + e\ :sup:`-`
     - \-
     - none
     - :c:data:`primordial_chemistry` >0                  
   * - H + e\ :sup:`-`
     - H\ :sup:`-` + :math:`{\gamma}`
     - k7
     - `Stancil et al. (1998) <https://ui.adsabs.harvard.edu/abs/1998ApJ...509....1S/abstract>`__
     - :c:data:`primordial_chemistry` >1                      
   * - H\ :sup:`-` + H
     - H\ :sub:`2` + e\ :sup:`-`
     - k49
     - `Kreckel et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010Sci...329...69K/abstract>`__
     - :c:data:`primordial_chemistry` >1                          
   * - H + H\ :sup:`+`
     - H\ :sub:`2`:sup:`+` + :math:`{\gamma}`
     - k9
     - `Latif et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.3163L/abstract>`__
     - :c:data:`primordial_chemistry` >1                     
   * - H\ :sub:`2`:sup:`+` + H
     - H\ :sub:`2` + H\ :sup:`+`
     - k10
     - `Karpas, Anicich & Huntress (1979) <https://ui.adsabs.harvard.edu/abs/1979JChPh..70.2877K/abstract>`__
     - :c:data:`primordial_chemistry` >1                 
   * - H\ :sub:`2` + H\ :sup:`+`
     - H\ :sub:`2`:sup:`+` + H
     - k11
     - Savin et al. (`2004a; <https://ui.adsabs.harvard.edu/abs/2004ApJ...606L.167S/abstract>`__  `b <https://ui.adsabs.harvard.edu/abs/2004ApJ...607L.147S/abstract>`__)
     - :c:data:`primordial_chemistry` >1           
   * - H\ :sub:`2` + e\ :sup:`-`
     - H + H + e\ :sup:`-`
     - k12
     - `Trevisan & Tennyson (2002) <https://ui.adsabs.harvard.edu/abs/2002PPCF...44.1263T/abstract>`__
     - :c:data:`primordial_chemistry` >1        
   * - H\ :sub:`2` + H
     - H + H + H
     - k13
     - `Martin, Schwarz & Mandy (1996) <https://ui.adsabs.harvard.edu/abs/1996ApJ...461..265M/abstract>`__
     - :c:data:`primordial_chemistry` >1                
   * - H\ :sup:`-` + e\ :sup:`-`
     - H + e\ :sup:`-` + e\ :sup:`-`
     - k14
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >1               
   * - H\ :sup:`-` + H
     - H + e\ :sup:`-` + H
     - k15
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >1                       
   * - H\ :sup:`-` + H\ :sup:`+`
     - H + H
     - k16
     - `Croft, Dickinson & Gadea (1999) <https://ui.adsabs.harvard.edu/abs/1999MNRAS.304..327C/abstract>`__
     - :c:data:`primordial_chemistry` >1
   * - H\ :sup:`-` + H\ :sup:`+`
     - H\ :sub:`2`:sup:`+` + e\ :sup:`-`
     - k17
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__  & `Shapiro & Kang. (1987) <https://ui.adsabs.harvard.edu/abs/1987ApJ...318...32S/abstract>`_
     - :c:data:`primordial_chemistry` >1                      
   * - H\ :sub:`2`:sup:`+` + e\ :sup:`-`
     - H + H
     - k18
     - `Abel et al. (1997) <https://ui.adsabs.harvard.edu/abs/1997NewA....2..181A/abstract>`__
     - :c:data:`primordial_chemistry` >1
   * - H\ :sub:`2`:sup:`+` + H\ :sup:`-`
     - H\ :sub:`2` + H
     - k19
     - `Dalgarno & Lepp (1985) <https://ui.adsabs.harvard.edu/abs/1987IAUS..120..109D/abstract>`__
     - :c:data:`primordial_chemistry` >1
   * - H + H + H
     - H\ :sub:`2` + H
     - k22
     - see table below
     - :c:data:`primordial_chemistry` >1
   * - H + H + H\ :sub:`2`
     - H\ :sub:`2`  + H\ :sub:`2`
     - k21
     - `Cohen & Westberg (1983) <https://ui.adsabs.harvard.edu/abs/1983JPCRD..12..531C/abstract>`__
     - :c:data:`primordial_chemistry` >1
   * - H\ :sup:`-` + :math:`{\gamma}`
     - H + e\ :sup:`-`
     - k27
     - none
     - :c:data:`primordial_chemistry` >1
   * - H\ :sub:`2`:sup:`+` + :math:`{\gamma}`
     - H + H\ :sup:`+`
     - k28
     - none
     - :c:data:`primordial_chemistry` >1
   * - H\ :sub:`2` + :math:`{\gamma}`
     - H\ :sub:`2`:sup:`+` + e\ :sup:`-`
     - k29
     - none
     - :c:data:`primordial_chemistry` >1
   * - H\ :sub:`2`:sup:`+` + :math:`{\gamma}`
     - H\ :sup:`+` +  H\ :sup:`+` + e\ :sup:`-`
     - k30
     - none
     - :c:data:`primordial_chemistry` >1
   * - H\ :sub:`2` + :math:`{\gamma}`
     - H + H
     - k31
     - none
     - :c:data:`primordial_chemistry` >1
   * - H + H + grain
     - H\ :sub:`2` + grain
     - k2dust
     - `Tielens & Hollenbach (1985) <https://ui.adsabs.harvard.edu/abs/1985ApJ...291..722T/abstract>`__
     - :c:data:`primordial_chemistry` >1 & :c:data:`h2_on_dust` ==1
   * - H\ :sup:`+` + D
     - H + D\ :sup:`+`
     - k50
     - `Savin (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...566..599S/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - D\ :sup:`+` + H
     - D + H\ :sup:`+`
     - k51
     - `Savin (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...566..599S/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - H\ :sub:`2` + D\ :sup:`+`
     - HD + H\ :sup:`+`
     - k52
     - `Galli & Palla (2002) <https://ui.adsabs.harvard.edu/abs/2002P%26SS...50.1197G/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - HD + H\ :sup:`+`
     - H\ :sub:`2` + D\ :sup:`+`
     - k53
     - `Galli & Palla (2002) <https://ui.adsabs.harvard.edu/abs/2002P%26SS...50.1197G/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - H\ :sub:`2` + D
     - HD + H
     - k54
     - `Clark et al. (2011) <https://ui.adsabs.harvard.edu/abs/2011ApJ...727..110C/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - HD + H
     - H\ :sub:`2` + D
     - k55
     - `Galli & Palla (2002) <https://ui.adsabs.harvard.edu/abs/2002P%26SS...50.1197G/abstract>`__  & `Ripamonti (2007) <https://ui.adsabs.harvard.edu/abs/2007MNRAS.376..709R/abstract>`__
     - :c:data:`primordial_chemistry` >2
   * - D + H\ :sup:`-`
     - HD + e\ :sup:`-`
     - k56
     - `Kreckel et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010Sci...329...69K/abstract>`__
     - :c:data:`primordial_chemistry` >2


.. note:: 
   For equations with :math:`{\gamma}` on the LHS, the rate for this equation come sfrom the choice of UV background models. See :c:data:`UVBackground` for more information.

This table below maps the chemical species used above to the relevent Grackle field pointers.

================
Chemical Species
================
====================== ========================
variable                Reaction Network
====================== ========================
H                      :c:data:`HI_density`
H\ :sup:`+`            :c:data:`HII_density`
H\ :sup:`-`            :c:data:`HM_density`
H\ :sub:`2`            :c:data:`H2I_density`
H\ :sub:`2`:sup:`+`    :c:data:`H2II_density`
He                     :c:data:`HeI_density`
He\ :sup:`+`           :c:data:`HeII_density`
He\ :sup:`+`\ :sup:`+` :c:data:`HeIII_density`
e\ :sup:`-`            :c:data:`e_density`
D                      :c:data:`DI_density`
D\ :sup:`+`            :c:data:`DII_density`
HD                     :c:data:`HDI_density`
====================== ========================

=============
K22 Citations
=============

============================================== =======================================================================================================
k22 :c:data:`three-body-rate<three_body_rate>` Citations
============================================== =======================================================================================================
0                                               `Abel, Bryan & Norman (2002) <https://ui.adsabs.harvard.edu/abs/2002Sci...295...93A/abstract>`_
1                                               `Palla, Salpeter & Stahler (1983) <https://ui.adsabs.harvard.edu/abs/1983ApJ...271..632P/abstract>`_
2                                               `Cohen & Westberg (1983) <https://ui.adsabs.harvard.edu/abs/1983JPCRD..12..531C/abstract>`_
3                                               `Flower & Harris (2007) <https://ui.adsabs.harvard.edu/abs/2007MNRAS.377..705F/abstract>`_
4                                               `Glover (2008) <https://ui.adsabs.harvard.edu/abs/2008IAUS..255....3G/abstract>`_   
5                                               `Forrey (2013) <https://ui.adsabs.harvard.edu/abs/2013ApJ...773L..25F/abstract>`_    
============================================== =======================================================================================================
