########################################################################
#
# Primordial chemistry and cooling equilibrium
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np

# Equilibrium abundances

def nHI(T, nH, rates='enzo'):
    return nH * alphaHII(T, rates=rates) / (alphaHII(T, rates=rates) +
                                            GammaeHI(T, rates=rates))

def nHII(T, nH, rates='enzo'):
    return nH - nHI(T, nH, rates=rates)

def nHeI(T, nH, Y=0.24, rates='enzo'):
    return (nHeII(T, nH, Y=Y, rates=rates) *
            (alphaHeII(T, rates=rates) + alphad(T, rates=rates)) /
            GammaeHeI(T, rates=rates))

def nHeII(T, nH, Y=0.24, rates='enzo'):
    y = Y / (4 - 4*Y)
    return y * nH / (1. + ((alphaHeII(T, rates=rates) +
                            alphad(T, rates=rates)) /
                           (GammaeHeI(T, rates=rates))) +
                     (GammaeHeII(T, rates=rates) /
                      alphaHeIII(T, rates=rates)))

def nHeIII(T, nH, Y=0.24, rates='enzo'):
    return nHeII(T, nH, Y=Y, rates=rates) * GammaeHeII(T, rates=rates) / \
      alphaHeIII(T, rates=rates)

def ne(T, nH, Y=0.24, rates='enzo'):
    return nHII(T, nH, rates=rates) + nHeII(T, nH, Y=Y, rates=rates) + \
      2 * nHeIII(T, nH, rates=rates)

### Ionization balance

# Recombination

def alphaHII(T, rates='enzo'):
    if rates == 'cen':
        return 8.4e-11 * np.power(T,-0.5) * np.power((T * 1e-3),-0.2) / \
            (1. + np.power((T * 1e-6),0.7))
    elif rates == 'enzo':
        T_eV = T/11605.e0
        log_T_eV = np.log(T_eV)
        alpha_rates = np.zeros_like(T)
        filter1 = T > 5500.e0
        alpha_rates[T > 5500.e0] = \
            np.exp(-28.61303380689232e0
                   - 0.7241125657826851e0*log_T_eV[filter1]
                   - 0.02026044731984691e0*log_T_eV[filter1]**2
                   - 0.002380861877349834e0*log_T_eV[filter1]**3
                   - 0.0003212605213188796e0*log_T_eV[filter1]**4
                   - 0.00001421502914054107e0*log_T_eV[filter1]**5
                   + 4.989108920299513e-6*log_T_eV[filter1]**6
                   + 5.755614137575758e-7*log_T_eV[filter1]**7
                   - 1.856767039775261e-8*log_T_eV[filter1]**8
                   - 3.071135243196595e-9*log_T_eV[filter1]**9)
        filter2 = T <= 5500.e0
        alpha_rates[filter2] = alphaHeII(T[filter2], rates=rates)
        return alpha_rates

def alphaHeII(T, rates='enzo'):
    if rates == 'cen':
        return 1.5e-10 * np.power(T,-0.6353)
    elif rates == 'enzo':
        T_eV = T/11605.e0
        return (1.54e-9*(1.e0+0.3e0/np.exp(8.099328789667e0/T_eV)) /
                (np.exp(40.49664394833662e0/T_eV)*T_eV**1.5e0) +
                3.92e-13/T_eV**0.6353e0)

def alphaHeIII(T, rates='enzo'):
    # same between cen and enzo
    return 3.36e-10 * np.power(T,-0.5) * np.power((T * 1e-3),-0.2) / \
      (1. + np.power((T * 1e-6),0.7))

# Dielectronic recombination

def alphad(T, rates='enzo'):
    if rates == 'cen':
        return 1.9e-3 * np.power(T,-1.5) * np.exp(-470000 / T) * \
            (1. + 0.3 * np.exp(-94000 / T))
    elif rates == 'enzo':
        return np.zeros_like(T)

# Collisional ionization

def GammaeHI(T, rates='enzo'):
    if rates == 'cen':
        return 5.85e-11 * np.power(T,0.5) * np.exp(-157809.1/T) / \
           (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        T_eV = T/11605.e0
        log_T_eV = np.log(T_eV)
        return np.exp(-32.71396786375e0
                      + 13.53655609057e0*log_T_eV
                      - 5.739328757388e0*log_T_eV**2
                      + 1.563154982022e0*log_T_eV**3
                      - 0.2877056004391e0*log_T_eV**4
                      + 0.03482559773736999e0*log_T_eV**5
                      - 0.00263197617559e0*log_T_eV**6
                      + 0.0001119543953861e0*log_T_eV**7
                      - 2.039149852002e-6*log_T_eV**8)

def GammaeHeI(T, rates='enzo'):
    if rates == 'cen':
        return 2.38e-11 * np.power(T,0.5) * np.exp(-285335.4/T) / \
            (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        T_eV = T/11605.e0
        log_T_eV = np.log(T_eV)
        return np.exp(-44.09864886561001e0
                      + 23.91596563469e0*log_T_eV
                      - 10.75323019821e0*log_T_eV**2
                      + 3.058038757198e0*log_T_eV**3
                      - 0.5685118909884001e0*log_T_eV**4
                      + 0.06795391233790001e0*log_T_eV**5
                      - 0.005009056101857001e0*log_T_eV**6
                      + 0.0002067236157507e0*log_T_eV**7
                      - 3.649161410833e-6*log_T_eV**8)

def GammaeHeII(T, rates='enzo'):
    if rates == 'cen':
        return 5.68e-12 * np.power(T,0.5) * np.exp(-631515.0/T) / \
            (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        T_eV = T/11605.e0
        log_T_eV = np.log(T_eV)
        return np.exp(-68.71040990212001e0
                      + 43.93347632635e0*log_T_eV
                      - 18.48066993568e0*log_T_eV**2
                      + 4.701626486759002e0*log_T_eV**3
                      - 0.7692466334492e0*log_T_eV**4
                      + 0.08113042097303e0*log_T_eV**5
                      - 0.005324020628287001e0*log_T_eV**6
                      + 0.0001975705312221e0*log_T_eV**7
                      - 3.165581065665e-6*log_T_eV**8)

### Cooling

# Collisional excitation cooling

def ceHI(T, nH, rates='enzo'):
    return 7.50e-19 * ne(T, nH, rates=rates) * nHI(T, nH, rates=rates) * \
      np.exp(-118348.0 / T) / (1. + np.power((T * 1e-5),0.5))

def ceHeII(T, nH, Y=0.24, rates='enzo'):
    return 5.54e-17 * ne(T, nH, rates=rates) * nHeII(T, nH, Y=Y, rates=rates) * \
      np.power(T,-0.397) * np.exp(-473638.0 / T) / (1. + np.power((T * 1e-5),0.5))

# Collisional ionization cooling

def ciHI(T, nH, rates='enzo'):
    if rates == 'cen':
        return 1.27e-21 * ne(T, nH, rates=rates) * nHI(T, nH, rates=rates) * \
          np.power(T,0.5) * np.exp(-157809.1 / T) / (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        return 2.18e-11 * GammaeHI(T, rates=rates) * \
          ne(T, nH, rates=rates) * nHI(T, nH, rates=rates)

def ciHeI(T, nH, rates='enzo'):
    if rates == 'cen':
        return 9.38e-22 * ne(T, nH, rates=rates) * nHeI(T, nH, rates=rates) * \
          np.power(T,0.5) * np.exp(-285335.4 / T) / (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        return 3.94e-11 * GammaeHeI(T, rates=rates) * \
        ne(T, nH, rates=rates) * nHeI(T, nH, rates=rates)

def ciHeII(T, nH, Y=0.24, rates='enzo'):
    if rates == 'cen':
        return 4.95e-22 * ne(T, nH, rates=rates) * nHeII(T, nH, Y=Y, rates=rates) * \
          np.power(T,0.5) * np.exp(-631515.0 / T) / (1. + np.power((T * 1e-5),0.5))
    elif rates == 'enzo':
        return 8.72e-11 * GammaeHeII(T, rates=rates) * \
          ne(T, nH, rates=rates) * nHeII(T, nH, Y=Y, rates=rates)

# Recombination cooling

def rHII(T, nH, rates='enzo'):
    return 8.70e-27 * ne(T, nH, rates=rates) * nHII(T, nH, rates=rates) * \
      np.power(T,0.5) * np.power((T * 1e-3),-0.2) / (1. + np.power((T * 1e-6),0.7))

def rHeII(T, nH, Y=0.24, rates='enzo'):
    return 1.55e-26 * ne(T, nH, rates=rates) * nHeII(T, nH, Y=Y, rates=rates) * \
      np.power(T,0.3647)

def rHeIII(T, nH, rates='enzo'):
    return 3.48e-26 * ne(T, nH, rates=rates) * nHeIII(T, nH, rates=rates) * np.power(T,0.5) * \
      np.power((T * 1e-3),-0.2) / (1. + np.power((T * 1e-6),0.7))

# Dielectronic recombination cooling

def drHeII(T, nH, Y=0.24, rates='enzo'):
    return 1.24e-13 * ne(T, nH, rates=rates) * nHeII(T, nH, Y=Y, rates=rates) * \
      np.power(T,-1.5) * np.exp(-470000.0 / T) * (1. + 0.3 * np.exp(-94000.0 / T))

# Free-free emission

def gff(T):
    return 1.1 + 0.34 * np.exp(-np.power((5.5 - np.log10(T)),2) / 3.0)

def freefree(T, nH, Y=0.24, rates='enzo'):
    return 1.42e-27 * gff(T) * np.power(T,0.5) * ne(T, nH, rates=rates) * \
      (nHII(T, nH, rates=rates) + nHeII(T, nH, Y=Y, rates=rates) + \
       4 * nHeIII(T, nH, rates=rates))

# Total cooling

def total_cooling(T, nH, rates='enzo'):
    return ceHI(T, nH, rates=rates) + ceHeII(T, nH, rates=rates) + \
      ciHI(T, nH, rates=rates) + ciHeI(T, nH, rates=rates) + \
      ciHeII(T, nH, rates=rates) + rHII(T, nH, rates=rates) + rHeII(T, nH, rates=rates) + \
      rHeIII(T, nH, rates=rates) + drHeII(T, nH, rates=rates) + freefree(T, nH, rates=rates)
