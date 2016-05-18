#! /usr/bin/env python
import matplotlib                as mpl
mpl.use("cairo")
import matplotlib.pyplot         as plt
import numpy                     as np
import sys
sys.path.insert(0, '/alphadata01/pfister/scripts/Bern3Ddiagtools/python')

def GillPressureRho(S,t,rho0,p):
    # given a density at surface pressure p=0 (rho0), calculate the density at pressure p
    
    Kw = 19652.21 + 148.4206 *t - 2.327105 *t**2 + 1.360477e-2 *t**3 - 5.155288e-5 *t**4
    K0 = Kw + S*(54.6746 - 0.603459 *t + 1.09987e-2 *t**2 - 6.1670e-5 *t**3)\
        + S**(3./2.)*(7.944e-2 + 1.6483e-2 *t - 5.3009e-4 *t**2)
    K  = K0 + p* (3.239908 + 1.43713e-3 *t + 1.16092e-4 *t**2 - 5.77905e-7 *t**3)\
        + p*S*   (2.2838e-3 - 1.0981e-5 *t -1.6078e-6 *t**2) + 1.91075e-4 *p*S**(3./2.)\
        + p**2*  (8.50935e-5 - 6.12293e-6 *t + 5.2787e-8 *t**2)\
        + p**2*S*(-9.9348e-7 + 2.0816e-8 *t + 9.1697e-10 *t**2)
    
    rho = rho0 / (1 - p/K)
    return rho

def GillRho(S,t,p,returnrho0=False):
    # function for calculating sea water density [kg m^-3] as function of:
    # pressure p [bar], temperature t [deg. Celsius], salinity [psu]
    # From Adrian E. Gill, Atmosphere and Ocean Dynamics, Appendix 3: Properties of Seawater

    rhow = 999.842594 + 6.793952e-2 *t - 9.095290e-3 *t**2 + 1.001685e-4 *t**3\
        - 1.120083e-6 *t**4 + 6.536332e-9 *t**5

    rho0 = rhow + S*(0.824493 - 4.0899e-3 *t + 7.6438e-5 *t**2 - 8.2467e-7 *t**3 + 5.3875e-9 *t**4)\
        + S**(3./2.)* (-5.72466e-3 + 1.0227e-4 *t - 1.6546e-6 *t**2) + 4.8314e-4 *S**2

    if returnrho0:
        return rho0
    else:
        rho = GillPressureRho(S,t,rho0,p)
        return rho

def TestRho():
    print "The below rho values should be near-identical: "
    print GillRho(0.,5.,0.), 999.96675
    print GillRho(35.,5.,0.), 1027.67547
    print GillRho(35.,25.,1000), 1062.53817

def TestRho_2():
    # function to show that rho_p decreases more rapidly than rho with increasing T 
    print "For Delta T = +2K: "
    print "Delta rho without p-correction: " 
    print GillRho(35.,6.,400,returnrho0=True)-GillRho(35.,4.,400.,returnrho0=True)
    print "Delta rho with p-correction: "
    print GillRho(35.,6.,400)-GillRho(35.,4.,400.)
    print "Therefore, in B3D diagnosis OHU(~rho) is slower and SSLR(~Delta rho) is faster with p-correction!"

def RefineRho(S, t, rho, z_w, z_t, epsilon = 0.001):
    # given S, T, rho (which is assumed to be only rho(S,T) ) and a vertical grid, 
    # this function calculates rho(S,T,p) according to function GillRho
    # 1. A hydrostatic pressure p is calculated from rho, 2. a new rho is calculated from p
    # 1. and 2. are iterated until the mean rho change in the deepest cell is < epsilon
    # (i.e., rho and p are in agreement)

    p = np.ma.masked_where(rho.mask, np.zeros(rho.shape))
    g = 9.81
    rhoz_old = rho.mean(axis=1).mean(axis=1)
    rhoz = rhoz_old - 2 * epsilon
    niterations = 0

    while (abs(rhoz_old[-1] - rhoz[-1]) > epsilon):
        p[~p.mask] = 0
        for i in np.arange(p.shape[2]):
            for j in np.arange(p.shape[1]):
                for k in np.arange(p.shape[0]):
                    if not rho.mask[k, j, i]:
                        if k == 0:
                            p[k, j, i] = (z_t[k] - z_w[k]) * rho[k, j, i] * g
                        else:
                            p[k, j, i] = p[k - 1, j, i] + ((z_w[k] - z_t[k - 1]) * rho[k - 1, j, i] + (z_t[k] - z_w[k]) * rho[k, j, i]) * g
                    else:
                        break

        p = p * 1e-05
        if niterations > 0:
            rhoz_old = rhoz
        rho = GillRho(S, t, p)
        rhoz = rho.mean(axis=1).mean(axis=1)
        niterations += 1

    return rho

# The following _simple functions are simplifications of their normal counterparts, that were used in the 
# Bern3D model for computational efficiency. Simplifications made:
# RefineRho: Only one iteration, rho=const=rho_ref is assumed in pressure integral
# GillPressureRho: Only terms linear in t, s and p are considered
# rho_B3D: much simpler incompressible equation of state, as used in B3D model
# with these simplifications, rho is very similar to full calculation (SSLR_Gill_simplified.py)

def RefineRho_simple(S, t, rho, rho_ref, z_w, z_t):
    p = np.ma.masked_where(rho.mask, np.zeros(rho.shape))
    g = 9.81

    p[~p.mask] = 0
    for i in np.arange(p.shape[2]):
        for j in np.arange(p.shape[1]):
            for k in np.arange(p.shape[0]):
                if not rho.mask[k, j, i]:
                    if k == 0:
                        p[k, j, i] = (z_t[k] - z_w[k]) * rho_ref * g
                    else:
                        p[k, j, i] = p[k-1,j,i] + (z_t[k]-z_t[k-1])*rho_ref*g
                else:
                    break

    p = p * 1e-05
    rho_new = rho_B3D(S, t, p)

    return rho_new

def GillPressureRho_simple(S,t,rho0,p):    
    K = 19652.21 + 148.4206 *t + 54.6746 *S + 3.239908 *p
    rho = rho0 / (1 - p/K)
    return rho

def rho_B3D(S,t,p):
    usc = 0.05 # velocity scale
    rsc = 6.37e6 # radius of the Earth             
    dsc = 5e3 # depth scale: 5000m
    fsc = 2*7.2921e-5 # Coriolis scale 2*Omega
    gsc = 9.81 # gravity acceleration (m/s2)
    rh0sc = 1e3 # density shift scale (1000 kg/m3)
    rhosc = rh0sc*fsc*usc*rsc/gsc/dsc # density scale
    saln0 = 34.78

    ec1 = - 0.0559 /rhosc
    ec2 = 0.7968   /rhosc
    ec3 = - 0.0063 /rhosc
    ec4 = 3.7315e-5/rhosc

    rho_mod = ec1*t + ec2*(S-saln0) + ec3*t**2 + ec4*t**3 # rho in model units
    rho0 = rho_mod*rhosc + 0.7968*saln0 + rh0sc # rho in SI units
    rho = GillPressureRho_simple(S,t,rho0,p) # include pressure dependence
    return rho

def SSLRcalc(rho0, rho, area, deltaz, totarea):
    # calculates global mean steric sea level change between rho0[k,j,i] (e.g. preindustrial) 
    # and rho[k,j,i] (e.g. future)
    # also needs the surface area (area[j,i]) and height (deltaz[k]) of all model cells
    SSLR = 0.0
    for i in np.arange(rho0.shape[2]):
        for j in np.arange(rho0.shape[1]):
            for k in np.arange(rho0.shape[0]):
                if not rho0.mask[k, j, i]:
                    SSLR += deltaz[k] * (rho0[k, j, i] / rho[k, j, i] - 1.0) * area[j, i]
                else:
                    break

    SSLR = SSLR / totarea
    return SSLR
