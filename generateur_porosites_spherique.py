# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 14:30:35 2025

@author: vincent.longchamp
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% FUNCTIONS

# Return random value, gaussian distribution centered on R, with 2*R*delta_R the width at half maximum
def random_radius_gauss(R,delta_R):
    r=-1
    r_min=0
    c = 2*R*delta_R/(2*np.power(2*np.log(2),1/2))
    while r < r_min :
        r = random.gauss(R, c)
        return r
    
# Random generation of pore (position and radius) 
def random_pore_coordinate(Lx, Ly, R, delta_R):
    if delta_R == 0 :
        r = R 
    elif random_radius == "gauss" :
        r = random_radius_gauss(R,delta_R)
    elif random_radius == "uni" :
        r = random.uniform(R*(1-delta_R), R*(1+delta_R))
    else :
        r = R 
    if allow_open_porosities :
        return [random.uniform(0,Lx) ,  random.uniform(0,Ly) , r]
    else :
        return [random.uniform(r+minimal_gap,Lx-r-minimal_gap) ,  random.uniform(r+minimal_gap,Ly-r-minimal_gap) , r]

# Compute the distance between two points of the plane
def compute_distance(x1,y1,x2,y2):
    return math.sqrt(pow(x1-x2,2)+pow(y1-y2,2))

# Check if pore intersect another pore in porosities
def check_intersection(pore, porosities):
    ok = True
    i = 0
    imax = len(porosities)
    while ok and i<imax :
        d = compute_distance(pore[0],pore[1],porosities[i][0],porosities[i][1])
        if d < pore[-1]+porosities[i][-1]+minimal_gap :
            ok = False
        i+=1
    return ok
    
# Ccompute the volume void fraction
def compute_vf(Lx,Ly,porosities):
    A0 = Lx*Ly
    Ap = 0
    for pore in porosities :
        Ap+=math.pi*pore[-1]*pore[-1]
    return (Ap/A0)

# Gaussian law
def gaussienne(x, offset=0, fwhm=1, puissance=2) :
    c = fwhm/(2*np.power(2*np.log(2),1/puissance))
    return np.exp(-1*np.power(x-offset,puissance)/(2*np.power(c,puissance)))
    
#%% GENERATION
Lx = 1      # length along X
Ly = 1      # length along Y
R = 0.025    # Pore radius
minimal_gap = 0.5*R # Minimal gap between pore 

random_radius = "none" # Choice of random type, possible value : gauss, uni, none
delta_R = 0.1  # Variation of pore radius : +-R*delta_R, only if random_radius != "none" 

void_fraction = 0.115 # Target void_fraction
allow_open_porosities = True # Allow pores to intersect a boundary (edge of the domain)

porosities = [] # List of generated porosities
current_vf = 0  # current void fraction
it = 0          # iteration number to stop the loop
while current_vf < void_fraction and it < 10000 :
    pore = random_pore_coordinate(Lx, Ly, R, delta_R)
    ok = check_intersection(pore,porosities)
    if ok :
        porosities.append(pore)
    current_vf = compute_vf(Lx,Ly,porosities)
    it+=1
print(it)


# R = 0.01    # Pore radius
# delta_R = 0.1  # Variation of pore radius : +-R*delta_R 
# minimal_gap = 0.1*R # Minimal gap between pore 
# void_fraction_2 = 0.15
# it= 0
# while current_vf < void_fraction_2 and it < 10000 :
#     pore = random_pore_coordinate(Lx, Ly, R, delta_R)
#     ok = check_intersection(pore,porosities)
#     if ok :
#         porosities.append(pore)
#     current_vf = compute_vf(Lx,Ly,porosities)
#     it+=1
# print(it)

    
#%% PLOT

# Pores positions
fig,ax =plt.subplots(1,1,figsize=(8,6))
ax.set_title("Void fraction : %.2f"%(current_vf*100)+" %")
plt.gca().set_aspect(aspect = 'equal')
ax.set_xlim(0,Lx)
ax.set_ylim(0,Ly)
ax.set_xlabel("X coordinate [m]")
ax.set_xlabel("Y coordinate [m]")
for pore in porosities : 
    circ = plt.Circle((pore[0],pore[1]), pore[-1],edgecolor="black",lw=1,alpha=0.7)
    ax.add_patch(circ)


# Histogram of pore radius
if random_radius != 'none':
    fig2, axh = plt.subplots(1,1,figsize=(8,6))
    axh.set_xlabel("Pore radius [m]") 
    axh.set_ylabel("Quantity") 
    Radius = np.array(porosities)[:,-1] 
    histo=axh.hist(Radius,histtype='bar',rwidth=0.9,bins=10)
    axh.set_xlim(xmin=0, xmax=np.max(Radius)*1.2)
    axh.set_title("Total pore number : %i "%(len(Radius)))
    
    if random_radius == "gauss" :
        x = np.linspace(0, np.max(Radius)*1.2, 1000)
        y = gaussienne(x, R, 2*R*delta_R) * np.max(histo[0])
        axh.plot(x,y)
    elif random_radius == "uni" :
        axh.plot([0,R*(1-delta_R),R*(1-delta_R),R*(1+delta_R),R*(1+delta_R),np.max(Radius)*1.2], [0,0,np.max(histo[0]),np.max(histo[0]),0,0])