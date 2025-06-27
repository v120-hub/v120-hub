# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 14:30:35 2025

@author: vincent.longchamp
"""

import numpy as np
import random
import math
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#%% FUNCTIONS

    
# Random generation of pore (position and radius) 
def random_pore_coordinate(Lx, Ly):
    a = Ra
    b = Rb
    rotation = angle + random.randint(-1*delta_angle,delta_angle)
    if allow_open_porosities :
        return [random.uniform(0,Lx) ,  random.uniform(0,Ly) , a,b, rotation ]
    else :
        return [random.uniform(a+minimal_gap,Lx-a-minimal_gap),  random.uniform(a+minimal_gap,Ly-a-minimal_gap) , a, b, rotation]

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
        if d < pore[2]+porosities[i][2]+minimal_gap :
            ok = False
        i+=1
    return ok
    
# Ccompute the volume void fraction
def compute_vf(Lx,Ly,porosities):
    A0 = Lx*Ly
    Ap = 0
    for pore in porosities :
        Ap+=math.pi*pore[2]*pore[3]
    return (Ap/A0)
    
#%% GENERATION
Lx = 1          # Length along X
Ly = 1          # Length along Y
Ra = 0.04       # Ellipse major axis
Rb = Ra*0.3     # Ellipse minor axis
angle = 0      # Mean angle
delta_angle = 20 # Angle variation around mean value (uniform distribution)
minimal_gap = 0.1*Ra 

void_fraction = 0.115 # Target void_fraction
allow_open_porosities = False # Allow pores to intersect a boundary (edge of the domain)

porosities = [] # List of generated porosities
current_vf = 0  # current void fraction
it = 0          # iteration number to stop the loop
while current_vf < void_fraction and it < 10000 :
    
    pore = random_pore_coordinate(Lx, Ly)
    ok = check_intersection(pore,porosities)
    if ok :
        porosities.append(pore)
    current_vf = compute_vf(Lx,Ly,porosities)
    it+=1
print(it)


    
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
    elli = ptc.Ellipse((pore[0],pore[1]), pore[2],pore[3], angle=pore[4], edgecolor="black",lw=1,alpha=0.7)
    ax.add_patch(elli)
