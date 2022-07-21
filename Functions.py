import ROOT
import numpy as np
import matplotlib.pyplot as plt
import Physics
from Physics import *
import ctypes
from ctypes import *
from array import array
nbin = int(128)
tpc_size = 250
MaxNtr = 20
max_depth=int(1)
def ToPixel(x,z):
    x+=250
    z+=250
    x_pix = int(x* nbin/tpc_size/2)
    z_pix = int(z* nbin/tpc_size/2)
    return x_pix,z_pix
def ToInt(y):
    y+=350
    y*=10
    return int(y)

def EventTag(tree):
    nk=0
    npi=0
    nprt=0
    particles=np.zeros(MaxNtr)
    nhittpc=tree.nhittpc
    ntrk = tree.ntrk
    idtpc=tree.idtpc
    for nh in range(0,nhittpc):
        particles[ntrk[nh]]=idtpc[nh]
    for particle in particles:
        if(abs(particle)==PionID):
            npi+=1
        if(abs(particle)==KaonID):
            nk+=1
        if(abs(particle)==ProtonID):
            nprt+=1
#    print("(npi,nk,np) is " ,npi,nk,nprt)
    if(npi==2 and nk==2 and nprt==1):
        return 1
    if(npi==1 and nk==2 and nprt==0):
        return 2
    if(npi==0 and nk==1 and nprt==0):
        return 3
    if(npi==4 and nk==0 and nprt==0):
        return 4
    else:    
        return 0

