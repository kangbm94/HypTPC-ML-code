import ROOT
import numpy as np
import matplotlib.pyplot as plt
import Physics
from Physics import *
from array import array
import sys
output_num = 9 
nbin = int(128)
tpc_size = 250
dl = 665
dl_offset=300
MaxNtr = 20
max_depth=int(1)
dedx_norm = 0.001


def ToPixel(x,z):
    x+=250
    z+=250
    x_pix = int(x* nbin/tpc_size/2)
    z_pix = int(z* nbin/tpc_size/2)
    if x_pix>nbin-1 or z_pix>nbin-1:
        print(x_pix,z_pix)
    return x_pix,z_pix
def ToInt(y):
    y+=dl_offset
    y_pix = int(y*nbin/dl)
    return y_pix

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
def CutMultiTrack(ntrk, cut):
    if(ntrk>cut-2):
        return cut-1
    else:
        return ntrk
def MaxCh(tags):
    a0 = 0
    cnt=0
    ch=0
    for a in tags:
        cnt+=1
        if(a>a0):
            a0 = a
            ch+=cnt
            cnt=0
    return ch-1
def Cut(ent,cut):
    if ent<cut:
        return ent
    else:
        return cut
def OneHotEncoding(a):
    label = np.zeros((a.size,output_num))
    for i in range(a.size):
        label[i][a[i]]=1
    return label
'''
def Sort(tags,cont,multi):
    lt = range(multi)
    ll = np.array([np.array(list) for _ in lt])
    listlength = tags.size
    for i in range(listlength):
        val = tags[i]
        for j in range(multi):
            ll[val].
    return list_ex
       ''' 
