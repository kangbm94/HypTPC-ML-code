import Functions
from Functions import *
def LoadEvent(tree,scale):
    TotEnt=int(tree.GetEntries()/scale)
    TPCEventTags = np.zeros(TotEnt)
    TPCNTrk = np.zeros(TotEnt)
    TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
    for i in range(TotEnt):
        tree.GetEntry(i)#We assume already mixed data
        if i%1000 ==0:
            print(i)
        ent=tree
        TPCEventTags[i]=ent.TPCEventTag
        tpcntrk = ent.ntrk 
        TPCNTrk[i] = CutMultiTrack(tpcntrk,output_num)
        nhit = ent.nhtpc
        for nh in range(0,nhit):
            x=ent.x[nh]
            y=ent.y[nh]
            z=ent.z[nh]
            x,z=ToPixel(x,z)
            y=ToInt(y)
            dedx=ent.dedx[nh]
            TPCEvents[i,x,z,0]=y/650
    return  TPCEvents,TPCEventTags,TPCNTrk 
