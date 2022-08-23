import Functions
from Functions import *
filename = "beamthroughG4_pi.root"
filename="kpxiG4_EM_HAD.root"
filename = "beamG4_EM_HAD.root"
filename = "unipiMTG4_EM_HAD.root"
filename = "Train_raw.root"
file = ROOT.TFile.Open(filename,"READ")
tree = file.Get("TPC_g")
TotEnt = tree.GetEntries()
print(TotEnt)
#TPCEvent.reshape((250,250))
cnt=0
#TPCEvents = np.zeros((TotEnt,nbin,nbin,5))
TPCEventTags = np.zeros(TotEnt)
TPCEventTag = 0
TPCEvent = np.zeros((nbin,nbin,max_depth),dtype=np.short)
outfile = ROOT.TFile("TrainData.root","recreate")
outtree = ROOT.TTree("tree","tree")
outtree.Branch("TPCEvent",TPCEvent,"TPCEvent[{}][{}][{}]/S".format(nbin,nbin,max_depth))
outtree.Branch("TPCEventTag",TPCEventTag,"TPCEventTag/I")
for ent in tree:
    if cnt%100 == 0:
        print("Event Process: ", cnt)
    if cnt>-100000:
        nhittpc=ent.nhittpc
        TPCEvent = np.zeros((nbin,nbin,max_depth),dtype=np.short)
        TPCEventTag=c_int(EventTag(ent))
        for nh in range(0,nhittpc):
            x=ent.xtpc[nh]
            y=ent.ytpc[nh]
            z=ent.ztpc[nh]
            xp,zp = ToPixel(x,z)
          #  y/=600
        #    print(xp,zp)
            for depth in (0,2):
                if(TPCEvent[xp,zp,depth]==0):
                    TPCEvent[xp,zp,depth] = ToInt(y)
                    break
             #   else:
                    #print(xp,zp,"Double Entries")
#        if TPCEventTags[cnt]==1 :
#            print("Event: ",cnt)
#            break
    outtree.Fill()
    cnt+=1
    #evt loop ends.
#    TPCEvents[cnt,:,:,:]=TPCEvent;
outtree.Write()
outfile.Close()
print("Drawing...")
plt.imshow(TPCEvent[:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
plt.show()
        
