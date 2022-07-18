import Functions
from Functions import *
outfile = ROOT.TFile("out_test.root","recreate")
outtree = ROOT.TTree("tree","tree")
c_short_p=POINTER(c_short)
event=np.zeros((nbin,nbin,max_depth))
#event=array(c_short,[nbin][nbin][max_depth])
print(type(event[0][0][0]))
#event=event.astype(np.short)
#event=event.ctypes.data_as(c_short_p)
outtree.Branch("event",event,"event[{}][{}][{}]/I".format(nbin,nbin,max_depth))
TPCEventTag=c_int(0);
print(type(TPCEventTag))
outtree.Branch("TPCEventTag",addressof(TPCEventTag),"TPCEventTag/I")
for i in range(0,10):
    event=np.zeros((nbin,nbin,max_depth),dtype=np.int16)
#    event=event.ctypes.data_as(c_short_p)
#    print(type(event[0][0][0]))
    TPCEventTag=0
#    event[i][i][0]=c_short(-202);
    print(i)
    event[int(i)][int(i)][0]=np.int16(-202);
    event = event.ctypes.data_as(c_short_p)
    TPCEventTag=c_int(i+5)
#    outtree.TPCEventTag=c_int(TPCEventTag)
#    event.ctypes.data_as(POINTER(c_short_p))
#    outtree.event=event
    outtree.Fill()
outtree.Write()
outfile.Close()
