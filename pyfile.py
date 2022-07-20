import Functions
from Functions import *
outfile = ROOT.TFile("out_test.root","recreate")
outtree = ROOT.TTree("tree","tree")
c_short_p=POINTER(c_short)
nbins=c_int(nbin)
#event=array(c_short,[nbin][nbin][max_depth])
#event=event.astype(np.short)
#event=event.ctypes.data_as(c_short_p)
TPCEventTag=c_int(0);
print(type(TPCEventTag))
event = array('f',[0]*nbin)
a = array('d',[0])
b = array('d',[0])
c = array('i',[0])
d = array('i',[0])
outtree.Branch("a",a,"a/D")
outtree.Branch("b",b,"b/D")
outtree.Branch("c",c,"c/I")
outtree.Branch("d",d,"d/I")
for i in range(0,10):
    '''
    a=c_double(i)
    b=c_double(i+10)
    c=c_double(i+20)
    d=c_double(i+30)
    '''
    
    ad=i+0.
    bd=i+10.
    cd=i+20.
    dd=i+30.
    a[0]=float(ad)
    b[0]=float(bd)
    c[0]=int(cd)
    d[0]=int(dd)
    outtree.Fill()
outtree.Write()
outfile.Close()
