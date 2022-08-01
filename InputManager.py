import Functions
from Functions import *
def LoadEventRDF(filename,scale):
    file = ROOT.TFile.Open(filename,"READ")
    tree = file.Get("tree")
    TotEnt=int(tree.GetEntries()/scale)
    df = ROOT.RDataFrame("tree",filename);
    npdf = df.AsNumpy(columns=["evnum","TPCEventTag","nhtpc","ntrk","x","y","z","dedx"]);
    Revnum = npdf.get("evnum");
    Rntrk = npdf.get("ntrk");
    RTPCEventTag = npdf.get("TPCEventTag");
    Rnhtpc = npdf.get("nhtpc");
    Rx = npdf.get("x");
    Ry = npdf.get("y");
    Rz = npdf.get("z");
    Rdedx = npdf.get("dedx");
    evnums = np.ndarray(shape=(TotEnt),dtype="int")
    TPCNTrk= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEventTags= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
    for i in range(TotEnt):
        evnums[i] = Revnum[i] 
        TPCEventTags[i]=RTPCEventTag[i]
        tpcntrk = Rntrk[i] 
        TPCNTrk[i] = CutMultiTrack(tpcntrk,output_num)
        nhit = Rnhtpc[i]
        if i%1000 ==0:
            print(i)
        for nh in range(0,nhit):
            x=Rx[i][nh]
            y=Ry[i][nh]
            z=Rz[i][nh]
            if np.isnan(x) or np.isnan(y) or np.isnan(z):
                print("NaN detected ",nh," : ",x,y,z)
                continue
            x,z=ToPixel(x,z)
            y=ToInt(y)
            dedx=Rdedx[i][nh]
            TPCEvents[i,x,z,0]=y/650
            TPCEvents[i,x,z,1]=dedx/dedx_norm
    return  TPCEvents,TPCEventTags,TPCNTrk,evnums,TotEnt 
def LoadEventXYZ(filename,scale):
    file = ROOT.TFile.Open(filename,"READ")
    tree = file.Get("tree")
    TotEnt=int(tree.GetEntries()/scale)
    df = ROOT.RDataFrame("tree",filename);
    npdf = df.AsNumpy(columns=["evnum","TPCEventTag","nhtpc","ntrk","x","y","z","dedx"]);
    Revnum = npdf.get("evnum");
    Rntrk = npdf.get("ntrk");
    RTPCEventTag = npdf.get("TPCEventTag");
    Rnhtpc = npdf.get("nhtpc");
    Rx = npdf.get("x");
    Ry = npdf.get("y");
    Rz = npdf.get("z");
    Rdedx = npdf.get("dedx");
    evnums = np.ndarray(shape=(TotEnt),dtype="int")
    TPCNTrk= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEventTags= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEventsXZ = np.zeros((TotEnt,nbin,nbin,max_depth))
    TPCEventsYZ = np.zeros((TotEnt,nbin,nbin,max_depth))
    TPCEventsXY = np.zeros((TotEnt,nbin,nbin,max_depth))
    for i in range(TotEnt):
        evnums[i] = Revnum[i] 
        TPCEventTags[i]=RTPCEventTag[i]
        tpcntrk = Rntrk[i] 
        TPCNTrk[i] = CutMultiTrack(tpcntrk,output_num)
        nhit = Rnhtpc[i]
        if i%1000 ==0:
            print(i)
        for nh in range(0,nhit):
            x=Rx[i][nh]
            y=Ry[i][nh]
            z=Rz[i][nh]
            if np.isnan(x) or np.isnan(y) or np.isnan(z):
                print("NaN detected ",nh," : ",x,y,z)
                continue
            x,z=ToPixel(x,z)
            y=ToInt(y)
            dedx=Rdedx[i][nh]/dedx_norm
            dedx=Cut(dedx,dedx_norm*3)
            TPCEventsXZ[i,x,z,0]=y
            TPCEventsYZ[i,y,z,0]=x
            TPCEventsXY[i,x,y,0]=z
            if max_depth>1:
               TPCEventsXZ[i,x,z,1]=dedx
               TPCEventsYZ[i,y,z,1]=dedx
               TPCEventsXY[i,x,y,1]=dedx
  #          TPCEventsXZ[i,x,z,0]=dedx
  #          TPCEventsYZ[i,y,z,0]=dedx
    return  TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCEventTags,TPCNTrk,evnums,TotEnt 
def LoadRealEventRDF(filename,scale):
    file = ROOT.TFile.Open(filename,"READ")
    tree = file.Get("tree")
    TotEnt=int(tree.GetEntries()/scale)
    df = ROOT.RDataFrame("tree",filename);
    npdf = df.AsNumpy(columns=["evnum","nhtpc","x","y","z","dedx"]);
    Revnum = npdf.get("evnum");
    Rnhtpc = npdf.get("nhtpc");
    Rx = npdf.get("x");
    Ry = npdf.get("y");
    Rz = npdf.get("z");
    Rdedx = npdf.get("dedx");
    evnums = np.ndarray(shape=(TotEnt),dtype="int")
    TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
    for i in range(TotEnt*frag,TotEnt*(frag+1)):
        if i%1000 ==0:
            print(i)
        evnums[i] = Revnum[i] 
        nhit = Rnhtpc[i]
        for nh in range(0,nhit):
            x=Rx[i][nh]
            y=Ry[i][nh]
            z=Rz[i][nh]
            x,z=ToPixel(x,z)
            y=ToInt(y)
            dedx=Rdedx[i][nh]
            TPCEvents[i,x,z,0]=y/650
            TPCEvents[i,x,z,1]=dedx/dedx_norm
    return  TPCEvents,evnums,TotEnt 
def LoadRealEventXYZ(filename,scale,frag):
    file = ROOT.TFile.Open(filename,"READ")
    tree = file.Get("tree")
    TotEnt=int(tree.GetEntries()/scale)
    df = ROOT.RDataFrame("tree",filename);
#    npdf = df.AsNumpy(columns=["evnum","TPCEventTag","nhtpc","ntrk","x","y","z","dedx"]);
    npdf = df.AsNumpy(columns=["evnum","nhtpc","x","y","z","dedx"]);
    Revnum = npdf.get("evnum");
#    Rntrk = npdf.get("ntrk");
#    RTPCEventTag = npdf.get("TPCEventTag");
    Rnhtpc = npdf.get("nhtpc");
    Rx = npdf.get("x");
    Ry = npdf.get("y");
    Rz = npdf.get("z");
    Rdedx = npdf.get("dedx");
    evnums = np.ndarray(shape=(TotEnt),dtype="int")
    TPCNTrk= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEventTags= np.ndarray(shape=(TotEnt),dtype="int")
    TPCEventsXZ = np.zeros((TotEnt,nbin,nbin,max_depth))
    TPCEventsYZ = np.zeros((TotEnt,nbin,nbin,max_depth))
    TPCEventsXY = np.zeros((TotEnt,nbin,nbin,max_depth))
    for i in range(TotEnt):
        event = i+TotEnt*frag;
        evnums[i] = Revnum[event] 
        nhit = Rnhtpc[event]
        if i%1000 ==0:
            print(event)
        for nh in range(0,nhit):
            x=Rx[event][nh]
            y=Ry[event][nh]
            z=Rz[event][nh]
            if np.isnan(x) or np.isnan(y) or np.isnan(z):
                print("NaN detected ",nh," : ",x,y,z)
                continue
            x,z=ToPixel(x,z)
            y=ToInt(y)
            dedx=Rdedx[event][nh]/dedx_norm
            dedx=Cut(dedx,dedx_norm*3)
            TPCEventsXZ[i,x,z,0]=y
            TPCEventsYZ[i,y,z,0]=x
            TPCEventsXY[i,x,y,0]=z
            if max_depth>1:
               TPCEventsXZ[i,x,z,1]=dedx
               TPCEventsYZ[i,y,z,1]=dedx
               TPCEventsXY[i,x,y,1]=dedx
  #          TPCEventsXZ[i,x,z,0]=dedx
  #          TPCEventsYZ[i,y,z,0]=dedx
    return  TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCEventTags,TPCNTrk,evnums,TotEnt 
