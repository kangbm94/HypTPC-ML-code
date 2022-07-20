import NeuralModel
from NeuralModel import *
import pandas as pd
from sklearn.utils import shuffle

print("Starting...")
filename = "TrainDataTagged3k.root"
#filename = "TrainDataTagged_wo_bg.root"
model = NeuralModel.model
model.load_weights("./Model_3")
cnt=0
print("Loading Data...")
file = ROOT.TFile.Open(filename,"READ")
tree = file.Get("tree")
scale = int(1)
TotEnt=int(tree.GetEntries()/scale)
TPCEventTags = np.zeros(TotEnt)
TPCEvnum = np.zeros(TotEnt)
TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
print("hello")
for i in range(0,TotEnt):
    tree.GetEntry(i)
#    tree.GetEntry(int(i*scale)+1)
    if i%1000 ==0:
        print(i)
    ent=tree
    TPCEvnum[i]=ent.evnum
    TPCEventTags[i]=ent.TPCEventTag
    nhit = ent.ntrk

    for nh in range(0,nhit):
        x=ent.x[nh]
        y=ent.y[nh]
        z=ent.z[nh]
        for dep in range(0,max_depth):
            if TPCEvents[i,x,z,dep]==0:
                TPCEvents[i,x,z,0]=y/650
                break


##########################Data Prepared###################
#cb=callbacks.ModelCheckpoint(filepath='./checkpoint',save_weights_only=True,save_freq=4)



TrainData=TPCEvents
predTrain=model.predict(TrainData)
print(predTrain.shape)

'''
	Else=0,
	L2PPi=1,
	L2NPi=2,
	KBeam=3
'''

fig = plt.figure()
fig,((a1t,a2t),(a3t,a4t)) = plt.subplots(2,2)
a1t.hist(predTrain[:,0],bins=20)
a1t.set_title("Background")
a2t.hist(predTrain[:,1],bins=20)
a2t.set_title("Lambda->P Pi")
a3t.hist(predTrain[:,2],bins=20)
a3t.set_title("Lambda->N Pi")
a4t.hist(predTrain[:,3],bins=20)
a4t.set_title("K Beam")

outfilename = "PredictedData.root"
outfile = ROOT.TFile.Open(outfilename,"recreate")
outtree = ROOT.TTree("tree","tree")

evnum = array('i',[0])
Background = array('d',[0])
L2PPi = array('d',[0])
L2NPi = array('d',[0])
KBeam = array('d',[0])

outtree.Branch("evnum",evnum,"evnum/I")
outtree.Branch("Background",Background,"Background/D")
outtree.Branch("L2PPi",L2PPi,"L2PPi/D")
outtree.Branch("L2NPi",L2NPi,"L2NPi/D")
outtree.Branch("KBeam",KBeam,"KBeam/D")
print("hel")
print(type(TPCEvnum[0]))
print(type(predTrain[0,0]))
for i in range(0,TotEnt):
    tree.GetEntry(i)
    evnum[0]=int(TPCEvnum[i])
    Background[0]=float(predTrain[i,0])
    L2PPi[0]=float(predTrain[i,1])
    L2NPi[0]=float(predTrain[i,2])
    KBeam[0]=float(predTrain[i,3])
    outtree.Fill()
outfile.Write()
outfile.Close()

plt.show()
'''
plt.imshow(TPCEvents[2500,:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
'''
