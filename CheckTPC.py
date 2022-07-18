import NeuralModel
from NeuralModel import *
import pandas as pd
from sklearn.utils import shuffle

print("Starting...")
filename = "TrainDataFairlyTagged.root"
#filename = "TrainDataTagged_wo_bg.root"
model = NeuralModel.model
model.load_weights("./Model_1")
cnt=0
print("Loading Data...")
file = ROOT.TFile.Open(filename,"READ")
tree = file.Get("tree")
scale = int(4)
TotEnt=int(tree.GetEntries()/scale)
TPCEventTags = np.zeros(TotEnt)
TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
for i in range(0,TotEnt-1):
    tree.GetEntry(int(i*scale)+1)
    if i%1000 ==0:
        print(i)
    ent=tree
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


fig = plt.figure()
fig,((a1t,a2t),(a3t,a4t)) = plt.subplots(2,2)
a1t.hist(predTrain[:,0],bins=20)
a1t.set_title("Training_Class1")
a2t.hist(predTrain[:,1],bins=20)
a2t.set_title("Training_Class2")
a3t.hist(predTrain[:,2],bins=20)
a3t.set_title("Training_Class3")
a4t.hist(predTrain[:,3],bins=20)
a4t.set_title("Training_Class4")

outfilename = "PredictedData.root"
outfile = ROOT.TFile.Open(outfilename,"recreate")
outtree = outfile.Get("tree")
TPCEventTag=c_int(0)
ntrk=c_int(0)
outtree.Branch("TPCEventTag",addressof(TPCEventTag),"TPCEventTag/I")
outtree.Branch("ntrk",addressof(ntrk),"ntrk/I")
outtree.Branch("x",tpc_x,"x[ntrk]/S")
outtree.Branch("y",tpc_y,"y[ntrk]/S")
outtree.Branch("z",tpc_z,"z[ntrk]/S")
for i in range(0,TotEnt-1):
plt.show()

'''
plt.imshow(TPCEvents[2500,:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
plt.show()
'''
