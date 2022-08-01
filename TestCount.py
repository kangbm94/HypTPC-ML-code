import ConcatNetwork as MODEL
from ConcatNetwork import *
import pandas as pd
from sklearn.utils import shuffle
import InputManager

#filename = "TaggedTrainDataMerged.root"
#filename = "TaggedTrainDataP500.root"
#filename = "TaggedTrainDataP700.root"
filename = "RealData05641.root"
#filename = "tagged_beam.root"
#filename = "tagged_train.root"
scale = int(20)
frag = 15
#TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCEventTags,TPCNTrk,evnums,TotEnt = InputManager.LoadEventXYZ(filename,scale)
TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCEventTags,TPCNTrk,evnums,TotEnt = InputManager.LoadRealEventXYZ(filename,scale,frag)
#TPCEvents,evnums,TotEnt = InputManager.LoadRealEventRDF(filename,scale)
model = MODEL.model
model.load_weights("./Model_5/Model_All")

predData=model.predict([TPCEventsXZ,TPCEventsYZ,TPCEventsXY])
print(predData.shape)
outfilename = "PredictedDataReal05641_"+str(frag)+".root"
outfile = ROOT.TFile.Open(outfilename,"recreate")
outtree = ROOT.TTree("tree","tree")
evnum = array('i',[0])
tag = array('i',[0])
pred = array('i',[0])
prb = np.ndarray(shape=(output_num),dtype="float64")
outtree.Branch("evnum",evnum,"evnum/I")
outtree.Branch("tag",tag,"tag/I")
outtree.Branch("pred",pred,"pred/I")
outtree.Branch("prb",prb,"prb[10]/D")
layer_outputs = [layer.output for layer in model.layers[:4]]
#active_model = models.Model(inputs=model.input, outputs = layer_outputs)
#activation = active_model.predict([TPCEventsXZ,TPCEventsYZ,TPCEventsXY])
#act= activation[3]
#print(act.shape)
#fig1=plt.figure()
#nev = 50
#plt.imshow(TPCEventsXZ[nev,:,:,0])
#fig2=plt.figure()
#plt.matshow(act[nev,:,:,0])
#plt.show()
for i in range(TotEnt):
    evnum[0]=evnums[i]
#    tag[0]=TPCNTrk[i]
    pred[0]=int(MaxCh(predData[i]))
    for j in range(output_num):
        prb[j]=predData[i][j]
    outtree.Fill()
outfile.Write()
