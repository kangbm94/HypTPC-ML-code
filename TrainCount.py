import CountingModel as MODEL
from CountingModel import *
import pandas as pd
from sklearn.utils import shuffle
from ROOT import TH2D
import InputManager
#filename = "Saved_codes/TrainDataFairlyTagged.root"
filename = "TaggedTrainDataP300.root"
#filename = "TaggedTrainDataKPXi.root"
#filename = "TaggedTrainDataKBeam.root"
#filename = "TrainDataTagged_wo_bg.root"
file = ROOT.TFile.Open(filename,"READ")
tree = file.Get("tree")
cnt=0
scale = int(1)
TotEnt=int(tree.GetEntries()/scale)
print("Starting...")
print("Loading Data...")
#fig1,((a1t,a2t,a3t,a4t),(a5t,a6t,a7t,a8t)) = plt.subplots(4,2)
#fig2,((a1v,a2v,a3v,a4v),(a5v,a6v,a7v,a8v)) = plt.subplots(4,2)
TPCEvents,TPCEventTags,TPCNTrk = InputManager.LoadEvent(tree,scale)
#TPCEvents,TPCEventTags = shuffle(TPCEvents,TPCEventTags)
TPCEvents,TPCNTrk = shuffle(TPCEvents,TPCNTrk)

print(TPCEvents.shape)
print(TPCNTrk.shape)
#TPCEventTags = pd.get_dummies(TPCEventTags)
print(TPCNTrk.shape)
ratio = 0.8
TrainEntity = int(TotEnt*ratio)
TrainData=TPCEvents[:TrainEntity]
TrainLabel=TPCNTrk[:TrainEntity]
ValiData=TPCEvents[TrainEntity:]
ValiLabel=TPCNTrk[TrainEntity:]
TL=TrainLabel;
VL=ValiLabel;
TrainLabel = pd.get_dummies(TrainLabel)
ValiLabel = pd.get_dummies(ValiLabel)
print(TrainData.shape)
print(TrainLabel.shape)
print(type(TrainLabel))
#print("Tag = ",TPCEventTags[5])
#print("Tag = ",TrainLabel[5])
print("Tag = ",TrainLabel)
fig3 = plt.figure()
plt.imshow(TPCEvents[0,:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
#plt.show()


##########################Data Prepared###################
model = MODEL.model
cb = [
    callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.00001),
    callbacks.ModelCheckpoint(filepath="./Model_3/Model",verbose=1,save_weights_only=True,monitor='val_accuracy',mode='max',save_best_only=True)
]
fit_result = model.fit(TrainData,TrainLabel,epochs=epoch,batch_size=batch,validation_data=(ValiData,ValiLabel),callbacks=[cb])

predTrain=model.predict(TrainData)
print(predTrain.shape)
predVali=model.predict(ValiData)
print(predTrain)
outfilename = "ValidationData.root"
outfile = ROOT.TFile.Open(outfilename,"recreate")
outtree = ROOT.TTree("tree","tree")
evnum = array('i',[0])
tag = array('i',[0])
pred = array('i',[0])
prb = np.ndarray(shape=(output_num),dtype="float64")
outtree.Branch("tag",tag,"tag/I")
outtree.Branch("pred",pred,"pred/I")
outtree.Branch("prb",prb,"prb[10]/D")
print(predVali.size)
print(predVali[:].size)
print(predVali[0][:].size)
for i in range(int((predVali.size)/(predVali[:][0].size))):
    tag[0]=int(VL[i])
    pred[0]=int(MaxCh(predVali[i]))
    for j in range(output_num):
        prb[j]=predVali[i][j]
    outtree.Fill()
outfile.Write()

outfilename2 = "TrainedData.root"
outfile2 = ROOT.TFile.Open(outfilename2,"recreate")
outtree2 = ROOT.TTree("tree","tree")
outtree2.Branch("tag",tag,"tag/I")
outtree2.Branch("pred",pred,"pred/I")
outtree2.Branch("prb",prb,"prb[10]/D")
for i in range(int((predTrain.size)/(predTrain[:][0].size))):
    tag[0]=int(TL[i])
    pred[0]=int(MaxCh(predTrain[i]))
    for j in range(output_num):
        prb[j]=predTrain[i][j]
    outtree2.Fill()
outfile2.Write()
'''
fig1,trainhist = plt.subplots(2,4)
fig2,validhist = plt.subplots(2,4)
trainhist=trainhist.flatten()
validhist=validhist.flatten()
'''
fit_history = fit_result.history
#fig1.title("Training_NTrk")
#fig2.title("Validation_NTrk")
loss=fit_history['loss']
val_loss=fit_history['val_loss']
val_acc=fit_history['val_accuracy']
acc=fit_history['accuracy']

fig4 = plt.figure()
epochs = range(1,len(loss)+1)
plt.subplot(2,1,1)
plt.plot(epochs,loss,'bo',label = 'Training Loss')
plt.plot(epochs,val_loss,'b',label = 'Validation Loss')
plt.xlabel('epochs')
plt.ylabel('loss')
plt.legend()


plt.subplot(2,1,2)
plt.plot(epochs,acc,'ro',label = 'Training Accuracy')
plt.plot(epochs,val_acc,'r',label = 'Validation Accuracy')
plt.xlabel('epochs')
plt.ylabel('accuracy')
plt.legend()
plt.show()

'''
plt.imshow(TPCEvents[2500,:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
plt.show()
'''
