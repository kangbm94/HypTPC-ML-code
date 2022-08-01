import ConcatNetwork as MODEL
from ConcatNetwork import *
import pandas as pd
from sklearn.utils import shuffle
import InputManager
#filename = "Saved_codes/TrainDataFairlyTagged.root"
#filename = "taggedTrainDataMerged.root"
#filename = "tagged_unipr_500.root"
#filename = "tagged_beam.root"
#filename = "tagged_beam.root"
filename = "tagged_train.root"
cnt=0
scale = int(1)
print("Starting...")
print("Loading Data...")
#TPCEvents,TPCEventTags,TPCNTrk,evnums,TotEnt = InputManager.LoadEvent(filename,scale)
TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCEventTags,TPCNTrk,evnums,TotEnt = InputManager.LoadEventXYZ(filename,scale)
TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCNTrk,evnums = shuffle(TPCEventsXZ,TPCEventsYZ,TPCEventsXY,TPCNTrk,evnums)
ratio = 0.8
TrainEntity = int(TotEnt*ratio)
TDXZ=TPCEventsXZ[:TrainEntity]
TDYZ=TPCEventsYZ[:TrainEntity]
TDXY=TPCEventsXY[:TrainEntity]
TrainLabel=TPCNTrk[:TrainEntity]
evnumsT=evnums[:TrainEntity]

VDXZ=TPCEventsXZ[TrainEntity:]
VDYZ=TPCEventsYZ[TrainEntity:]
VDXY=TPCEventsXY[TrainEntity:]

ValiLabel=TPCNTrk[TrainEntity:]
evnumsV=evnums[TrainEntity:]
TL=TrainLabel;
VL=ValiLabel;
TrainLabel = OneHotEncoding(TrainLabel)
ValiLabel = OneHotEncoding(ValiLabel)
#TrainLabel = pd.get_dummies(TrainLabel)
#ValiLabel = pd.get_dummies(ValiLabel)
fig1 = plt.figure()
plt.imshow(TPCEventsXZ[0,:,:,0],interpolation='nearest')
fig2 = plt.figure()
plt.imshow(TPCEventsYZ[0,:,:,0],interpolation='nearest')
#plt.gca().invert_yaxis()


##########################Data Prepared###################
model = MODEL.model
cb = [
    callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.00001),
    callbacks.ModelCheckpoint(filepath="./Model_5/Model_All",verbose=1,save_weights_only=True,monitor='val_accuracy',mode='max',save_best_only=True)
]
fit_result = model.fit([TDXZ,TDYZ,TDXY],TrainLabel,epochs=epoch,batch_size=batch,validation_data=([VDXZ,VDYZ,VDXY],ValiLabel),callbacks=[cb])
#fit_result = model.fit([TDXZ,TDYZ],TrainLabel,epochs=epoch,batch_size=batch,validation_data=([VDXZ,VDYZ],ValiLabel),callbacks=[cb])
#fit_result = model.fit(TDXZ,TrainLabel,epochs=epoch,batch_size=batch,validation_data=(VDXZ,ValiLabel),callbacks=[cb])
predTrain=model.predict([TDXZ,TDYZ,TDXY])
print(predTrain.shape)
predVali=model.predict([VDXZ,VDYZ,VDXY])
print(predTrain)
outfilename = "ValidationData.root"
outfile = ROOT.TFile.Open(outfilename,"recreate")
outtree = ROOT.TTree("tree","tree")
evnum = array('i',[0])
tag = array('i',[0])
pred = array('i',[0])
prb = np.ndarray(shape=(output_num),dtype="float64")
outtree.Branch("tag",tag,"tag/I")
outtree.Branch("evnum",evnum,"evnum/I")
outtree.Branch("pred",pred,"pred/I")
outtree.Branch("prb",prb,"prb[10]/D")
print(predVali.size)
print(predVali[:].size)
print(predVali[0][:].size)
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

for i in range(int((predVali.size)/(predVali[:][0].size))):
    tag[0]=int(VL[i])
    evnum[0]=evnumsV[i]
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
    evnum[0]=evnumsT[i]
    for j in range(output_num):
        prb[j]=predTrain[i][j]
    outtree2.Fill()
outfile2.Write()
'''
plt.imshow(TPCEvents[2500,:,:,0],interpolation='nearest')
plt.gca().invert_yaxis()
plt.show()
'''
