import NeuralModel
from NeuralModel import *
import pandas as pd
from sklearn.utils import shuffle


#filename = "Saved_codes/TrainDataFairlyTagged.root"
filename = "TrainDataFairlyTaggedS.root"
#filename = "TrainDataTagged_wo_bg.root"
cb=callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.00001)
cb=callbacks.ModelCheckpoint(filepath="./Model_3/Model",save_weights_only=True,monitor='val_accuracy',mode='max',save_best_only=True)
file = ROOT.TFile.Open(filename,"READ")
tree = file.Get("tree")
cnt=0
scale = int(4)
TotEnt=int(tree.GetEntries()/scale)
print("Starting...")
TPCEventTags = np.zeros(TotEnt)
print("Loading Data...")
TPCEvents = np.zeros((TotEnt,nbin,nbin,max_depth))
for i in range(0,TotEnt):
    tree.GetEntry(int(i*scale))
    if i%1000 ==0:
        print(i)
    ent=tree
    TPCEventTags[i]=ent.TPCEventTag
    nhit = ent.nhtpc
    for nh in range(0,nhit):
        x=ent.x[nh]
        y=ent.y[nh]
        z=ent.z[nh]
        x,z=ToPixel(x,z)
        y=ToInt(y)
        for dep in range(0,max_depth):
            if TPCEvents[i,x,z,dep]==0:
                TPCEvents[i,x,z,0]=y/650
                break

TPCEvents,TPCEventTags = shuffle(TPCEvents,TPCEventTags)

print(TPCEvents.shape)
print(TPCEventTags.shape)
TPCEventTags = pd.get_dummies(TPCEventTags)
print(TPCEventTags.shape)
ratio = 0.8
TrainEntity = int(TotEnt*ratio)
TrainData=TPCEvents[:TrainEntity]
TrainLabel=TPCEventTags[:TrainEntity]
ValiData=TPCEvents[TrainEntity:]
ValiLabel=TPCEventTags[TrainEntity:]
print(TrainData.shape)
print(TrainLabel.shape)
print(type(TrainLabel))
#print("Tag = ",TPCEventTags[5])
#print("Tag = ",TrainLabel[5])
print("Tag = ",TrainLabel)
#plt.imshow(TPCEvents[0,:,:,0],interpolation='nearest')
#plt.gca().invert_yaxis()
#plt.show()


##########################Data Prepared###################
model = NeuralModel.model
fit_result = model.fit(TrainData,TrainLabel,epochs=epoch,batch_size=batch,validation_data=(ValiData,ValiLabel),callbacks=[cb])
cb = [
callbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.00001),
callbacks.ModelCheckpoint(filepath="./Model_3/Model",verbose=1,save_weights_only=True,monitor='val_accuracy',mode='max',save_best_only=True)
]

predTrain=model.predict(TrainData)
print(predTrain.shape)
predVali=model.predict(ValiData)


fit_history = fit_result.history
fig = plt.figure()
fig,((a1t,a2t),(a3t,a4t)) = plt.subplots(2,2)
a1t.hist(predVali[:,0],bins=20)
a1t.set_title("Validation_Class1")
a2t.hist(predVali[:,1],bins=20)
a2t.set_title("Validation_Class2")
a3t.hist(predVali[:,2],bins=20)
a3t.set_title("Validation_Class3")
a4t.hist(predVali[:,3],bins=20)
a4t.set_title("Validation_Class4")

fig2 = plt.figure()
fig2,((a1v,a2v),(a3v,a4v)) = plt.subplots(2,2)
a1v.hist(predTrain[:,0],bins=20)
a1v.set_title("Traning_Class1")
a2v.hist(predTrain[:,1],bins=20)
a2v.set_title("Traning_Class2")
a3v.hist(predTrain[:,2],bins=20)
a3v.set_title("Traning_Class3")
a4v.hist(predTrain[:,3],bins=20)
a4v.set_title("Traning_Class4")

loss=fit_history['loss']
val_loss=fit_history['val_loss']
val_acc=fit_history['val_accuracy']
acc=fit_history['accuracy']

fig3 = plt.figure()
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
