import CountingModel as MODEL
from CountingModel import *
import pandas as pd
from sklearn.utils import shuffle
import InputManager

#filename = "TaggedTrainDataMerged.root"
#filename = "TaggedTrainDataP500.root"
#filename = "TaggedTrainDataP700.root"
#filename = "RealData05641.root"
#filename = "tagged_beam.root"
filename = "tagged_unipr_500.root"
scale = int(10)
model.load_weights("./Model_5/Model_All")
#TPCEvents,evnums,TotEnt = InputManager.LoadRealEventRDF(filename,scale)
model = MODEL.model

layer_outputs = [layer.output for layer in model.layers[:]]
active_model = models.Model(inputs=model.input, outputs = layer_outputs)

nev = 50
TPCEvents,TPCEventTags,TPCNTrk,evnums,TotEnt = InputManager.LoadEventRDF(filename,scale)
activation = active_model.predict(TPCEvents)
act= activation[10]
print(act.shape)
fig1=plt.figure()
plt.imshow(TPCEvents[nev,:,:,0])
fig2=plt.figure()
plt.matshow(act[nev,:,:,10])
plt.show()
for i in range(TotEnt):
    evnum[0]=evnums[i]
    tag[0]=TPCNTrk[i]
    pred[0]=int(MaxCh(predData[i]))
    for j in range(output_num):
        prb[j]=predData[i][j]
    outtree.Fill()
outfile.Write()
