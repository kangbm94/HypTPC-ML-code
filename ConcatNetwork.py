import Functions
from Functions import *
import tensorflow as tf
import os
from tensorflow.keras.models import Model
from tensorflow.keras import datasets, layers, models,optimizers,callbacks
from keras.layers import *
lr=.1
batch = 64 
epoch = 100 
depth = 8

def LowCNN(X):
    H= layers.Conv2D(depth,kernel_size=3,activation='relu')(X)
    H= layers.Conv2D(depth,kernel_size=3,activation='relu')(H)
    H= layers.MaxPooling2D()(H)
    H= layers.Conv2D(depth*2,kernel_size=3,activation='relu')(H)
    H= layers.Conv2D(depth*2,kernel_size=3,activation='relu')(H)
    H= layers.MaxPooling2D()(H)
    H= layers.Conv2D(depth*4,kernel_size=3,activation='relu')(H)
    H= layers.Conv2D(depth*4,kernel_size=3,activation='relu')(H)
    H= layers.MaxPooling2D()(H)
    H= layers.Conv2D(depth*8,kernel_size=3,activation='relu')(H)
    H= layers.Conv2D(depth*8,kernel_size=3,activation='relu')(H)
    H= layers.Dropout(0.5)(H)
    H= layers.Flatten()(H)
#    H= layers.Dense(output_num,activation="softmax")(H)
    return H


'''
Y=LowCNN(X)
model = Model(X,Y)
model.summary()

#DL = layers.Dense(output_num,activation='softmax')(XYZ)
#model=models.Model(inputs=[XY,YZ,ZX],Y)
#model=models.Model(inputs=[XY,YZ,XY],outputs=DL)
'''
XZ= layers.Input(shape=[nbin,nbin,max_depth])
YZ= layers.Input(shape=[nbin,nbin,max_depth])
XY= layers.Input(shape=[nbin,nbin,max_depth])
B_XZ= LowCNN(XZ)
B_YZ= LowCNN(YZ)
B_XY= LowCNN(XY)
XYZ = layers.Concatenate()([B_XZ,B_YZ,B_XY])
#XYZ = layers.Concatenate()([B_XZ,B_YZ])
#XYZ = layers.Dropout(0.2)(XYZ)
DL= layers.Dense(output_num,activation='softmax')(XYZ)
#DL= layers.Dense(output_num,activation='softmax')(B_XZ)
#DL= layers.Dense(output_num,activation='softmax')(DL)
#model=models.Model(inputs=[XZ,YZ,XY],outputs=DL);
model=models.Model(inputs=[XZ,YZ,XY],outputs=DL);
#model=models.Model(inputs=XZ,outputs=DL);
#model.build(input_shape=[nbin,nbin,max_depth])
model.summary()
model.compile(loss='categorical_crossentropy',optimizer=optimizers.Adadelta(learning_rate=lr),metrics='accuracy')
#model.compile(loss='categorical_crossentropy',optimizer=Adam(learn_rate),metrics='accuracy',ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.004))

