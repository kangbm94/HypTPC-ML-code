import Functions
from Functions import *
import tensorflow as tf
import os
from tensorflow.keras import datasets, layers, models,optimizers,callbacks
from keras.layers.convolutional import Conv2D
lr=0.1
batch = 64 
epoch = 40 
depth = 32
X= layers.Input(shape=[nbin,nbin,max_depth])
#X=layers.BatchNormalization()(X)
H= layers.Conv2D(depth,kernel_size=3,activation='relu')(X)
H= layers.Conv2D(depth,kernel_size=3,activation='relu')(X)
H= layers.AveragePooling2D()(H)
#H= layers.MaxPooling2D()(H)
H= layers.Conv2D(depth*2,kernel_size=3,activation='relu')(H)
H= layers.Conv2D(depth*2,kernel_size=3,activation='relu')(H)
H= layers.AveragePooling2D()(H)
#H= layers.MaxPooling2D()(H)
H= layers.Conv2D(depth*4,kernel_size=3,activation='relu')(H)
H= layers.Conv2D(depth*4,kernel_size=3,activation='relu')(H)
H= layers.AveragePooling2D()(H)
#H= layers.MaxPooling2D()(H)
H= layers.Conv2D(depth*8,kernel_size=3,activation='relu')(H)
H= layers.Conv2D(depth*8,kernel_size=3,activation='relu')(H)
#H=layers.BatchNormalization()(H)
H= layers.Dropout(0.5)(H)

H= layers.Flatten()(H)
Y= layers.Dense(output_num,activation='softmax')(H)
model=models.Model(X,Y);
model.summary()
model.compile(loss='categorical_crossentropy',optimizer=optimizers.Adadelta(learning_rate=lr),metrics='accuracy')
#model.compile(loss='categorical_crossentropy',optimizer=Adam(learn_rate),metrics='accuracy',ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.004))


