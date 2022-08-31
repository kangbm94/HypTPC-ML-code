import Functions
from Functions import *
import tensorflow as tf
import os
from tensorflow.keras import datasets, layers, models,optimizers,callbacks
from keras.layers.convolutional import Conv2D
lr=1.
batch = 30 
epoch = 10 
depth = 16
X= layers.Input(shape=[nbin,nbin,max_depth])
X=layers.BatchNormalization()(X)
H= layers.Conv2D(depth,kernel_size=5,activation='relu',padding = 'same')(X)
H= layers.Conv2D(depth,kernel_size=3,activation='relu',padding = 'same')(H)
H= layers.AveragePooling2D()(H)
H= layers.Conv2D(depth*2,kernel_size=3,activation='relu',padding = 'same')(H)
H= layers.Conv2D(depth*2,kernel_size=3,activation='relu',padding = 'same')(H)
H= layers.AveragePooling2D()(H)
H= layers.Dropout(0.5)(H)
H= layers.Conv2D(depth*4,kernel_size=3,activation='relu',padding = 'same')(H)
H= layers.Conv2D(depth*4,kernel_size=3,activation='relu',padding = 'same')(H)

#X = tf.keras.applications.VGG16(include_top=False, input_shape=(nbin,nbin,max_depth),weights='imagenet')
H= layers.Flatten()(H)
#H= layers.Dense(500,activation='relu')(H)
#H= layers.Dense(32,activation='relu')(H)
Y= layers.Dense(4,activation='softmax')(H)
model=models.Model(X,Y);
model.summary()
model.compile(loss='categorical_crossentropy',optimizer=optimizers.Adadelta(learning_rate=lr),metrics='accuracy')
#model.compile(loss='categorical_crossentropy',optimizer=Adam(learn_rate),metrics='accuracy',ReduceLROnPlateau(monitor='val_loss',factor=0.2,patience=3,min_lr=0.004))


