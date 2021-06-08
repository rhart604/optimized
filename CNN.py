#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

from Bio import SeqIO
from Bio.SeqIO import SeqRecord
import tensorflow

from keras.utils import np_utils
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Embedding, Conv1D, AveragePooling1D,MaxPooling1D,GlobalMaxPooling1D, Flatten,Dense
from tensorflow.keras.layers import Dropout, Activation
from tensorflow.keras.layers import Convolution2D, MaxPooling2D
from tensorflow.keras.preprocessing import sequence


# In[ ]:


## select input files for both human and mouse
hmfile = "SRR7661022_hm.fastq"
msfile = "SRR9637128_ms.fastq"


# In[ ]:


record_iterator = SeqIO.parse(hmfile,"fastq")
first_record = next(record_iterator)


# In[ ]:


padded_file_size = 1000 ## change to the idea size of dataset


def getRecordGen(filename,file_format):
    with open(filename) as handle:
        for record in SeqIO.parse(handle,file_format):
            yield record
            
def computeMaxSeqLength(original_fn,original_format):
    '''
    define the maximal length of the reads for padding
    '''
    
    max_seq_len = 0
    for record in getRecordGen(original_fn,original_format):
        if len(record.seq) > max_seq_len:
            max_seq_len = len(record.seq)
    return max_seq_len  
    
def cleanFastQFile(original_fn,output_fn_prefix,per_file_record_limit=50000,max_seq_len=76,original_format='fastq'):
    train_output_fn = output_fn_prefix+"_train.fasta"
    val_output_fn = output_fn_prefix+"_val.fasta"
    train_handle =  open(train_output_fn,"w")
    val_handle = open(val_output_fn,"w")
    handle = train_handle
    record_count = 0
    for record in getRecordGen(original_fn,original_format):
        padding = "A"*(max_seq_len-len(record.seq))
        padded_seq = (record.seq+padding)
        description = "Padded With A bp. Original length: " + str(len(record.seq)) 
        fasta_record = SeqRecord(padded_seq,record.id,description=description)
        
        SeqIO.write(fasta_record,handle,"fasta")
        record_count += 1
        
        #Start writting over validation file after reaching target for training file.
        if record_count >= per_file_record_limit:
            handle = val_handle
        
        #Stop iterating over source file
        if record_count >= per_file_record_limit*2:
            break
        
    return train_output_fn,val_output_fn
                  


# In[ ]:


hmpad = "hmpadded"
mspad = "mspadded"
hmpad_train,hmpad_val = cleanFastQFile(hmfile,hmpad)
mspad_train,mspad_val = cleanFastQFile(msfile,mspad)


# In[ ]:


def plotNullBbDistribution(filename,file_formalt,species):
    
    '''
    Get the number of 'N' in the input sequence files.
    input: sequence files, specific farmat, 'Mouse' or 'Human'
    output: histrogram of number of 'N' in inputfiles.
    '''
    
    null_values = []
    for record in getRecordGen(filename,"fasta"):
        if(record.seq.count('N') > 0):
            null_values.append(record.seq.count('N'))

            
    plt.figure(figsize=(8,6))
    plt.hist(null_values,bins=50,density=True)
    plt.xlabel('Number of "N" in read')
    plt.ylabel('ratio in dataset')
    plt.title(species)
    return null_values


# In[ ]:


sample = plotNullBbDistribution(mspad_sample,"fasta",'Mouse')
sample2 = plotNullBbDistribution(hmpad_sample,"fasta",'Human')


# In[ ]:


def getLabelDict(posIndexList,negIndexList):
    labels = {}
    for posIndex in posIndexList:
        labels[posIndex] = 1
    for negIndex in negIndexList:
        labels[negIndex] = 0
    return labels

def getIndex(fasta_file_path):
    with open(fasta_file_path) as handle:
        return SeqIO.index(fasta_file_path,"fasta")
    
def getIndexList(index):
    return list(index)

def getValueFromIndexes(index,pos_index,neg_index):
    if index in pos_index:
        return pos_index[index]
    if index in neg_index:
        return neg_index[index]

def createSamplesWithLabels(pos_fn,neg_fn):
    pos_index = getIndex(pos_fn)
    neg_index = getIndex(neg_fn)
    pos_index_list = getIndexList(pos_index)
    neg_index_list = getIndexList(neg_index)
    labels = getLabelDict(pos_index_list,neg_index_list)
    return (pos_index_list+neg_index_list,pos_index,neg_index,labels)
    
indexes_train,pos_index_train,neg_index_train,labels_train = createSamplesWithLabels(hmpad_train,mspad_train)
indexes_val,pos_index_val,neg_index_val,labels_val = createSamplesWithLabels(hmpad_val,mspad_val)


# In[ ]:


def seqToNumArray(seq):
    
    ## ACGTN in order
    alpha = "ACGTN" 
    seq_array = np.zeros(len(seq))
    for i in range(len(seq)):
        seq_array[i] = alpha.find(seq[i])
    return seq_array


# In[ ]:


class SplicingDataGenerator(tensorflow.keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self, samples,labels,pos_index,neg_index,seq_length=76,batch_size=128,shuffle=True):
        'Initialization'
        self.batch_size = batch_size
        self.labels = labels
        self.samples = samples
        self.seq_length = seq_length
        self.pos_index = pos_index
        self.neg_index = neg_index
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.samples) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        tmp_IDs = [self.samples[i] for  i in indexes]

        # Generate data
        X, y = self.__data_generation(tmp_IDs)

        return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.samples))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, tmp_IDs):
        'Generates data containing batch_size samples'
        # Initialization
        X = np.empty((self.batch_size,self.seq_length),dtype=int)
        y = np.empty((self.batch_size), dtype=int)

        # Generate data
        for i, ID in enumerate(tmp_IDs):
            # Store sample
            record = getValueFromIndexes(ID,self.pos_index,self.neg_index)
            X[i,] = seqToNumArray(record.seq)

            # Store class
            y[i] = self.labels[ID]

        return X,y
        


# In[ ]:


training_generator = SplicingDataGenerator(indexes_train,labels_train,pos_index_train,neg_index_train,shuffle=True)
val_generator = SplicingDataGenerator(indexes_val,labels_val,pos_index_val,neg_index_val,shuffle=True)


# In[ ]:


def createCNNCategoricalEmbedding(seq_length,kernel_size,kernel_number,fst_layer_avg_pooling,snd_layer_global_pooling,hidden_layer_units):
    model = tensorflow.keras.Sequential()
    model.add(Embedding(5,5,input_length=seq_length))
    model.add(Conv1D(kernel_number[0],kernel_size,padding="same",activation="relu"))
    if fst_layer_avg_pooling:
        model.add(AveragePooling1D())
    else:
        model.add(MaxPooling1D())
#     model.add(Dropout(0.2)) # if dropout is used
#     model.add(Conv1D(kernel_number[1],kernel_size,padding="valid",activation="relu"))
#     if snd_layer_global_pooling:
#         model.add(GlobalMaxPooling1D())
#     else:
#         model.add(MaxPooling1D())
#     model.add(Dropout(0.2)
    model.add(Flatten())
    model.add(Dense(units=76,activation='relu'))
    model.add(Dense(1,activation='sigmoid'))
    model.compile(loss="binary_crossentropy",optimizer='adam',metrics=['accuracy'])
    
    return model

    


# In[ ]:


def trainCNNCategoricalEmbedding(seq_length,kernel_size,kernel_number,fst_layer_pooling,snd_layer_pooling,hidden_layer_units):
    training_generator = SplicingDataGenerator(indexes_train,labels_train,pos_index_train,neg_index_train,shuffle=False)
    val_generator = SplicingDataGenerator(indexes_val,labels_val,pos_index_val,neg_index_val,shuffle=True)
    model = createCNNCategoricalEmbedding(seq_length,kernel_size,kernel_number,fst_layer_pooling,snd_layer_pooling,hidden_layer_units)
    print(model.summary())
    hist = model.fit(training_generator,validation_data=val_generator,epochs=50, verbose=1)
#     parameter_summary = "kernel_size: {kernel_size}".format(kernel_size)
#     print(parameter_summary)
    plt.plot(hist.history['accuracy'])
    plt.plot(hist.history['val_accuracy'])
    plt.title('model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.show()
    return model


# In[ ]:


model = trainCNNCategoricalEmbedding(76,3,(20,30),False,False,76)

