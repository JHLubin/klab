#Title: Neural Network on prediction of cleavage sites
#Author: Changpeng Lu
# TensorFlow and tf.keras
import tensorflow as tf
from tensorflow import keras
from keras import regularizers
#Helper Libraries
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, auc
import random

data_train = pd.read_csv("new_X_train_30000.csv", index_col = 0)
y_train = pd.read_csv("new_y_train_30000.csv",index_col=0)
data_test = pd.read_csv("X_test.csv", index_col = 0)
remain = np.genfromtxt("test_index_30000",dtype='str')
data_test = data_test.loc[remain,:]
#result = y_train.copy().values

class_names = ["CLEAVED","MIDDLE","UNCLEAVED"]
#class_names = ["CLEAVED","UNCLEAVED"]
#cn = sum(data_train['result'] == 'CLEAVED')
#un = sum(data_train['result'] == 'UNCLEAVED')
#mn = sum(data_train['result'] == 'MIDDLE')
#rcn = np.random.randint(0,cn,size = un-cn)
#rmn = np.random.randint(0,mn,size = un-mn)
m = data_train.shape[0]
#data_train.index = range(0,data_train.shape[0])
#new_data_train = data_train.copy()
#data_train_c = data_train[data_train['result'] == 'CLEAVED'].copy()
#for num in rcn:
#	new_data_train.loc[len(result)+1] = data_train_c.iloc[num,:]
#	m=m+1
#data_train_m = data_train[data_train['result'] == 'MIDDLE'].copy()
#for nn in rmn:
#	new_data_train.loc[m+1] = data_train_m.iloc[nn,:]
#	m=m+1
#print(new_data_train.shape[0] == data_train[data_train['result']=='UNCLEAVED'].shape[0]*3)

result_n = y_train['result'].copy().values
newre = np.zeros(len(result_n))
for i in range(0,len(result_n)):
    if result_n[i] == 'CLEAVED':
        newre[i] = 0
    elif result_n[i] == 'MIDDLE':
        newre[i] = 1.0
    else:
        newre[i] = 2.0
newre = newre.astype(int)
model = keras.Sequential([keras.layers.Dense(1024,activation=tf.nn.relu),keras.layers.Dense(3, activation=tf.nn.softmax)])

model.compile(optimizer=tf.train.AdamOptimizer(),loss='sparse_categorical_crossentropy',metrics=['accuracy'])

model.fit(data_train.values, newre, epochs=5)

result_test = data_test.copy().values
newre_test = np.zeros(len(result_test))
for i in range(0,len(result_test)):
    if result_test[i] == 'CLEAVED':
        newre_test[i] = 0
    elif result_test[i] == 'MIDDLE':
        newre_test[i] = 1
    else:
        newre_test[i] = 2
newre_test = newre_test.astype(int)

test_loss, test_acc = model.evaluate(data_test.values, newre_test)

print('Test accuracy:', test_acc)

predictions = model.predict(data_test.values)

label_pre = []
for i in range(0,len(predictions)):
        label_pre.append(np.argmax(predictions[i]))
#label_df = pd.DataFrame({'result' : label_pre},index=data_test.index.values)
#df1 = pd.DataFrame({'score': predictions[:,0]},index=data_test.index.values)
#df2 = pd.DataFrame({'score' : predictions[:,1]},index=data_test.index.values)
#df3 = pd.DataFrame({'score' : predictions[:,2]},index=data_test.index.values)
#result_cleaved = pd.concat([data_test,df1,label_df],axis=1)
#result_middle =pd.concat([data_test,df2,label_df],axis=1)
#result_uncleaved = pd.concat([data_test,df3,label_df],axis=1)
#result_cleaved.to_csv("result_cleaved_ann")
#result_middle.to_csv("result_middle_ann")
#result_uncleaved.to_csv("result_uncleaved_ann")
np.savetxt("label_result_ann_resample_30000",label_pre,delimiter=",")
np.savetxt("result_ann_resample_30000",predictions,delimiter=",")
print(data_train.shape[0],data_test.shape[0])

#def plotROC(preds, truth, classification, name):
#    fpr, tpr, thresholds = roc_curve(truth, preds , pos_label = classification)
#    roc_auc = auc(fpr, tpr)
    # chooses a random color for plotting
#    c = (np.random.rand(), np.random.rand(), np.random.rand())
    #create the plot
#    plt.plot(fpr, tpr, color = c, label = name + ' (AUC = %0.3f)' % roc_auc)
#    plt.plot([0, 1], [0, 1], 'k--')
#    plt.xlim([0.0, 1.0])
#    plt.ylim([0.0, 1.0])
#    plt.xlabel('FPR')
#    plt.ylabel('TPR')
#    plt.title('ROC')
#    plt.legend(loc="lower right")
    
#    return roc_auc


#mu_data_test = pd.read_csv("mu_data_test.csv",index_col=0)
#mc_data_test = pd.read_csv("mc_data_test.csv",index_col=0)
#mm_data_test = pd.read_csv("mm_data_test.csv",index_col=0)
# Plot ROC for both models
#fig = plt.figure(figsize = (16, 12))
#plt.subplot(221)
#plotROC(predictions[:,0], mu_data_test['result'] ,'CLEAVED','CLEAVED versus REST')
#plt.subplot(222)
#plotROC(predictions[:,2], mc_data_test['result'],'UNCLEAVED','UNCLEAVED versus REST')
#plt.subplot(223)
#plotROC(predictions[:,1], mm_data_test['result'], 'MIDDLE', 'MIDDLE versus REST')
#plt.show()
