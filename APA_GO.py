
# Name:Python code AffinityPropagation clustering for GO data
# Author: Songtao Wei

import numpy
import pandas
import sklearn
import matplotlib.pyplot
import sklearn.cluster
from sklearn.datasets import load_digits
import umap
from matplotlib import pyplot,colors
from sklearn.cluster import AffinityPropagation, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

dataFrame = pandas.read_excel('GO_data_afterpreprocessing.xlsx', sheet_name=0)
Samplename=list(dataFrame)
Samplename=Samplename[1:1809]
Data= dataFrame
Number_miss=Data.isnull().sum()
Data=Data.dropna(axis=0,subset=Samplename)
Genename=Data.to_numpy()[:,0]
data_np=numpy.transpose(Data.to_numpy()[:,1:1809])

pca=PCA(n_components=30, random_state=20191)
pca.fit(data_np)
Gene_data_reduction=pca.transform(data_np)
#print(pca.explained_variance_ratio_)
#print(pca.explained_variance_)
#print(numpy.sum(pca.explained_variance_ratio_))
model = AffinityPropagation(damping=0.6, max_iter=5000, preference=-12500, random_state=20191)
model.fit(Gene_data_reduction)
labels=model.labels_
Center=model.cluster_centers_
#print('silhouette_score: \n', silhouette_score(Gene_data_reduction,labels,metric='euclidean'))
r1=pandas.Series(model.labels_).value_counts()
r2=pandas.Series(model.labels_)
Title="Clustering result for GO dataset"
pyplot.figure()
r1.plot(kind='bar',grid=True,title=Title)

reducer=umap.UMAP(random_state=20191)
embedding=reducer.fit_transform(Gene_data_reduction)
print(embedding.shape)
pyplot.figure()
pyplot.rcParams['font.sans-serif'] = ['SimHei']
pyplot.rcParams['axes.unicode_minus'] = False
pyplot.scatter(embedding[:, 0], embedding[:, 1], c=labels, s=5, cmap='Spectral')
pyplot.colorbar()
pyplot.title('UMAP figure for GO dataset')
pyplot.xlabel('UMAP_1')
pyplot.ylabel('UMAP_2')
pyplot.show()

Table_AR = pandas.read_excel('Age_Region_data.xlsx', sheet_name=0)
AR=Table_AR.to_numpy()
AR_P7=AR[AR[:,1]=="P7",:]
AR_P7_CB=AR_P7[AR_P7[:,2]=="CB",:]
Index_P7_CB=AR_P7_CB[:,0]
labels_P7_CB=labels[Index_P7_CB.astype(int)]
R_P7_CB=pandas.Series(labels_P7_CB).value_counts()
Title="Clustering result for P7~CB"
pyplot.figure()
R_P7_CB.plot(kind='bar',grid=True,title=Title)
AR_P7_CP=AR_P7[AR_P7[:,2]=="CP",:]
Index_P7_CP=AR_P7_CP[:,0]
labels_P7_CP=labels[Index_P7_CP.astype(int)]
R_P7_CP=pandas.Series(labels_P7_CP).value_counts()
Title="Clustering result for P7~CP"
pyplot.figure()
R_P7_CP.plot(kind='bar',grid=True,title=Title)
AR_P7_CTX=AR_P7[AR_P7[:,2]=="CTX",:]
Index_P7_CTX=AR_P7_CTX[:,0]
labels_P7_CTX=labels[Index_P7_CTX.astype(int)]
R_P7_CTX=pandas.Series(labels_P7_CTX).value_counts()
Title="Clustering result for P7~CTX"
pyplot.figure()
R_P7_CTX.plot(kind='bar',grid=True,title=Title)
AR_P7_HIP=AR_P7[AR_P7[:,2]=="HIP",:]
Index_P7_HIP=AR_P7_HIP[:,0]
labels_P7_HIP=labels[Index_P7_HIP.astype(int)]
R_P7_HIP=pandas.Series(labels_P7_HIP).value_counts()
Title="Clustering result for P7~HIP"
pyplot.figure()
R_P7_HIP.plot(kind='bar',grid=True,title=Title)
AR_P7_OB=AR_P7[AR_P7[:,2]=="OB",:]
Index_P7_OB=AR_P7_OB[:,0]
labels_P7_OB=labels[Index_P7_OB.astype(int)]
R_P7_OB=pandas.Series(labels_P7_OB).value_counts()
Title="Clustering result for P7~OB"
pyplot.figure()
R_P7_OB.plot(kind='bar',grid=True,title=Title)
AR_P7_STR=AR_P7[AR_P7[:,2]=="STR",:]
Index_P7_STR=AR_P7_STR[:,0]
labels_P7_STR=labels[Index_P7_STR.astype(int)]
R_P7_STR=pandas.Series(labels_P7_STR).value_counts()
Title="Clustering result for P7~STR"
pyplot.figure()
R_P7_STR.plot(kind='bar',grid=True,title=Title)

AR_P60=AR[AR[:,1]=="P60",:]
AR_P60_CB=AR_P60[AR_P60[:,2]=="CB",:]
Index_P60_CB=AR_P60_CB[:,0]
labels_P60_CB=labels[Index_P60_CB.astype(int)]
R_P60_CB=pandas.Series(labels_P60_CB).value_counts()
Title="Clustering result for P60~CB"
pyplot.figure()
R_P60_CB.plot(kind='bar',grid=True,title=Title)
AR_P60_CP=AR_P60[AR_P60[:,2]=="CP",:]
Index_P60_CP=AR_P60_CP[:,0]
labels_P60_CP=labels[Index_P60_CP.astype(int)]
R_P60_CP=pandas.Series(labels_P60_CP).value_counts()
Title="Clustering result for P60~CP"
pyplot.figure()
R_P60_CP.plot(kind='bar',grid=True,title=Title)
AR_P60_CTX=AR_P60[AR_P60[:,2]=="CTX",:]
Index_P60_CTX=AR_P60_CTX[:,0]
labels_P60_CTX=labels[Index_P60_CTX.astype(int)]
R_P60_CTX=pandas.Series(labels_P60_CTX).value_counts()
Title="Clustering result for P60~CTX"
pyplot.figure()
R_P60_CTX.plot(kind='bar',grid=True,title=Title)
AR_P60_HIP=AR_P60[AR_P60[:,2]=="HIP",:]
Index_P60_HIP=AR_P60_HIP[:,0]
labels_P60_HIP=labels[Index_P60_HIP.astype(int)]
R_P60_HIP=pandas.Series(labels_P60_HIP).value_counts()
Title="Clustering result for P60~HIP"
pyplot.figure()
R_P60_HIP.plot(kind='bar',grid=True,title=Title)
AR_P60_OB=AR_P60[AR_P60[:,2]=="OB",:]
Index_P60_OB=AR_P60_OB[:,0]
labels_P60_OB=labels[Index_P60_OB.astype(int)]
R_P60_OB=pandas.Series(labels_P60_OB).value_counts()
Title="Clustering result for P60~OB"
pyplot.figure()
R_P60_OB.plot(kind='bar',grid=True,title=Title)
AR_P60_STR=AR_P60[AR_P60[:,2]=="STR",:]
Index_P60_STR=AR_P60_STR[:,0]
labels_P60_STR=labels[Index_P60_STR.astype(int)]
R_P60_STR=pandas.Series(labels_P60_STR).value_counts()
Title="Clustering result for P60~STR"
pyplot.figure()
R_P60_STR.plot(kind='bar',grid=True,title=Title)

