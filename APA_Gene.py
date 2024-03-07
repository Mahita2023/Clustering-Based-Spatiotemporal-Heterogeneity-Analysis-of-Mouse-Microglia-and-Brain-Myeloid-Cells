
# Name:Python code AffinityPropagation clustering for Gene data
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

dataFrame = pandas.read_excel('Gene_data_afterpreprocessing.xlsx', sheet_name=0)
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
model = AffinityPropagation(damping=0.6, max_iter=5000, preference=-4000, random_state=20191)
model.fit(Gene_data_reduction)
labels=model.labels_
Center=model.cluster_centers_
#print('silhouette_score: \n', silhouette_score(Gene_data_reduction,labels,metric='euclidean'))
r1=pandas.Series(model.labels_).value_counts()
r2=pandas.Series(model.labels_)
Title="Clustering result for Gene dataset"
pyplot.figure()
r1.plot(kind='bar',grid=True,title=Title)

reducer=umap.UMAP(random_state=20191)
embedding=reducer.fit_transform(Gene_data_reduction)
pyplot.figure()
pyplot.rcParams['font.sans-serif'] = ['SimHei']
pyplot.rcParams['axes.unicode_minus'] = False
pyplot.scatter(embedding[:, 0], embedding[:, 1], c=labels, s=5, cmap='Spectral')
pyplot.colorbar()
pyplot.title('UMAP figure for Gene dataset')
pyplot.xlabel('UMAP_1')
pyplot.ylabel('UMAP_2')
pyplot.show()