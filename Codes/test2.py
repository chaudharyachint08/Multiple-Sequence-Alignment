import numpy as np
import sklearn.cluster
import distance

import time

init = time.clock()

words = open('myio.py').read().split('\n') #Replace this line

words = np.asarray(words) #So that indexing with a list will work

lev_similarity = -1*np.array( [ [distance.levenshtein(w1,w2) for w1 in words] for w2 in words ] )


affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=0.5)
affprop.fit(lev_similarity)

count = 0
for cluster_id in np.unique(affprop.labels_):
    count+=1
    head = words[affprop.cluster_centers_indices_[cluster_id]]
    cluster = np.unique(words[np.nonzero(affprop.labels_==cluster_id)])
    cluster_str = ", ".join(cluster)
    #print("[%s] --> %s" % (head, cluster_str))

print('Total Clusters Formed',count)
print(time.clock()-init)



'''
kmns = sklearn.cluster.KMeans(n_clusters=5)
kmns.fit(lev_similarity)

count = 0
for cluster_id in np.unique(kmns.labels_):
    count+=1
    head = words[kmns.cluster_centers_indices_[cluster_id]]
    cluster = np.unique(words[np.nonzero(kmns.labels_==cluster_id)])
    cluster_str = ", ".join(cluster)
    #print("[%s] --> %s" % (head, cluster_str))

print('Total Clusters Formed',count)
print(time.clock()-init)
'''
