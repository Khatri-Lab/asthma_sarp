###############################################################
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier
###############################################################


###############################################################
out_folder = '/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/'

GSE115770_expr = pd.read_csv(out_folder + 'GSE115770_baseline_expr.csv', index_col = 0)
raw_sarp = pd.read_csv(out_folder + 'SARP_cluster_pheno_expr.csv')
GSE115770_pheno = pd.read_csv(out_folder + 'GSE115770_baseline_pheno.csv')

train_x = np.log2(raw_sarp[GSE115770_expr.columns] + 1)
train_y = raw_sarp.Cluster
test_x = np.log2(GSE115770_expr + 1)
###############################################################


###############################################################
kNN_cluster = KNeighborsClassifier(n_neighbors = 5)
kNN_cluster.fit(train_x, train_y)
test_y_baseline = kNN_cluster.predict(test_x)

GSE115770_pheno.index = GSE115770_pheno['library.sampleId']
test_y = pd.Series(test_y_baseline, name = 'pred_cluster_raw_sarp')
test_y.index = GSE115770_expr.index
GSE115770_pheno = GSE115770_pheno.merge(test_y, how = 'left', left_index = True, right_index = True)
GSE115770_pheno.to_csv(out_folder + 'GSE115770_baseline_pheno_clusters.csv')
###############################################################