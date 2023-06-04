import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve
import numpy as np
from sklearn import metrics
import pdb
pos_scores = pd.read_csv('../../data/dev/blackbox_cracking/all_pos_scores.csv', header=None)
neg_scores = pd.read_csv('../../data/dev/blackbox_cracking/all_neg_scores.csv', header=None)
#Check overlap:
print('Overlap:', neg_scores[neg_scores[0].isin(pos_scores[0])])
#plt.violinplot([pos_scores[3], neg_scores[3]])
#plt.show()
print('Positive set:')
print(pos_scores.describe())
print('Negative set:')
print(neg_scores.describe())

y_true = np.zeros(len(pos_scores)+len(neg_scores))
y_true[:len(pos_scores)]=1
probas_pred = np.concatenate([pos_scores[3], neg_scores[3]])
precision, recall, thresholds = precision_recall_curve(y_true, probas_pred)
AUC = metrics.auc(recall, precision)
print('PR AUC:', AUC)
# plt.plot(recall, precision, label='AUC:'+str(np.round(AUC,3)))
# plt.legend()
# plt.xlabel('Recall')
# plt.ylabel('Precision')
# plt.ylim([0.5,1])
# plt.show()
#
# #ROC
# fpr, tpr, thresholds = metrics.roc_curve(y_true, probas_pred)
# AUC = metrics.auc(fpr, tpr)
# plt.plot(fpr, tpr, label='AUC:'+str(np.round(AUC,3)))
# plt.legend()
# plt.xlabel('FPR')
# plt.ylabel('TPR')
# plt.show()

#y pred
for t in np.arange(0,1,0.01):
    y_pred = np.zeros(len(pos_scores)+len(neg_scores))
    y_pred[probas_pred>=t]=1
    print('Accuracy:',t, metrics.accuracy_score(y_true, y_pred))
