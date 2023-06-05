import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve
import numpy as np
from sklearn import metrics
import pdb
pos_scores = pd.read_csv('../../data/dev/blackbox_cracking/all_pos_scores.csv', header=None)
neg_scores = pd.read_csv('../../data/dev/blackbox_cracking/all_neg_scores.csv', header=None)
outdir='../../data/plots/'
#Check overlap:
print('Overlap:', neg_scores[neg_scores[0].isin(pos_scores[0])])
print('Positive set:')
print(pos_scores.describe())
print('Negative set:')
print(neg_scores.describe())

#PR curve
y_true = np.zeros(len(pos_scores)+len(neg_scores))
y_true[:len(pos_scores)]=1
probas_pred = np.concatenate([pos_scores[3], neg_scores[3]])
precision, recall, thresholds = precision_recall_curve(y_true, probas_pred)
PR_AUC = metrics.auc(recall, precision)
print('PR AUC:', PR_AUC)
fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
plt.plot(recall, precision, label='AUC:'+str(np.round(PR_AUC,3)))
plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.5,1])
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(outdir+'pr_curve.png', format='png', dpi=300)
plt.close()

#ROC curve
fpr, tpr, thresholds = metrics.roc_curve(y_true, probas_pred)
AUC = metrics.auc(fpr, tpr)
fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
plt.plot(fpr, tpr, label='AUC:'+str(np.round(AUC,3)))
plt.legend()
plt.xlabel('FPR')
plt.ylabel('TPR')
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(outdir+'roc_curve.png', format='png', dpi=300)
plt.close()

#y pred
x, acc = [], []
for t in np.arange(0,max(probas_pred),0.1):
    y_pred = np.zeros(len(pos_scores)+len(neg_scores))
    y_pred[probas_pred>=t]=1
    #print('Accuracy:',t, metrics.accuracy_score(y_true, y_pred))
    x.append(t)
    acc.append(metrics.accuracy_score(y_true, y_pred))

fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
plt.plot(x, acc)
plt.xlabel('pDockQ threshold')
plt.ylabel('Accuracy')
plt.tight_layout()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(outdir+'acc_vs_pdockqs.png', format='png', dpi=300)
plt.close()

#Plot the comparison methods
#Accuracy
accuracies = [0.50, 0.50, 0.51, 0.51, 0.52, max(acc)]
labels = [ 'Richoux-LSTM', 'DeepFE', 'Richoux-FC', 'SPRINT', 'PIPR', 'AlphaFold']
fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
plt.scatter(range(len(accuracies)), accuracies)
plt.xticks(ticks=range(len(accuracies)), labels=labels, rotation=45)
plt.axhline(y=0.5, xmin=0, xmax=5, linestyle='--', color='grey', label='random')
plt.ylabel('Accuracy')
plt.title('Accuracies')
plt.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(outdir+'accuracies.png', format='png', dpi=300)
plt.close()
print('Accuracy:', max(acc))

#PR AUC
PR_AUCs = [0.50, 0.50, 0.51, 0.51, 0.52, PR_AUC]
labels = [ 'Richoux-LSTM', 'DeepFE', 'Richoux-FC', 'SPRINT', 'PIPR', 'AlphaFold']
fig, ax = plt.subplots(figsize=(9/2.54, 9/2.54))
plt.scatter(range(len(PR_AUCs)), PR_AUCs)
plt.xticks(ticks=range(len(PR_AUCs)), labels=labels, rotation=45)
plt.axhline(y=0.5, xmin=0, xmax=5, linestyle='--', color='grey', label='random')
plt.ylabel('PR AUC')
plt.title('PR AUCs')
plt.legend()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(outdir+'PR_AUCs.png', format='png', dpi=300)
plt.close()
