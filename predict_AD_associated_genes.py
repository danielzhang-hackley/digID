import numpy as np
import pandas as pd
import os
from tensorflow import keras
from ladder_net import get_ladder_network_fc
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from statistics import stdev

# run start
data_dir = 'data directory'
output = 'result directory'

features = 'path/to/your/patient/and/knowledge/data.csv'
class_label = 'AD_labels.csv'

# get the dataset
print('Input data...')
x = pd.read_csv(os.path.join(data_dir, features), index_col="entrezId", sep=",", header=0, na_values=["?"])

# input class labels
classes = pd.read_csv(os.path.join(data_dir, class_label), index_col="entrezId", sep=",", header=0, na_values=["?"])
x = x.join(classes, how='left', lsuffix="_lsuf")
x = x.fillna(0)
Y = x['class_label'].to_numpy()
x.drop('class_label', axis=1, inplace=True)

keras.utils.set_random_seed(1337)
np.random.seed(2023)

inp_size = x.shape[1]
n_classes = 2
gene = x.index
x = x.to_numpy()
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)

AUC, ACC, specificity, sensitivity = [], [], [], []
count = 1
for i, (train, test) in enumerate(kf.split(x, Y)):
    print('############')
    print(f'#  fold {i+1}  #')
    print('############')
    sc = StandardScaler()
    X_train = sc.fit_transform(x[train])
    X_test = sc.transform(x[test])
    X_test = X_test[np.where(Y[test] > 0)]
    y_test = Y[test][np.where(Y[test] > 0)]
    X_train_unlabeled = X_train

    X_train_P = X_train[np.where(Y[train] == 1)]
    temp = np.where(Y[train] == 2, 0, Y[train])
    X_train_U = X_train[np.where(temp == 0)]
    y_train_P = Y[train][np.where(Y[train] == 1)]
    y_train_U = Y[train][np.where(temp == 0)]

    idxs_annot_n = np.random.choice(X_train_U.shape[0], 95, replace=False)
    X_train_labeled = np.concatenate((X_train_P, X_train_U[idxs_annot_n]))
    y_train_labeled = np.concatenate((y_train_P, y_train_U[idxs_annot_n]), axis=None)

    y_train_labeled = np.where(y_train_labeled == 2, 0, y_train_labeled)
    y_train_labeled = keras.utils.to_categorical(y_train_labeled, n_classes)
    y_test = np.where(y_test == 2, 0, y_test)
    y_test = keras.utils.to_categorical(y_test, n_classes)

    n_rep = X_train_unlabeled.shape[0] // X_train_labeled.shape[0]
    X_train_labeled_rep = np.concatenate([X_train_labeled]*n_rep)
    addition = range(X_train_unlabeled.shape[0] - X_train_labeled_rep.shape[0])

    X_train_labeled_rep = np.concatenate((X_train_labeled_rep, X_train_labeled[addition]))
    y_train_labeled_rep = np.concatenate([y_train_labeled]*n_rep)
    y_train_labeled_rep = np.concatenate((y_train_labeled_rep, y_train_labeled[addition]))

    # initialize the model
    model = get_ladder_network_fc(layer_sizes=[inp_size, 1000, 500, 250, 250, 250, n_classes])
    model.fit([X_train_labeled_rep, X_train_unlabeled], y_train_labeled_rep, epochs=36, verbose=2, batch_size=100)

    y_test_pr = model.test_model.predict(X_test, batch_size=100)
    test_auc = roc_auc_score(y_test, y_test_pr, average='micro')
    test_acc = accuracy_score(y_test.argmax(-1), y_test_pr.argmax(-1))
    tn, fp, fn, tp = confusion_matrix(y_test.argmax(-1), y_test_pr.argmax(-1)).ravel()
    test_sen = tp / (tp+fn)
    test_spe = tn / (tn+fp)
    print("Test accuracy: %f" % test_acc, "AUC: %f" % test_auc, "specificity: %f" % test_spe, "sensitivity: %f" % test_sen)

    y_all_pr = model.test_model.predict(sc.transform(x), batch_size=100)
    y_all_pr = np.where(y_all_pr[:, 1] <= 0.75, 0, y_all_pr[:, 1])
    y_all_pr = np.where(y_all_pr > 0.75, 1, y_all_pr)
    if i == 0:
        All_prob = pd.DataFrame(gene)
        All_prob = All_prob.set_index('entrezId')
        All_prob[f"{count}"] = y_all_pr
    else:
        All_prob[f"{count}"] = y_all_pr
    count = count + 1

    ACC.append(test_acc)
    AUC.append(test_auc)
    specificity.append(test_spe)
    sensitivity.append(test_sen)

AD_genes = pd.DataFrame.sum(All_prob, axis=1)
AD_genes = pd.DataFrame(AD_genes, columns=['sum'])
AD_genes = AD_genes[AD_genes['sum'] >= 4]

print('Performance summary:')
print("Average accuracy: %f " % np.average(ACC), "Standard deviation: %f" % stdev(ACC))
print("Average AUC: %f " % np.average(AUC), "Standard deviation: %f" % stdev(AUC))
print("Average specificity: %f " % np.average(specificity), "Standard deviation: %f" % stdev(specificity))
print("Average sensitivity: %f " % np.average(sensitivity), "Standard deviation: %f" % stdev(sensitivity))

AD_genes.to_csv(output + 'predicted_AD_genes.csv')
