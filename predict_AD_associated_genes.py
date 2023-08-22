import numpy as np
import pandas as pd
from tensorflow import keras
from ladder_net import get_ladder_network_fc
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from scipy.stats import sem

# run start
class_label = r'./data/AD_labels.csv'
official_name = r'./data/gene_ID_name.csv'

# get the dataset
print('Input data...')
x = pd.read_csv(r'./data/GO_PATHWAY_EXP_PPI_combined_data.csv', index_col="entrezId", sep=",", header=0, na_values=["?"])

# input class labels
classes = pd.read_csv(class_label, index_col="entrezId", sep=",", header=0, na_values=["?"])
x = x.join(classes, how='left', lsuffix="_lsuf")
x = x.fillna(0)
Y = x['class_label'].to_numpy()
x.drop('class_label', axis=1, inplace=True)

inp_size = x.shape[1]
n_classes = 2
gene = x.index
x = x.to_numpy()
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
AUC, ACC, AUC_samples, f1 = np.zeros(25), np.zeros(25), np.zeros(25), np.zeros(25)
training = 45
count = 0

for i, (train, test) in enumerate(kf.split(x, Y)):
    print('############')
    print(f'# fold {i+1}   #')
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

    np.random.seed(2023)
    idxs_annot_n = np.random.choice(X_train_U.shape[0], 150, replace=False)
    X_train_labeled = np.concatenate((X_train_P, X_train_U[idxs_annot_n]))
    y_train_labeled = np.concatenate((y_train_P, y_train_U[idxs_annot_n]), axis=None)

    y_train_labeled = np.where(y_train_labeled == 2, 0, y_train_labeled)
    y_train_labeled = keras.utils.to_categorical(y_train_labeled, n_classes)
    y_test = np.where(y_test == 2, 0, y_test)
    y_TEST = keras.utils.to_categorical(y_test, n_classes)

    n_rep = X_train_unlabeled.shape[0] // X_train_labeled.shape[0]
    X_train_labeled_rep = np.concatenate([X_train_labeled]*n_rep)
    addition = range(X_train_unlabeled.shape[0] - X_train_labeled_rep.shape[0])

    X_train_labeled_rep = np.concatenate((X_train_labeled_rep, X_train_labeled[addition]))
    y_train_labeled_rep = np.concatenate([y_train_labeled]*n_rep)
    y_train_labeled_rep = np.concatenate((y_train_labeled_rep, y_train_labeled[addition]))

    for run in range(5):
        # initialize the model
        model = get_ladder_network_fc(layer_sizes=[inp_size, 1000, 500, 250, 250, 250, n_classes])

        # train the model
        model.fit([X_train_labeled_rep, X_train_unlabeled], y_train_labeled_rep, epochs=training, verbose=2)

        for k in range(2):
            hist = model.fit([X_train_labeled_rep, X_train_unlabeled], y_train_labeled_rep, epochs=1, verbose=0)
            train_acc = hist.history['accuracy'][0]

            # retrieve testing performance
            y_test_pr = model.test_model.predict(X_test, batch_size=100)
            test_auc = roc_auc_score(y_TEST, y_test_pr, average='micro')
            test_acc = accuracy_score(y_TEST.argmax(-1), y_test_pr.argmax(-1))
            test_f1 = f1_score(y_TEST.argmax(-1), y_test_pr.argmax(-1), average='weighted')
            if k == 0:
                old_train_acc = train_acc
                old_auc = test_auc
                old_acc = test_acc
                old_f1 = test_f1
                y_all_pr = model.test_model.predict(sc.transform(x), batch_size=100)
            elif train_acc >= old_train_acc:
                old_train_acc = train_acc
                old_auc = test_auc
                old_acc = test_acc
                old_f1 = test_f1
                y_all_pr = model.test_model.predict(sc.transform(x), batch_size=100)
        print("Test accuracy: %f" % old_acc, "AUC: %f" % old_auc, "f1: %f" % old_f1)
        print("########################################################")
        print("########################################################")       
        ACC[count] = old_acc
        AUC[count] = old_auc
        f1[count] = old_f1
        if i == 0 and run == 0:
            All_probability = pd.DataFrame(gene)
            All_probability = All_probability.set_index('entrezId')
            All_probability[f"{count}"] = y_all_pr[:, 1]
        else:
            All_probability[f"{count}"] = y_all_pr[:, 1]
        count = count + 1

AD_genes = pd.DataFrame.median(All_probability, axis=1)
AD_genes = pd.Series.to_frame(AD_genes)
AD_genes = AD_genes[AD_genes[0] > 0.77]
geneName = pd.read_csv(official_name, index_col="entrezId", sep=",", header=0)
AD_genes = AD_genes.join(geneName, how='left', lsuffix="_lsuf").sort_index().drop(columns=0)

print('Performance summary:')
print("Average accuracy: %f " % np.average(ACC), "Standard Error: %f" % sem(ACC))
print("Average AUC: %f " % np.average(AUC), "Standard Error: %f" % sem(AUC))
print("Average f1: %f " % np.average(f1), "Standard Error: %f" % sem(f1))
AD_genes.to_csv(r"./results/Predicted_AD_associated_genes.csv")


