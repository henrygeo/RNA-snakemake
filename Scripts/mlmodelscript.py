import sklearn as skl
from sklearn import linear_model as lm
from sklearn import model_selection as ms
#from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
#from sklearn.feature_selection import VarianceThreshold
import pandas as pd
import matplotlib.pyplot as plt
#from sklearn.neural_network import MLPRegressor
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier
#from sklearn.ensemble import GradientBoostingRegressor
import numpy as np
import joblib

#from sklearn.decomposition import PCA
cpmdata = pd.read_csv("/scratch/henrygeo/RNAproject/cpmTMM.csv",  index_col=0)
#fielddata = pd.read_csv("/scratch/henrygeo/RNAproject/datasheets/rnaseqDatasheet.csv",  index_col=8)
#DEGdata = pd.read_csv("/scratch/henrygeo/RNAproject/DEGcountsfitness.csv",  index_col=0)
#fitness = pd.read_csv("/scratch/henrygeo/RNAproject/datasheets/fitness.csv",  index_col=0)

df = pd.DataFrame(DEGdata)
names = df.index
features = df.columns[5:]
data = df.values

#for cpmdata
cpmdata = pd.read_csv("/scratch/henrygeo/RNAproject/cpmTMM.csv",  index_col=0)
df = pd.DataFrame(cpmdata)
names = df.index

f = []
for i in range(len(df.Seednum)):
     if df.Seednum[i] > 0:
         f.append("Set")
     else:
         f.append("NS")

df['SeedSuc'] = f
features = df.columns[:-3]

dff = df.drop(df[df.Seednum>50].index)
data = dff.values


###model time
X, y = data[:,5:],data[:,4] #Seednum
X, y = data[:,5:],data[:,1] #Success
X, y = data[:,5:],data[:,0] #Group
X, y = data[:,5:],data[:,2] #Biomass

#for cpm data
 #Seednum
X, y = data[:,:-3],data[:,-2] #t_SN
X, y = data[:,:-3],data[:,-1] #Seed Success
X, y = data[:,:-3],data[:,-3] #Seednum

np.random.seed(525)
X_train, X_test, y_train, y_test = ms.train_test_split(X, y, test_size=0.40)
set = StandardScaler()
set.fit(X_train)
X_train = set.transform(X_train)
X_test = set.transform(X_test)


set.fit(y_train.reshape(-1,1))
y_train = set.transform(y_train.reshape(-1,1))
y_test = set.transform(y_test.reshape(-1,1))
y_train = y_train.ravel()
y_test = y_test.ravel()
 
param_grid = {
    'hidden_layer_sizes': [(150,100,50), (120,80,40), (100,50,20), (50,20,10)],
    'max_iter': [300, 600],
    'activation': [ 'logistic', 'relu'],
    'solver': ['sgd',  'lbfgs', 'adam'],
    'alpha': np.arange(0.1,4,0.1),
    'learning_rate': ['adaptive'],
}

mlp_reg = MLPRegressor()

search = ms. GridSearchCV(mlp_reg, param_grid, n_jobs= -1, cv=5, verbose = 3, scoring = 'neg_mean_squared')
search.fit(X_train, y_train)
print(search.best_params_) 
#{'activation': 'logistic', 'alpha': 0.5, 'hidden_layer_sizes': (150, 100, 50), 'learning_rate': 'constant', 'max_iter': 400, 'solver': 'sgd'} dtf, not a good fit
#{'activation': 'logistic', 'alpha': 0.5, 'hidden_layer_sizes': (100, 50, 20), 'learning_rate': 'adaptive', 'max_iter': 200, 'solver': 'adam'} #for biomass DEG
#{'activation': 'logistic', 'alpha': 1.0, 'hidden_layer_sizes': (50, 20, 10), 'learning_rate': 'adaptive', 'max_iter': 600, 'solver': 'adam'} #for seed number DEG

mlp_reg = MLPRegressor(activation ='logistic', alpha= 4.0, hidden_layer_sizes = (80, 50, 10), learning_rate= 'adaptive', max_iter= 600, solver= 'adam')
mlp_reg.fit(X_train, y_train)
y_pred = mlp_reg.predict(X_test)
skl.metrics.r2_score(y_test, y_pred)
skl.metrics.explained_variance_score(y_test, y_pred)

#r2 0.15, variance explained 0.18


#{'activation': 'logistic', 'alpha': 0.9, 'hidden_layer_sizes': (120, 80, 40), 'learning_rate': 'adaptive', 'max_iter': 600, 'solver': 'sgd'} set seed?
mlp_reg = MLPRegressor(activation ='logistic', alpha= 0.9, hidden_layer_sizes = (120, 80, 40), learning_rate= 'adaptive', max_iter= 600, solver= 'sgd')
mlp_reg.fit(X_train, y_train)
y_pred = mlp_reg.predict(X_test)
skl.metrics.r2_score(y_test, y_pred)
skl.metrics.explained_variance_score(y_test, y_pred)

df_temp = pd.DataFrame({'Actual': y_test, 'Predicted': y_pred})
df_temp.head()

#classifier

mlp_clf = MLPClassifier()
clfsearch = ms.GridSearchCV(mlp_clf, param_grid, n_jobs= -1, cv=5, verbose = 3)
clfsearch.fit(X_train, y_train)
print(clfsearch.best_params_) 
#{'activation': 'relu', 'alpha': 0.9, 'hidden_layer_sizes': (150, 100, 50), 'learning_rate': 'constant','max_iter': 400, 'solver': 'sgd'}
#{'activation': 'relu', 'alpha': 0.2, 'hidden_layer_sizes': (100, 50, 20), 'learning_rate': 'constant', 'max_iter': 600, 'solver': 'sgd'} for seedset
#{'activation': 'relu', 'alpha': 0.30000000000000004, 'hidden_layer_sizes': (50, 20, 10), 'learning_rate': 'adaptive', 'max_iter': 600, 'solver': 'sgd'} #for seed set
#{'activation': 'relu', 'alpha': 0.30000000000000004, 'hidden_layer_sizes': (120, 80, 40), 'learning_rate': 'adaptive', 'max_iter': 600, 'solver': 'sgd'} #second run
mlp_clf = MLPClassifier(activation ='logistic', alpha= 0.8, hidden_layer_sizes = (50, 10), learning_rate= 'adaptive', max_iter= 1000, solver= 'sgd')
mlp_clf.fit(X_train, y_train)
y_pred = mlp_clf.predict(X_test)

skl.metrics.accuracy_score(y_test, y_pred)
skl.metrics.precision_score(y_pred, y_test, pos_label = 'Set') #for identifying accurately if a plant will set seed
##0.67, with 0.81 for identifying if a plant will not set seed


display = skl.metrics.PrecisionRecallDisplay.from_estimator(
    mlp_clf, X_test, y_test, name="MLPClassifier"
)
_ = display.ax_.set_title("2-class Precision-Recall curve")
plt.show()

skl.metrics.ConfusionMatrixDisplay.from_estimator(mlp_clf, X_test, y_test)
plt.show()
#save model

joblib_mlp = "MLP_clf.pkl"
joblib.dump(mlp_clf, joblib_mlp)


##Code for function to get feature importance from MLP classifier

def get_feature_importance(j, n):
  s = skl.metrics.accuracy_score(y_test, y_pred) # baseline score
  total = 0.0
  for i in range(n):
    perm = np.random.permutation(range(X_test.shape[0]))
    X_test_ = X_test.copy()
    X_test_[:, j] = X_test[perm, j]
    y_pred_ = mlp_clf.predict(X_test_)
    s_ij = skl.metrics.accuracy_score(y_test, y_pred_)
    total += s_ij
  return s - total / n

f = []
for j in range(X_test.shape[1]):
  f_j = get_feature_importance(j, 200)
  f.append(f_j)

 
importance={'importance':f, 'gene':features}
importance=pd.DataFrame(importance)
importance=importance.sort_values(by='importance', ascending = False)
topimp=importance[importance.importance>0.07]

plt.bar(range(topimp.shape[0]), topimp.importance, alpha = 0.7, color = "r")
plt.xticks(ticks=range(topimp.shape[0]), labels = topimp.gene, rotation=45, ha="right")
plt.xlabel("Feature")
plt.ylabel("Importance")
plt.title("Feature importances")
plt.show()

importance.to_csv("importantgenesMLPclf.txt")

coefficients = mlp_clf.coefs_
importance = np.abs(coefficients)
importantgenes.txt = np.array(features)[importance > 0]



#### Classifier with SGD

param_grid = {
     'loss': ['squared_hinge'],
     'penalty':['l2', 'elasticnet'],
     'l1_ratio': np.arange(0,1,0.1),
     'alpha': np.arange(1,3,0.2),
     'learning_rate': ['optimal'],
     'eta0': np.arange(2,10,0.5)
}


sgd_clf = lm.SGDClassifier()
SGDsearch = ms.GridSearchCV(sgd_clf, param_grid, scoring = 'balanced_accuracy', n_jobs= -1, cv=5, verbose = 3)
SGDsearch.fit(X_train, y_train)
print(SGDsearch.best_params_) 

#{'alpha': 2.0, 'eta0': 9, 'l1_ratio': 0.5, 'learning_rate': 'constant', 'loss': 'squared_hinge', 'penalty': 'elasticnet'}
#{'alpha': 1.2, 'eta0': 2.0, 'l1_ratio': 0.6000000000000001, 'learning_rate': 'optimal', 'loss': 'squared_hinge', 'penalty': 'l2'}


#{'alpha': 0.4, 'eta0': 2, 'l1_ratio': 0.30000000000000004, 'learning_rate': 'constant', 'loss': 'squared_hinge', 'penalty': 'l2'} #total genes filtered/Success

#{'alpha': 2.2, 'eta0': 8, 'l1_ratio': 0.7000000000000001, 'learning_rate': 'optimal', 'loss': 'squared_hinge', 'penalty': 'elasticnet'}

sgd_clf = lm.SGDClassifier(alpha= 2, eta0= 9, l1_ratio= 0.5, learning_rate= 'optimal', loss= 'squared_hinge', penalty= 'elasticnet')
sgd_clf.fit(X_train, y_train)
y_pred = sgd_clf.predict(X_test)

sgd_clf.score(X_test, y_test)

skl.metrics.accuracy_score(y_test, y_pred)
skl.metrics.precision_score(y_test, y_pred, pos_label = 'Set')

skl.metrics.ConfusionMatrixDisplay.from_estimator(sgd_clf, X_test, y_test)


def get_feature_importance(j, n):
  s = skl.metrics.accuracy_score(y_test, y_pred) # baseline score
  total = 0.0
  for i in range(n):
    perm = np.random.permutation(range(X_test.shape[0]))
    X_test_ = X_test.copy()
    X_test_[:, j] = X_test[perm, j]
    y_pred_ = sgd_clf.predict(X_test_)
    s_ij = skl.metrics.precision_score(y_test, y_pred_, pos_label = 'Set')
    total += s_ij
  return s - total / n

f = []
for j in range(X_test.shape[1]):
  f_j = get_feature_importance(j, 200)
  f.append(f_j)

 
importance={'importance':f, 'gene':features}
importance=pd.DataFrame(importance)
importance=importance.sort_values(by='importance', ascending = False)
topimp=importance[importance.importance>0.248]

plt.bar(range(topimp.shape[0]), topimp.importance, alpha = 0.7, color = "r")
plt.xticks(ticks=range(topimp.shape[0]), labels = topimp.gene, rotation=45, ha="right")
plt.xlabel("Feature")
plt.ylabel("Importance")
plt.title("Feature importances")
plt.show()

importance.to_csv("importantgenesSGDclfALL.txt")

skl.metrics.ConfusionMatrixDisplay.from_estimator(sgd_clf, X_test, y_test)
plt.show()

#GB classifier
param_grid = {
     'loss': ['log_loss', 'deviance', 'exponential'],
     'n_estimators': np.arange(100,1000,100),
     'criterion': ['friedman_mse'],
     'max_depth': np.arange(2,10,1),
     'learning_rate': np.arange(0,2,0.2),
}


GB_clf = GradientBoostingClassifier()
GBsearch = ms.GridSearchCV(GB_clf, param_grid, n_jobs= -1, cv=5, verbose = 3)
GBsearch.fit(X_train, y_train)
print(GBsearch.best_params_) 
#{'criterion': 'friedman_mse', 'learning_rate': 0.6000000000000001, 'loss': 'log_loss', 'max_depth': 2, 'n_estimators': 500}



GB_clf = GradientBoostingClassifier(alpha= 2, eta0= 9, l1_ratio= 0.5, learning_rate= 'optimal', loss= 'squared_hinge', penalty= 'elasticnet')
GB_clf.fit(X_train, y_train)
y_pred = GB_clf.predict(X_test)

sgd_clf.score(X_test, y_test)

skl.metrics.accuracy_score(y_test, y_pred)
skl.metrics.precision_score(y_test, y_pred, pos_label = 'Set')


#ELastic net
param_grid = {
    'alpha': np.arange(0.1,10,0.1),
    'l1_ratio': np.arange(0,1,0.01),
}

enet =  lm.ElasticNet()


#Use cv to find alpha and l1 ratio for EN
search = ms.GridSearchCV(enet, param_grid, n_jobs= -1, cv = 10, verbose=3)

search.fit(X_train,y_train)

search.best_params_
search.best_score_

#alpha = 9.9

alpha, l1rat = 9.9, 0

enetfit = lm.ElasticNet( alpha = 9.9, l1_ratio =0)

enetfit.fit(X_train, y_train)

y_pred  = enetfit.predict(X_test)
skl.metrics.r2_score(y_test, y_pred)
skl.metrics.explained_variance_score(y_test, y_pred)


def get_feature_importance(j, n):
  s = skl.metrics.explained_variance_score(y_test, y_pred) # baseline score
  total = 0.0
  for i in range(n):
    perm = np.random.permutation(range(X_test.shape[0]))
    X_test_ = X_test.copy()
    X_test_[:, j] = X_test[perm, j]
    y_pred_ = enetfit.predict(X_test_)
    s_ij = skl.metrics.explained_variance_score(y_test, y_pred_)
    total += s_ij
  return s - total / n

f = []
for j in range(X_test.shape[1]):
  f_j = get_feature_importance(j, 200)
  f.append(f_j)

importance={'importance':f, 'gene':features}
importance=pd.DataFrame(importance)
importance=importance.sort_values(by='importance', ascending = False)
topimp=importance[importance.importance>0.0]

plt.bar(range(topimp.shape[0]), topimp.importance, alpha = 0.7, color = "r")
plt.xticks(ticks=range(topimp.shape[0]), labels = topimp.gene, rotation=45, ha="right")
plt.xlabel("Feature")
plt.ylabel("Importance")
plt.title("Feature importances")
plt.show()


importance.to_csv("importantgenesenet.csv")

### Gradient Boosting regression

param_grid = {
     'loss': ['squared_error', 'absolute_error'],
     'min_samples_split':[2, 4],
     'max_depth': np.arange(1,8,1),
     'learning_rate': np.arange(0.01,1,0.2),
     'n_estimators': [100,200,500],
     'random_state': [525]
}


GBreg = GradientBoostingRegressor()
GBsearch = ms.GridSearchCV(GBreg, param_grid, scoring = 'neg_mean_squared_error', n_jobs= -1, cv=5, verbose = 3)
GBsearch.fit(X_train, y_train)

#{'alpha': 0.2, 'learning_rate': 0.2, 'loss': 'huber', 'max_depth': 1, 'min_samples_split': 2, 'n_estimators': 100, 'random_state': 525}

GBreg = GradientBoostingRegressor(random_state = 525, n_estimators = 100)
GBreg.fit(X_train, y_train)
y_pred = GBreg.predict(X_test)
skl.metrics.r2_score(y_test, y_pred)
skl.metrics.explained_variance_score(y_test, y_pred)

importance={'importance':GBreg.feature_importances_, 'gene':features}
importance=pd.DataFrame(importance)
importance=importance.sort_values(by='importance', ascending = False)
topimp=importance[importance.importance>0.01]
plt.bar(range(topimp.shape[0]), topimp.importance, alpha = 0.7, color = "r")
plt.xticks(ticks=range(topimp.shape[0]), labels = topimp.gene, rotation=45, ha="right")
plt.xlabel("Feature")
plt.ylabel("Importance")
plt.title("Feature importances")
plt.show()

importance.to_csv("importantgenesgradboostALL.csv")

importance.quantile(q=0.80, numeric_only = True)

### command line to get ref gene ids
grep -E -f "/scratch/henrygeo/RNAproject/Results/importantTrans.txt" "/scratch/henrygeo/RNAproject/Results/stringz/gtfs/merged.gtf"| awk -F" " '!_[$1]++' >refmatch.txt
 seqtk subseq "/scratch/henrygeo/RNAproject/HederaceaGenes.fasta" "/scratch/henrygeo/TopPCAregGenes.txt" > PCRgenes.fa
## Next is subsetting data with just the informative genes and rerunning
## Also maybe using transcript abundance instead of gene abundance
