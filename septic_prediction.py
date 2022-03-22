import os
os.environ["CUDA_VISIBLE_DEVICES"]="1"

import sys
import math
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from sklearn.metrics import precision_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import recall_score
from sklearn import metrics
from sklearn.model_selection import GridSearchCV
from collections import Counter

def create_missing_indicator(org, data, columns):
    for col in columns:
        selected = pd.isnull(org[col])
        data.loc[selected, col + "_missing"] = 1
        data.loc[~selected, col + "_missing"] = 0
    return data.loc[:,[c for c in data.columns.values if 'missing' in c]]

def create_model():
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.LSTM(128, input_dim=X_train.shape[2]))
    model.add(tf.keras.layers.Dense(128, activation='relu'))
    model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate,
        decay_steps=decay_steps,
        decay_rate=decay_rate,
        staircase=True)
    optimizer = tf.keras.optimizers.Adam(learning_rate = lr_schedule)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    early_stopping = tf.keras.callbacks.EarlyStopping(monitor='loss', verbose=1, mode='auto')
    return model


if __name__ == '__main__':
	dataset = sys.argv[1]
	missing_indicator = int(sys.argv[2])

	# path storing complete data
	path = './result/' + dataset + '_24h_imputation_result/'
	# path storing data with native missing values for MI
	org_path = './data/'+dataset+ '/data_with_missing/'
	# path storing prediction results
	save = './result/' + dataset + '_24h_prediction/'

	if dataset == 'cchs':
	    NUM = 2842
	elif dataset == 'mayo':
	    NUM = 6836
	elif dataset == 'mimic':
	    NUM = 772

	max_length = 70 # maximum sequence length

	# parameters of LSTM
	initial_learning_rate = 0.005
	decay_steps = 1000
	decay_rate = 0.96
	nodeN = 128

	df = pd.read_csv(path + str(1) + '.csv')
	org_df = pd.read_csv(org_path + str(1) + '.csv')
	visitIdx = [1 for d in range(len(df))]
	df['VisitIdx'] = visitIdx

	for i in range(2, NUM+1):
	    tmp_df = pd.read_csv(path + str(i) + '.csv')
	    tmp_org_df = pd.read_csv(org_path + str(i) + '.csv')
	    visitIdx = [i for d in range(len(tmp_df))]
	    tmp_df['VisitIdx'] = visitIdx
	    df = pd.concat([df, tmp_df], ignore_index=True)
	    org_df = pd.concat([org_df, tmp_org_df], ignore_index=True)

	numeric_columns = [c for c in df.columns.values if 'X' in c]
	mi_columns = create_missing_indicator(org_df, df, numeric_columns)
	print('Imputed data contains any NaN?', df[numeric_columns].isnull().values.any())
	new_columns = ['VisitIdx'] + numeric_columns
	if missing_indicator:
	    new_columns += [c for c in mi_columns.columns.values if 'missing' in c]
	df[numeric_columns] = (df[numeric_columns] - df[numeric_columns].mean())/df[numeric_columns].std()

	X = []
	y = []
	lengths = []

	# Separate input sequence X and target y
	targets = {}
	selected = (df.VisitIdx <= int(NUM/2))
	for v in df.loc[selected, :].VisitIdx.unique().tolist():
	    targets[str(int(v))] = 1
	selected = (df.VisitIdx > int(NUM/2))
	for v in df.loc[selected, :].VisitIdx.unique().tolist():
	    targets[str(int(v))] = 0

	sequences = {}
	for row in df.loc[:, new_columns].values:
	    if str(int(row[0])) not in sequences:
	        sequences[str(int(row[0]))] = []
	    sequences[str(int(row[0]))].append(row[1:])

	for key in sequences.keys():
	    X.append(sequences[key])
	    y.append(targets[key])
	    lengths.append(len(sequences[key]))

	print("mean length of seqs:", np.array(lengths).mean())
	print ('Number of septic shock/ non septic shock patients')
	print(Counter(y))
	X = np.array(X)
	y = np.array(y)


	X = tf.compat.v1.keras.preprocessing.sequence.pad_sequences(X, maxlen=max_length, padding = 'pre', truncating='pre',dtype = 'float64')
	X = np.reshape(X, (X.shape[0], max_length, X.shape[2]))
	print('X dim after padding:', X.shape)


	tf.compat.v1.reset_default_graph()
	tf.compat.v1.disable_eager_execution()
	np.random.seed(1000)

	cv_folds = 5

	kf = KFold(n_splits=cv_folds, shuffle=True)

	performances = pd.DataFrame(columns=['Accuracy', 'std' ,'Precision', 'std', 'Recall', 'std', 'F-measure', 'std', 'AUC', 'std'])
	accuracy, precision, recall, f_measure, auc_roc = [0] * cv_folds, [0] * cv_folds, [0] * cv_folds, [0] * cv_folds, [0] * cv_folds
	fold = 0
	for train_index, test_index in kf.split(X):
	    print('#'*15, "Cross-Validation fold", fold)

	    X_train, X_test = X[train_index], X[test_index]
	    y_train, y_test = y[train_index], y[test_index]

	    model = tf.keras.wrappers.scikit_learn.KerasClassifier(build_fn=create_model, verbose=1)
	    epochs = [25]
	    batches = [64]
	    param_grid = dict(epochs=epochs, batch_size=batches)
	    grid = GridSearchCV(estimator=model, param_grid=param_grid, cv=2)
	    # with tf.device('/gpu:2'):
	    grid_result = grid.fit(X_train, y_train)

	    print(grid_result.best_params_)

	    predict = grid_result.predict(X_test)
	    predict_prob = grid_result.predict_proba(X_test)
	    test_pred, test_pred_prob = [], []
	    for each in predict:
	        test_pred.append(each[0])
	    for each in predict_prob:
	        test_pred_prob.append(each[0])

	    accuracy[fold] = accuracy_score(y_test, test_pred)
	    recall[fold] =  recall_score(y_test, test_pred)
	    f_measure[fold] = f1_score(y_test, test_pred)
	    precision[fold] = precision_score(y_test, test_pred)
	    fpr, tpr, thresholds = metrics.roc_curve(y_test, test_pred_prob, pos_label=0)
	    auc_roc[fold] = metrics.auc(fpr, tpr)
	    print ("Confusion Matrix:")
	    print (confusion_matrix(y_test, test_pred))
	    fold += 1

	performances.loc[0] = [np.mean(accuracy), np.std(accuracy), np.mean(precision), np.std(precision),\
	                         np.mean(recall), np.std(recall), np.mean(f_measure), np.std(recall), np.mean(auc_roc), np.std(auc_roc)]
	performances.to_csv(save+'lstm_'+dataset + '_MI' + str(missing_indicator)+'_lr_' + str(initial_learning_rate) + "_" + str(decay_steps) + "_" + str(decay_rate) + '.csv', sep='\t')