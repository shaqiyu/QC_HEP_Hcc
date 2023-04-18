import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sklearn
from sklearn.svm import SVC
from sklearn.metrics import normalized_mutual_info_score, roc_curve, auc, confusion_matrix, accuracy_score, roc_auc_score
from sklearn.preprocessing import label_binarize

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, IBMQ, execute, transpile, Aer, assemble
import statistics
from sklearn.svm import SVC, LinearSVC
from sklearn.cluster import SpectralClustering
from sklearn.model_selection import GridSearchCV, cross_val_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.model_selection import train_test_split
from qiskit import Aer, BasicAer
from qiskit.circuit.library import ZZFeatureMap, TwoLocal
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit_machine_learning.algorithms import QSVC, PegasosQSVC
from qiskit_machine_learning.kernels import QuantumKernel
from qiskit.aqua.utils.dataset_helper import get_feature_dimension
from Qhep_Modules_new import Plotter, Utilities, customised_feature_maps
from sklearn.preprocessing import StandardScaler

rng = np.random.RandomState(0)

def get_kernel(i,j,Type):
  if Type == 'test':
    fn = 'results_test/result_%d_%d.txt' % (i,j)
  elif Type == 'train':
    fn = 'results_train/result_%d_%d.txt' % (i,j)
  else:
    return 0.0

  f = open(fn,'r')
  value = float(f.readline())
  f.close()
  return value

def eval_model(c):

  score = []
  for i in range(1):

    data = pd.read_csv("sample_%d.csv" % i)
    train = data[0:100]
    test = data[100:200]
    train_label = train.pop('tag')
    test_label = test.pop('tag')
    train_label_oh = label_binarize(train_label, classes=[1,-1])
    test_label_oh = label_binarize(test_label, classes=[1,-1])

    test_kernel = []
    for i in range(100):
      test_kernel_line = []
      for j in range(100):
        value = get_kernel(i,j,'test')
        test_kernel_line.append(value)
      test_kernel.append(test_kernel_line)
    test_kernel_array = np.array(test_kernel, np.float32)

    train_kernel = []
    for i in range(100):
      train_kernel_line = []
      for j in range(100):
        if i == j:
          train_kernel_line.append(1.0)
        else:
          if i>j:
            # the matrix is symmetric
            # switch i and j as j is always greater than i in data
            value = get_kernel(j,i,'train')
            train_kernel_line.append(value)
          else:
            value = get_kernel(i,j,'train')
            train_kernel_line.append(value)
      train_kernel.append(train_kernel_line)
    train_kernel_array = np.array(train_kernel, np.float32)
    print("train_kernel_array:")
    print(train_kernel_array)
    print("test_kernel_array:")
    print(test_kernel_array)

#Classical SVM(SVC)
    train_data = train.to_numpy()
    test_data = test.to_numpy()
    train_array = np.array(train_data, np.float32)
    test_array = np.array(test_data, np.float32)
    X = np.concatenate([train_array, test_array], axis=0)
    y = np.concatenate([train_label, test_label])
    X_train = train_data
    X_test = test_data
    scaler = StandardScaler()
    X_train_svm = scaler.fit_transform(X_train)
    X_test_svm = scaler.transform(X_test)
    y_train = train_label
    y_test = test_label

    svc_svm = SVC(kernel='rbf',probability=True, C=10, gamma=0.1)   
    clf = svc_svm.fit(X_train_svm, y_train)
    y_score = clf.decision_function(X_test_svm)
    fpr_svc, tpr_svc, _ = sklearn.metrics.roc_curve(y_test, y_score)
    
#simulate
    seed=2022
    backend = BasicAer.get_backend("statevector_simulator")
    feature_map_cus = customised_feature_maps.FeatureMap(num_qubits=6, depth=1, degree=1, entanglement='full', inverse=False)
    #print("Custom feature map:\n", feature_map_cus)
    quantum_instance = QuantumInstance(backend, shots=10, seed_simulator=1000, seed_transpiler=seed)
    qkernel = QuantumKernel(feature_map=feature_map_cus, quantum_instance=quantum_instance)
    qsvm_kernel_matrix_train = qkernel.evaluate(x_vec=X_train)
    qsvm_kernel_matrix_test = qkernel.evaluate(x_vec=X_test,y_vec=X_train)
    kernel_train_IBM = np.asmatrix(qsvm_kernel_matrix_train)
    kernel_test_IBM = np.asmatrix(qsvm_kernel_matrix_test)
    kernel_train_IBM_array = np.array(qsvm_kernel_matrix_train, np.float32)
    kernel_test_IBM_array = np.array(qsvm_kernel_matrix_test, np.float32)
    #print("kernel_train_IBM:")
    #print(kernel_train_IBM)
    #print("kernel_test_IBM:")
    #print(kernel_test_IBM)
    
    #qsvc = QSVC(quantum_kernel=qkernel, probability=True, C=5)
    #qsvc = PegasosQSVC(probability=True, C=5, precomputed=True)
    qsvc = SVC(C=30, probability=True, kernel="precomputed")

    clfq = qsvc.fit(kernel_train_IBM_array, train_label)
    predictions_IBM = clfq.predict_proba(qsvm_kernel_matrix_test)
    fpr_qsvc, tpr_qsvc, _ = sklearn.metrics.roc_curve(test_label_oh, predictions_IBM[:, 0] )

#IBM hardware
    kernel_train_ibm_array = np.loadtxt('hardware_ibm_train_tran.txt', dtype=np.float32)
    qsvm_kernel_matrix_test_ibm = np.loadtxt('hardware_ibm_test_tran.txt', dtype=np.float32)

    qsvc_ibm = SVC(kernel='precomputed', probability=True, C=1)

    clfq_ibm = qsvc_ibm.fit(kernel_train_ibm_array, train_label)
    predictions_ibm = clfq_ibm.predict_proba(qsvm_kernel_matrix_test_ibm)
    fpr_qsvc_ibm, tpr_qsvc_ibm, _ = sklearn.metrics.roc_curve(test_label_oh, predictions_ibm[:, 0] )

#QSVM
    qsvc_origin = SVC(kernel="precomputed", C=10, probability=True )
    clfq_origin = qsvc_origin.fit(train_kernel_array, train_label)
    predictions_origin = clfq_origin.predict_proba(test_kernel)
    fpr_origin, tpr_origin, _ = sklearn.metrics.roc_curve(test_label_oh, predictions_origin[:, 0])

    
    roc_auc_origin = sklearn.metrics.auc(fpr_origin, tpr_origin) #Origin hardware
    roc_auc_svc = sklearn.metrics.auc(fpr_svc, tpr_svc)  
    roc_auc_qsvc = sklearn.metrics.auc(fpr_qsvc, tpr_qsvc) #simulator
    roc_auc_qsvc_ibm = sklearn.metrics.auc(fpr_qsvc_ibm, tpr_qsvc_ibm)
    std_qsvc =0
    std_svc =0


    print("SVM:")
    print(roc_auc_svc)
    print("origin hardware:")
    print(roc_auc_origin)
    print("IBM hardware:")
    print(roc_auc_qsvc_ibm)
    print("simulator:")
    print(roc_auc_qsvc)



    Plotter.plotROCcurve(tpr_qsvc, fpr_qsvc, roc_auc_qsvc, std_qsvc, tpr_origin, fpr_origin, roc_auc_origin, std_svc, tpr_qsvc_ibm, fpr_qsvc_ibm, roc_auc_qsvc_ibm, "ROC_scan")

  print("C: %d" %c)


if __name__ == "__main__":
    eval_model(10)
