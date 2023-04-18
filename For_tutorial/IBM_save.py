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
    print("Custom feature map:\n", feature_map_cus)
    quantum_instance = QuantumInstance(backend, shots=10, seed_simulator=1000, seed_transpiler=seed)
    qkernel = QuantumKernel(feature_map=feature_map_cus, quantum_instance=quantum_instance)
    qsvm_kernel_matrix_train = qkernel.evaluate(x_vec=X_train)
    qsvm_kernel_matrix_test = qkernel.evaluate(x_vec=X_test,y_vec=X_train)
    kernel_train_IBM = np.asmatrix(qsvm_kernel_matrix_train)
    kernel_test_IBM = np.asmatrix(qsvm_kernel_matrix_test)
    kernel_train_IBM_array = np.array(qsvm_kernel_matrix_train, np.float32)
    kernel_test_IBM_array = np.array(qsvm_kernel_matrix_test, np.float32)
    print("kernel_train_IBM:")
    print(kernel_train_IBM)
    print("kernel_test_IBM:")
    print(kernel_test_IBM)
    
    #qsvc = QSVC(quantum_kernel=qkernel, probability=True, C=5)
    #qsvc = PegasosQSVC(probability=True, C=5, precomputed=True)
    qsvc = SVC(C=30, probability=True, kernel="precomputed")

    clfq = qsvc.fit(kernel_train_IBM_array, train_label)
    predictions_IBM = clfq.predict_proba(qsvm_kernel_matrix_test)
    fpr_qsvc, tpr_qsvc, _ = sklearn.metrics.roc_curve(test_label_oh, predictions_IBM[:, 0] )

#IBM hardware
    IBMQ.save_account("75685791c48132479d6c609db91a6f636051da392dba2f4dde1190c8b1f222aead93df9c00bb30fb1df78a7789d59b52844fefba66d0d3681683dcad95080e3f", overwrite=True)
    provider = IBMQ.load_account()

    provider = IBMQ.get_provider(hub='ibm-q', group='open', project='main')

    feature_map_cus_ibm = customised_feature_maps.FeatureMap(
        num_qubits=6, depth=1, degree=1, entanglement='full', inverse=False)

    print("Custom feature map for hardware:\n", feature_map_cus_ibm)
    hard_name = 'ibm_oslo' # 5 qubits: 'ibmq_manila', 'ibmq_santiago', 'ibmq_belem', 'ibmq_lima', 'ibmq_quito'; 7 qubits: 'ibm_oslo', 'ibm_nairobi'
    backend_ibm_hardware = provider.get_backend(hard_name) 

    mapped_circuit = transpile(feature_map_cus_ibm, backend=backend_ibm_hardware)
    quantum_instance_ibm = QuantumInstance(backend_ibm_hardware, shots=1024) 
    qkernel_ibm = QuantumKernel(feature_map=mapped_circuit, quantum_instance=quantum_instance_ibm) 
    qsvm_kernel_matrix_train_ibm = qkernel_ibm.evaluate(x_vec=X_train)
    qsvm_kernel_matrix_test_ibm = qkernel_ibm.evaluate(x_vec=X_test,y_vec=X_train)
    kernel_train_ibm = np.asmatrix(qsvm_kernel_matrix_train_ibm)
    kernel_test_ibm = np.asmatrix(qsvm_kernel_matrix_test_ibm)
    kernel_train_ibm_array = np.array(qsvm_kernel_matrix_train_ibm, np.float32)
    kernel_test_ibm_array = np.array(qsvm_kernel_matrix_test_ibm, np.float32)
    print("kernel_train_IBM_hardware:")
    print(kernel_train_ibm)
    print("kernel_test_IBM_hardware:")
    print(kernel_test_ibm)
    
    np.set_printoptions(threshold=np.inf)
    log_file_train = open("./hardware_ibm_train.txt","a")
    print("qsvm_kernel_train_IBM:",qsvm_kernel_matrix_train_ibm, file=log_file_train)
    log_file_train1 = open("./hardware_ibm_train_as.txt","a")
    print("kernel_train_IBM:",kernel_train_ibm, file=log_file_train1)
    log_file_test = open("./hardware_ibm_test.txt","a")
    print("qsvm_kernel_test_IBM:",qsvm_kernel_matrix_test_ibm, file=log_file_test)
    log_file_test1 = open("./hardware_ibm_test_as.txt","a")
    print("kernel_test_IBM:",kernel_test_ibm, file=log_file_test1)

    qsvc_ibm = SVC(kernel='precomputed', probability=True, C=30)

    clfq_ibm = qsvc.fit(kernel_train_ibm_array, train_label)
    predictions_ibm = clfq_ibm.predict_proba(qsvm_kernel_matrix_test_ibm)
    fpr_qsvc_ibm, tpr_qsvc_ibm, _ = sklearn.metrics.roc_curve(test_label_oh, predictions_ibm[:, 0] )

#QSVM
    svc = SVC(C=30, probability=True, kernel="precomputed")
    #svc = QSVC(C=5, probability=True, quantum_kernel="precomputed")
    csvc = svc.fit(train_kernel_array, train_label)
    predictions = csvc.predict_proba(test_kernel)
    fpr, tpr, _ = sklearn.metrics.roc_curve(test_label_oh, predictions[:, 0])

    f = open("roc_1.txt", 'w')
    j = 0
    while j < len(fpr):
      f.write(str(fpr[j]) + ' ' + str(tpr[j]) + '\n')
      j = j + 1
    f.close()

    
    roc_auc = sklearn.metrics.auc(fpr, tpr)
    roc_auc_svc = sklearn.metrics.auc(fpr_svc, tpr_svc)
    roc_auc_qsvc = sklearn.metrics.auc(fpr_qsvc, tpr_qsvc)
    roc_auc_qsvc_ibm = sklearn.metrics.auc(fpr_qsvc_ibm, tpr_qsvc_ibm)

    score.append(roc_auc)
    print(roc_auc)
    print(roc_auc_svc)
    print(roc_auc_qsvc)
    print(roc_auc_qsvc_ibm)

    plt.figure(figsize=(5, 5))
    ax = plt.gca()
    plt.plot([0, 1], [0, 1], 'r--')
    plt.plot(tpr,1 - fpr, 'b', label='Wuyuan Hardware(AUC = %0.3f)' %(roc_auc))
    plt.plot(tpr_qsvc,1 - fpr_qsvc, 'r--', label='IBM Simulator(AUC = %0.3f)' %(roc_auc_qsvc))
    plt.plot(tpr_qsvc_ibm,1 - fpr_qsvc_ibm, 'y', label='IBM Hardware(AUC = %0.3f)' %(roc_auc_qsvc_ibm))
    plt.plot(tpr_svc,1 - fpr_svc, 'k', label='Classical SVM(AUC = %0.3f)' %(roc_auc_svc))
    plt.xlabel('Signal efficiency')
    plt.ylabel('Background rejection')
    #plt.title('ROC curve')
    ax.set_xlabel('Signal efficiency')
    ax.tick_params(axis="x", direction='in', length=6)
    ax.tick_params(axis="y", direction='in', length=6)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.legend(loc='best')
    plt.savefig('ROC.pdf')



  print("C: %d" %c)
  print(score)

if __name__ == "__main__":
    eval_model(10)
