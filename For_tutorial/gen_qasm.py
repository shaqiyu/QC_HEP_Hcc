import numpy as np

# Import Qiskit
from qiskit import QuantumCircuit
from qiskit import Aer, transpile
from qiskit.tools.visualization import plot_histogram, plot_state_city
import qiskit.quantum_info as qi
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from qiskit import Aer
from qiskit.circuit.library import PauliFeatureMap
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit_machine_learning.kernels import QuantumKernel
from qiskit.circuit import QuantumCircuit, ParameterVector
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from qiskit import QuantumRegister
from qiskit import Aer, BasicAer
import sklearn
from sklearn.metrics import normalized_mutual_info_score, roc_curve, auc, confusion_matrix, accuracy_score, roc_auc_score
from sklearn.preprocessing import label_binarize
import statistics
from sklearn.svm import SVC, LinearSVC
from sklearn.cluster import SpectralClustering
from sklearn.metrics import normalized_mutual_info_score
from sklearn.model_selection import train_test_split
from qiskit import Aer, BasicAer
from qiskit.circuit.library import ZZFeatureMap, TwoLocal
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit_machine_learning.algorithms import QSVC
from qiskit_machine_learning.kernels import QuantumKernel
from qiskit.aqua.utils.dataset_helper import get_feature_dimension
from Qhep_Modules_new import Plotter, Utilities, customised_feature_maps

import statistics
import math
import sys
rng = np.random.RandomState(0)
def gen_qasm():

  seed=2022
  #backend = BasicAer.get_backend("aer_simulator_statevector")
  #backend = Aer.get_backend("aer_simulator_statevector")
  backend = BasicAer.get_backend("statevector_simulator")
  feature_map_cus = customised_feature_maps.FeatureMap(num_qubits=6, depth=1, degree=1, entanglement='full', inverse=False)

  print("Custom feature map:\n", feature_map_cus)

  for i in range(1):
    seed = 2022
    algorithm_globals.random_seed = seed

    data = pd.read_csv("sample_%d.csv" % i)
    train = data[0:100]
    test = data[100:200]
    train_label = train.pop('tag')
    test_label = test.pop('tag')
    print (train)
    print (test)
    print (train_label)
    print (test_label)
    
    train_data = train.to_numpy()
    test_data = test.to_numpy()
    X_train = train_data
    X_test = test_data
   # test_label_oh = label_binarize(test_label, classes=[1,-1])
   # print (test_label_oh)

    #backend_options = {'max_parallel_threads':8}
    q_backend = QuantumInstance(backend, shots=10, seed_simulator=None, seed_transpiler=None)
    q_kernel = QuantumKernel(feature_map=feature_map_cus, quantum_instance=q_backend)
    
    qsvm_kernel_matrix_train = q_kernel.evaluate(x_vec=X_train)
    qsvm_kernel_matrix_test = q_kernel.evaluate(x_vec=X_test,y_vec=X_train)
    kernel_train_IBM = np.asmatrix(qsvm_kernel_matrix_train)
    kernel_test_IBM = np.asmatrix(qsvm_kernel_matrix_test)
    print("kernel_train_IBM:")
    print(kernel_train_IBM)
    print("kernel_test_IBM:")
    print(kernel_test_IBM)

    x = ParameterVector('x', length=6)
    y = ParameterVector('y', length=6)
    circuit = q_kernel.construct_circuit(x,y)
    circuit = q_backend.transpile(circuit)[0]


    for i in range(99):
      for j in range(i+1,100):
        f = open("all_qasm_train/qasm_%d_%d.txt"%(i,j),'w')
        cirtemp = circuit.assign_parameters({x:X_train[i], y:X_train[j]}, inplace=False)
        f.write(cirtemp.qasm())
        f.close()

if __name__ == "__main__":
  gen_qasm()
