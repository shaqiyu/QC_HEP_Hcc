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
  nqubits=6
  backend = BasicAer.get_backend("statevector_simulator")
  #backend = Aer.get_backend("aer_simulator_statevector")
  feature_map_cus = customised_feature_maps.FeatureMap(num_qubits=6, depth=1, degree=1, entanglement='full', inverse=False)

  print(feature_map_cus)

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
   # test_label_oh = label_binarize(test_label, classes=[1,-1])
   # print (test_label_oh)

    train_data = train.to_numpy()
    test_data = test.to_numpy()
    X_train = train_data
    X_test = test_data
    #backend_options = {'max_parallel_threads':8}
    
    q_backend = QuantumInstance(backend, shots=10, seed_simulator=seed, seed_transpiler=seed)
    q_kernel = QuantumKernel(feature_map=feature_map_cus, quantum_instance=q_backend)

    x = ParameterVector('x', length=nqubits)
    y = ParameterVector('y', length=nqubits)
    circuit = q_kernel.construct_circuit(x,y)
    circuit = q_backend.transpile(circuit)[0]


    for i in range(100):
      for j in range(100):
        f = open("all_qasm_test/qasm_%d_%d.txt"%(i,j),'w')
        cirtemp = circuit.assign_parameters({x:test_data[i], y:train_data[j]}, inplace=False)
        f.write(cirtemp.qasm())
        f.close()

if __name__ == "__main__":
  gen_qasm()
