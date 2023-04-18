from pyqpanda import *

def run(i,j):
    QCM = QCloud()

    QCM.init_qvm("6425065D78304CE2A9EA7EF5252C18E5")
#    QCM.init_qvm("4EDE33E2E2F448CA8EC66A03885BA408")
#    QCM.init_qvm("F4B3CDC6400E43DFA0DED791D9545BFE")
#   QCM.init_qvm("BAFB044038EE4521A7815A2563ADE198")
#    QCM.init_qvm("4D077A6E21774F99A1E6415EF4D8AF85")

    qlist = QCM.qAlloc_many(5)
    clist = QCM.cAlloc_many(5)

    qvm = init_quantum_machine(QMachineType.CPU)
    qvm.init_qvm()

    prog_trans, qv, cv = convert_qasm_to_qprog("all_qasm_test/qasm_%d_%d.txt"%(i,j), qvm)
    try:
      result = QCM.real_chip_measure(prog_trans, 10000, real_chip_type.origin_wuyuan_d5)
    except:
      print("job %d %d failed" % (i,j))
      return
    value = result['00000']
    f = open("results_test/result_%d_%d.txt" %(i,j),'w')
    f.write(str(value))
    f.close()

if __name__=="__main__":

  for i in range(60,70):
    for j in range(100):
      run(i,j)
