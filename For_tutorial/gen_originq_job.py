from pyqpanda import *

def run(i,j):
    QCM = QCloud()
#    QCM.init_qvm("4DCD448B223A488E878E9654FEA0C619")   
    QCM.init_qvm("58DCD70A14814A4DB73ACF5C8F854FA2")
#    QCM.init_qvm("DD62C839BDE1406C81D4337039494970")   #shaqiyu20 e-mail
#    QCM.init_qvm("AEF7CD7BE2744F9D90D87A57AB09D68A")   #180 e-mail
#    QCM.init_qvm("645B9D74C11B4F999D1C2F9D12471E33")   #163 e-mail

    qlist = QCM.qAlloc_many(5)
    clist = QCM.cAlloc_many(5)

    qvm = init_quantum_machine(QMachineType.CPU)
    qvm.init_qvm()

    prog_trans, qv, cv = convert_qasm_to_qprog("all_qasm_train/qasm_%d_%d.txt"%(i,j), qvm)
    try:
      #result = QCM.real_chip_measure(prog_trans, 10000, real_chip_type.origin_wuyuan_d5)
      result = QCM.real_chip_measure(prog_trans, 10000, real_chip_type.origin_wuyuan_d5)
    except:
      print("job %d %d failed" % (i,j))
      return
    value = result['00000']
    f = open("results_train_try/result_%d_%d.txt" %(i,j),'w')
    f.write(str(value))
    f.close()

if __name__=="__main__":

  for i in range(19,20):
    for j in range(i+1,100):
      run(i,j)
