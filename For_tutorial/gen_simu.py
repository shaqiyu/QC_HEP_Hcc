from pyqpanda import *
from numpy import pi


def run(i,j):

    qvm = MPSQVM()
    qvm.init_qvm()
    qvm.set_configure(64, 64)

    prog_trans, qv, cv = convert_qasm_to_qprog("all_qasm_train/qasm_%d_%d.txt"%(i,j), qvm)
#    print(prog_trans)
    try:

        result1 = qvm.prob_run_dict(prog_trans, qv, -1)
    except:
      print("job %d %d failed" % (i,j))
      return
    value = result1['000000']
    f = open("results_train/result_%d_%d.txt" %(i,j),'w')
    f.write(str(value))
    f.close()

if __name__=="__main__":

  for i in range(100):
    for j in range(i+1,100):
      run(i,j)
