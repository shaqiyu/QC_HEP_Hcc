#python3 bin/classification.py -n 100 -g 1 -b sim 
#python3 bin/classification.py -n 200 -g 1 -b sim 
#python3 bin/classification.py -n 10000 -g 1 -b sim 
#python3 bin/classification.py -n 18000 -g 1 -b sim 
#python3 bin/classification.py -n 20000 -g 1 -b sim 
#python3 bin/classification.py -n 600 -g 1 -b sim 
#python3 bin/classification.py -n 800 -g 1 -b sim 
#python3 bin/classification.py -n 1000 -g 1 -b sim
#python3 bin/classification.py -n 2000 -g 1 -b sim
#python3 bin/classification.py -n 4000 -g 1 -b sim
#python3 bin/classification.py -n 6000 -g 1 -b sim
#python3 bin/classification.py -n 24000 -b sim 
python3 bin/classification.py -n 2000 -b sim |tee ./log/Qsvm_2000_loop
python3 bin/classification.py -n 4000 -b sim |tee ./log/Qsvm_4000_loop
python3 bin/classification.py -n 6000 -b sim |tee ./log/Qsvm_6000_loop
python3 bin/classification.py -n 8000 -b sim |tee ./log/Qsvm_8000_loop
#python3 bin/classification.py -n 10000 -b sim |tee ./log/Qsvm_10000_loop
#python3 bin/classification.py -n 12000 -b sim |tee ./log/svm_12000_loop
#python3 bin/classification.py -n 14000 -b sim |tee ./log/svm_14000_loop
#python3 bin/classification.py -n 16000 -b sim |tee ./log/svm_16000_loop
#python3 bin/classification.py -n 18000 -b sim |tee ./log/svm_18000_loop
#python3 bin/classification.py -n 20000 -b sim |tee ./log/svm_20000_loop
#python3 bin/classification.py -n 2000 -b noQSVM |tee ./log/svm_2000_loop
#python3 bin/classification.py -n 4000 -b noQSVM |tee ./log/svm_4000_loop
#python3 bin/classification.py -n 6000 -b noQSVM |tee ./log/svm_6000_loop
#python3 bin/classification.py -n 8000 -b noQSVM |tee ./log/svm_8000_loop
#python3 bin/classification.py -n 10000 -b noQSVM |tee ./log/svm_10000_loop
#python3 bin/classification.py -n 12000 -b noQSVM |tee ./log/svm_12000_loop
#python3 bin/classification.py -n 14000 -b noQSVM |tee ./log/svm_14000_loop
#python3 bin/classification.py -n 16000 -b noQSVM |tee ./log/svm_16000_loop
#python3 bin/classification.py -n 18000 -b noQSVM |tee ./log/svm_18000_loop
#python3 bin/classification.py -n 20000 -b noQSVM |tee ./log/svm_20000_loop