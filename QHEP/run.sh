#-------
#Test
#--------
#python bin/Optclassification.py  -e Parforward  -n 2000 -d 1
#python bin/Eventclassification.py Parforward -d 1
#-------------------------------------------------------
#python bin/cv.py -e Parforward  -n 200 -d 1

#----------------
# Hyper-parameters tuning:
#------------------
#python3 bin/hyperParam.py -e Parforward  -n 100 -d 1 -b sim


#python3 bin/classification.py -e Full -n 100 -d 1 -b sim #ibmq #sim
#python3 bin/classification2.py -n 10 -d 1 -b both -c 10
python3 bin/classification.py -n 4000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 7000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 8000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 9000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 10000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 15000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 20000 -d 1 -b sim #ibmq #sim
#python3 bin/classification.py -n 25000 -d 1 -b sim #ibmq #sim
#python bin/classification.py -e Parbackward -n 2000 -d 10
#----
#python bin/classification.py -e Secondforward -n  2000 -d 10
#python bin/classification.py -e Secondbackward -n 2000 -d 10


#python bin/classification.py -e Fullforward -n 2000 -d 1
#python bin/classification.py -e Fullbackward -n 2000 -d 10

#IBMQ.get_provider(hub='ibm-q', group='open', project='main')

#python bin/param.py -e Parforward  -n 1000 -d 1
