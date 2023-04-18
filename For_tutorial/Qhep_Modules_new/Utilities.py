import ROOT

# import root_numpy
from root_numpy import root2array, rec2array, array2hist, tree2array
import logging
from array import array
import numpy as np
import sys, re
import os

cwd = os.getcwd()
sys.path.insert(0, "{}/modules/".format(cwd))
import json
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import pandas as pd


def readConfig(configFile="config_qqyy_8.json"):

    with open("config/" + configFile, "r") as f_obj:
        config = json.load(f_obj)
        # print(config)
        logging.info("Signal sample(s) %s", config["signal_file_list"])
        logging.info("Background sample(s) %s", config["bkg_file_list"])
        logging.info("Signal(s) tree name %s", config["signal_tree_name"])
        logging.info("Background(s) tree name %s", config["bkg_tree_name"])
        logging.info("Signal(s) weight variable %s", config["signal_weight_name"])
        logging.info("Background(s) weight variable %s", config["bkg_weight_name"])
        logging.info("Signal(s) variables %s", config["var_signal"])
        logging.info("Background(s) variables %s", config["var_bkg"])

    return config


#def Trandform(dataArray, rangConst="_0_1"):
#    dataConstrain = []

    # mini = np.quantile(dataArray, q=0.05, axis=0)
    # maxi = np.quantile(dataArray, q=0.95, axis=0)

#    if rangConst == "_0_1":  # [0, 1]
#        dataConstrain = (dataArray - np.min(dataArray, axis=0)) / np.ptp(
#            dataArray, axis=0
#        )
#    elif rangConst == "_n1_1":  # [-1,1]
        
#        mini = 0.0 #[min(idx) for idx in zip(*dataArray)]
#        maxi = 200.0 #[max(idx) for idx in zip(*dataArray)]

        #nominator   =  np.subtract(dataArray,dataMin) #[a - b for a, b in zip(data, dataMin)] #data - dataMin[:,None]
        #denomenator =  np.subtract(dataMax,dataMin)
        #ratio = 2*np.divide(nominator, denomenator) - 1

        #dataConstrain = ratio
        
        #dataConstrain = 2.*(dataArray - np.min(dataArray))/np.ptp(dataArray)-1 #T1
        #dataConstrain = ( 2.0 * (dataArray - np.min(dataArray, axis=0)) / np.ptp(dataArray, axis=0)) - 1  # T2
        #dataConstrain = ( 2.0 * (dataArray - np.min(dataArray, axis=0)) / np.ptp(dataArray, axis=0)) - 1  # T2
        #dataConstrain = ( 2.*( dataArray - mini)/(maxi - mini) ) - 1
        # dataConstrain = 2.*(np.clip(dataArray, a_min=np.min(dataArray,axis=0), a_max=np.max(dataArray,axis=0)) - np.min(dataArray,axis=0))/(np.max(dataArray,axis=0) - np.min(dataArray,axis=0))-1 #T3

        #print("list:\n",dataArray)
        #print("dat-min\n:", nomin)
        #print("max-min\n:", demon)
        #print("ratio:\n", dataConstrain)

        #print("list:\n",dataArray)
        #print("list after:\n",dataConstrain)
        #print("mim:\n", np.min(dataArray, axis=0))
        #print("max:\n", np.max(dataArray, axis=0))
        #print("max-min\n:", np.max(dataArray, axis=1) - np.min(dataArray, axis=1))
        #print("ptp(max-min)\n:", np.ptp(dataArray, axis=0))

def Trandform(dataArray, dataArray2, rangConst="_0_1"):
    dataConstrain = []

    if rangConst == "_0_1":  # [0, 1]
        dataConstrain = (dataArray - np.min(dataArray, axis=0)) / np.ptp(
            dataArray, axis=0
        )
    elif rangConst == "_n1_1":  # [-1,1]
         mini1 = np.min(dataArray, axis=0)
         mini2 = np.min(dataArray2, axis=0)
         maxi1 = np.max(dataArray, axis=0)
         maxi2 = np.max(dataArray2, axis=0)
         mini = np.minimum(mini1, mini2)
         maxi = np.maximum(maxi1, maxi2)
         dataConstrain = (2.*(dataArray - mini)/(maxi-mini))-1 #T1
        #dataConstrain = (
            #2.0 * (dataArray - 0)/ 100) - 1  # T2
        # dataConstrain = ( 2.*(np.clip(dataArray, a_min=mini, a_max=maxi) - mini)/(maxi - mini) ) - 1
        # dataConstrain = 2.*(np.clip(dataArray, a_min=np.min(dataArray,axis=0), a_max=np.max(dataArray,axis=0)) - np.min(dataArray,axis=0))/(np.max(dataArray,axis=0) - np.min(dataArray,axis=0))-1 #T

         print("min:\n", mini)
         print("max:\n", maxi)
         print("min1:\n", mini1)
         print("max1:\n", maxi1)
         print("min2:\n", mini2)
         print("max2:\n", maxi2)

    elif rangConst == "_0_2pi":  # [0, 2pi]
        dataConstrain = (
            2
            * np.pi
            * (dataArray - np.min(dataArray, axis=0))
            / np.ptp(dataArray, axis=0)
        )
    elif rangConst == "_npi_pi":  # [-pi, pi]
        dataConstrain = -np.pi + np.pi * (
            dataArray - np.min(dataArray, axis=0)
        ) / np.ptp(dataArray, axis=0)

    return dataConstrain


def preparingData(
    confiFile="config_qqyy_8.json", prossEvent=10, removeNegativeWeight=True, fraction=0.5, seed=None, dataType="Classical"
):

    config = readConfig(confiFile)
    #root2array(filenames, treename=None, branches=None, selection=None, object_selection=None, start=None, stop=None, step=None, include_weight=False, weight_name='weight', cache_size=-1, warn_missing_tree=False)
    signal_dataset = root2array(
        filenames=config["signal_file_list"],
        treename=config["signal_tree_name"],
        branches=config["var_signal"],
        selection=config["signal_selection"],
        include_weight=False,
        weight_name="weight",
        stop=prossEvent,
    )
    #signal_dataset = signal_dataset.astype(np.float32)

    signal_weights = root2array(
        filenames=config["signal_file_list"],
        treename=config["signal_tree_name"],
        branches=config["signal_weight_name"],
        include_weight=False,
        stop=prossEvent,
    )

    bkg_dataset = root2array(
        filenames=config["bkg_file_list"],
        treename=config["bkg_tree_name"],
        branches=config["var_bkg"],
        selection=config["bkg_selection"],
        include_weight=False,
        weight_name="weight",
        stop=prossEvent,
    )

    bkg_weights = root2array(
        filenames=config["bkg_file_list"],
        treename=config["bkg_tree_name"],
        branches=config["bkg_weight_name"],
        include_weight=False,
        stop=prossEvent,
    )
 

    #balance the weight. It's not important for now since we don't use the weight
    #signal_weights = signal_weights*(bkg_weights/signal_weights)
    signal_dataset1 = rec2array(signal_dataset)
    bkg_dataset1 = rec2array(bkg_dataset)
    signal_dataset = rec2array(signal_dataset)
    bkg_dataset = rec2array(bkg_dataset)

    #print("check weight",signal_weights[0], bkg_weights[0])

    # Transformation applies for both Classical and Quantum
    #signal_dataset = Trandform(signal_dataset, "_n1_1")
    #bkg_dataset = Trandform(bkg_dataset, "_n1_1")
    signal_dataset = Trandform(signal_dataset, bkg_dataset1, "_n1_1")
    bkg_dataset = Trandform(bkg_dataset, signal_dataset1, "_n1_1")

    print("number of signal:",len(signal_dataset),"number of background:",len(bkg_dataset))
    # remove negative weights and then normalise each events by the it's weight
    # signal weight
    if removeNegativeWeight == True:
        for i in range(len(signal_weights)):
            if signal_weights[i] < 0:
                signal_weights[i] = signal_weights[i] * -1

        for i in range(len(bkg_weights)):
            if bkg_weights[i] < 0:
                bkg_weights[i] = bkg_weights[i] * -1

    train_size = int(len(signal_dataset) * fraction)
    test_size = int(len(signal_dataset) * fraction)

    logging.info("Total number of signal : %s", str(len(signal_dataset)))
    logging.info("Total number of backgrouns : %s", str(len(bkg_dataset)))
    logging.info("Train size : %s", str(train_size))
    logging.info("Testing size : %s", str(test_size))

    X_signal = signal_dataset
    X_signal_weights = signal_weights
    X_background = bkg_dataset
    X_bkg_weights = bkg_weights

    y_signal = np.ones(X_signal.shape[0])
    y_background = np.ones(X_background.shape[0])
    y_background = -1 * y_background

    X = np.concatenate([X_signal, X_background], axis=0)
    y = np.concatenate([y_signal, y_background])
    XW = np.concatenate([X_signal_weights, X_bkg_weights], axis=0)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, train_size=train_size, test_size=test_size, random_state=seed
    )
    XW_train, XW_test = train_test_split(XW, train_size=train_size, test_size=test_size)
    #else:
    #   X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size, test_size=test_size, shuffle=False) #random_state=None)
    #   XW_train, XW_test = train_test_split(XW, train_size=train_size, test_size=test_size, random_state=None)

    # TODO: Use standarised only in the classical case
    if dataType=="Classical":
       print("The data will be prepared for the classical case so transformation is needed.")
       scaler = StandardScaler()
       X_train = scaler.fit_transform(X_train)
       X_test = scaler.transform(X_test)
    else:
       # Do nothing here
       print("The data will be prepared for the quantum case so no transformation is needed.")
     

    labels = {1: "S", -1: "B"}

    train_dataset = {
        labels[1]: X_train[y_train == 1],
        labels[-1]: X_train[y_train == 1],
    }
    test_dataset = {labels[1]: X_test[y_test == 1], labels[-1]: X_test[y_test == -1]}

    # print(signal_dataset)
    # return signal_dataset,signal_weights,bkg_dataset,bkg_weights
    return (
        X_train,
        X_test,
        y_train,
        y_test,
        XW_train,
        XW_test,
        train_dataset,
        test_dataset,
        signal_dataset,
        bkg_dataset,
        signal_weights,
        bkg_weights,
        X,
        y,
    )


def getSigbkgProb(label, probs):
    sig = []
    backG = []
    for i in range(len(probs)):
        if label[i] == 1:
            sig.append(probs[i])

        elif label[i] == -1:
            backG.append(probs[i])

    return sig, backG
