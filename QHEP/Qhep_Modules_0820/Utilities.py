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

def readConfig(configFile="config_qqyy.json"):

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

def Trandform(sigArray, bkgArray, rangConst="_0_1"):
    sigConstrain = []
    bkgConstrain = []

    mini1 = np.min(sigArray, axis=0)
    mini2 = np.min(bkgArray, axis=0)
    maxi1 = np.max(sigArray, axis=0)
    maxi2 = np.max(bkgArray, axis=0)
    mini = np.minimum(mini1, mini2)
    maxi = np.maximum(maxi1, maxi2)


    if rangConst == "_0_1":  # [0, 1]
        sigConstrain = (sigArray - np.min(sigArray, axis=0)) / np.ptp(
            sigArray, axis=0
        )
    elif rangConst == "_n1_1":  # [-1,1]

         sigConstrain = (2.*(sigArray - mini)/(maxi-mini))-1 #T1
         bkgConstrain = (2.*(bkgArray - mini)/(maxi-mini))-1 #T1


    elif rangConst == "_0_2pi":  # [0, 2pi]
        sigConstrain = (
            2
            * np.pi
            * (sigArray - np.min(sigArray, axis=0))
            / np.ptp(sigArray, axis=0)
        )
    elif rangConst == "_npi_pi":  # [-pi, pi]
        sigConstrain = -np.pi + np.pi * (
            sigArray - np.min(sigArray, axis=0)
        ) / np.ptp(sigArray, axis=0)

    return sigConstrain, bkgConstrain

def preparingData(confiFile="config_qqyy_8.json", prossEvent=100, fraction=1, seed=None, dataType="Classical"):

    config = readConfig(confiFile)

    signal_dataset = root2array(
        filenames=config["signal_file_list"],
        treename=config["signal_tree_name"],
        branches=config["var_signal"],
        selection=config["signal_selection"],
        include_weight=False,
        #include_weight=True,
        weight_name=config["signal_weight_name"],
        stop=prossEvent,
    )

    bkg_dataset = root2array(
        filenames=config["bkg_file_list"],
        treename=config["bkg_tree_name"],
        branches=config["var_bkg"],
        selection=config["bkg_selection"],
        include_weight=False,
        #include_weight=True,
        weight_name=config["bkg_weight_name"],
        stop=prossEvent,
    )
 
    signal_dataset = rec2array(signal_dataset)
    bkg_dataset = rec2array(bkg_dataset)
    
    signal_dataset, bkg_dataset= Trandform(signal_dataset, bkg_dataset, "_n1_1")

    print("number of signal:",len(signal_dataset),"number of background:",len(bkg_dataset))

    train_size = int(len(signal_dataset) * fraction)
    test_size = int(len(signal_dataset) * fraction)

    logging.info("Total number of signal : %s", str(len(signal_dataset)))
    logging.info("Total number of backgrouns : %s", str(len(bkg_dataset)))
    logging.info("Train size : %s", str(train_size))
    logging.info("Testing size : %s", str(test_size))

    X_signal = signal_dataset
    X_background = bkg_dataset

    y_signal = np.ones(X_signal.shape[0])
    y_background = np.ones(X_background.shape[0])
    y_background = -1 * y_background

    X = np.concatenate([X_signal, X_background], axis=0)
    y = np.concatenate([y_signal, y_background])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, train_size=train_size, test_size=test_size, random_state=seed
    )

    # Use standarised only in the classical case
    if dataType=="Classical":
       print("The data will be prepared for the classical case so transformation is needed.")
       scaler = StandardScaler()
       X_train = scaler.fit_transform(X_train)
       X_test = scaler.transform(X_test)
       X = scaler.fit_transform(X)
    else:
       # Do nothing here
       print("The data will be prepared for the quantum case so no transformation is needed.")
     
    labels = {1: "S", -1: "B"}

    train_dataset = {
        labels[1]: X_train[y_train == 1],
        labels[-1]: X_train[y_train == -1],
    }
    test_dataset = {labels[1]: X_test[y_test == 1], labels[-1]: X_test[y_test == -1]}

    return (
        X_train,
        X_test,
        y_train,
        y_test,
        signal_dataset,
        bkg_dataset,
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
