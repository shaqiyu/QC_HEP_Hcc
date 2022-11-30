import ROOT as r

def WeightSum(tree_sig, tree_bkg):
    weightSum_sig = 0.
    weightSum_bkg = 0.
    print("unweighted signal events is",tree_sig.GetEntries(),"and bkg is",tree_bkg.GetEntries())
    i=0 
    for event in tree_sig:
        weighsig = event.weight
        if i <= 2000:
            weightSum_sig += event.weight
        i+=1
    print("signal weight per event:",weighsig)
    j=0 
    for event in tree_bkg:
        weightbkg = event.weight
        if j<= 2000:
            weightSum_bkg += event.weight
           
            
        j+=1
    print("bkg weight per events:",weightbkg)
    return weightSum_sig, weightSum_bkg

file_sig = r.TFile('/media/abdualazem/Work/WorkSpace/QC/code/Qiskit/QHEP/samples/qqyy/Results_qqsig.root')
file_bkg = r.TFile('/media/abdualazem/Work/WorkSpace/QC/code/Qiskit/QHEP/samples/qqyy/Results_qqbkg.root')

tree_sig = file_sig.Get('BDTtrain')
tree_bkg = file_bkg.Get('BDTtrain')

weight_signal_events, weight_bkg_events = WeightSum(tree_sig, tree_bkg)

print("Total weighted events for the signal is ",weight_signal_events,"and bkg is",weight_bkg_events)

weight_signal_events = weight_signal_events*(weight_bkg_events/weight_signal_events)

print("new signal weight", weight_signal_events)