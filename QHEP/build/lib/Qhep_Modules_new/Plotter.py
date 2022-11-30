import ROOT
import logging
from array import array
import sys, re
import os
import matplotlib.pyplot as plt
import matplotlib

cwd = os.getcwd()
sys.path.insert(0, "{}/modules/".format(cwd))
import numpy as np


def getHistograms(config, dataset, weights, j, min, max, sample="signal"):
    # h=ROOT.TH1F(sample+"_"+config["var_signal"][j],sample+"_"+config["var_signal"][j],100,min,int(max/10+12))
    h = ROOT.TH1F(
        sample + "_" + config["var_signal"][j],
        sample + "_" + config["var_signal"][j],
        100,
        min,
        max,
    )
    for i in range(len(dataset)):
        # if i%10000==0 :
        #    logging("processing variable "+str(j)+" : "+str(i)+"/"+str(len(dataset)))
        # if i%10000==0 :
        #   print(weights[i])
        h.Fill(dataset[i][j], weights[i])
        # print("Column:",dataset[i][j], "weight:",weights[i])
        # h.Fill(dataset[i][j])
    h.GetXaxis().SetTitle(config["var_name"][j])
    # print (config["var_name"][j])
    return h


def plotVars(config, signal_dataset, bkg_dataset, outName):
    j = 0
    for variable in config["var_signal"]:
        max = 1.0
        min =-1.0
        for i in range(len(signal_dataset)):
           max=np.maximum(signal_dataset[i][j],max)
           min=np.minimum(signal_dataset[i][j],min)
        #for i in range(len(bkg_dataset)):
        #   max=np.maximum(bkg_dataset[i][j],max)
        #    min=np.minimum(bkg_dataset[i][j],min)
        logging.info("range of " + variable + " : [" + str(min) + ", " + str(max) + "]")
        logging.info("plotting " + variable)
        hist_signal = getHistograms(
            config, signal_dataset, j, min, max, sample="signal"
        )
        hist_bkg = getHistograms(
            config, bkg_dataset, j, min, max, sample="bkg"
        )
        canvas = ROOT.TCanvas("canvas", "canvas", 600, 600)
        canvas.SetTickx()
        canvas.SetTicky()
        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gStyle.SetOptStat(0)
        canvas.SetName(config["var_signal"][j])
        if hist_bkg.Integral() != 0 and hist_bkg.Integral() != 0:
            hist_signal.Scale(1 / float(hist_signal.Integral()))
            hist_bkg.Scale(1 / float(hist_bkg.Integral()))
        hist_signal.SetLineColor(2)
        hist_bkg.SetLineColor(4)
        hist_signal.SetMarkerColor(2)
        hist_signal.GetYaxis().SetTitle("Arbitrary units")
        hist_signal.GetYaxis().SetTitleOffset(1.7)
        hist_signal.GetXaxis().SetTitleOffset(1.2)
        hist_bkg.SetMarkerColor(4)
        hist_signal.SetFillColor(2)
        hist_bkg.SetFillColor(4)
        hist_signal.SetFillStyle(3004)
        hist_signal.SetTitle("")
        hist_bkg.SetFillStyle(3005)
        hist_bkg.SetTitle("")
        ymax = np.maximum(hist_signal.GetMaximum(), hist_bkg.GetMaximum())
        hist_signal.GetYaxis().SetRangeUser(0, 1.5 * ymax)
        #hist_signal.GetXaxis().SetRangeUser(-1, 1)

        hist_signal.Draw("hist")
        hist_bkg.Draw("hist same")
        legend = ROOT.TLegend(0.60, 0.74, 0.85, 0.86)
        legend.SetBorderSize(0)
        legend.SetTextAlign(12)
        legend.SetTextFont(42)
        legend.SetTextSize(0.035)
        legend.SetLineColor(0)
        legend.SetLineStyle(0)
        legend.SetLineWidth(0)
        legend.SetFillColor(0)
        legend.SetFillStyle(1001)

        legend.AddEntry(hist_signal, "Signal", "F")
        legend.AddEntry(hist_bkg, "Background", "F")
        legend.Draw()
        canvas.Draw()
        # WorkDir = config["WorkDir"]
        # print(WorkDir+"plots/"+config["var_signal"][j]+".png")
        canvas.Print("./plots/qqyy/Variables/" + config["var_signal"][j] +"_"+ outName + ".pdf")
        # print("plots/"+config["var_signal"][j]+".png")
        # canvas.SaveAs("plot/"+config["var_signal"][j]+".png")
        #%jsroot on
        # canvas.Write()
        # hist_signal.Write()
        # hist_bkg.Write()
        j = j + 1
        # print('proccessed', j)


def plot_the_response(
    scores_train_signal,
    scores_test_signal,
    weight_train_signal,
    weight_test_signal,
    scores_train_bkg,
    scores_test_bkg,
    weight_train_bkg,
    weight_test_bkg,
    method,
):

    c = ROOT.TCanvas("c", "c", 600, 600)
    c.SetTickx()
    c.SetTicky()
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gStyle.SetOptStat(0)

    max = np.maximum(scores_train_signal.max(), scores_test_signal.max())
    max = np.maximum(0, max)
    min = np.minimum(scores_train_signal.max(), scores_test_signal.max())
    min = np.minimum(0, min)

    hist_train_signal = ROOT.TH1F(
        "train_signal" + method, "train_signal" + method, 20, -max, max
    )
    hist_test_signal = ROOT.TH1F(
        "test_signal" + method, "test_signal" + method, 20, -max, max
    )
    hist_train_bkg = ROOT.TH1F(
        "train_bkg" + method, "train_bkg" + method, 20, -max, max
    )
    hist_test_bkg = ROOT.TH1F("test_bkg" + method, "test_bkg" + method, 20, -max, max)

    for i in range(len(scores_train_signal)):
        hist_train_signal.Fill(scores_train_signal[i], weight_train_signal[i])
        hist_test_signal.Fill(scores_test_signal[i], weight_test_signal[i])
        hist_train_bkg.Fill(scores_train_bkg[i], weight_train_bkg[i])
        hist_test_bkg.Fill(scores_test_bkg[i], weight_test_bkg[i])

    hist_train_signal.Draw("hist")
    hist_train_signal.SetTitle("")
    hist_train_signal.GetXaxis().SetTitle(method + " Response")
    hist_train_signal.GetXaxis().SetTitleOffset(1.2)
    hist_train_signal.GetYaxis().SetTitle("Arbitrary units")
    hist_train_signal.SetMarkerColor(2)
    hist_train_signal.SetFillColor(2)
    hist_train_signal.SetLineColor(2)
    hist_train_signal.SetFillStyle(3004)

    # hist_test_signal.Draw("same P")
    hist_test_signal.SetTitle("")
    hist_test_signal.SetMarkerStyle(20)
    hist_test_signal.SetMarkerColor(2)
    hist_test_signal.SetLineColor(2)
    hist_test_signal.SetMarkerSize(1.2)

    hist_train_bkg.Draw("Hist same")
    hist_train_bkg.SetTitle("")
    hist_train_bkg.SetMarkerColor(4)
    hist_train_bkg.SetLineColor(4)
    hist_train_bkg.SetFillColor(4)
    hist_train_bkg.SetFillStyle(3005)

    # hist_test_bkg.Draw("same P")
    hist_test_bkg.SetTitle("")
    hist_test_bkg.SetMarkerStyle(20)
    hist_test_bkg.SetMarkerColor(4)
    hist_test_bkg.SetLineColor(4)
    hist_test_bkg.SetMarkerSize(1.2)

    hist_train_signal.Scale(1 / hist_train_signal.Integral())
    hist_test_signal.Scale(1 / hist_test_signal.Integral())
    hist_test_bkg.Scale(1 / hist_test_bkg.Integral())

    max = (
        np.maximum(
            np.maximum(hist_train_signal.GetMaximum(), hist_train_bkg.GetMaximum()),
            np.maximum(hist_test_signal.GetMaximum(), hist_test_bkg.GetMaximum()),
        )
        * 1.5
    )

    hist_train_signal.GetYaxis().SetRangeUser(0, max)

    # legend = ROOT.TLegend(0.58,0.71,0.74,0.87)
    legend = ROOT.TLegend(0.58, 0.73, 0.78, 0.87)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineColor(0)
    legend.SetLineStyle(0)
    legend.SetLineWidth(0)
    legend.SetFillColor(0)
    legend.SetFillStyle(1001)

    legend.AddEntry(hist_train_signal, "Signal", "F")
    legend.AddEntry(hist_train_bkg, "Background", "F")
    # legend.AddEntry(hist_train_signal,"Signal (train)","F")
    # legend.AddEntry(hist_train_bkg,"Background (train)","F")
    # legend.AddEntry(hist_test_signal,"Signal (test)","ep")
    # legend.AddEntry(hist_test_bkg,"Background (test)","ep")
    legend.Draw()

    c.SaveAs("./plots/qqyy/Responses/Response_" + method + ".pdf")


def plotOvertraining(
    signal_train_prob,
    bkg_train_prob,
    signal_test_prob,
    bkg_test_prob,
    weight_train,
    weight_test,
    method,
):

    c = ROOT.TCanvas("c", "c", 600, 600)
    c.SetTickx()
    c.SetTicky()
    ROOT.gPad.SetLeftMargin(0.13)

    max = np.maximum(signal_train_prob.max(), signal_test_prob.max())
    max = np.maximum(1, max)
    min = np.minimum(bkg_train_prob.max(), bkg_test_prob.max())
    min = np.minimum(-1, min)
    # max_b=np.maximum(bkg_train_prob.max(),bkg_test_prob.max())
    # max=np.maximum(1,max)
    # min_b=np.minimum(bkg_train_prob.max(),bkg_test_prob.max())
    # min=np.minimum(-1,min)
    # min =-1
    # max = 1.0
    h_s_train = ROOT.TH1F("TrainS_" + method, "TrainS_" + method, 100, min, max)
    h_s_test = ROOT.TH1F("TestS_" + method, "TestS_" + method, 100, min, max)
    h_b_train = ROOT.TH1F("TrainB_" + method, "TrainS_" + method, 100, min, max)
    h_b_test = ROOT.TH1F("TestB_" + method, "TestS_" + method, 100, min, max)

    for i in range(len(signal_train_prob)):
        h_s_train.Fill(signal_train_prob[i], weight_train[i])
        # print("signal train:", signal_train_prob[i], "weight train:", weight_train[i])

        h_s_test.Fill(signal_test_prob[i], weight_test[i])
        # print("signal test:", signal_test_prob[i], "weight test:", weight_test[i])

    for i in range(len(signal_train_prob)):
        h_b_train.Fill(bkg_train_prob[i], weight_train[i])
        # print("bkg train:", bkg_train_prob[i], "weight train:", weight_train[i])

        h_b_test.Fill(bkg_test_prob[i], weight_test[i])
        # print("bkg train:", bkg_test_prob[i], "weight test:", weight_test[i])
    #
    h_s_train.SetTitle("")
    h_s_train.GetXaxis().SetTitle("Classification " + method + " Output")
    h_s_train.GetXaxis().SetTitleOffset(1.2)
    # h_s_train.GetYaxis().SetTitle("Arbitrary units")
    h_s_train.GetYaxis().SetTitle("Events")
    h_s_train.GetXaxis().SetRangeUser(bkg_train_prob.min(), signal_train_prob.max())
    h_s_train.SetMarkerColor(2)
    h_s_train.SetFillColor(2)
    h_s_train.SetLineColor(2)
    h_s_train.SetFillStyle(3004)

    h_b_train.SetTitle("")
    h_b_train.SetMarkerColor(4)
    h_b_train.SetLineColor(4)
    h_b_train.SetFillColor(4)
    h_b_train.SetFillStyle(3005)

    h_s_test.SetTitle("")
    h_s_test.SetMarkerStyle(20)
    h_s_test.SetMarkerColor(2)
    h_s_test.SetLineColor(2)
    h_s_test.SetMarkerSize(1.2)

    h_b_test.SetTitle("")
    h_b_test.SetMarkerStyle(20)
    h_b_test.SetMarkerColor(4)
    h_b_test.SetLineColor(4)
    h_b_test.SetMarkerSize(1.2)

    # h_s_train.Scale(1/h_s_train.Integral())
    # h_s_test.Scale(1/h_s_test.Integral())
    # h_b_train.Scale(1/h_b_train.Integral())
    # h_b_test.Scale(1/h_b_test.Integral())

    max = (
        np.maximum(
            np.maximum(h_s_train.GetMaximum(), h_b_train.GetMaximum()),
            np.maximum(h_s_test.GetMaximum(), h_b_test.GetMaximum()),
        )
        * 1.5
    )

    h_s_train.GetYaxis().SetRangeUser(0, max)

    legend = ROOT.TLegend(0.15, 0.78, 0.80, 0.89)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineColor(0)
    legend.SetLineStyle(0)
    legend.SetLineWidth(0)
    legend.SetFillColor(0)
    # legend.SetFillStyle(1001)

    legend.AddEntry(h_s_train, "Signal (Train)", "F")
    legend.AddEntry(h_s_test, "Signal (Test)", "ep")
    legend.AddEntry(h_b_train, "Background (Train)", "F")
    legend.AddEntry(h_b_test, "Background (Test)", "ep")

    h_s_train.Draw("HIST")
    h_b_train.Draw("same HIST")
    h_s_test.Draw("same P")
    h_b_test.Draw("same P")

    legend.Draw()

    c.SaveAs("./plots/qqyy/Responses/OverTraining_" + method + ".pdf")


#def plotROCcurve(tpr_QSVM, fpr_QSVM, aruc_QSVM, tpr_SVM, fpr_SVM, aruc_SVM, name):
def plotROCcurve(tpr_QSVM, fpr_QSVM, aruc_QSVM, std_qauc, tpr_SVM, fpr_SVM, aruc_SVM, std_auc, name):
    c2 = ROOT.TCanvas( 'c2', '', 600, 600 )
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.05)
    c2.SetTickx()
    c2.SetTicky()
    c2.SetGridx();
    c2.SetGridy();
    
    MultiGraph  = ROOT.TMultiGraph()
    
    graph_QSVM = ROOT.TGraph( len(fpr_QSVM), tpr_QSVM,  1 - fpr_QSVM )
    graph_QSVM.SetLineColor( ROOT.kBlue )
    graph_QSVM.SetLineWidth( 3 )
    graph_QSVM.SetMarkerStyle( 21 )
    MultiGraph.Add(graph_QSVM,"A")
    
    graph_SVM  = ROOT.TGraph( len(fpr_SVM), tpr_SVM,  1 - fpr_SVM )
    graph_SVM.SetLineColor( ROOT.kRed )
    graph_SVM.SetLineStyle(ROOT.kDashed)
    graph_SVM.SetLineWidth( 3 )
    graph_SVM.SetTitle( '' )
    graph_SVM.SetMarkerStyle( 21 )
    MultiGraph.Add(graph_SVM,"A")
    
    MultiGraph.Draw("A")
    MultiGraph.GetXaxis().SetTitle( 'Signal efficiency' )
    MultiGraph.GetYaxis().SetTitle( 'Background rejection' )
    MultiGraph.GetXaxis().SetTitleSize(0.04)
    MultiGraph.GetXaxis().SetLabelSize(0.04)
    MultiGraph.GetXaxis().SetTitleOffset(1.1)
    MultiGraph.GetYaxis().SetTitleSize(0.04)
    MultiGraph.GetYaxis().SetLabelSize(0.04)
    MultiGraph.GetYaxis().SetTitleOffset(1.5)
    MultiGraph.GetXaxis().SetTickLength(0.018);
    MultiGraph.GetYaxis().SetTickLength(0.018);
    MultiGraph.GetXaxis().SetNdivisions(515)
    MultiGraph.GetYaxis().SetNdivisions(515)
    MultiGraph.GetYaxis().SetRangeUser(0,1.01)
    MultiGraph.GetXaxis().SetRangeUser(0,1.01)
    
    legend = ROOT.TLegend(0.16, 0.15, 0.72, 0.28)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.SetTextSize(0.032)
    legend.SetLineColor(0)
    legend.SetLineStyle(0)
    #legend.SetLineWidth(3)
    legend.SetFillColor(0)
    
    #legend.AddEntry(graph_QSVM, "QSVM (AUC = %0.3f)" % (aruc_QSVM), "l")
    #legend.AddEntry(graph_SVM, "SVM (AUC = %0.3f)" % (aruc_SVM), "l")
    legend.AddEntry(graph_QSVM, "QSVM (AUC = %0.3f \pm %0.3f)" % (aruc_QSVM, std_qauc), "l")
    legend.AddEntry(graph_SVM, "SVM (AUC = %0.3f \pm %0.3f)" % (aruc_SVM, std_auc), "l")

    legend.Draw()
    
    c2.Update()
    c2.Modified()
    c2.Update()
    c2.Print("./plots/qqyy/ROCs/"+ name +".pdf")

def plotROC_IBM(tpr_QSVM, fpr_QSVM, aruc_QSVM, std_qauc, tpr_QSVM_ibm, fpr_QSVM_ibm, aruc_QSVM_ibm, std_qauc_ibm, tpr_SVM, fpr_SVM, aruc_SVM, std_auc, name):
    c2 = ROOT.TCanvas( 'c2', '', 600, 600 )
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetRightMargin(0.05)
    c2.SetTickx()
    c2.SetTicky()
    c2.SetGridx();
    c2.SetGridy();

    MultiGraph  = ROOT.TMultiGraph()

    graph_QSVM = ROOT.TGraph( len(fpr_QSVM), tpr_QSVM,  1 - fpr_QSVM )
    graph_QSVM.SetLineColor( ROOT.kBlack )
    graph_QSVM.SetLineWidth( 3 )
    graph_QSVM.SetMarkerStyle( 21 )
    MultiGraph.Add(graph_QSVM,"A")

    graph_QSVM_ibm = ROOT.TGraph( len(fpr_QSVM_ibm), tpr_QSVM_ibm,  1 - fpr_QSVM_ibm )
    graph_QSVM_ibm.SetLineColor( ROOT.kBlue )
    graph_QSVM_ibm.SetLineWidth( 3 )
    graph_QSVM_ibm.SetMarkerStyle( 21 )
    MultiGraph.Add(graph_QSVM_ibm,"A")


    graph_SVM  = ROOT.TGraph( len(fpr_SVM), tpr_SVM,  1 - fpr_SVM )
    graph_SVM.SetLineColor( ROOT.kRed )
    graph_SVM.SetLineStyle(ROOT.kDashed)
    graph_SVM.SetLineWidth( 3 )
    graph_SVM.SetTitle( '' )
    graph_SVM.SetMarkerStyle( 21 )
    MultiGraph.Add(graph_SVM,"A")

    MultiGraph.Draw("A")
    MultiGraph.GetXaxis().SetTitle( 'Signal efficiency' )
    MultiGraph.GetYaxis().SetTitle( 'Background rejection' )
    MultiGraph.GetXaxis().SetTitleSize(0.04)
    MultiGraph.GetXaxis().SetLabelSize(0.04)
    MultiGraph.GetXaxis().SetTitleOffset(1.1)
    MultiGraph.GetYaxis().SetTitleSize(0.04)
    MultiGraph.GetYaxis().SetLabelSize(0.04)
    MultiGraph.GetYaxis().SetTitleOffset(1.5)
    MultiGraph.GetXaxis().SetTickLength(0.018);
    MultiGraph.GetYaxis().SetTickLength(0.018);
    MultiGraph.GetXaxis().SetNdivisions(515)
    MultiGraph.GetYaxis().SetNdivisions(515)
    MultiGraph.GetYaxis().SetRangeUser(0,1.01)
    MultiGraph.GetXaxis().SetRangeUser(0,1.01)

    legend = ROOT.TLegend(0.16, 0.15, 0.72, 0.28)
    legend.SetBorderSize(0)
    legend.SetTextAlign(12)
    legend.SetTextFont(42)
    legend.SetTextSize(0.032)
    legend.SetLineColor(0)
    legend.SetLineStyle(0)
    #legend.SetLineWidth(3)
    legend.SetFillColor(0)

    legend.AddEntry(graph_QSVM_ibm, "IBM (AUC = %0.3f \pm %0.3f)" % (aruc_QSVM_ibm, std_qauc_ibm), "l")
    legend.AddEntry(graph_QSVM, "SVM (AUC = %0.3f \pm %0.3f)" % (aruc_QSVM, std_qauc), "l")
    legend.AddEntry(graph_SVM, "SVM (AUC = %0.3f \pm %0.3f)" % (aruc_SVM, std_auc), "l")


    legend.Draw()

    c2.Update()
    c2.Modified()
    c2.Update()
    c2.Print("./plots/qqyy/ROCs/"+ name +"_comp.pdf")


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="left",
             rotation_mode="anchor")

    plt.xlabel("$\gamma$", labelpad=0.05)
    plt.ylabel("C", labelpad=0.05)
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    #ax.ticklabel_format(style='plain')

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
