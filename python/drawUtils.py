#! /usr/bin/env python

import os, multiprocessing
import copy
import math
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory, gPad
from ROOT import TFile, TChain, TTree, TCut, TH1, TH1F, TH2F, THStack, TGraph, TGraphAsymmErrors
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine, TBox, TGaxis

#### IMPORT SAMPLES AND VARIABLES DICTIONARIES ####

from Analyzer.LLPonAOD.samples import sample, samples
from Analyzer.LLPonAOD.variables import *
#from Analyzer.LLPonAOD.skimmed_variables import *
from Analyzer.LLPonAOD.selections import *

##################
#    PROJECT     #
##################

def project(var, cut, cut_s, weight, samplelist, pd, ntupledir, treename="ntuple/tree", formula=""):
#def project(var, cut, cut_s, weight, samplelist, pd, ntupledir, treename="trigger/tree"):
    # Create dict
    file = {}
    tree = {}
    chain = {}
    hist = {}
    
    ### Create and fill MC histograms ###
    for i, s in enumerate(samplelist):
        if "HIST" in cut: # Histogram written to file
            for j, ss in enumerate(samples[s]['files']):
                file[ss] = TFile(ntupledir + ss + ".root", "READ")
                hist[s] = file[ss].Get("ntuple/" + histDir[var[0:2]] + "/" + var) if not s in hist else hist[s].Add(file[ss].Get("ntuple/" + histDir[var[0:2]] + "/" + var))
        else: # Project from tree
            chain[s] = TChain(treename)
            for j, ss in enumerate(samples[s]['files']):
                if not 'data' in s or ('data' in s and ss in pd):
                    chain[s].Add(ntupledir + ss + ".root")
                    #print "Sample: ", ss
                    #print "current weight of chain: ", chain[s].GetWeight()
                    #filename = TFile(ntupledir + ss +'.root')
                    #tree = filename.Get("ntuple/tree") 
                    #print "real weight: ", tree.GetWeight()
                    #chain[s].SetWeight(tree.GetWeight(),"global")
                    #print "forcing weight of chain: ", chain[s].GetWeight()
            if variable[var]['nbins']>0: hist[s] = TH1F(s, ";"+variable[var]['title'], variable[var]['nbins'], variable[var]['min'], variable[var]['max']) # Init histogram
            else: hist[s] = TH1F(s, ";"+variable[var]['title'], len(variable[var]['bins'])-1, array('f', variable[var]['bins']))
            hist[s].Sumw2()
            tmpcut = cut
            tmpcut_s = cut_s
            if not 'data' in s:
                if s.endswith('_0b'): tmpcut += " && nBJets==0"
                elif s.endswith('_1b'): tmpcut += " && nBJets==1"
                elif s.endswith('_2b'): tmpcut += " && nBJets>=2"
                if s.endswith('_0l'): tmpcut += " && genNl==0"
                elif s.endswith('_1l'): tmpcut += " && genNl==1"
                elif s.endswith('_2l'): tmpcut += " && genNl>=2"
            cutstring = "("+weight+")" + ("*("+tmpcut+")" if len(tmpcut)>0 else "")
            cutstring_s = "("+weight+")" + ("*("+tmpcut_s+")" if len(tmpcut_s)>0 else "")
            if "(" in formula:
                var_updated = formula+"("+var+"))"
            else:
                var_updated = formula+"("+var+")"
            if "VBFH_M" in s:#important bugfix! Not applying jet matching to signal!
                chain[s].Project(s, var, cutstring_s) if formula=="" else chain[s].Project(s, var_updated, cutstring_s)
            elif "ggH_M" in s:#important bugfix! Not applying jet matching to signal!
                chain[s].Project(s, var, cutstring_s) if formula=="" else chain[s].Project(s, var_updated, cutstring_s)
            elif "ZH_M" in s:#important bugfix! Not applying jet matching to signal!
                chain[s].Project(s, var, cutstring_s) if formula=="" else chain[s].Project(s, var_updated, cutstring_s)
            else:
                chain[s].Project(s, var, cutstring) if formula=="" else chain[s].Project(s, var_updated, cutstring_s)
            hist[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
            hist[s].Scale(samples[s]['weight'] if hist[s].Integral() >= 0 else 0)
            #if s in sign:
                #print "Is it empty?"
                #print s, hist[s].Integral()

        hist[s].SetFillColor(samples[s]['fillcolor'])
        hist[s].SetFillStyle(samples[s]['fillstyle'])
        hist[s].SetLineColor(samples[s]['linecolor'])
        hist[s].SetLineStyle(samples[s]['linestyle'])
    
    if "HIST" in cut: hist["files"] = file
    return hist


##################
#      DRAW      #
##################

def draw(hist, data, back, sign, snorm=1, ratio=0, poisson=False, log=False):
    # If not present, create BkgSum
    if not 'BkgSum' in hist.keys():
        hist['BkgSum'] = hist['data_obs'].Clone("BkgSum") if 'data_obs' in hist else hist[back[0]].Clone("BkgSum")
        hist['BkgSum'].Reset("MICES")
        for i, s in enumerate(back): hist['BkgSum'].Add(hist[s])
    hist['BkgSum'].SetMarkerStyle(0)
    
    # Some style
    for i, s in enumerate(data):
        hist[s].SetMarkerStyle(21)
        hist[s].SetMarkerSize(1.25)
    for i, s in enumerate(sign):
        hist[s].SetLineWidth(3)
        
    for i, s in enumerate(data+back+sign+['BkgSum']):
        addOverflow(hist[s], False) # Add overflow
    
    # Set Poisson error bars
    #if len(data) > 0: hist['data_obs'].SetBinErrorOption(1) # doesn't work
    
    # Poisson error bars for data
    if poisson:
        alpha = 1 - 0.6827
        hist['data_obs'].SetBinErrorOption(TH1.kPoisson)
        data_graph = TGraphAsymmErrors(hist['data_obs'].GetNbinsX())
        data_graph.SetMarkerStyle(hist['data_obs'].GetMarkerStyle())
        data_graph.SetMarkerSize(hist['data_obs'].GetMarkerSize())
        res_graph = data_graph.Clone()
        for i in range(hist['data_obs'].GetNbinsX()):
            N = hist['data_obs'].GetBinContent(i+1)
            B = hist['BkgSum'].GetBinContent(i+1)
            L =  0 if N==0 else ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            data_graph.SetPoint(i, hist['data_obs'].GetXaxis().GetBinCenter(i+1), N if not N==0 else -1.e99)
            data_graph.SetPointError(i, hist['data_obs'].GetXaxis().GetBinWidth(i+1)/2., hist['data_obs'].GetXaxis().GetBinWidth(i+1)/2., N-L, U-N)
            res_graph.SetPoint(i, hist['data_obs'].GetXaxis().GetBinCenter(i+1), N/B if not B==0 and not N==0 else -1.e99)
            res_graph.SetPointError(i, hist['data_obs'].GetXaxis().GetBinWidth(i+1)/2., hist['data_obs'].GetXaxis().GetBinWidth(i+1)/2., (N-L)/B if not B==0 else -1.e99, (U-N)/B if not B==0 else -1.e99)
    
    
    # Create stack
    bkg = THStack("Bkg", ";"+hist['BkgSum'].GetXaxis().GetTitle()+";Events")
    for i, s in enumerate(back): bkg.Add(hist[s])
    
    # Legend
    n = len([x for x in data+back+['BkgSum']+sign if samples[x]['plot']])
    for i, s in enumerate(sign):
        if 'sublabel' in samples[s]: n+=1
        if 'subsublabel' in samples[s]: n+=1
    #leg = TLegend(0.68, 0.9-0.05*n, 0.93, 0.9)
    leg = TLegend(0.68-0.05, 0.9-0.05*n, 0.93, 0.9)#DCMS
    leg.SetTextSize(0.03)#DCMS
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader("Signal x-sec=%.0f pb"%(1*snorm))
    if len(data) > 0:
        leg.AddEntry(hist[data[0]], samples[data[0]]['label'], "ple1")
    for i, s in reversed(list(enumerate(['BkgSum']+back))):
        leg.AddEntry(hist[s], samples[s]['label'], "f")    
    for i, s in enumerate(sign):
        leg.AddEntry(hist[s], samples[s]['label'], "f")

    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[-1].GetXaxis().GetTitle(), 1000, 800 if ratio else 700)
    
    if ratio:
        c1.Divide(1, 2)
        setTopPad(c1.GetPad(1), ratio)
        setBotPad(c1.GetPad(2), ratio)
    c1.cd(1)
    c1.GetPad(bool(ratio)).SetTopMargin(0.06)
    c1.GetPad(bool(ratio)).SetRightMargin(0.05)
    c1.GetPad(bool(ratio)).SetTicks(1, 1)
    if log:
        c1.GetPad(bool(ratio)).SetLogy()
        #c1.GetPad(bool(ratio)).SetLogx()
        
    # Draw
    bkg.Draw("HIST") # stack
    hist['BkgSum'].Draw("SAME, E2") # sum of bkg
    if poisson: data_graph.Draw("SAME, PE")
    elif len(data) > 0: hist['data_obs'].Draw("SAME, PE")
    for i, s in enumerate(sign):
        if samples[s]['plot']:
            hist[s].DrawNormalized("SAME, HIST", hist[s].Integral()*snorm) # signals

    bkg.GetYaxis().SetTitleOffset(bkg.GetYaxis().GetTitleOffset()*1.075)

    # Determine range
    if 'data_obs' in hist:
        bkg.SetMaximum((2.5 if log else 1.2)*max(bkg.GetMaximum(), hist['data_obs'].GetBinContent(hist['data_obs'].GetMaximumBin())+hist['data_obs'].GetBinError(hist['data_obs'].GetMaximumBin())))
        bkg.SetMinimum(max(min(hist['BkgSum'].GetBinContent(hist['BkgSum'].GetMinimumBin()), hist['data_obs'].GetMinimum()), 5.e-1)  if log else 0.)
    else:
        bkg.SetMaximum(bkg.GetMaximum()*(2.5 if log else 1.2))
        bkg.SetMinimum(5.e-1 if log else 0.)
    if log:
        bkg.GetYaxis().SetNoExponent(bkg.GetMaximum() < 1.e4)
        bkg.GetYaxis().SetMoreLogLabels(True)
    
    leg.Draw()
    #drawCMS(LUMI, "Preliminary")
    #drawRegion(channel)
    #drawAnalysis("LL")
    
    setHistStyle(bkg, 1.2 if ratio else 1.1)
    setHistStyle(hist['BkgSum'], 1.2 if ratio else 1.1)

    if ratio:
        c1.cd(2)
        err = hist['BkgSum'].Clone("BkgErr;")
        err.SetTitle("")
        err.GetYaxis().SetTitle("Data / Bkg")
        for i in range(1, err.GetNbinsX()+1):
            err.SetBinContent(i, 1)
            if hist['BkgSum'].GetBinContent(i) > 0:
                err.SetBinError(i, hist['BkgSum'].GetBinError(i)/hist['BkgSum'].GetBinContent(i))
        setBotStyle(err)
        errLine = err.Clone("errLine")
        errLine.SetLineWidth(2)
        errLine.SetFillStyle(0)
        errLine.SetLineColor(2)#L#
        errLine.SetLineStyle(2)#L#
        #err.GetXaxis().SetLabelOffset(err.GetXaxis().GetLabelOffset()*5)
        #err.GetXaxis().SetTitleOffset(err.GetXaxis().GetTitleOffset()*2)
        err.Draw("E2")
        errLine.Draw("SAME, HIST")
        if 'data_obs' in hist:
            res = hist['data_obs'].Clone("Residues")
            for i in range(0, res.GetNbinsX()+1):
                if hist['BkgSum'].GetBinContent(i) > 0: 
                    res.SetBinContent(i, res.GetBinContent(i)/hist['BkgSum'].GetBinContent(i))
                    res.SetBinError(i, res.GetBinError(i)/hist['BkgSum'].GetBinContent(i))
            setBotStyle(res)
            if poisson: res_graph.Draw("SAME, PE0")
            else: res.Draw("SAME, PE0")
            if len(err.GetXaxis().GetBinLabel(1))==0: # Bin labels: not a ordinary plot
                drawRatio(hist['data_obs'], hist['BkgSum'])
                drawKolmogorov(hist['data_obs'], hist['BkgSum'])
        else: res = None
    c1.Update()
    
    # return list of objects created by the draw() function
    return [c1, bkg, leg, err if ratio else None, errLine if ratio else None, res if ratio else None, data_graph if poisson else None, res_graph if poisson else None]




def drawSignal(hist, sign, log=False):
    
    # Legend
    n = len(sign)
    leg = TLegend(0.7, 0.9-0.05*n, 0.95, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    for i, s in enumerate(sign): leg.AddEntry(hist[s], samples[s]['label'], "fl")
    
    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[-1].GetXaxis().GetTitle(), 800, 600)
    
    c1.cd(1)
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.05)
    c1.GetPad(0).SetTicks(1, 1)
    if log:
        c1.GetPad(0).SetLogy()
        
    # Draw
    max_val = 0
    for i, s in enumerate(sign): 
        hist[s].SetLineWidth(3)
        hist[s].Draw("SAME, HIST" if i>0 else "HIST") # signals
        max_val = max(max_val,hist[s].GetMaximum())
        addOverflow(hist[s], True) # Add overflow

    
    #?#hist[sign[0]].GetXaxis().SetRangeUser(0., 1500)
    #?hist[sign[0]].GetYaxis().SetTitleOffset(hist[sign[-1]].GetYaxis().GetTitleOffset()*1.075)
    #?hist[sign[0]].SetMaximum(max(hist[sign[0]].GetMaximum(), hist[sign[-1]].GetMaximum())*1.25)
    #?hist[sign[0]].SetMinimum(0.)

    hist[sign[0]].GetYaxis().SetTitleOffset(hist[sign[-1]].GetYaxis().GetTitleOffset()*1.075)
    hist[sign[0]].SetMaximum(max_val*1.25)
    #hist[sign[0]].SetMinimum(0.)
    
    if log:
        hist[sign[0]].GetYaxis().SetNoExponent(hist[sign[0]].GetMaximum() < 1.e4)
        hist[sign[0]].GetYaxis().SetMoreLogLabels(True)

    if log:
        c1.GetPad(0).SetLogy()

    
    leg.Draw()
    #drawCMS(LUMI, "Preliminary")
    
    c1.Update()
    
    # return list of objects created by the draw() function
    return [c1, leg]

def drawRatio(data, bkg):
    errData = array('d', [1.0])
    errBkg = array('d', [1.0])
    intData = data.IntegralAndError(1, data.GetNbinsX(), errData)
    intBkg = bkg.IntegralAndError(1, bkg.GetNbinsX(), errBkg)
    ratio = intData / intBkg if intBkg!=0 else 0.
    error = math.hypot(errData[0]*ratio/intData,  errBkg[0]*ratio/intBkg) if intData>0 and intBkg>0 else 0
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextFont(62)
    latex.SetTextSize(0.085)
    latex.DrawLatex(0.25, 0.85, "Data/Bkg = %.3f #pm %.3f" % (ratio, error))
    print "  Ratio:\t%.3f +- %.3f" % (ratio, error)
    #return [ratio, error]

def drawKolmogorov(data, bkg, fontsize=0.085):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextFont(62)
    latex.SetTextSize(fontsize)
    latex.DrawLatex(0.55, 0.85, "#chi^{2}/ndf = %.2f,   K-S = %.3f" % (data.Chi2Test(bkg, "CHI2/NDF"), data.KolmogorovTest(bkg)))

def printTable(hist, sign=[], SIGNAL=1):
    samplelist = [x for x in hist.keys() if not 'data' in x and not 'BkgSum' in x and not x in sign and not x=="files"]
    print "Sample                  Events          Entries         %"
    print "-"*80
    for i, s in enumerate(['data_obs']+samplelist+['BkgSum']):
        if i==1 or i==len(samplelist)+1: print "-"*80
        print "%-20s" % s, "\t%-10.2f" % hist[s].Integral(), "\t%-10.0f" % (hist[s].GetEntries()-2), "\t%-10.2f" % (100.*hist[s].Integral()/hist['BkgSum'].Integral()) if hist['BkgSum'].Integral() > 0 else 0, "%"
    print "-"*80
    #for i, s in enumerate(sign):
    for s in sorted(sign):
        if not samples[s]['plot']: continue
        print "%-20s" % s, "\t%-10.2f" % (hist[s].Integral()*SIGNAL), "\t%-10.0f" % (hist[s].GetEntries()-2), "\t%-10.2f" % (100.*hist[s].GetEntries()/float(hist[s].GetOption())) if float(hist[s].GetOption()) > 0 else 0, "%"    
    print "-"*80




##################
#     OTHERS     #
##################

def getPrimaryDataset(cut):
    pd = []
#    if 'HLT_PFMET' in cut: pd += [x for x in samples['data_obs']['files'] if "MET" in x]
#    if 'HLT_' in cut: pd += [x for x in samples['data_obs']['files'] if "MET" in x]
    pd += [x for x in samples['data_obs']['files'] if ("MET" in x or "DisplacedJet" in x or "SingleMuon" in x or "JetHT" in x)]
    return pd


def addOverflow(hist, addUnder=True):
    n = hist.GetNbinsX()
    hist.SetBinContent(n, hist.GetBinContent(n) + hist.GetBinContent(n+1))
    hist.SetBinError(n, math.sqrt( hist.GetBinError(n)**2 + hist.GetBinError(n+1)**2 ) )
    hist.SetBinContent(n+1, 0.)
    hist.SetBinError(n+1, 0.)
    if addUnder:
        hist.SetBinContent(1, hist.GetBinContent(0) + hist.GetBinContent(1))
        hist.SetBinError(1, math.sqrt( hist.GetBinError(0)**2 + hist.GetBinError(1)**2 ) )
        hist.SetBinContent(0, 0.)
        hist.SetBinError(0, 0.)

def setTopPad(TopPad, r=4):
    TopPad.SetPad("TopPad", "", 0., 1./r, 1.0, 1.0, 0, -1, 0)
    TopPad.SetTopMargin(0.24/r)
    TopPad.SetBottomMargin(0.04/r)
    TopPad.SetRightMargin(0.05)
    TopPad.SetTicks(1, 1)

def setBotPad(BotPad, r=4, forcetop=0):
    BotPad.SetPad("BotPad", "", 0., 0., 1.0, 1./r, 0, -1, 0)
    if forcetop==0:
        forcetop = r/100
    BotPad.SetTopMargin(forcetop)
    BotPad.SetBottomMargin(r/10.)
    BotPad.SetRightMargin(0.05)
    BotPad.SetTicks(1, 1)

def setHistStyle(hist, r=1.1):
    hist.GetXaxis().SetTitleSize(hist.GetXaxis().GetTitleSize()*r*r)
    hist.GetYaxis().SetTitleSize(hist.GetYaxis().GetTitleSize()*r*r)
    hist.GetXaxis().SetLabelSize(hist.GetXaxis().GetLabelSize()*r)
    hist.GetYaxis().SetLabelSize(hist.GetYaxis().GetLabelSize()*r)
    hist.GetXaxis().SetLabelOffset(hist.GetXaxis().GetLabelOffset()*r*r*r*r)
    hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset()*r)
    hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset())
    if hist.GetXaxis().GetTitle().find("GeV") != -1: # and not hist.GetXaxis().IsVariableBinSize()
        div = (hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) / hist.GetXaxis().GetNbins()
        hist.GetYaxis().SetTitle("Events / %.1f GeV" % div)

def setBotStyle(h, r=4, fixRange=True, miny=0., maxy=2.):
    h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize()*(r-1));
    h.GetXaxis().SetLabelOffset(h.GetXaxis().GetLabelOffset()*(r-1));
    h.GetXaxis().SetTitleSize(h.GetXaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize()*(r-1));
    h.GetYaxis().SetNdivisions(505);
    h.GetYaxis().SetTitleSize(h.GetYaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()/(r-1));
    if fixRange:
        h.GetYaxis().SetRangeUser(miny, maxy)
        for i in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(i)<1.e-6:
                h.SetBinContent(i, -1.e-6)

##################
### DRAW UTILS ###
##################

def drawCMS(LUMI, text, onTop=False, left_marg_CMS=0.15,data_obs=[]):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)
    if (type(LUMI) is float or type(LUMI) is int) and float(LUMI) > 0: latex.DrawLatex(0.95, 0.985, "%.1f fb^{-1}  (13 TeV)" % (float(LUMI)/1000.))
    elif type(LUMI) is str: latex.DrawLatex(0.95, 0.985, "%s fb^{-1}  (13 TeV)" % LUMI)
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.05 if len(text)>0 else 0.06)
    if not onTop: latex.DrawLatex(left_marg_CMS, 0.87 if len(text)>0 else 0.84, "CMS")
    else:
        latex.DrawLatex(0.20, 0.9, "CMS")#DCMS
    latex.SetTextSize(0.045)
    latex.SetTextFont(52)
    if not onTop:
        latex.DrawLatex(left_marg_CMS, 0.83, text)
    else:
        #latex.DrawLatex(0.40, 0.98, text)
        latex.DrawLatex(0.35, 0.89, text)#DCMS
    dat = ""
    if len(data_obs)>0:
        print samples[data_obs[0]]['files'][0]
        if "SingleMuon" in (samples[data_obs[0]]['files'][0]):
            dat = "SingleLepton dataset"
        elif "SingleElectron" in (samples[data_obs[0]]['files'][0]):
            dat = "SingleLepton dataset"
        elif "DisplacedJet" in (samples[data_obs[0]]['files'][0]):
            dat = "DisplacedJet dataset"
        elif "MET" in (samples[data_obs[0]]['files'][0]):
            dat = "MET dataset"
        print "dat: ", dat
        latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextFont(72) #52
        latex2.SetTextSize(0.04)
        latex2.SetTextAlign(10)
        latex2.DrawLatex(0.45, 0.95, dat)

def drawAnalysis(s, center=False):
    analyses = {
        "LL" : "VBF H #rightarrow #pi #pi #rightarrow b#bar{b} b#bar{b}",
        "LLZH" : "ZH #rightarrow #pi #pi #rightarrow b#bar{b} b#bar{b}",
        "LLVBF" : "VBF H #rightarrow #pi #pi #rightarrow b#bar{b} b#bar{b}",
        "LLggH" : "ggH #rightarrow #pi #pi #rightarrow b#bar{b} b#bar{b}",
        }
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    #latex.SetTextAlign(33)
    latex.DrawLatex(0.15 if not center else 0.3, 0.95, s if not s in analyses else analyses[s])

def drawRegion(channel, left=False, left_marg_CMS=0.15, top=0.75):
    region = {
        "VBFtrigger": "VBF triggers",
        "VBF": "VBF",
        "DisplacedJets" : "Displaced Jets phase-space",
        "MET" : "MET phase-space",
        "BTagCSV" : "BTagCSV phase-space",
        "VBFplusDisplacedTrigger" : "VBF + Displaced Jet trigger",
        "VBFplusDisplacedHadronicTrigger" : "VBF + Displaced hadronic Jet trigger",
        "VBFplusMETTrigger" : "VBF + MET trigger",
        "VBFplusVBFTrigger" : "VBF + VBF trigger",
        "ZtoMMCR": "Z #rightarrow #mu#mu CR",
        "ZtoEECR": "Z #rightarrow ee CR",
        "ZtoMMVBFCR": "VBF + Z #rightarrow #mu#mu CR",
        "ZtoEEVBFCR": "VBF + Z #rightarrow ee CR",
        "WtoMNCR": "W #rightarrow #mu#nu CR",
        "TopEMCR": "t #rightarrow #mue+X CR",
        "DisplacedZtoMMCR": "VBF + Displaced jets + Z #rightarrow #mu#mu CR",
        "DisplacedWtoMNCR": "VBF + Displaced jets + W #rightarrow #mu#nu CR",
        "DisplacedTopEMCR": "VBF + Displaced jets + 1 #mu 1 e",
        "DisplacedCR0Tag" : "VBF + Displaced jets CR (0 calo tag)",
        "DisplacedHadronicCR0Tag" : "VBF + Displaced hadronic jets CR (0 calo tag)",
        "TSG" : "IsoMu24 + 1#mu",
        "L1seed" : "IsoMu24 + 1#mu + L1seed",
        "hltTripleJet50" : "IsoMu24 + 1#mu + L1seed + TripleJet50",
        "VBFplusPFMETNoMuTrigger" : "VBF + PFMETNoMu120 trigger",
        "VBFplusDisplacedHadronicTrigger" : "VBF + Displaced jet trigger",
        "ZtoMM": "ZH #rightarrow #mu#mu H",
        "ZtoEE": "ZH #rightarrow ee H",
        "ZHMM": "ZH #rightarrow #mu#mu H",
        "ZHEE": "ZH #rightarrow ee H",
        "ZH": "ZH #rightarrow ll H",
        "METMuCR": "MEt.pt>250 & HT>200 & 1 muon CR",
        "METCR": "MEt.pt>250 & HT>200 & veto leptons",
        "METHT": "E_{T}^{miss}>200 GeV & H_{T}>200 GeV",
        "METHTVeto": "E_{T}^{miss}>200 GeV & H_{T}>200 GeV & veto #l, #gamma",
        "METHTNoVeto": "MEt.pt>200 & HT>100, no veto",
        "METPreSel": "E_{T}^{miss}>200 GeV & H_{T}>100 GeV & veto #l, #gamma",
        }
    
    text = ""
    if channel in region:
        text = region[channel]
    else:
        text = ""
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(72) #52
    latex.SetTextSize(0.035)
    if left: latex.DrawLatex(left_marg_CMS, top, text)
    else:
        latex.SetTextAlign(10)
        #latex.DrawLatex(0.12, 0.75, text)
        latex.DrawLatex(0.15, top, text)#DCMS

def drawTagVar(tagvar, left=False, left_marg_CMS=0.15):
    tagvarlist = ["nCaloTagJets","nLooseCaloTagJets","nCaloTagJetsRebuilt","nHardCaloTagJets","nLeadingCaloTagJets","nGenMatchedJets","1Loose1Tight"]
    
    text = ""
    if tagvar in tagvarlist:
        text = tagvar
    else:
        text = ""
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(62) #52
    latex.SetTextSize(0.035)
    if left: latex.DrawLatex(left_marg_CMS, 0.6, text)
    else:
        latex.SetTextAlign(10)
        latex.DrawLatex(0.12, 0.7, text)


def drawBox(x1, y1, x2, y2, t="", fillstyle=3005):
    box = TBox(x1, y1, x2, y2)
    box.SetFillColor(1)
    box.SetFillStyle(fillstyle)
    box.Draw()
    if not t=="":
        text = TLatex()
        text.SetTextColor(1)
        text.SetTextFont(42)
        text.SetTextAlign(23)
        text.SetTextSize(0.04)
        text.DrawLatex((x1+x2)/2., y2/1.15, t)
        text.Draw()
    return box

def drawLine(x1, y1, x2, y2,color=1):
    line = TLine(x1, y1, x2, y2)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.SetLineColor(color)
    line.Draw()
    return line

def drawText(x, y, t, col=1):
    text = TLatex()
    text.SetTextColor(col)
    text.SetTextFont(42)
    text.SetTextAlign(23)
    text.SetTextSize(0.04)
    text.DrawLatex(x, y, t)
    text.Draw()
    return text



