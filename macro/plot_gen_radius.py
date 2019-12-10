#! /usr/bin/env python

import os, multiprocessing
import copy
import math
import numpy as np
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TH1, TH1F, TH2F, THStack, TGraph, TMultiGraph, TGraphAsymmErrors, TSpline, TSpline3
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine, TBox


#from Analyzer.LLPonAOD.samples_v3 import sample, samples
from Analyzer.LLPonAOD.samples import sample, samples
from Analyzer.LLPonAOD.selections import selection
from Analyzer.LLPonAOD.variables import *

gStyle.SetOptStat(0)
NTUPLEDIR   = "/nfs/dust/cms/group/cms-llp/v0_calo_AOD/"
sign = ['VBFH_M40_ctau10000']

def plot_2D(var,nbins=50,minimum=0,maximum=2000,filename=""):
    chain = {}
    hist = {}
    r_ecal = 129
    r_hcal = 179
    r_magnet = 295
    r_mb1 = 402

    z_ecal = 300
    z_hcal = 376
    z_magnet = 0
    z_mb1 = 560

    if var=="radius2D":
        v_ecal = TLine(r_ecal,minimum,r_ecal,maximum)
        v_hcal = TLine(r_hcal,minimum,r_hcal,maximum)
        v_magnet = TLine(r_magnet,minimum,r_magnet,maximum)
        v_mb1 = TLine(r_mb1,minimum,r_mb1,maximum)
        h_ecal = TLine(minimum,r_ecal,maximum,r_ecal)
        h_hcal = TLine(minimum,r_hcal,maximum,r_hcal)
        h_magnet = TLine(minimum,r_magnet,maximum,r_magnet)
        h_mb1 = TLine(minimum,r_mb1,maximum,r_mb1)
    elif var=="z":
        v_ecal = TLine(z_ecal,minimum,z_ecal,maximum)
        v_hcal = TLine(z_hcal,minimum,z_hcal,maximum)
        v_magnet = TLine(z_magnet,minimum,z_magnet,maximum)
        v_mb1 = TLine(z_mb1,minimum,z_mb1,maximum)
        h_ecal = TLine(minimum,z_ecal,maximum,z_ecal)
        h_hcal = TLine(minimum,z_hcal,maximum,z_hcal)
        h_magnet = TLine(minimum,z_magnet,maximum,z_magnet)
        h_mb1 = TLine(minimum,z_mb1,maximum,z_mb1)
    else:
        v_ecal = TLine(r_ecal,minimum,r_ecal,maximum)
        v_hcal = TLine(r_hcal,minimum,r_hcal,maximum)
        v_magnet = TLine(r_magnet,minimum,r_magnet,maximum)
        v_mb1 = TLine(r_mb1,minimum,r_mb1,maximum)
        h_ecal = TLine(minimum,r_ecal,maximum,r_ecal)
        h_hcal = TLine(minimum,r_hcal,maximum,r_hcal)
        h_magnet = TLine(minimum,r_magnet,maximum,r_magnet)
        h_mb1 = TLine(minimum,r_mb1,maximum,r_mb1)

    v_ecal.SetLineColor(2)
    h_ecal.SetLineColor(2)
    v_hcal.SetLineColor(881)
    h_hcal.SetLineColor(881)
    v_magnet.SetLineColor(1)
    h_magnet.SetLineColor(1)
    v_mb1.SetLineColor(801)
    h_mb1.SetLineColor(801)

    v_ecal.SetLineWidth(4)
    h_ecal.SetLineWidth(4)
    v_hcal.SetLineWidth(4)
    h_hcal.SetLineWidth(4)
    v_magnet.SetLineWidth(4)
    h_magnet.SetLineWidth(4)
    v_mb1.SetLineWidth(4)
    h_mb1.SetLineWidth(4)

    v_ecal.SetLineStyle(3)
    h_ecal.SetLineStyle(3)
    v_hcal.SetLineStyle(2)
    h_hcal.SetLineStyle(2)
    v_magnet.SetLineStyle(4)
    h_magnet.SetLineStyle(4)
    v_mb1.SetLineStyle(8)
    h_mb1.SetLineStyle(8)

    leg = TLegend(0.75, 0.75, 0.9, 0.9)
    leg.AddEntry(v_ecal,"ECAL","L")
    leg.AddEntry(v_hcal,"HCAL","L")
    leg.AddEntry(v_magnet,"solenoid","L")
    leg.AddEntry(v_mb1,"MB1","L")

    #pal= 68 #kAvocado
    #pal= 64 #kAquamarine, very readable
    #pal= 75 #kCherry 75, awful
    #pal= 85 #kIsland 85, not beautiful but readable
    #pal= 86 #kLake 86, too violet
    #pal= 87 #kLightTemperature 87, used for trigger
    #pal= 91 #kPastel 91, too purple
    #pal= 100 #kSolar 100, very red and orange
    pal= 98 #kSandyTerrain 98, quite fine
    #pal= 99 #kSienna 99, a bit hard to read
    gStyle.SetPalette(pal)
    gStyle.SetPaintTextFormat(".0f")
    for i, s in enumerate(sign):
        chain[s] = TChain("ntuple/tree")
        if filename=="":
            for p, ss in enumerate(samples[s]['files']):
                chain[s].Add(NTUPLEDIR + ss + ".root")
        else:
            chain[s].Add(filename+".root")
        #filename[s] = TFile("VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000.root", "READ")
        hist[s] = TH2F(s, "", nbins, minimum, maximum, nbins, minimum, maximum)
        hist[s].Sumw2()
        cutstring = "(fabs(GenBquarks[0].eta)<2.4 && fabs(GenBquarks[3].eta)<2.4)*(EventWeight)"
        cutstring = "(fabs(GenBquarks[0].eta)<2.4 && fabs(GenBquarks[3].eta)<2.4 && HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v)*(EventWeight)"
        #cutstring = "(fabs(GenBquarks[0].eta)<2.4 && fabs(GenBquarks[3].eta)<2.4 && HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v && MEt.pt>200)*(EventWeight)"
        #cutstring = "(fabs(GenBquarks[0].eta)<2.4 && fabs(GenBquarks[3].eta)<2.4 && HLT_VBF_DisplacedJet40_VTightID_Hadronic_v)*(EventWeight)"
        #cutstring = "(fabs(GenBquarks[0].eta)<2.4 && fabs(GenBquarks[3].eta)<2.4 && HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack_v)*(EventWeight)"
        if var=="z":
            #chain[s].Project(s, "sqrt(pow(GenBquarks[0].radius,2) - pow(GenBquarks[0].radius2D,2)) * GenBquarks[0].eta/abs(GenBquarks[0].eta):sqrt(pow(GenBquarks[3].radius,2) - pow(GenBquarks[3].radius2D,2)) * GenBquarks[0].eta/abs(GenBquarks[0].eta)", cutstring)
            #sign of eta for getting the right z value!
            chain[s].Project(s, "sqrt(pow(GenBquarks[0].radius,2) - pow(GenBquarks[0].radius2D,2)):sqrt(pow(GenBquarks[3].radius,2) - pow(GenBquarks[3].radius2D,2))", cutstring)
        else:
            chain[s].Project(s, "GenBquarks[0]."+var+":GenBquarks[3]."+var+"", cutstring)
        hist[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
        c1 = TCanvas("c1", "c1", 1000, 1000)
        c1.cd()
        c1.SetGrid()
        #c1.SetLogx()
        #c1.SetLogy()
        hist[s].GetYaxis().SetTitle("GenBquarks[0] "+var+" (cm)")
        hist[s].GetYaxis().SetTitleOffset(1.4)
        hist[s].GetXaxis().SetTitle("GenBquarks[3] "+var+" (cm)")
        hist[s].SetTitle(samples[s]['label'] if filename=="" else filename)
        hist[s].SetMarkerColor(0)
        hist[s].Draw("colztext")
        v_ecal.Draw("sames")
        h_ecal.Draw("sames")
        v_hcal.Draw("sames")
        h_hcal.Draw("sames")
        v_magnet.Draw("sames")
        h_magnet.Draw("sames")
        v_mb1.Draw("sames")
        h_mb1.Draw("sames")
        leg.Draw("sames")
        c1.Print("macro/2D_gen_b_quark_"+var+"_"+(s if filename=="" else filename)+".png")
        c1.Print("macro/2D_gen_b_quark_"+var+"_"+(s if filename=="" else filename)+".pdf")

        raw_input("Press Enter to continue...")
        c1.Close()

#plot_2D("radius",nbins=20,minimum=0,maximum=1000,filename="")
plot_2D("radius2D",nbins=20,minimum=0,maximum=1000,filename="VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000")
#plot_2D("z",nbins=10,minimum=0,maximum=1000,filename="VBFH_HToSSTobbbb_MH-125_MS-30_ctauS-1000")
