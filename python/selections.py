#! /usr/bin/env python

long_string = "("
long_string += "Jets.Jets[0].isGenMatched" #new version from v3
long_string += ")"

selection = {
    "none" : "",
    "isMC" : "isMC",
    "VBF" : "isVBF",
    ##Comment: including only triggers from BTagCSV, DisplacedJet and MET datasets, as per: https://docs.google.com/spreadsheets/d/1oBxzCCM1XP_dfezelrlamR6sfuAdnWm3cHTcaKbt1xA/edit?usp=sharing
    "METfilters" : "(isMC?Flag_eeBadScFilter:1) && (Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_globalTightHalo2016Filter && Flag_goodVertices && Flag_BadPFMuon && Flag_BadChCand)",
    "PFMETNoMuTrigger" : "(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v) && (isMC?1:Flag_EcalDeadCellTriggerPrimitiveFilter) && (isMC?1:Flag_HBHENoiseFilter) && (isMC?1:Flag_HBHENoiseIsoFilter) && (isMC?1:Flag_globalTightHalo2016Filter) && (isMC?1:Flag_goodVertices) && Flag_BadPFMuon && Flag_BadChCand",
    "VBFDisplHadTrigger" : "(HLT_VBF_DisplacedJet40_VTightID_Hadronic_v || HLT_VBF_DisplacedJet40_VVTightID_Hadronic_v) && (isMC?1:Flag_EcalDeadCellTriggerPrimitiveFilter) && (isMC?1:Flag_HBHENoiseFilter) && (isMC?1:Flag_HBHENoiseIsoFilter) && (isMC?1:Flag_globalTightHalo2016Filter) && (isMC?1:Flag_goodVertices) && Flag_BadPFMuon && Flag_BadChCand",
    "VetoLeptons" : "nMuons==0 && nElectrons==0 && nPhotons==0 && nTaus==0",#enriches in QCD
    "L" : "(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v) && Flag_HBHENoiseIsoFilter",
}

selection["VBFplusPFMETNoMuTrigger"] = selection["VBF"] + " && " + selection["PFMETNoMuTrigger"] + " && HT>100"
selection["VBFplusPFMETNoMuTriggerPlateau"] = selection["VBF"] + " && " + selection["PFMETNoMuTrigger"] + " && MEt.pt>250 && HT>100"

selection["METCR"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && MEt.pt>250  && HT>200"
selection["METNoVBFCR"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && !isVBF && MEt.pt>250  && HT>200"
selection["METVBFCR"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && isVBF && MEt.pt>250  && HT>200"
selection["METMuCR"] = selection["PFMETNoMuTrigger"] + " && MEt.pt>250  && HT>200 && nTightMuons==1"

selection["DisplHadPreSel"] = selection["VBFDisplHadTrigger"] + " && " + selection["VetoLeptons"] + " && HT>100 && isVBF"

selection["METPreSel"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && HT>100 && isVBF"#No Met cuts
selection["METPreSel200"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && HT>100 && isVBF && MEt.pt>200"
selection["METPreSel120"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && HT>100 && isVBF && MEt.pt>120"
selection["METPreSel120QCDKiller"] = selection["PFMETNoMuTrigger"] + " && " + selection["VetoLeptons"] + " && HT>100 && isVBF && MEt.pt>120 && MinJetMetDPhi>0.5"


selection["METHTLowPt"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>200"
selection["METHTLowMet"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>120"
selection["METHTLowMetHT"] = selection["PFMETNoMuTrigger"] + " && HT>100 && MEt.pt>120"
selection["METHT"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>200 && CHSJets.pt>15"
selection["METHTNoVBF"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>250 && !isVBF"

selection["METHTv0"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>200 && Jets.pt>15"

selection["METHTv0miniAOD"] = selection["PFMETNoMuTrigger"] + " && HT>200 && MEt.pt>200 && (isMC?CHSJets.pt>15:1)"
