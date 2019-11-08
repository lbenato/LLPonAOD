
#include "ObjectsFormat.h"

//*******************//
//  Leptons (e+mu)   //
//*******************//

void ObjectsFormat::FillElectronType(LeptonType& I, const pat::Electron* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.inTrkPt     = R->hasUserFloat("inTrkPt") ? R->userFloat("inTrkPt") : -1.;
    I.pfIso03     = R->hasUserFloat("pfIso03") ? R->userFloat("pfIso03") : -1.;
    I.pfIso04     = R->hasUserFloat("pfIso04") ? R->userFloat("pfIso04") : -1.;
    I.trkIso      = R->hasUserFloat("trkIso") ? R->userFloat("trkIso") : R->pfIsolationVariables().sumChargedHadronPt;
    I.miniIso     = R->hasUserFloat("miniIso") ? R->userFloat("miniIso") : -1.;
    I.dxy         = R->hasUserFloat("dxy") ? R->userFloat("dxy") : R->dB();
    I.dz          = R->hasUserFloat("dz") ? R->userFloat("dz") : 0.;
    I.ip3d        = R->dB(pat::Electron::PV3D);
    I.sip3d       = R->dB(pat::Electron::PV3D)/R->edB(pat::Electron::PV3D);
    I.nPixelHits  = R->gsfTrack()->hitPattern().numberOfValidPixelHits();
    I.dPhi_met    = R->hasUserFloat("dPhi_met") ? R->userFloat("dPhi_met") : -1.;
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.isElectron  = true;
    I.isMuon      = false;
    I.isVeto      = R->hasUserInt("isVeto") ? R->userInt("isVeto") : false;
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = R->hasUserInt("isMedium") ? R->userInt("isMedium") : false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isHighPt    = R->hasUserInt("isHEEP") ? R->userInt("isHEEP") : false;

    I.SSscale               = R->hasUserFloat("SSscale") ? R->userFloat("SSscale") : -1.;
    I.SSsigma               = R->hasUserFloat("SSsigma") ? R->userFloat("SSsigma") : -1.;
    I.SSscaleUnc            = R->hasUserFloat("SSscaleUnc") ? R->userFloat("SSscaleUnc") : -1.;
    I.SSsigmaUncUp          = R->hasUserFloat("SSsigmaUncUp") ? R->userFloat("SSsigmaUncUp") : -1.;
    I.SSsigmaUncDown        = R->hasUserFloat("SSsigmaUncDown") ? R->userFloat("SSsigmaUncDown") : -1.;
    I.SScorr                = R->hasUserFloat("SScorr") ? R->userFloat("SScorr") : -1.;
    I.energySScorr          = R->hasUserFloat("energySScorr") ? R->userFloat("energySScorr") : -1.;
    I.energySScorrUncUp     = R->hasUserFloat("energySScorrUncUp") ? R->userFloat("energySScorrUncUp") : -1.;
    I.energySScorrUncDown   = R->hasUserFloat("energySScorrUncDown") ? R->userFloat("energySScorrUncDown") : -1.;
    I.ptSScorr              = R->hasUserFloat("ptSScorr") ? R->userFloat("ptSScorr") : -1.;
    I.ptSScorrUncUp         = R->hasUserFloat("ptSScorrUncUp") ? R->userFloat("ptSScorrUncUp") : -1.;
    I.ptSScorrUncDown       = R->hasUserFloat("ptSScorrUncDown") ? R->userFloat("ptSScorrUncDown") : -1.;

    if(isMC && R->genLepton()) I.isMatched = false;//(Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genLepton()))==23);
//    I.isMVANonTrigMedium      = R->hasUserInt("isMVANonTrigMedium") ? R->userInt("isMVANonTrigMedium") : false;
//    I.isMVANonTrigTight      = R->hasUserInt("isMVANonTrigTight") ? R->userInt("isMVANonTrigTight") : false;
//    I.isMVATrigMedium      = R->hasUserInt("isMVATrigMedium") ? R->userInt("isMVATrigMedium") : false;
//    I.isMVATrigTight      = R->hasUserInt("isMVATrigTight") ? R->userInt("isMVATrigTight") : false;
}

void ObjectsFormat::FillMuonType(LeptonType& I, const pat::Muon* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.inTrkPt     = R->hasUserFloat("inTrkPt") ? R->userFloat("inTrkPt") : -1.;
    I.pfIso03     = R->hasUserFloat("pfIso03") ? R->userFloat("pfIso03") : -1.;
    I.pfIso04     = R->hasUserFloat("pfIso04") ? R->userFloat("pfIso04") : -1.;
    I.trkIso      = R->hasUserFloat("trkIso") ? R->userFloat("trkIso") : R->trackIso();
    I.miniIso     = R->hasUserFloat("miniIso") ? R->userFloat("miniIso") : -1.;
    I.dxy         = R->hasUserFloat("dxy") ? R->userFloat("dxy") : R->dB();
    I.dz          = R->hasUserFloat("dz") ? R->userFloat("dz") : 0.;
    I.ip3d        = R->dB(pat::Muon::PV3D);
    I.sip3d       = R->dB(pat::Muon::PV3D)/R->edB(pat::Muon::PV3D);
    I.nPixelHits  = R->innerTrack()->hitPattern().numberOfValidPixelHits();
    I.dPhi_met    = R->hasUserFloat("dPhi_met") ? R->userFloat("dPhi_met") : -1.;
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.isElectron  = false;
    I.isMuon      = true;
    I.isVeto      = R->isPFMuon();
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = R->hasUserInt("isMedium") ? R->userInt("isMedium") : false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isHighPt    = R->hasUserInt("isHighPt") ? R->userInt("isHighPt") : false;
    I.isTrackerHighPt = R->hasUserInt("isTrackerHighPt") ? R->userInt("isTrackerHighPt") : false;

    if(isMC && R->genLepton()) I.isMatched = false;//(Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genLepton()))==23);
}

void ObjectsFormat::ResetLeptonType(LeptonType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.inTrkPt     = -1.;
    I.pfIso03     = -1.;
    I.pfIso04     = -1.;
    I.trkIso      = -1.;
    I.miniIso     = -1.;
    I.dxy         = -99.;
    I.dz          = -99.;
    I.ip3d        = -99.;
    I.sip3d       = -99.;
    I.nPixelHits  = -1.;
    I.dPhi_met    = -1.;
    I.charge      = 0;
    I.pdgId       = 0;
    I.isElectron  = false;
    I.isMuon      = false;
    I.isVeto      = false;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isHighPt    = false;
    I.isTrackerHighPt = false;
//    I.isMVANonTrigMedium = false;
//    I.isMVANonTrigTight = false;
//    I.isMVATrigMedium = false;
//    I.isMVATrigTight = false;

    I.SSscale               = -1.;
    I.SSsigma               = -1.;
    I.SSscaleUnc            = -1.;
    I.SSsigmaUncUp          = -1.;
    I.SSsigmaUncDown        = -1.;
    I.SScorr                = -1.;
    I.energySScorr          = -1.;
    I.energySScorrUncUp     = -1.;
    I.energySScorrUncDown   = -1.;
    I.ptSScorr              = -1.;
    I.ptSScorrUncUp         = -1.;
    I.ptSScorrUncDown       = -1.;

    I.isMatched   = false;
}

std::string ObjectsFormat::ListLeptonType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:inTrkPt/F:pfIso03/F:pfIso04/F:trkIso/F:miniIso/F:dxy/F:dz/F:ip3d/F:sip3d/F:nPixelHits/F:dPhi_met/F:charge/I:pdgId/I:isElectron/O:isMuon/O:isVeto/O:isLoose/O:isMedium/O:isTight/O:isHighPt/O:isTrackerHighPt/O:SSscale/F:SSsigma/F:SSscaleUnc/F:SSsigmaUncUp/F:SSsigmaUncDown/F:SScorr/F:energySScorr/F:energySScorrUncUp/F:energySScorrUncDown/F:ptSScorr/F:ptSScorrUncUp/F:ptSScorrUncDown/F:isMatched/O";} // isHEEP/O:isMVANonTrigMedium/O:isMVANonTrigTight/O:isMVATrigMedium/O:isMVATrigTight/O:

//********************//
//    Photons         //
//********************//

void ObjectsFormat::FillPhotonType(PhotonType& I, const pat::Photon* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.pfIso       = R->hasUserFloat("pfIso") ? R->userFloat("pfIso") : -1.;
    I.dz          = R->hasUserFloat("dz") ? R->userFloat("dz") : 0.;
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = R->hasUserInt("isMedium") ? R->userInt("isMedium") : false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isMVANonTrigMedium      = R->hasUserInt("isMVANonTrigMedium") ? R->userInt("isMVANonTrigMedium") : false;
    if(isMC && R->genPhoton()) I.isMatched = false;//(Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genLepton()))==23);
}


void ObjectsFormat::ResetPhotonType(PhotonType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.pfIso       = -1.;
    I.dz          = -99.;
    I.charge      = 0;
    I.pdgId       = 0;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isMVANonTrigMedium = false;
    I.isMatched   = false;
}

std::string ObjectsFormat::ListPhotonType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:pfIso/F:dz/F:charge/I:pdgId/I:isLoose/O:isMedium/O:isTight/O:isMVANonTrigMedium/O:isMatched/O";}


//********************//
//       Taus         //
//********************//

void ObjectsFormat::FillTauType(TauType& I, const pat::Tau* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.pfIso       = R->hasUserFloat("pfIso") ? R->userFloat("pfIso") : -1.;
    I.dz          = R->hasUserFloat("dz") ? R->userFloat("dz") : 0.;
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = R->hasUserInt("isMedium") ? R->userInt("isMedium") : false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    if(isMC) I.isMatched = false;//(Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genLepton()))==23);
}


void ObjectsFormat::ResetTauType(TauType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.pfIso       = -1.;
    I.dz          = -99.;
    I.charge      = 0;
    I.pdgId       = 0;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isMatched   = false;
}

std::string ObjectsFormat::ListTauType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:pfIso/F:dz/F:charge/I:pdgId/I:isLoose/O:isMedium/O:isTight/O:isMatched/O";}


//*******************//
//        Jets       //
//*******************//

void ObjectsFormat::FillJetType(JetType& I, const pat::Jet* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    if(isMC && R->genJet()) {
      I.ptGenJ    = R->genJet()->pt();
      I.etaGenJ   = R->genJet()->eta();
      I.phiGenJ   = R->genJet()->phi();
      I.massGenJ  = R->genJet()->mass();
    }
    if(isMC && R->genParton()) {
      I.ptGen     = R->genParton()->pt();
      I.etaGen    = R->genParton()->eta();
      I.phiGen    = R->genParton()->phi();
      I.massGen   = R->genParton()->mass();
      I.pdgIdGen   = R->genParton()->pdgId();
    }
    //if(isMC && R->genParton()) {
    //I.ptLhe     = R->userFloat("ptLhe");
    //I.etaLhe    = R->userFloat("etaLhe");
    //I.phiLhe    = R->userFloat("phiLhe");
    //}
    I.energyRaw   = R->correctedJet(0).energy();
    I.ptRaw       = R->correctedJet(0).pt();
    I.ptUnc       = R->hasUserFloat("JESUncertainty") ? R->userFloat("JESUncertainty") : -1;
    //I.dPhi_met    = R->hasUserFloat("dPhi_met") ? R->userFloat("dPhi_met") : -1.;
    //I.dPhi_Jet1   = R->hasUserFloat("dPhi_Jet1") ? R->userFloat("dPhi_Jet1") : -1.;
    //I.puId        = R->hasUserFloat("pileupJetId:fullDiscriminant")? R->userFloat("pileupJetId:fullDiscriminant") : -1.; //-1.; // FIXME
    I.CSV         = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVR        = R->hasUserFloat("ReshapedDiscriminator") ? R->userFloat("ReshapedDiscriminator") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVRUp      = R->hasUserFloat("ReshapedDiscriminatorUp") ? R->userFloat("ReshapedDiscriminatorUp") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVRDown    = R->hasUserFloat("ReshapedDiscriminatorDown") ? R->userFloat("ReshapedDiscriminatorDown") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.CMVA        = R->bDiscriminator("pfCombinedMVAV2BJetTags");
    //I.CMVAR       = R->hasUserFloat("CMVAR") ? R->userFloat("CMVAR") : -99.;
    //I.CMVARUp     = R->hasUserFloat("CMVARUp") ? R->userFloat("CMVARUp") : -99.;
    //I.CMVARDown   = R->hasUserFloat("CMVARDown") ? R->userFloat("CMVARDown") : -99.;
  //  I.CSVV1       = R->bDiscriminator("combinedSecondaryVertexV1BJetTags");
  //  I.CSVSL       = R->bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags");
  //  I.JPro        = R->bDiscriminator("jetProbabilityBJetTags"); // jetBProbabilityBJetTags
    //I.QGLikelihood = R->hasUserFloat("QGLikelihood") ? R->userFloat("QGLikelihood") : -1.;
  //  if(isMC) {
  //    if(abs(R->partonFlavour())==5 || abs(R->partonFlavour())==4) {
  //      I.LooseSF      = theBTagWeight->ReturnScaleFactor(1, R->pt(), R->eta());
  //      I.MediumSF     = theBTagWeight->ReturnScaleFactor(2, R->pt(), R->eta());
  //      I.TightSF      = theBTagWeight->ReturnScaleFactor(3, R->pt(), R->eta());
  //      I.LooseSFa     = (1.-0.7*theBTagWeight->ReturnScaleFactor(1, R->pt(), R->eta()))/(1.-0.7);
  //      I.MediumSFa    = (1.-0.6*theBTagWeight->ReturnScaleFactor(2, R->pt(), R->eta()))/(1.-0.6);
  //      I.TightSFa     = (1.-0.5*theBTagWeight->ReturnScaleFactor(3, R->pt(), R->eta()))/(1.-0.5);
  //    }
  //    else {
  //      I.LooseSF      = theBTagWeight->ReturnScaleFactorMistag(1, R->pt(), R->eta());
  //      I.MediumSF     = theBTagWeight->ReturnScaleFactorMistag(2, R->pt(), R->eta());
  //      I.TightSF      = theBTagWeight->ReturnScaleFactorMistag(3, R->pt(), R->eta());
  //      I.LooseSFa     = (1.-1.e-1*theBTagWeight->ReturnScaleFactorMistag(1, R->pt(), R->eta()))/(1.-1.e-1);
  //      I.MediumSFa    = (1.-1.e-2*theBTagWeight->ReturnScaleFactorMistag(2, R->pt(), R->eta()))/(1.-1.e-2);
  //      I.TightSFa     = (1.-1.e-3*theBTagWeight->ReturnScaleFactorMistag(3, R->pt(), R->eta()))/(1.-1.e-3);
  //    }
  //  }

    I.cHadE       = R->chargedHadronEnergy();
    I.nHadE       = R->neutralHadronEnergy();
    I.eleE        = R->electronEnergy();
    I.photonE     = R->photonEnergy();
    I.muE         = R->muonEnergy();
    I.nEmE        = R->neutralEmEnergy();
    I.cEmE        = R->chargedEmEnergy();
    I.cmuE        = R->chargedMuEnergy();
    I.cHadEFrac   = R->hasUserFloat("cHadEFrac") ? R->userFloat("cHadEFrac") : -1.;//R->chargedHadronEnergyFraction();
    I.nHadEFrac   = R->hasUserFloat("nHadEFrac") ? R->userFloat("nHadEFrac") : -1.;//R->neutralHadronEnergyFraction();
    I.eleEFrac    = R->hasUserFloat("eleEFrac") ? R->userFloat("eleEFrac") : -1.;//R->electronEnergyFraction();
    I.photonEFrac = R->hasUserFloat("photonEFrac") ? R->userFloat("photonEFrac") : -1.;//R->photonEnergyFraction();
    I.muEFrac     = R->hasUserFloat("muEFrac") ? R->userFloat("muEFrac") : -1.;//R->muonEnergyFraction();
    I.nEmEFrac    = R->hasUserFloat("nEmEFrac") ? R->userFloat("nEmEFrac") : -1.;//R->neutralEmEnergyFraction();
    I.cEmEFrac    = R->hasUserFloat("cEmEFrac") ? R->userFloat("cEmEFrac") : -1.;//R->chargedEmEnergyFraction();
    I.cmuEFrac    = R->hasUserFloat("cmuEFrac") ? R->userFloat("cmuEFrac") : -1.;//R->chargedMuEnergyFraction();
    I.cHadMulti   = R->chargedHadronMultiplicity();
    I.nHadMulti   = R->neutralHadronMultiplicity();
    I.eleMulti    = R->electronMultiplicity();
    I.photonMulti = R->photonMultiplicity();
    I.muMulti     = R->muonMultiplicity();
    I.cMulti      = R->chargedMultiplicity();
    I.nMulti      = R->neutralMultiplicity();
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.cHadMultiFrac   = float(R->chargedHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nHadMultiFrac   = float(R->neutralHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.eleMultiFrac    = float(R->electronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.photonMultiFrac = float(R->photonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.muMultiFrac     = float(R->muonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.cMultiFrac      = float(R->chargedMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nMultiFrac      = float(R->neutralMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    //I.emEFrac     = -1.;//R->emEnergyFraction();
    //I.emEinEB     = -1.;//R->emEnergyInEB();
    //I.emEinEE     = -1.;//R->emEnergyInEE();
    //I.emEinHF     = -1.;//R->emEnergyInHF();
    //I.EFracHad    = -1.;//R->energyFractionHadronic();
    //I.hadEinHB    = -1.;//R->hadEnergyInHB();
    //I.hadEinHE    = -1.;//R->hadEnergyInHE();
    //I.hadEinHF    = -1.;//R->hadEnergyInHF();
    //I.hadEinHO    = -1.;//R->hadEnergyInHO();
    I.partonFlavour     = R->partonFlavour();
    I.hadronFlavour     = R->hadronFlavour();
    if(isMC && R->genParton()) I.mother = Utilities::FindMotherId(R->genParton());
    //I.isMatched   = (I.mother==25);
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isTightLepVeto     = R->hasUserInt("isTightLepVeto") ? R->userInt("isTightLepVeto") : false;
    I.isCSVL      = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.5426 ? true : false;
    I.isCSVM      = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.8484 ? true : false;
    I.isCSVT      = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>0.9535 ? true : false;
    //I.isMatched   = false;
    //I.dR_q1       = R->hasUserFloat("dR_q1") ? R->userFloat("dR_q1") : 1000;
    //I.dR_q2       = R->hasUserFloat("dR_q2") ? R->userFloat("dR_q2") : 1000;
    //I.dR_q3       = R->hasUserFloat("dR_q3") ? R->userFloat("dR_q3") : 1000;
    //I.dR_q4       = R->hasUserFloat("dR_q4") ? R->userFloat("dR_q4") : 1000;
    //I.m_q1        = R->hasUserFloat("dR_q1") ? (R->userFloat("dR_q1")<0.4 ? true : false) : false;
    //I.m_q2        = R->hasUserFloat("dR_q2") ? (R->userFloat("dR_q2")<0.4 ? true : false) : false;
    //I.m_q3        = R->hasUserFloat("dR_q3") ? (R->userFloat("dR_q3")<0.4 ? true : false) : false;
    //I.m_q4        = R->hasUserFloat("dR_q4") ? (R->userFloat("dR_q4")<0.4 ? true : false) : false;
    //I.dR_pi1      = R->hasUserFloat("dR_pi1") ? R->userFloat("dR_pi1") : 1000;
    //I.dR_pi2      = R->hasUserFloat("dR_pi2") ? R->userFloat("dR_pi2") : 1000;
    I.matchBquark = R->hasUserInt("hasMatchedBquarks") ? R->userInt("hasMatchedBquarks") : -1;
    I.matchLL     = R->hasUserInt("hasMatchedLL") ? R->userInt("hasMatchedLL") : -1;
    //I.original_jet_index     = R->hasUserInt("original_jet_index") ? R->userInt("original_jet_index") : -1;
    I.isGenMatched = R->hasUserInt("isGenMatched") ? R->userInt("isGenMatched") : 0;
    I.isVBFGenMatched = R->hasUserInt("isVBFGenMatched") ? R->userInt("isVBFGenMatched") : 0;
    I.alphaMax    = R->hasUserFloat("alphaMax") ? R->userFloat("alphaMax") : -1000;
    I.sigIP2DMedian = R->hasUserFloat("sigIP2DMedian") ? R->userFloat("sigIP2DMedian") : -1000;
    I.theta2DMedian = R->hasUserFloat("theta2DMedian") ? R->userFloat("theta2DMedian") : -1000;
    I.POCA_theta2DMedian = R->hasUserFloat("POCA_theta2DMedian") ? R->userFloat("POCA_theta2DMedian") : -1000;
    I.nPixelHitsMedian = R->hasUserFloat("nPixelHitsMedian") ? R->userFloat("nPixelHitsMedian") : -1.0;
    I.nHitsMedian = R->hasUserFloat("nHitsMedian") ? R->userFloat("nHitsMedian") : -1.0;
    I.hcalE       = R->chargedHadronEnergy() + R->neutralHadronEnergy();
    I.ecalE       = R->chargedEmEnergy() + R->neutralEmEnergy();
    I.FracCal     = (R->chargedEmEnergy() + R->neutralEmEnergy())/(R->chargedHadronEnergy() + R->neutralHadronEnergy());
    if (R->tagInfoLabels().size() > 0){
      //      if(R->hasTagInfo("pfImpactParameter")){
	//      const reco::BaseTagInfo* IPTagInfo = R->tagInfo("pfImpactParameterTagInfos");
	//	std::cout << "ip2d " << R->tagInfoCandIP("pfImpactParameter")->impactParameterData()[0].ip2d.value() << std::endl;
      //      }
      if(R->hasTagInfo("pfSecondaryVertex")){
	I.nSV = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertices();
	I.nSVCand = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexCandidates();
	I.nVertexTracks = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexTracks();
	I.nSelectedTracks = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nSelectedTracks();
	if (R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->taggingVariables().getList(reco::btau::vertexJetDeltaR,false).size() > 0) I.dRSVJet =  R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->taggingVariables().getList(reco::btau::vertexJetDeltaR,false)[0];
	if(R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertices() > 0){
	  I.flightDist2d  = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0, true).value(); //TODO: here is only first vertex used, should be done like:? (https://github.com/cms-btv-pog/RecoBTag-PerformanceMeasurements/blob/b6d8bf40bb046bfa177e6e9b8c305256c3983e4c/plugins/BTagAnalyzer.cc#L3288)
	  I.flightDist2dError = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0, true).error(); //TODO: see above
	  I.flightDist3d  = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0, false).value(); //TODO: see above
	  I.flightDist3dError = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0, false).error(); //TODO: see above
	  I.SV_x = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vx(); //TODO: see above
	  I.SV_y = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vy(); //TODO: see above
	  I.SV_z = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vz(); //TODO: see above
	  I.SV_dx = -1.;//R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vx(); //TODO: error on vx not available
	  I.SV_dy = -1.;//R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vy(); //TODO: error on vy not available
	  I.SV_dz = -1.;//R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).vz(); //TODO: error on vz not available
	  I.nTracksSV = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexTracks(0);
	  I.SV_mass = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).mass();
	}
      }
    }
    I.isCaloTag   = R->hasUserInt("isCaloTag") ? R->userInt("isCaloTag") : 0;
    //I.VBF_DisplacedJet40_VTightID_Hadronic_match = R->hasUserInt("VBF_DisplacedJet40_VTightID_Hadronic_match") ? R->userInt("VBF_DisplacedJet40_VTightID_Hadronic_match") : 0;
    //I.VBF_DisplacedJet40_VVTightID_Hadronic_match = R->hasUserInt("VBF_DisplacedJet40_VVTightID_Hadronic_match") ? R->userInt("VBF_DisplacedJet40_VVTightID_Hadronic_match") : 0;
    I.ptJESUp = R->hasUserFloat("ptJESUp") ? R->userFloat("ptJESUp") : -1.;
    I.ptJESDown = R->hasUserFloat("ptJESDown") ? R->userFloat("ptJESDown") : -1.;
    I.ptJERUp = R->hasUserFloat("ptJERUp") ? R->userFloat("ptJERUp") : -1.;
    I.ptJERDown = R->hasUserFloat("ptJERDown") ? R->userFloat("ptJERDown") : -1.;
    I.tau1 = R->hasUserFloat("tau1") ? R->userFloat("tau1") : -1.;
    I.tau2 = R->hasUserFloat("tau2") ? R->userFloat("tau2") : -1.;
    I.tau3 = R->hasUserFloat("tau3") ? R->userFloat("tau3") : -1.;
    I.nSubJets = R->hasUserFloat("nSubJets") ? R->userFloat("nSubJets") : -1;
    if (R->hasUserFloat("tau1") && R->hasUserFloat("tau2") && R->hasUserFloat("tau3")){
      I.tau21 = R->userFloat("tau2")/R->userFloat("tau1");
      I.tau31 = R->userFloat("tau3")/R->userFloat("tau1");
      I.tau32 = R->userFloat("tau3")/R->userFloat("tau2");
    }
    I.tau1_neutral = R->hasUserFloat("tau1_neutral") ? R->userFloat("tau1_neutral") : -1.;
    I.tau2_neutral = R->hasUserFloat("tau2_neutral") ? R->userFloat("tau2_neutral") : -1.;
    I.tau1_charged = R->hasUserFloat("tau1_charged") ? R->userFloat("tau1_charged") : -1.;
    I.tau2_charged = R->hasUserFloat("tau2_charged") ? R->userFloat("tau2_charged") : -1.;
    if (R->hasUserFloat("tau1_neutral") && R->hasUserFloat("tau2_neutral"))
      I.tau21_neutral = R->userFloat("tau2_neutral")/R->userFloat("tau1_neutral");
    if (R->hasUserFloat("tau1_charged") && R->hasUserFloat("tau2_charged"))
      I.tau21_charged =R->userFloat("tau2_charged")/R->userFloat("tau1_charged");
    I.TriggerMatched_VBFJet = R->hasUserInt("TriggerMatched_VBFJet") ? R->userInt("TriggerMatched_VBFJet") : 0;
    I.TriggerMatched_DisplacedJet = R->hasUserInt("TriggerMatched_DisplacedJet") ? R->userInt("TriggerMatched_DisplacedJet") : 0;
    I.TriggerMatched_TripleJet50 = R->hasUserInt("TriggerMatched_TripleJet50") ? R->userInt("TriggerMatched_TripleJet50") : 0;
    I.nConstituents = R->hasUserInt("nConstituents") ? R->userInt("nConstituents") : -1;
    I.nTrackConstituents = R->hasUserInt("nTrackConstituents") ? R->userInt("nTrackConstituents") : -1;
    I.nTracks0PixelHits = R->hasUserInt("nTracks0PixelHits") ? R->userInt("nTracks0PixelHits") : -1;
    I.nTracks1PixelHit  = R->hasUserInt("nTracks1PixelHit")  ? R->userInt("nTracks1PixelHit")  : -1;
    I.nTracks2PixelHits = R->hasUserInt("nTracks2PixelHits") ? R->userInt("nTracks2PixelHits") : -1;
    I.nTracks3PixelHits = R->hasUserInt("nTracks3PixelHits") ? R->userInt("nTracks3PixelHits") : -1;
    I.nTracks4PixelHits = R->hasUserInt("nTracks4PixelHits") ? R->userInt("nTracks4PixelHits") : -1;
    I.nTracks5PixelHits = R->hasUserInt("nTracks5PixelHits") ? R->userInt("nTracks5PixelHits") : -1;
    I.nTracksAtLeast6PixelHits = R->hasUserInt("nTracksAtLeast6PixelHits") ? R->userInt("nTracksAtLeast6PixelHits") : -1;
    I.nTracksValidHitInBPix1 = R->hasUserInt("nTracksValidHitInBPix1") ? R->userInt("nTracksValidHitInBPix1") : -1;
    I.nTracks0LostInnerHits = R->hasUserInt("nTracks0LostInnerHits") ? R->userInt("nTracks0LostInnerHits") : -1;
    I.nTracks1LostInnerHit = R->hasUserInt("nTracks1LostInnerHit") ? R->userInt("nTracks1LostInnerHit") : -1;
    I.nTracksAtLeast2LostInnerHits = R->hasUserInt("nTracksAtLeast2LostInnerHits") ? R->userInt("nTracksAtLeast2LostInnerHits") : -1;
    I.nMatchedGenBquarks = R->hasUserInt("nMatchedGenBquarks") ? R->userInt("nMatchedGenBquarks") : -1;
}

void ObjectsFormat::ResetJetType(JetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.energyRaw   = -1.;
    I.ptRaw       = -1.;
    I.ptUnc       = -1.;
    //I.dPhi_met    = -1.;
    //I.dPhi_Jet1   = -1.;
    //I.puId        = -1.;
    I.CSV         = -99.;
    //I.CSVR        = -99.;
    //I.CSVRUp      = -99.;
    //I.CSVRDown    = -99.;
    I.CMVA        = -99.;
    //I.CMVAR       = -99.;
    //I.CMVARUp     = -99.;
    //I.CMVARDown   = -99.;
    //I.QGLikelihood = -1.;
    // I.chf         = -1.;
    // I.nhf         = -1.;
    // I.phf         = -1.;
    // I.elf         = -1.;
    // I.muf         = -1.;
    I.cHadE       = -1.;
    I.nHadE       = -1.;
    I.eleE        = -1.;
    I.photonE     = -1.;
    I.muE         = -1.;
    I.nEmE        = -1.;
    I.cEmE        = -1.;
    I.cmuE        = -1.;
    I.cHadEFrac   = -1.;
    I.nHadEFrac   = -1.;
    I.eleEFrac    = -1.;
    I.photonEFrac = -1.;
    I.muEFrac     = -1.;
    I.nEmEFrac    = -1.;
    I.cEmEFrac    = -1.;
    I.cmuEFrac    = -1.;
    I.cHadMulti   = -1.;
    I.nHadMulti   = -1.;
    I.eleMulti    = -1.;
    I.photonMulti = -1.;
    I.muMulti     = -1.;
    I.cMulti      = -1.;
    I.nMulti      = -1.;
    I.npr         = -1.;
    I.cHadMultiFrac   = -1.;
    I.nHadMultiFrac   = -1.;
    I.eleMultiFrac    = -1.;
    I.photonMultiFrac = -1.;
    I.muMultiFrac     = -1.;
    I.cMultiFrac      = -1.;
    I.nMultiFrac      = -1.;

    //I.emEFrac       = -1.;
    //I.emEinEB       = -1.;
    //I.emEinEE       = -1.;
    //I.emEinHF       = -1.;
    //I.EFracHad      = -1.;
    //I.hadEinHB      = -1.;
    //I.hadEinHE      = -1.;
    //I.hadEinHF      = -1.;
    //I.hadEinHO      = -1.;

    I.ptGenJ      = -10.;
    I.etaGenJ     = -4.;
    I.phiGenJ     = -4.;
    I.massGenJ    = -10.;
    I.ptGen       = -10.;
    I.etaGen      = -4.;
    I.phiGen      = -4.;
    I.massGen     = -10.;
    I.pdgIdGen     = 0.;
//    I.ptLhe       = -10.;
//    I.etaLhe      = -4.;
//    I.phiLhe      = -4.;
//    I.chm         = -1;
//    I.npr         = -1;
//    I.cm          = -1;
//   I.nm          = -1;
    I.partonFlavour     = 0;
    I.hadronFlavour     = 0;
    I.mother      = false;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isTightLepVeto     = false;
    I.isCSVL      = false;
    I.isCSVM      = false;
    I.isCSVT      = false;
    //I.isMatched   = false;
    //I.dR_q1       = 1000.;
    //I.dR_q2       = 1000.;
    //I.dR_q3       = 1000.;
    //I.dR_q4       = 1000.;
    //I.m_q1        = false;
    //I.m_q2        = false;
    //I.m_q3        = false;
    //I.m_q4        = false;
    //I.dR_pi1      = 1000.;
    //I.dR_pi2      = 1000.;
    I.matchBquark = -1;
    I.matchLL     = -1;
    //I.original_jet_index = -1;
    I.isGenMatched = 0;
    I.isVBFGenMatched = 0;
    I.alphaMax = -100.;
    I.sigIP2DMedian = -100.;
    I.theta2DMedian = -100.;
    I.POCA_theta2DMedian = -100.;
    I.nPixelHitsMedian = -1.0;
    I.nHitsMedian = -1.0;
    I.hcalE       = -100.;
    I.ecalE       = -100.;
    I.FracCal     = -100.;
    I.nSV = 0;
    I.nSVCand = 0;
    I.nVertexTracks = 0;
    I.nSelectedTracks = 0;
    I.dRSVJet =  -100.;
    I.flightDist2d  = -100.;
    I.flightDist2dError = -100.;
    I.flightDist3d  = -100.;
    I.flightDist3dError = -100.;
    I.SV_x = -100.;
    I.SV_y =-100.;
    I.SV_z = -100.;
    I.SV_dx = -100.;
    I.SV_dy = -100.;
    I.SV_dz = -100.;
    I.nTracksSV = 0;
    I.SV_mass = -100.;
    I.isCaloTag = 0;
    //I.VBF_DisplacedJet40_VTightID_Hadronic_match = 0;
    //I.VBF_DisplacedJet40_VVTightID_Hadronic_match = 0;
    I.ptJESUp = -1.;
    I.ptJESDown = -1.;
    I.ptJERUp = -1.;
    I.ptJERDown = -1.;
    I.TriggerMatched_VBFJet = 0;
    I.TriggerMatched_DisplacedJet = 0;
    I.TriggerMatched_TripleJet50 = 0;
    I.tau1 = -1.;
    I.tau2 = -1.;
    I.tau3 = -1.;
    I.nSubJets = -1;
    I.tau21 = -1.;
    I.tau31 = -1.;
    I.tau32 = -1.;
    I.tau1_neutral = -1.;
    I.tau2_neutral = -1.;
    I.tau21_neutral = -1.;
    I.tau1_charged = -1.;
    I.tau2_charged = -1.;
    I.tau21_charged = -1.;
    I.nConstituents = -1;
    I.nTrackConstituents = -1;
    I.nTracks0PixelHits = -1;
    I.nTracks1PixelHit = -1;
    I.nTracks2PixelHits = -1;
    I.nTracks3PixelHits = -1;
    I.nTracks4PixelHits = -1;
    I.nTracks5PixelHits = -1;
    I.nTracksAtLeast6PixelHits = -1;
    I.nTracksValidHitInBPix1 = -1;
    I.nTracks0LostInnerHits = -1;
    I.nTracks1LostInnerHit = -1;
    I.nTracksAtLeast2LostInnerHits = -1;
    I.nMatchedGenBquarks = -1;
}

std::string ObjectsFormat::ListJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:energyRaw/F:ptRaw/F:ptUnc/F:CSV/F:CMVA/F:cHadE/F:nHadE/F:eleE/F:photonE/F:muE/F:nEmE/F:cEmE/F:cmuE/F:cHadEFrac/F:nHadEFrac/F:eleEFrac/F:photonEFrac/F:muEFrac/F:nEmEFrac/F:cEmEFrac/F:cmuEFrac/F:cHadMulti/F:nHadMulti/F:eleMulti/F:photonMulti/F:muMulti/F:cMulti/F:nMulti/F:npr/F:cHadMultiFrac/F:nHadMultiFrac/F:eleMultiFrac/F:photonMultiFrac/F:muMultiFrac/F:cMultiFrac/F:nMultiFrac/F:ptGenJ/F:etaGenJ/F:phiGenJ/F:massGenJ/F:ptGen/F:etaGen/F:phiGen/F:massGen/F:pdgIdGen/I:partonFlavour/I:hadronFlavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:isCSVL/O:isCSVM/O:isCSVT/O:matchBquark/I:matchLL/I:isGenMatched/I:isVBFGenMatched/I:alphaMax/F:sigIP2DMedian/F:theta2DMedian/F:POCA_theta2DMedian/F:nPixelHitsMedian/F:nHitsMedian/F:hcalE/F:ecalE/F:FracCal/F:flightDist2d/F:flightDist2dError/F:flightDist3d/F:flightDist3dError/F:nSV/I:nSVCand/I:nVertexTracks/I:nSelectedTracks/I:dRSVJet/F:SV_x/F:SV_y/F:SV_z/F:SV_dx/F:SV_dy/F:SV_dz/F:nTracksSV/I:SV_mass/F:isCaloTag/I:ptJESUp/F:ptJESDown/F:ptJERUp/F:ptJERDown/F:TriggerMatched_VBFJet/I:TriggerMatched_DisplacedJet/I:TriggerMatched_TripleJet50/I:nSubJets/I:tau1/F:tau2/F:tau3/F:tau21/F:tau31/F:tau23/F:tau1_neutral/F:tau2_neutral/F:tau21_neutral/F:tau1_charged/F:tau2_charged/F:tau21_charged/F:nConstituents/I:nTrackConstituents/I:nTracks0PixelHits/I:nTracks1PixelHit/I:nTracks2PixelHits/I:nTracks3PixelHits/I:nTracks4PixelHits/I:nTracks5PixelHits/I:nTracksAtLeast6PixelHits/I:nTracksValidHitInBPix1/I:nTracks0LostInnerHits/I:nTracks1LostInnerHit/I:nTracksAtLeast2LostInnerHits/I:nMatchedGenBquarks/I";}

//"pt/F:eta/F:phi/F:mass/F:energy/F:ptRaw/F:ptUnc/F:dPhi_met/F:dPhi_Jet1/F:puId/F:CSV/F:CSVR/F:CSVRUp/F:CSVRDown/F:CMVA/F:CMVAR/F:CMVARUp/F:CMVARDown/F:QGLikelihood/F:cHadE/F:nHadE/F:cHadEFrac/F:nHadEFrac/F:nEmE/F:nEmEFrac/F:cEmE/F:cEmEFrac/F:cmuE/F:cmuEFrac/F:muE/F:muEFrac/F:eleE/F:eleEFrac/F:eleMulti/F:photonE/F:photonEFrac/F:photonMulti/F:cHadMulti/F:nHadMulti/F:npr/F:cMulti/F:nMulti/F:emEFrac/F:emEinEB/F:emEinEE/F:emEinHF/F:EFracHad/F:hadEinHB/F:hadEinHE/F:hadEinHF/F:hadEinHO/F:ptGenJ/F:etaGenJ/F:phiGenJ/F:massGenJ/F:ptGen/F:etaGen/F:phiGen/F:massGen/F:pdgIdGen/I:partonFlavour/I:hadronFlavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:isCSVL/O:isCSVM/O:isCSVT/O:isMatched/O:dR_q1/F:dR_q2/F:dR_q3/F:dR_q4/F:m_q1/O:m_q2/O:m_q3/O:m_q4/O:dR_pi1/F:dR_pi2/F:matchBquark/I:matchLL/I:original_jet_index/I:isGenMatched/I:alphaMax/F:sigIP2DMedian/F:theta2DMedian/F:hcalE/F:ecalE/F:FracCal/F:flightDist2d/F:flightDist2dSig/F:flightDist3d/F:flightDist3dSig/F:nSV/I:nSVCand/I:nVertexTracks/I:nSelectedTracks/I:dRSVJet/F:SV_x/F:SV_y/F:SV_z/F:SV_dx/F:SV_dy/F:SV_dz/F:nTracksSV/I:SV_mass/F:indexSV/I:isCaloTag/I:VBF_DisplacedJet40_VTightID_Hadronic_match/I:VBF_DisplacedJet40_VVTightID_Hadronic_match/I:ptJESUp/F:ptJESDown/F:ptJERUp/F:ptJERDown/F:nConstituents/I:nTrackConstituents/I:TriggerMatched_VBFJet/I:TriggerMatched_DisplacedJet/I";}
////std::string ObjectsFormat::ListJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:ptRaw/F:ptUnc/F:dPhi_met/F:dPhi_Jet1/F:puId/F:CSV/F:CSVR/F:CSVRUp/F:CSVRDown/F:CMVA/F:CMVAR/F:CMVARUp/F:CMVARDown/F:QGLikelihood/F:chf/F:nhf/F:phf/F:elf/F:muf/F:ptGenJ/F:etaGenJ/F:phiGenJ/F:massGenJ/F:ptGen/F:etaGen/F:phiGen/F:massGen/F:pdgIdGen/I:ptLhe/F:etaLhe/F:phiLhe/I:chm/I:npr/I:cm/I:nm/I:partonFlavour/I:hadronFlavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:isCSVL/O:isCSVM/O:isCSVT/O:isMatched/O:dR_q1/F:dR_q2/F:dR_q3/F:dR_q4/F:m_q1/O:m_q2/O:m_q3/O:m_q4/O:dR_pi1/F:dR_pi2/F:matchBquark/I:matchLL/I:original_jet_index/I:isGenMatched/I:isMatchedToMatchedCHSJet/I";}


//*******************//
//     Fat Jet       //
//*******************//

void ObjectsFormat::FillFatJetType(FatJetType& I, const pat::Jet* R, std::string SoftdropPuppiMassString, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
  //  if(isMC && R->genJet()) {
  //    I.ptGenJ    = R->genJet()->pt();
  //    I.etaGenJ   = R->genJet()->eta();
  //    I.phiGenJ   = R->genJet()->phi();
  //    I.massGenJ  = R->genJet()->mass();
  //  }
  //  if(isMC && R->genParton()) {
  //    I.ptGen     = R->genParton()->pt();
  //    I.etaGen    = R->genParton()->eta();
  //    I.phiGen    = R->genParton()->phi();
  //    I.massGen   = R->genParton()->mass();
  //  }
  //  if(isMC && R->genParton()) {
  //    I.ptLhe     = R->userFloat("ptLhe");
  //    I.etaLhe    = R->userFloat("etaLhe");
  //    I.phiLhe    = R->userFloat("phiLhe");
  //  }
    I.ptRaw       = R->correctedJet(0).pt();
    I.ptUnc       = R->hasUserFloat("JESUncertainty") ? R->userFloat("JESUncertainty") : -1;
    I.dPhi_met    = R->hasUserFloat("dPhi_met") ? R->userFloat("dPhi_met") : -1.;
    I.dPhi_Jet1   = R->hasUserFloat("dPhi_Jet1") ? R->userFloat("dPhi_Jet1") : -1.;
    I.puId        = -1.; //R->userFloat("pileupJetId:fullDiscriminant");
    I.CSV         = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.CSVR        = R->hasUserFloat("ReshapedDiscriminator") ? R->userFloat("ReshapedDiscriminator") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.CSVRUp      = R->hasUserFloat("ReshapedDiscriminatorUp") ? R->userFloat("ReshapedDiscriminatorUp") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.CSVRDown    = R->hasUserFloat("ReshapedDiscriminatorDown") ? R->userFloat("ReshapedDiscriminatorDown") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.pfBoostedDoubleSVAK8 = R->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
  //  I.CSVV1       = R->bDiscriminator("combinedSecondaryVertexV1BJetTags");
  //  I.CSVSL       = R->bDiscriminator("combinedSecondaryVertexSoftPFLeptonV1BJetTags");
  //  I.JPro        = R->bDiscriminator("jetProbabilityBJetTags"); // jetBProbabilityBJetTags
  //  if(isMC) {
  //    if(abs(R->partonFlavour())==5 || abs(R->partonFlavour())==4) {
  //      I.LooseSF      = theBTagWeight->ReturnScaleFactor(1, R->pt(), R->eta());
  //      I.MediumSF     = theBTagWeight->ReturnScaleFactor(2, R->pt(), R->eta());
  //      I.TightSF      = theBTagWeight->ReturnScaleFactor(3, R->pt(), R->eta());
  //      I.LooseSFa     = (1.-0.7*theBTagWeight->ReturnScaleFactor(1, R->pt(), R->eta()))/(1.-0.7);
  //      I.MediumSFa    = (1.-0.6*theBTagWeight->ReturnScaleFactor(2, R->pt(), R->eta()))/(1.-0.6);
  //      I.TightSFa     = (1.-0.5*theBTagWeight->ReturnScaleFactor(3, R->pt(), R->eta()))/(1.-0.5);
  //    }
  //    else {
  //      I.LooseSF      = theBTagWeight->ReturnScaleFactorMistag(1, R->pt(), R->eta());
  //      I.MediumSF     = theBTagWeight->ReturnScaleFactorMistag(2, R->pt(), R->eta());
  //      I.TightSF      = theBTagWeight->ReturnScaleFactorMistag(3, R->pt(), R->eta());
  //      I.LooseSFa     = (1.-1.e-1*theBTagWeight->ReturnScaleFactorMistag(1, R->pt(), R->eta()))/(1.-1.e-1);
  //      I.MediumSFa    = (1.-1.e-2*theBTagWeight->ReturnScaleFactorMistag(2, R->pt(), R->eta()))/(1.-1.e-2);
  //      I.TightSFa     = (1.-1.e-3*theBTagWeight->ReturnScaleFactorMistag(3, R->pt(), R->eta()))/(1.-1.e-3);
  //    }
  //  }
    I.CHSprunedMass            = R->hasUserFloat("ak8PFJetsCHSPrunedMass") ? R->userFloat("ak8PFJetsCHSPrunedMass") : -1.;
    I.CHSsoftdropMass          = R->hasUserFloat("ak8PFJetsCHSSoftDropMass") ? R->userFloat("ak8PFJetsCHSSoftDropMass") : -1.;
    I.softdropPuppiMass        = R->hasUserFloat(SoftdropPuppiMassString) ? R->userFloat(SoftdropPuppiMassString) : -1.;
    I.CHSprunedMassCorr        = R->hasUserFloat("ak8PFJetsCHSPrunedMassCorr") ? R->userFloat("ak8PFJetsCHSPrunedMassCorr") : -1.;
    I.CHSsoftdropMassCorr      = R->hasUserFloat("ak8PFJetsCHSSoftDropMassCorr") ? R->userFloat("ak8PFJetsCHSSoftDropMassCorr") : -1.;
    if(!isMC) I.softdropPuppiMassCorr = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorr") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorr") : -1.;
    if(isMC) I.softdropPuppiMassCorr = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") : -1.;//smeared softdrop puppi mass for MC
    I.softdropPuppiMassCorrNotSmeared = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrNotSmeared") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrNotSmeared") : -1.;
    I.nSoftDropSubJets    = (R->hasSubjets("SoftDrop") && R->subjets("SoftDrop").size() > 0) ? R->subjets("SoftDrop").size() : -1;//calculated from SoftDrop sub-jets
    I.pt1         = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->pt() : -1.) : -1.;
    I.eta1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->eta() : -9.) : -9.;
    I.phi1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->phi() : -9.) : -9.;
    I.mass1       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->mass() : -1.) : -1.;
    I.CSV1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    //I.CSVR1       = R->hasUserFloat("ReshapedDiscriminator1") ? R->userFloat("ReshapedDiscriminator1") : -99.;
    //I.CSVR1Up     = R->hasUserFloat("ReshapedDiscriminatorUp1") ? R->userFloat("ReshapedDiscriminatorUp1") : -99.;
    //I.CSVR1Down   = R->hasUserFloat("ReshapedDiscriminatorDown1") ? R->userFloat("ReshapedDiscriminatorDown1") : -99.;
    I.CMVA1       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    //I.CMVAR1      = R->hasUserFloat("CMVAR1") ? R->userFloat("CMVAR1") : -99.;
    //I.CMVAR1Up    = R->hasUserFloat("CMVAR1Up") ? R->userFloat("CMVAR1Up") : -99.;
    //I.CMVAR1Down  = R->hasUserFloat("CMVAR1Down") ? R->userFloat("CMVAR1Down") : -99.;
    I.flavour1    = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->hadronFlavour() : -1.) : -1.;
    if (R->hasSubjets("SoftDrop")){
      if( R->subjets("SoftDrop").size() > 0){
	if (R->subjets("SoftDrop")[0]->hasTagInfo("pfSecondaryVertex")){
	  I.nSV1 = R->subjets("SoftDrop")[0]->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertices();
	  I.nVertexTracks1 = R->subjets("SoftDrop")[0]->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexTracks();
	}
      }
    }
    I.pt2         = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->pt() : -1.) : -1.;
    I.eta2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->eta() : -9.) : -9.;
    I.phi2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->phi() : -9.) : -9.;
    I.mass2       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->mass() : -1.) : -1.;
    I.CSV2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    //I.CSVR2       = R->hasUserFloat("ReshapedDiscriminator2") ? R->userFloat("ReshapedDiscriminator2") : -99.;
    //I.CSVR2Up     = R->hasUserFloat("ReshapedDiscriminatorUp2") ? R->userFloat("ReshapedDiscriminatorUp2") : -99.;
    //I.CSVR2Down   = R->hasUserFloat("ReshapedDiscriminatorDown2") ? R->userFloat("ReshapedDiscriminatorDown2") : -99.;
    I.CMVA2       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    //I.CMVAR2      = R->hasUserFloat("CMVAR2") ? R->userFloat("CMVAR2") : -99.;
    //I.CMVAR2Up    = R->hasUserFloat("CMVAR2Up") ? R->userFloat("CMVAR2Up") : -99.;
    //I.CMVAR2Down  = R->hasUserFloat("CMVAR2Down") ? R->userFloat("CMVAR2Down") : -99.;
    I.flavour2    = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->hadronFlavour() : -1.) : -1.;
    if (R->hasSubjets("SoftDrop")){
      if( R->subjets("SoftDrop").size() > 1){
	if (R->subjets("SoftDrop")[1]->hasTagInfo("pfSecondaryVertex")){
	  I.nSV2 = R->subjets("SoftDrop")[1]->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertices();
	  I.nVertexTracks2 = R->subjets("SoftDrop")[1]->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nVertexTracks();
	}
      }
    }
    I.dR          = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? deltaR(*R->subjets("SoftDrop")[0], *R->subjets("SoftDrop")[1]) : -1.) : -1.;
    I.chsTau21    = R->hasUserFloat("NjettinessAK8:tau1") ? (R->userFloat("NjettinessAK8:tau1") != 0 ? R->userFloat("NjettinessAK8:tau2")/R->userFloat("NjettinessAK8:tau1") : -1.) : -1.;
    I.puppiTau21  = R->hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ? (R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") != 0 ? R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2")/R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") : -1.) : -1.;
    I.ddtTau21    = R->hasUserFloat("ddtTau21") ? R->userFloat("ddtTau21") : -1.;
    I.BDSV        = R->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
    I.cHadE       = R->chargedHadronEnergy();
    I.nHadE       = R->neutralHadronEnergy();
    I.eleE        = R->electronEnergy();
    I.photonE     = R->photonEnergy();
    I.muE         = R->muonEnergy();
    I.nEmE        = R->neutralEmEnergy();
    I.cEmE        = R->chargedEmEnergy();
    I.cmuE        = R->chargedMuEnergy();
    I.cHadEFrac   = R->hasUserFloat("cHadEFrac") ? R->userFloat("cHadEFrac") : -1.;//R->chargedHadronEnergyFraction();
    I.nHadEFrac   = R->hasUserFloat("nHadEFrac") ? R->userFloat("nHadEFrac") : -1.;//R->neutralHadronEnergyFraction();
    I.eleEFrac    = R->hasUserFloat("eleEFrac") ? R->userFloat("eleEFrac") : -1.;//R->electronEnergyFraction();
    I.photonEFrac = R->hasUserFloat("photonEFrac") ? R->userFloat("photonEFrac") : -1.;//R->photonEnergyFraction();
    I.muEFrac     = R->hasUserFloat("muEFrac") ? R->userFloat("muEFrac") : -1.;//R->muonEnergyFraction();
    I.nEmEFrac    = R->hasUserFloat("nEmEFrac") ? R->userFloat("nEmEFrac") : -1.;//R->neutralEmEnergyFraction();
    I.cEmEFrac    = R->hasUserFloat("cEmEFrac") ? R->userFloat("cEmEFrac") : -1.;//R->chargedEmEnergyFraction();
    I.cmuEFrac    = R->hasUserFloat("cmuEFrac") ? R->userFloat("cmuEFrac") : -1.;//R->chargedMuEnergyFraction();
    I.cHadMulti   = R->chargedHadronMultiplicity();
    I.nHadMulti   = R->neutralHadronMultiplicity();
    I.eleMulti    = R->electronMultiplicity();
    I.photonMulti = R->photonMultiplicity();
    I.muMulti     = R->muonMultiplicity();
    I.cMulti      = R->chargedMultiplicity();
    I.nMulti      = R->neutralMultiplicity();
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.cHadMultiFrac   = float(R->chargedHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nHadMultiFrac   = float(R->neutralHadronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.eleMultiFrac    = float(R->electronMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.photonMultiFrac = float(R->photonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.muMultiFrac     = float(R->muonMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.cMultiFrac      = float(R->chargedMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.nMultiFrac      = float(R->neutralMultiplicity())/float(R->chargedMultiplicity() + R->neutralMultiplicity());
    I.partonFlavour     = R->partonFlavour();
    I.hadronFlavour     = R->hadronFlavour();
    if(isMC && R->genParton()) I.mother = Utilities::FindMotherId(R->genParton());//verify
    //if(isMC && R->genParton()) I.mother = false;//Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genParton()));
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isTightLepVeto     = R->hasUserInt("isTightLepVeto") ? R->userInt("isTightLepVeto") : false;
    I.isMatched   = (I.mother==25);
    I.JESUnc      = R->hasUserFloat("JESUncertainty") ? R->userFloat("JESUncertainty") : -1.;
    I.ptJERUp     = R->hasUserFloat("ptJERUp") ? R->userFloat("ptJERUp") : -1.;
    I.etaJERUp    = R->hasUserFloat("etaJERUp") ? R->userFloat("etaJERUp") : -1.;
    I.phiJERUp    = R->hasUserFloat("phiJERUp") ? R->userFloat("phiJERUp") : -9.;
    I.energyJERUp = R->hasUserFloat("energyJERUp") ? R->userFloat("energyJERUp") : -1.;
    I.ptJERDown   = R->hasUserFloat("ptJERDown") ? R->userFloat("ptJERDown") : -1.;
    I.etaJERDown  = R->hasUserFloat("etaJERDown") ? R->userFloat("etaJERDown") : -1.;
    I.phiJERDown  = R->hasUserFloat("phiJERDown") ? R->userFloat("phiJERDown") : -9.;
    I.energyJERDown = R->hasUserFloat("energyJERDown") ? R->userFloat("energyJERDown") : -1.;
    I.smearFact   = R->hasUserFloat("smearFactor") ? R->userFloat("smearFactor") : -1.;
    I.smearFactUp   = R->hasUserFloat("smearFactorUp") ? R->userFloat("smearFactorUp") : -1.;
    I.smearFactDown = R->hasUserFloat("smearFactorDown") ? R->userFloat("smearFactorDown") : -1.;
    I.softdropPuppiMassCorrJMS = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMS") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMS") : -1.;
    I.softdropPuppiMassCorrJMSUp = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSUp") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMSUp") : -1.;
    I.softdropPuppiMassCorrJMSDown = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSDown") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMSDown") : -1.;
    I.softdropPuppiMassCorrJMR = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") : -1.;
    I.softdropPuppiMassCorrJMRUp = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRUp") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMRUp") : -1.;
    I.softdropPuppiMassCorrJMRDown = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRDown") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMRDown") : -1.;
    I.dR_q1       = R->hasUserFloat("dR_q1") ? R->userFloat("dR_q1") : 1000;
    I.dR_q2       = R->hasUserFloat("dR_q2") ? R->userFloat("dR_q2") : 1000;
    I.dR_q3       = R->hasUserFloat("dR_q3") ? R->userFloat("dR_q3") : 1000;
    I.dR_q4       = R->hasUserFloat("dR_q4") ? R->userFloat("dR_q4") : 1000;
    I.m_q1        = R->hasUserFloat("dR_q1") ? (R->userFloat("dR_q1")<0.8 ? true : false) : false;
    I.m_q2        = R->hasUserFloat("dR_q2") ? (R->userFloat("dR_q2")<0.8 ? true : false) : false;
    I.m_q3        = R->hasUserFloat("dR_q3") ? (R->userFloat("dR_q3")<0.8 ? true : false) : false;
    I.m_q4        = R->hasUserFloat("dR_q4") ? (R->userFloat("dR_q4")<0.8 ? true : false) : false;
    I.dR_pi1      = R->hasUserFloat("dR_pi1") ? R->userFloat("dR_pi1") : 1000;
    I.dR_pi2      = R->hasUserFloat("dR_pi2") ? R->userFloat("dR_pi2") : 1000;
    I.matchBquark = R->hasUserInt("hasMatchedBquarks") ? R->userInt("hasMatchedBquarks") : -1;
    I.matchLL     = R->hasUserInt("hasMatchedLL") ? R->userInt("hasMatchedLL") : -1;
    I.isGenMatched = R->hasUserInt("isGenMatched") ? R->userInt("isGenMatched") : 0;
    I.alphaMax    = R->hasUserFloat("alphaMax") ? R->userFloat("alphaMax") : -1000;
    I.sigIP2DMedian = R->hasUserFloat("sigIP2DMedian") ? R->userFloat("sigIP2DMedian") : -1000;
    I.theta2DMedian = R->hasUserFloat("theta2DMedian") ? R->userFloat("theta2DMedian") : -1000;
    I.POCA_theta2DMedian = R->hasUserFloat("POCA_theta2DMedian") ? R->userFloat("POCA_theta2DMedian") : -1000;
    I.nPixelHitsMedian = R->hasUserFloat("nPixelHitsMedian") ? R->userFloat("nPixelHitsMedian") : -1.0;
    I.nHitsMedian = R->hasUserFloat("nHitsMedian") ? R->userFloat("nHitsMedian") : -1.0;
    I.nConstituents = R->hasUserInt("nConstituents") ? R->userInt("nConstituents") : -1;
    I.nTrackConstituents = R->hasUserInt("nTrackConstituents") ? R->userInt("nTrackConstituents") : -1;
    I.nTracks0PixelHits = R->hasUserInt("nTracks0PixelHits") ? R->userInt("nTracks0PixelHits") : -1;
    I.nTracks1PixelHit  = R->hasUserInt("nTracks1PixelHit")  ? R->userInt("nTracks1PixelHit")  : -1;
    I.nTracks2PixelHits = R->hasUserInt("nTracks2PixelHits") ? R->userInt("nTracks2PixelHits") : -1;
    I.nTracks3PixelHits = R->hasUserInt("nTracks3PixelHits") ? R->userInt("nTracks3PixelHits") : -1;
    I.nTracks4PixelHits = R->hasUserInt("nTracks4PixelHits") ? R->userInt("nTracks4PixelHits") : -1;
    I.nTracks5PixelHits = R->hasUserInt("nTracks5PixelHits") ? R->userInt("nTracks5PixelHits") : -1;
    I.nTracksAtLeast6PixelHits = R->hasUserInt("nTracksAtLeast6PixelHits") ? R->userInt("nTracksAtLeast6PixelHits") : -1;
    I.nTracksValidHitInBPix1 = R->hasUserInt("nTracksValidHitInBPix1") ? R->userInt("nTracksValidHitInBPix1") : -1;
    I.nTracks0LostInnerHits = R->hasUserInt("nTracks0LostInnerHits") ? R->userInt("nTracks0LostInnerHits") : -1;
    I.nTracks1LostInnerHit = R->hasUserInt("nTracks1LostInnerHit") ? R->userInt("nTracks1LostInnerHit") : -1;
    I.nTracksAtLeast2LostInnerHits = R->hasUserInt("nTracksAtLeast2LostInnerHits") ? R->userInt("nTracksAtLeast2LostInnerHits") : -1;
    I.nMatchedGenBquarks = R->hasUserInt("nMatchedGenBquarks") ? R->userInt("nMatchedGenBquarks") : -1;
}

void ObjectsFormat::ResetFatJetType(FatJetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.ptRaw       = -1.;
    I.ptUnc       = -1.;
    I.dPhi_met    = -1.;
    I.dPhi_Jet1   = -1.;
    I.puId        = -1.;
    I.CSV         = -99.;
    I.CSVR        = -99.;
    I.CSVRUp      = -99.;
    I.CSVRDown    = -99.;
    I.pfBoostedDoubleSVAK8 = -99.;
    I.CHSprunedMass            = -1.;
    I.CHSsoftdropMass          = -1.;
    I.softdropPuppiMass     = -1.;
    I.CHSprunedMassCorr        = -1.;
    I.CHSsoftdropMassCorr      = -1.;
    I.softdropPuppiMassCorr = -1.;
    I.softdropPuppiMassCorrNotSmeared = -1.;
    I.nSoftDropSubJets = -1;
    I.pt1         = -1.;
    I.eta1        = -9.;
    I.phi1        = -9.;
    I.mass1       = -1.;
    I.CSV1        = -99.;
    //I.CSVR1       = -99.;
    //I.CSVR1Up     = -99.;
    //I.CSVR1Down   = -99.;
    I.CMVA1       = -99.;
    //I.CMVAR1      = -99.;
    //I.CMVAR1Up    = -99.;
    //I.CMVAR1Down  = -99.;
    I.flavour1    = -1.;
    I.nSV1  = -1;
    I.nVertexTracks1  = -1;
    I.pt2         = -1.;
    I.eta2        = -9.;
    I.phi2        = -9.;
    I.mass2       = -1.;
    I.CSV2        = -99.;
    //I.CSVR2       = -99.;
    //I.CSVR2Up     = -99.;
    //I.CSVR2Down   = -99.;
    I.CMVA2       = -99.;
    //I.CMVAR2      = -99.;
    //I.CMVAR2Up    = -99.;
    //I.CMVAR2Down  = -99.;
    I.flavour2    = -1.;
    I.nSV2  = -1;
    I.nVertexTracks2  = -1;
    I.dR          = -1.;
    I.chsTau21    = -1.;
    I.puppiTau21  = -1.;
    I.ddtTau21    = -1.;
    I.BDSV        = -99.;
    I.cHadE       = -1.;
    I.nHadE       = -1.;
    I.eleE        = -1.;
    I.photonE     = -1.;
    I.muE         = -1.;
    I.nEmE        = -1.;
    I.cEmE        = -1.;
    I.cmuE        = -1.;
    I.cHadEFrac   = -1.;
    I.nHadEFrac   = -1.;
    I.eleEFrac    = -1.;
    I.photonEFrac = -1.;
    I.muEFrac     = -1.;
    I.nEmEFrac    = -1.;
    I.cEmEFrac    = -1.;
    I.cmuEFrac    = -1.;
    I.cHadMulti   = -1.;
    I.nHadMulti   = -1.;
    I.eleMulti    = -1.;
    I.photonMulti = -1.;
    I.muMulti     = -1.;
    I.cMulti      = -1.;
    I.nMulti      = -1.;
    I.npr         = -1.;
    I.cHadMultiFrac   = -1.;
    I.nHadMultiFrac   = -1.;
    I.eleMultiFrac    = -1.;
    I.photonMultiFrac = -1.;
    I.muMultiFrac     = -1.;
    I.cMultiFrac      = -1.;
    I.nMultiFrac      = -1.;
    I.partonFlavour     = 0;
    I.hadronFlavour     = 0;
    I.mother      = false;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isTightLepVeto     = false;
    I.isMatched   = false;
    I.JESUnc      = -1.;
    I.ptJERUp     = -1.;
    I.etaJERUp    = -1.;
    I.phiJERUp    = -9.;
    I.energyJERUp = -1.;
    I.ptJERDown   = -1.;
    I.etaJERDown  = -1.;
    I.phiJERDown  = -9.;
    I.energyJERDown = -1.;
    I.smearFact   = -1.;
    I.smearFactUp   = -1.;
    I.smearFactDown = -1.;
    I.softdropPuppiMassCorrJMS = -1.;
    I.softdropPuppiMassCorrJMSUp = -1.;
    I.softdropPuppiMassCorrJMSDown = -1.;
    I.softdropPuppiMassCorrJMR = -1.;
    I.softdropPuppiMassCorrJMRUp = -1.;
    I.softdropPuppiMassCorrJMRDown = -1.;
    I.dR_q1       = 1000.;
    I.dR_q2       = 1000.;
    I.dR_q3       = 1000.;
    I.dR_q4       = 1000.;
    I.m_q1        = false;
    I.m_q2        = false;
    I.m_q3        = false;
    I.m_q4        = false;
    I.dR_pi1      = 1000.;
    I.dR_pi2      = 1000.;
    I.matchBquark = -1;
    I.matchLL     = -1;
    I.isGenMatched = false;
    I.alphaMax         = -100.;
    I.sigIP2DMedian    = -100.;
    I.theta2DMedian    = -100.;
    I.POCA_theta2DMedian = -100.;
    I.nPixelHitsMedian = -1.;
    I.nHitsMedian      = -1.;
    I.nConstituents = -1;
    I.nTrackConstituents = -1;
    I.nTracks0PixelHits = -1;
    I.nTracks1PixelHit  = -1;
    I.nTracks2PixelHits = -1;
    I.nTracks3PixelHits = -1;
    I.nTracks4PixelHits = -1;
    I.nTracks5PixelHits = -1;
    I.nTracksAtLeast6PixelHits = -1;
    I.nTracksValidHitInBPix1 = -1;
    I.nTracks0LostInnerHits = -1;
    I.nTracks1LostInnerHit = -1;
    I.nTracksAtLeast2LostInnerHits = -1;
    I.nMatchedGenBquarks = -1;
}

//to be removed: CSVR1, etc
std::string ObjectsFormat::ListFatJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:ptRaw/F:ptUnc/F:dPhi_met/F:dPhi_Jet1/F:puId/F:CSV/F:CSVR/F:CSVRUp/F:CSVRDown/F:pfBoostedDoubleSVAK8/F:CHSprunedMass/F:CHSsoftdropMass/F:softdropPuppiMass/F:CHSprunedMassCorr/F:CHSsoftdropMassCorr/F:softdropPuppiMassCorr/F:softdropPuppiMassCorrNotSmeared/F:nSoftDropSubJets/I:pt1/F:eta1/F:phi1/F:mass1/F:CSV1/F:CSVR1/F:CSVR1Up/F:CSVR1Down/F:CMVA1/F:CMVAR1/F:CMVAR1Up/F:CMVAR1Down/F:flavour1/F:nSV1/I:nVertexTracks1/I:pt2/F:eta2/F:phi2/F:mass2/F:CSV2/F:CSVR2/F:CSVR2Up/F:CSVR2Down/F:CMVA2/F:CMVAR2/F:CMVAR2Up/F:CMVAR2Down/F:flavour2/F:nSV2/I:nVertexTracks2/I:dR/F:chsTau21/F:puppiTau21/F:ddtTau21/F:BDSV/F:cHadE/F:nHadE/F:eleE/F:photonE/F:muE/F:nEmE/F:cEmE/F:cmuE/F:cHadEFrac/F:nHadEFrac/F:eleEFrac/F:photonEFrac/F:muEFrac/F:nEmEFrac/F:cEmEFrac/F:cmuEFrac/F:cHadMulti/F:nHadMulti/F:eleMulti/F:photonMulti/F:muMulti/F:cMulti/F:nMulti/F:npr/F:cHadMultiFrac/F:nHadMultiFrac/F:eleMultiFrac/F:photonMultiFrac/F:muMultiFrac/F:cMultiFrac/F:nMultiFrac/F:partonFlavour/I:hadronFlavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:isCSVL/O:isCSVM/O:isCSVT/O:isMatched/O:JESUnc/F:ptJERUp/F:etaJERUp/F:phiJERUp/F:energyJERUp/F:ptJERDown/F:etaJERDown/F:phiJERDown/F:energyJERDown/F:smearFact/F:smearFactUp/F:smearFactDown/F:softdropPuppiMassCorrJMS/F:softdropPuppiMassCorrJMSUp/F:softdropPuppiMassCorrJMSDown/F:softdropPuppiMassCorrJMR/F:softdropPuppiMassCorrJMRUp/F:softdropPuppiMassCorrJMRDown/F:dR_q1/F:dR_q2/F:dR_q3/F:dR_q4/F:m_q1/O:m_q2/O:m_q3/O:m_q4/O:dR_pi1/F:dR_pi2/F:matchBquark/I:matchLL/I:isGenMatched/O:alphaMax/F:sigIP2DMedian/F:theta2DMedian/F:POCA_theta2DMedian/F:nPixelHitsMedian/F:nHitsMedian/F:nConstituents/I:nTrackConstituents/I:nTracks0PixelHits/I:nTracks1PixelHit/I:nTracks2PixelHits/I:nTracks3PixelHits/I:nTracks4PixelHits/I:nTracks5PixelHits/I:nTracksAtLeast6PixelHits/I:nTracksValidHitInBPix1/I:nTracks0LostInnerHits/I:nTracks1LostInnerHit/I:nTracksAtLeast2LostInnerHits/I:nMatchedGenBquarks/I";}


//*******************//
// Custom Fat Jet    //
//*******************//

void ObjectsFormat::FillCustomFatJetType(CustomFatJetType& I, const pat::Jet* R, std::string SoftdropPuppiMassString, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.ptRaw       = R->correctedJet(0).pt();
    I.ptUnc       = R->hasUserFloat("JESUncertainty") ? R->userFloat("JESUncertainty") : -1;
    I.dPhi_met    = R->hasUserFloat("dPhi_met") ? R->userFloat("dPhi_met") : -1.;
    I.dPhi_Jet1   = R->hasUserFloat("dPhi_Jet1") ? R->userFloat("dPhi_Jet1") : -1.;
    I.puId        = -1.; //R->userFloat("pileupJetId:fullDiscriminant");
    I.CSV         = R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVR        = R->hasUserFloat("ReshapedDiscriminator") ? R->userFloat("ReshapedDiscriminator") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVRUp      = R->hasUserFloat("ReshapedDiscriminatorUp") ? R->userFloat("ReshapedDiscriminatorUp") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //I.CSVRDown    = R->hasUserFloat("ReshapedDiscriminatorDown") ? R->userFloat("ReshapedDiscriminatorDown") : R->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    I.CHSprunedMass            = R->hasUserFloat("ak8PFJetsCHSPrunedMass") ? R->userFloat("ak8PFJetsCHSPrunedMass") : -1.;
    I.CHSsoftdropMass          = R->hasUserFloat("ak8PFJetsCHSSoftDropMass") ? R->userFloat("ak8PFJetsCHSSoftDropMass") : -1.;
    I.prunedMass            = R->hasUserFloat("ak8PFJetsPrunedMass") ? R->userFloat("ak8PFJetsPrunedMass") : -1.;
    I.softdropMass          = R->hasUserFloat("ak8PFJetsSoftDropMass") ? R->userFloat("ak8PFJetsSoftDropMass") : -1.;
    I.softdropPuppiMass        = R->hasUserFloat(SoftdropPuppiMassString) ? R->userFloat(SoftdropPuppiMassString) : -1.;
    I.CHSprunedMassCorr        = R->hasUserFloat("ak8PFJetsCHSPrunedMassCorr") ? R->userFloat("ak8PFJetsCHSPrunedMassCorr") : -1.;
    I.CHSsoftdropMassCorr      = R->hasUserFloat("ak8PFJetsCHSSoftDropMassCorr") ? R->userFloat("ak8PFJetsCHSSoftDropMassCorr") : -1.;
    I.prunedMassCorr        = R->hasUserFloat("ak8PFJetsPrunedMassCorr") ? R->userFloat("ak8PFJetsPrunedMassCorr") : -1.;
    I.softdropMassCorr      = R->hasUserFloat("ak8PFJetsSoftDropMassCorr") ? R->userFloat("ak8PFJetsSoftDropMassCorr") : -1.;
    if(!isMC) I.softdropPuppiMassCorr = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorr") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorr") : -1.;
    if(isMC) I.softdropPuppiMassCorr = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") : -1.;//smeared softdrop puppi mass for MC
    I.softdropPuppiMassCorrNotSmeared = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrNotSmeared") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrNotSmeared") : -1.;
    //////////////////////
    //subjets SoftDrop
    //////////////////////
    I.pt1         = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->pt() : -1.) : -1.;
    I.eta1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->eta() : -9.) : -9.;
    I.phi1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->phi() : -9.) : -9.;
    I.mass1       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->mass() : -1.) : -1.;
    I.CSV1        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    //I.CSVR1       = R->hasUserFloat("ReshapedDiscriminator1") ? R->userFloat("ReshapedDiscriminator1") : -99.;
    //I.CSVR1Up     = R->hasUserFloat("ReshapedDiscriminatorUp1") ? R->userFloat("ReshapedDiscriminatorUp1") : -99.;
    //I.CSVR1Down   = R->hasUserFloat("ReshapedDiscriminatorDown1") ? R->userFloat("ReshapedDiscriminatorDown1") : -99.;
    I.CMVA1       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    //I.CMVAR1      = R->hasUserFloat("CMVAR1") ? R->userFloat("CMVAR1") : -99.;
    //I.CMVAR1Up    = R->hasUserFloat("CMVAR1Up") ? R->userFloat("CMVAR1Up") : -99.;
    //I.CMVAR1Down  = R->hasUserFloat("CMVAR1Down") ? R->userFloat("CMVAR1Down") : -99.;
    I.flavour1    = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->hadronFlavour() : -1.) : -1.;
    I.pt2         = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->pt() : -1.) : -1.;
    I.eta2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->eta() : -9.) : -9.;
    I.phi2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->phi() : -9.) : -9.;
    I.mass2       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->mass() : -1.) : -1.;
    I.CSV2        = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    //I.CSVR2       = R->hasUserFloat("ReshapedDiscriminator2") ? R->userFloat("ReshapedDiscriminator2") : -99.;
    //I.CSVR2Up     = R->hasUserFloat("ReshapedDiscriminatorUp2") ? R->userFloat("ReshapedDiscriminatorUp2") : -99.;
    //I.CSVR2Down   = R->hasUserFloat("ReshapedDiscriminatorDown2") ? R->userFloat("ReshapedDiscriminatorDown2") : -99.;
    I.CMVA2       = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    //I.CMVAR2      = R->hasUserFloat("CMVAR2") ? R->userFloat("CMVAR2") : -99.;
    //I.CMVAR2Up    = R->hasUserFloat("CMVAR2Up") ? R->userFloat("CMVAR2Up") : -99.;
    //I.CMVAR2Down  = R->hasUserFloat("CMVAR2Down") ? R->userFloat("CMVAR2Down") : -99.;
    I.flavour2    = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->hadronFlavour() : -1.) : -1.;
    //////////////////////
    //subjets SoftDropPuppi
    //////////////////////
    I.pt1SDP      = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->pt() : -1.) : -1.;
    I.eta1SDP     = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->eta() : -9.) : -9.;
    I.phi1SDP     = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->phi() : -9.) : -9.;
    I.mass1SDP    = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->mass() : -1.) : -1.;
    I.CSV1SDP      = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    I.CMVA1SDP    = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    I.flavour1SDP = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 0 ? R->subjets("SoftDrop")[0]->hadronFlavour() : -1.) : -1.;
    I.pt2SDP      = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->pt() : -1.) : -1.;
    I.eta2SDP     = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->eta() : -9.) : -9.;
    I.phi2SDP     = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->phi() : -9.) : -9.;
    I.mass2SDP    = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->mass() : -1.) : -1.;
    I.CSV2SDP     = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99.) : -99.;
    I.CMVA2SDP    = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->bDiscriminator("pfCombinedMVAV2BJetTags") : -99.) : -99.;
    I.flavour2SDP = R->hasSubjets("SoftDropPuppi") ? (R->subjets("SoftDrop").size() > 1 ? R->subjets("SoftDrop")[1]->hadronFlavour() : -1.) : -1.;
    I.dR          = R->hasSubjets("SoftDrop") ? (R->subjets("SoftDrop").size() > 1 ? deltaR(*R->subjets("SoftDrop")[0], *R->subjets("SoftDrop")[1]) : -1.) : -1.;
    I.Tau21       = R->hasUserFloat("NjettinessAK8:tau1") ? (R->userFloat("NjettinessAK8:tau1") != 0 ? R->userFloat("NjettinessAK8:tau2")/R->userFloat("NjettinessAK8:tau1") : -1.) : -1.;
    I.puppiTau21  = R->hasUserFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") ? (R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") != 0 ? R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2")/R->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1") : -1.) : -1.;
    I.ddtTau21    = R->hasUserFloat("ddtTau21") ? R->userFloat("ddtTau21") : -1.;
    I.BDSV        = R->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
    I.chf         = R->hasUserFloat("chf") ? R->userFloat("chf") : -1.;//R->chargedHadronEnergyFraction();
    I.nhf         = R->hasUserFloat("nhf") ? R->userFloat("nhf") : -1.;//R->neutralHadronEnergyFraction();
    I.phf         = R->hasUserFloat("phf") ? R->userFloat("phf") : -1.;//R->neutralEmEnergyFraction();
    I.elf         = R->hasUserFloat("elf") ? R->userFloat("elf") : -1.;//R->chargedEmEnergyFraction();
    I.muf         = R->hasUserFloat("muf") ? R->userFloat("muf") : -1.;//R->muonEnergyFraction();
    I.chm         = R->chargedHadronMultiplicity();
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.partonFlavour     = R->partonFlavour();
    I.hadronFlavour     = R->hadronFlavour();
    if(isMC && R->genParton()) I.mother = false;//Utilities::FindMotherId(dynamic_cast<const reco::Candidate*>(R->genParton()));
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isTightLepVeto     = R->hasUserInt("isTightLepVeto") ? R->userInt("isTightLepVeto") : false;
    //I.isMatched   = (I.mother==25);
    //I.JESUnc      = R->hasUserFloat("JESUncertainty") ? R->userFloat("JESUncertainty") : -1.;
    //I.ptJERUp     = R->hasUserFloat("ptJERUp") ? R->userFloat("ptJERUp") : -1.;
    //I.etaJERUp    = R->hasUserFloat("etaJERUp") ? R->userFloat("etaJERUp") : -1.;
    //I.phiJERUp    = R->hasUserFloat("phiJERUp") ? R->userFloat("phiJERUp") : -9.;
    //I.energyJERUp = R->hasUserFloat("energyJERUp") ? R->userFloat("energyJERUp") : -1.;
    //I.ptJERDown   = R->hasUserFloat("ptJERDown") ? R->userFloat("ptJERDown") : -1.;
    //I.etaJERDown  = R->hasUserFloat("etaJERDown") ? R->userFloat("etaJERDown") : -1.;
    //I.phiJERDown  = R->hasUserFloat("phiJERDown") ? R->userFloat("phiJERDown") : -9.;
    //I.energyJERDown = R->hasUserFloat("energyJERDown") ? R->userFloat("energyJERDown") : -1.;
    //I.smearFact   = R->hasUserFloat("smearFactor") ? R->userFloat("smearFactor") : -1.;
    //I.smearFactUp   = R->hasUserFloat("smearFactorUp") ? R->userFloat("smearFactorUp") : -1.;
    //I.smearFactDown = R->hasUserFloat("smearFactorDown") ? R->userFloat("smearFactorDown") : -1.;
    //I.softdropPuppiMassCorrJMS = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMS") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMS") : -1.;
    //I.softdropPuppiMassCorrJMSUp = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSUp") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMSUp") : -1.;
    //I.softdropPuppiMassCorrJMSDown = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSDown") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMSDown") : -1.;
    //I.softdropPuppiMassCorrJMR = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMR") : -1.;
    //I.softdropPuppiMassCorrJMRUp = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRUp") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMRUp") : -1.;
    //I.softdropPuppiMassCorrJMRDown = R->hasUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRDown") ? R->userFloat("ak8PFJetsPuppiSoftDropMassCorrJMRDown") : -1.;
    I.dR_q1       = R->hasUserFloat("dR_q1") ? R->userFloat("dR_q1") : 1000;
    I.dR_q2       = R->hasUserFloat("dR_q2") ? R->userFloat("dR_q2") : 1000;
    I.dR_q3       = R->hasUserFloat("dR_q3") ? R->userFloat("dR_q3") : 1000;
    I.dR_q4       = R->hasUserFloat("dR_q4") ? R->userFloat("dR_q4") : 1000;
    I.m_q1        = R->hasUserFloat("dR_q1") ? (R->userFloat("dR_q1")<0.8 ? true : false) : false;
    I.m_q2        = R->hasUserFloat("dR_q2") ? (R->userFloat("dR_q2")<0.8 ? true : false) : false;
    I.m_q3        = R->hasUserFloat("dR_q3") ? (R->userFloat("dR_q3")<0.8 ? true : false) : false;
    I.m_q4        = R->hasUserFloat("dR_q4") ? (R->userFloat("dR_q4")<0.8 ? true : false) : false;
    I.dR_pi1      = R->hasUserFloat("dR_pi1") ? R->userFloat("dR_pi1") : 1000;
    I.dR_pi2      = R->hasUserFloat("dR_pi2") ? R->userFloat("dR_pi2") : 1000;
    I.matchBquark = R->hasUserInt("hasMatchedBquarks") ? R->userInt("hasMatchedBquarks") : -1;
    I.matchLL     = R->hasUserInt("hasMatchedLL") ? R->userInt("hasMatchedLL") : -1;
    I.dR_q1_sj1   = R->hasUserFloat("dR_q1_sj1") ? R->userFloat("dR_q1_sj1") : 1000;
    I.dR_q2_sj1   = R->hasUserFloat("dR_q2_sj1") ? R->userFloat("dR_q2_sj1") : 1000;
    I.dR_q3_sj1   = R->hasUserFloat("dR_q3_sj1") ? R->userFloat("dR_q3_sj1") : 1000;
    I.dR_q4_sj1   = R->hasUserFloat("dR_q4_sj1") ? R->userFloat("dR_q4_sj1") : 1000;
    I.dR_pi1_sj1  = R->hasUserFloat("dR_pi1_sj1") ? R->userFloat("dR_pi1_sj1") : 1000;
    I.dR_pi2_sj1  = R->hasUserFloat("dR_pi2_sj1") ? R->userFloat("dR_pi2_sj1") : 1000;
    I.dR_q1_sj2   = R->hasUserFloat("dR_q1_sj2") ? R->userFloat("dR_q1_sj2") : 1000;
    I.dR_q2_sj2   = R->hasUserFloat("dR_q2_sj2") ? R->userFloat("dR_q2_sj2") : 1000;
    I.dR_q3_sj2   = R->hasUserFloat("dR_q3_sj2") ? R->userFloat("dR_q3_sj2") : 1000;
    I.dR_q4_sj2   = R->hasUserFloat("dR_q4_sj2") ? R->userFloat("dR_q4_sj2") : 1000;
    I.dR_pi1_sj2  = R->hasUserFloat("dR_pi1_sj2") ? R->userFloat("dR_pi1_sj2") : 1000;
    I.dR_pi2_sj2  = R->hasUserFloat("dR_pi2_sj2") ? R->userFloat("dR_pi2_sj2") : 1000;
}

void ObjectsFormat::ResetCustomFatJetType(CustomFatJetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.ptRaw       = -1.;
    I.ptUnc       = -1.;
    I.dPhi_met    = -1.;
    I.dPhi_Jet1   = -1.;
    I.puId        = -1.;
    I.CSV         = -99.;
    //I.CSVR        = -99.;
    //I.CSVRUp      = -99.;
    //I.CSVRDown    = -99.;
    I.CHSprunedMass  = -1.;
    I.CHSsoftdropMass = -1.;
    I.prunedMass  = -1.;
    I.softdropMass = -1.;
    I.softdropPuppiMass = -1.;
    I.CHSprunedMassCorr = -1.;
    I.CHSsoftdropMassCorr = -1.;
    I.prunedMassCorr = -1.;
    I.softdropMassCorr = -1.;
    I.softdropPuppiMassCorr = -1.;
    I.softdropPuppiMassCorrNotSmeared = -1.;
    I.pt1         = -1.;
    I.eta1        = -9.;
    I.phi1        = -9.;
    I.mass1       = -1.;
    I.CSV1        = -99.;
    //I.CSVR1       = -99.;
    //I.CSVR1Up     = -99.;
    //I.CSVR1Down   = -99.;
    I.CMVA1       = -99.;
    //I.CMVAR1      = -99.;
    //I.CMVAR1Up    = -99.;
    //I.CMVAR1Down  = -99.;
    I.flavour1    = -1.;
    I.pt2         = -1.;
    I.eta2        = -9.;
    I.phi2        = -9.;
    I.mass2       = -1.;
    I.CSV2        = -99.;
    //I.CSVR2       = -99.;
    //I.CSVR2Up     = -99.;
    //I.CSVR2Down   = -99.;
    I.CMVA2       = -99.;
    //I.CMVAR2      = -99.;
    //I.CMVAR2Up    = -99.;
    //I.CMVAR2Down  = -99.;
    I.flavour2    = -1.;

    I.pt1SDP      = -1.;
    I.eta1SDP     = -9.;
    I.phi1SDP     = -9.;
    I.mass1SDP    = -1.;
    I.CSV1SDP     = -99.;
    I.CMVA1SDP    = -99.;
    I.flavour1SDP = -1.;
    I.pt2SDP      = -1.;
    I.eta2SDP     = -9.;
    I.phi2SDP     = -9.;
    I.mass2SDP    = -1.;
    I.CSV2SDP     = -99.;
    I.CMVA2SDP    = -99.;
    I.flavour2SDP = -1.;

    I.dR          = -1.;
    I.Tau21       = -1.;
    I.puppiTau21  = -1.;
    I.ddtTau21    = -1.;
    I.BDSV        = -99.;
    I.chf         = -1.;
    I.nhf         = -1.;
    I.phf         = -1.;
    I.elf         = -1.;
    I.muf         = -1.;
    I.chm         = -1.;
    I.npr         = -1.;
    I.partonFlavour     = 0;
    I.hadronFlavour     = 0;
    I.mother      = false;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isTightLepVeto     = false;
    //I.isMatched   = false;
    //I.JESUnc      = -1.;
    //I.ptJERUp     = -1.;
    //I.etaJERUp    = -1.;
    //I.phiJERUp    = -9.;
    //I.energyJERUp = -1.;
    //I.ptJERDown   = -1.;
    //I.etaJERDown  = -1.;
    //I.phiJERDown  = -9.;
    //I.energyJERDown = -1.;
    //I.smearFact   = -1.;
    //I.smearFactUp   = -1.;
    //I.smearFactDown = -1.;
    //I.softdropPuppiMassCorrJMS = -1.;
    //I.softdropPuppiMassCorrJMSUp = -1.;
    //I.softdropPuppiMassCorrJMSDown = -1.;
    //I.softdropPuppiMassCorrJMR = -1.;
    //I.softdropPuppiMassCorrJMRUp = -1.;
    //I.softdropPuppiMassCorrJMRDown = -1.;
    I.dR_q1       = 1000.;
    I.dR_q2       = 1000.;
    I.dR_q3       = 1000.;
    I.dR_q4       = 1000.;
    I.m_q1        = false;
    I.m_q2        = false;
    I.m_q3        = false;
    I.m_q4        = false;
    I.dR_pi1      = 1000.;
    I.dR_pi2      = 1000.;
    I.matchBquark = -1;
    I.matchLL     = -1;
    I.dR_q1_sj1   = 1000;
    I.dR_q2_sj1   = 1000;
    I.dR_q3_sj1   = 1000;
    I.dR_q4_sj1   = 1000;
    I.dR_pi1_sj1  = 1000;
    I.dR_pi2_sj1  = 1000;
    I.dR_q1_sj2   = 1000;
    I.dR_q2_sj2   = 1000;
    I.dR_q3_sj2   = 1000;
    I.dR_q4_sj2   = 1000;
    I.dR_pi1_sj2  = 1000;
    I.dR_pi2_sj2  = 1000;
}

//std::string ObjectsFormat::ListCustomFatJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:ptRaw/F:ptUnc/F:dPhi_met/F:dPhi_Jet1/F:puId/F:CSV/F:CSVR/F:CSVRUp/F:CSVRDown/F:CHSprunedMass/F:CHSsoftdropMass/F:prunedMass/F:softdropMass/F:softdropPuppiMass/F:CHSprunedMassCorr/F:CHSsoftdropMassCorr/F:prunedMassCorr/F:softdropMassCorr/F:softdropPuppiMassCorr/F:softdropPuppiMassCorrNotSmeared/F:pt1/F:eta1/F:phi1/F:mass1/F:CSV1/F:CSVR1/F:CSVR1Up/F:CSVR1Down/F:CMVA1/F:CMVAR1/F:CMVAR1Up/F:CMVAR1Down/F:flavour1/F:pt2/F:eta2/F:phi2/F:mass2/F:CSV2/F:CSVR2/F:CSVR2Up/F:CSVR2Down/F:CMVA2/F:CMVAR2/F:CMVAR2Up/F:CMVAR2Down/F:flavour2/F:dR/F:Tau21/F:puppiTau21/F:ddtTau21/F:BDSV/F:chf/F:nhf/F:phf/F:elf/F:muf/F:chm/I:npr/I:flavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:isCSVL/O:isCSVM/O:isCSVT/O:isMatched/O:JESUnc/F:ptJERUp/F:etaJERUp/F:phiJERUp/F:energyJERUp/F:ptJERDown/F:etaJERDown/F:phiJERDown/F:energyJERDown/F:smearFact/F:smearFactUp/F:smearFactDown/F:softdropPuppiMassCorrJMS/F:softdropPuppiMassCorrJMSUp/F:softdropPuppiMassCorrJMSDown/F:softdropPuppiMassCorrJMR/F:softdropPuppiMassCorrJMRUp/F:softdropPuppiMassCorrJMRDown/F:dR_q1/F:dR_q2/F:dR_q3/F:dR_q4/F:m_q1/O:m_q2/O:m_q3/O:m_q4/O:dR_pi1/F:dR_pi2/F:matchBquark/I:matchLL/I";}

//simplified version:
std::string ObjectsFormat::ListCustomFatJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:ptRaw/F:ptUnc/F:dPhi_met/F:dPhi_Jet1/F:puId/F:CSV/F:CHSprunedMass/F:CHSsoftdropMass/F:prunedMass/F:softdropMass/F:softdropPuppiMass/F:CHSprunedMassCorr/F:CHSsoftdropMassCorr/F:prunedMassCorr/F:softdropMassCorr/F:softdropPuppiMassCorr/F:softdropPuppiMassCorrNotSmeared/F:pt1/F:eta1/F:phi1/F:mass1/F:CSV1/F:CMVA1/F:flavour1/F:pt2/F:eta2/F:phi2/F:mass2/F:CSV2/F:CMVA2/F:flavour2/F:pt1SDP/F:eta1SDP/F:phi1SDP/F:mass1SDP/F:CSV1SDP/F:CMVA1SDP/F:flavour1SDP/F:pt2SDP/F:eta2SDP/F:phi2SDP/F:mass2SDP/F:CSV2SDP/F:CMVA2SDP/F:flavour2SDP/F:dR/F:Tau21/F:puppiTau21/F:ddtTau21/F:BDSV/F:chf/F:nhf/F:phf/F:elf/F:muf/F:chm/I:npr/I:partonFlavour/I:hadronFlavour/I:mother/I:isLoose/O:isMedium/O:isTight/O:isTightLepVeto/O:dR_q1/F:dR_q2/F:dR_q3/F:dR_q4/F:m_q1/O:m_q2/O:m_q3/O:m_q4/O:dR_pi1/F:dR_pi2/F:matchBquark/I:matchLL/I:dR_q1_sj1/F:dR_q2_sj1/F:dR_q3_sj1/F:dR_q4_sj1/F:dR_pi1_sj1/F:dR_pi2_sj1/F:dR_q1_sj2/F:dR_q2_sj2/F:dR_q3_sj2/F:dR_q4_sj2/F:dR_pi1_sj2/F:dR_pi2_sj2/F";}



//*******************//
//  Missing energy   //
//*******************//

//void ObjectsFormat::FillMEtType(MEtType& I, const pat::MET* R, bool isMC) {
//    I.pt          = R->pt();
//    I.eta         = R->eta();
//    I.phi         = R->phi();
//    I.sign        = R->metSignificance();
//}

//void ObjectsFormat::ResetMEtType(MEtType& I) {
//    I.pt          = -1.;
//    I.eta         = -9.;
//    I.phi         = -9.;
//    I.sign        = -1.;
//}

//std::string ObjectsFormat::ListMEtType() {return "pt/F:eta/F:phi/F:sign/F";}


void ObjectsFormat::FillMEtType(MEtType& I, const pat::MET* R, bool isMC) {
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.sign        = R->metSignificance();
    I.ptShiftJetResUp = R->hasUserFloat("ptShiftJetResUp") ? R->userFloat("ptShiftJetResUp") : -1;
    I.ptShiftJetResDown = R->hasUserFloat("ptShiftJetResDown") ? R->userFloat("ptShiftJetResDown") : -1;
    I.ptShiftJetEnUp = R->hasUserFloat("ptShiftJetEnUp") ? R->userFloat("ptShiftJetEnUp") : -1;
    I.ptShiftJetEnDown = R->hasUserFloat("ptShiftJetEnDown") ? R->userFloat("ptShiftJetEnDown") : -1;
    I.ptShiftUnclusteredEnUp = R->hasUserFloat("ptShiftUnclusteredEnUp") ? R->userFloat("ptShiftUnclusteredEnUp") : -1;
    I.ptShiftUnclusteredEnDown = R->hasUserFloat("ptShiftUnclusteredEnDown") ? R->userFloat("ptShiftUnclusteredEnDown") : -1;
    I.ptShiftJetResUpSmear = R->hasUserFloat("ptShiftJetResUpSmear") ? R->userFloat("ptShiftJetResUpSmear") : -1.;
    I.ptShiftJetResDownSmear = R->hasUserFloat("ptShiftJetResDownSmear") ? R->userFloat("ptShiftJetResDownSmear") : -1.;
    I.ptRaw       = R->hasUserFloat("ptRaw") ? R->userFloat("ptRaw") : -1.;
    I.phiRaw      = R->hasUserFloat("phiRaw") ? R->userFloat("phiRaw") : -9.;
    //I.ptType1     = R->hasUserFloat("ptType1") ? R->userFloat("ptType1") : -1.;
    //I.phiType1    = R->hasUserFloat("phiType1") ? R->userFloat("phiType1") : -9.;
    if(isMC && R->genMET()) {I.ptGen       = R->genMET()->pt();}
    if(isMC && R->genMET()) {I.phiGen      = R->genMET()->phi();}
    //I.ptScaleUp   = R->hasUserFloat("ptScaleUp") ? R->userFloat("ptScaleUp") : -1.;
    //I.ptScaleDown = R->hasUserFloat("ptScaleDown") ? R->userFloat("ptScaleDown") : -1.;
    //I.ptResUp     = R->hasUserFloat("ptResUp") ? R->userFloat("ptResUp") : -1.;
    //I.ptResDown   = R->hasUserFloat("ptResDown") ? R->userFloat("ptResDown") : -1.;
    I.ptCalo      = R->caloMETPt();
}

void ObjectsFormat::ResetMEtType(MEtType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.sign        = -1.;
    I.ptShiftJetResUp = -1.;
    I.ptShiftJetResDown = -1.;
    I.ptShiftJetEnUp = -1.;
    I.ptShiftJetEnDown = -1.;
    I.ptShiftUnclusteredEnUp = -1.;
    I.ptShiftUnclusteredEnDown = -1.;
    I.ptShiftJetResUpSmear = -1.;
    I.ptShiftJetResDownSmear = -1.;
    I.ptRaw       = -1.;
    I.phiRaw      = -9.;
    //I.ptType1     = -1.;
    //I.phiType1    = -9.;
    I.ptGen       = -1.;
    I.phiGen      = -9.;
    //I.ptScaleUp   = -1.;
    //I.ptScaleDown = -1.;
    //I.ptResUp     = -1.;
    //I.ptResDown   = -1.;
    I.ptCalo      = -1.;
}

//std::string ObjectsFormat::ListMEtType() {return "pt/F:eta/F:phi/F:sign/F:ptRaw/F:phiRaw/F:ptGen/F:phiGen/F:ptCalo/F";}
std::string ObjectsFormat::ListMEtType() {return "pt/F:eta/F:phi/F:sign/F:ptShiftJetResUp/F:ptShiftJetResDown/F:ptShiftJetEnUp/F:ptShiftJetEnDown/F:ptShiftUnclusteredEnUp/F:ptShiftUnclusteredEnDown/F:ptShiftJetResUpSmear/F:ptShiftJetResDownSmear/F:ptRaw/F:phiRaw/F:ptGen/F:phiGen/F:ptCalo/F";}



void ObjectsFormat::FillMEtFullType(MEtFullType& I, const pat::MET* R, bool isMC) {
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.sign        = R->metSignificance();
    I.ptRaw       = R->uncorPt();
    I.phiRaw      = R->uncorPhi();
    if(isMC && R->genMET()) {I.ptGen       = R->genMET()->pt();}
    if(isMC && R->genMET()) {I.phiGen      = R->genMET()->phi();}
    I.ptJERUp     = R->shiftedPt(pat::MET::METUncertainty::JetResUp);
    I.ptJERDown   = R->shiftedPt(pat::MET::METUncertainty::JetResDown);
    I.ptJERUpSmear     = R->shiftedPt(pat::MET::METUncertainty::JetResUpSmear);
    I.ptJERDownSmear   = R->shiftedPt(pat::MET::METUncertainty::JetResDownSmear);
    I.ptJESUp     = R->shiftedPt(pat::MET::METUncertainty::JetEnUp);
    I.ptJESDown   = R->shiftedPt(pat::MET::METUncertainty::JetEnDown);
    I.ptMUSUp     = R->shiftedPt(pat::MET::METUncertainty::MuonEnUp);
    I.ptMUSDown   = R->shiftedPt(pat::MET::METUncertainty::MuonEnDown);
    I.ptELSUp     = R->shiftedPt(pat::MET::METUncertainty::ElectronEnUp);
    I.ptELSDown   = R->shiftedPt(pat::MET::METUncertainty::ElectronEnDown);
    I.ptTAUUp     = R->shiftedPt(pat::MET::METUncertainty::TauEnUp);
    I.ptTAUDown   = R->shiftedPt(pat::MET::METUncertainty::TauEnDown);
    I.ptUNCUp     = R->shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp);
    I.ptUNCDown   = R->shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown);
    I.ptPHOUp     = R->shiftedPt(pat::MET::METUncertainty::PhotonEnUp);
    I.ptPHODown   = R->shiftedPt(pat::MET::METUncertainty::PhotonEnDown);
    I.phf         = R->NeutralEMFraction();
    I.nhf         = R->NeutralHadEtFraction();
    I.elf         = R->ChargedEMEtFraction();
    I.chf         = R->ChargedHadEtFraction();
    I.muf         = R->MuonEtFraction();
}

void ObjectsFormat::ResetMEtFullType(MEtFullType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.sign        = -1.;
    I.ptRaw       = -1.;
    I.phiRaw      = -9.;
    I.ptGen       = -1.;
    I.phiGen      = -9.;
    I.ptJERUp     = -1.;
    I.ptJERDown   = -1.;
    I.ptJERUpSmear     = -1.;
    I.ptJERDownSmear   = -1.;
    I.ptJESUp     = -1.;
    I.ptJESDown   = -1.;
    I.ptMUSUp     = -1.;
    I.ptMUSDown   = -1.;
    I.ptELSUp     = -1.;
    I.ptELSDown   = -1.;
    I.ptTAUUp     = -1.;
    I.ptTAUDown   = -1.;
    I.ptUNCUp     = -1.;
    I.ptUNCDown   = -1.;
    I.ptPHOUp     = -1.;
    I.ptPHODown   = -1.;
    I.phf         = -1.;
    I.nhf         = -1.;
    I.elf         = -1.;
    I.chf         = -1.;
    I.muf         = -1.;
}

std::string ObjectsFormat::ListMEtFullType() {return "pt/F:eta/F:phi/F:sign/F:ptRaw/F:phiRaw/F:ptGen/F:phiGen/F:ptJERUp/F:ptJERDown/F:ptJERUpSmear/F:ptJERDownSmear/F:ptJESUp/F:ptJESDown/F:ptMUSUp/F:ptMUSDown/F:ptELSUp/F:ptELSDown/F:ptTAUUp/F:ptTAUDown/F:ptUNCUp/F:ptUNCDown/F:ptPHOUp/F:ptPHODown/F:phf/F:nhf/F:elf/F:chf/F:muf/F";}


void ObjectsFormat::FillCandidateType(CandidateType& I, pat::CompositeCandidate* R, bool isMC) {
  if(!R) return;
  if(R->numberOfDaughters() == 0) return;
  I.pt          = R->pt();
  I.eta         = R->eta();
  I.phi         = R->phi();
  I.mass        = R->mass();
//  I.tmass       = R->numberOfDaughters()>1 ? sqrt( 2.*R->daughter(0)->pt()*R->daughter(1)->pt()*(1.-cos(deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi())) ) ) : -1.;
//Lisa
  I.tmass       = R->numberOfDaughters()>1 ? sqrt( 2.*R->daughter(0)->pt()*R->daughter(1)->pt()*(1.-cos(deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi())) ) ) : (R->numberOfDaughters()==1 ? R->mt() : -1.);
//  I.softdropMass= R->hasUserFloat("softdropMass") ? R->userFloat("softdropMass") : -1.;
  I.dR          = R->numberOfDaughters()>1 ? deltaR(*R->daughter(0), *R->daughter(1)) : -1.;
  I.dEta        = R->numberOfDaughters()>1 ? fabs( R->daughter(0)->eta() - R->daughter(1)->eta() ) : -1.;
  I.dPhi        = R->numberOfDaughters()>1 ? fabs( deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi()) ) : -1.;
  I.twist       = R->numberOfDaughters()>1 ? fabs(atan( deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi())/fabs(R->daughter(0)->eta()-R->daughter(1)->eta()) )) : -1.;
//  I.angle       = (R->daughter(0)->momentum().unit()).Dot(R->daughter(1)->momentum().unit()); // acos()
//  I.ptBalance   = ((R->daughter(0)->pt()-R->daughter(1)->pt()) / (R->daughter(0)->pt()+R->daughter(1)->pt())/2.);
//  I.centrality  = (R->daughter(0)->pt()+R->daughter(0)->pt()) / (R->daughter(0)->p()+R->daughter(0)->p());// /R->mass();
//  I.charge      = R->charge();
}

void ObjectsFormat::ResetCandidateType(CandidateType& I) {
  I.pt          = -1.;
  I.eta         = -9.;
  I.phi         = -9.;
  I.mass        = -1.;
  I.tmass       = -1.;
//  I.softdropMass = -1.;
  I.dR          = -1.;
  I.dEta        = -1.;
  I.dPhi        = -1.;
  I.twist       = -1.;
//  I.angle       = -1.;
//  I.ptBalance   = -1.;
//  I.centrality  = -1.;
//  I.charge      = -1.;
}

std::string ObjectsFormat::ListCandidateType() {return "pt/F:eta/F:phi/F:mass/F:tmass/F:dR/F:dEta/F:dPhi/F:twist/F";}

/*
void ObjectsFormat::FillCandidateType(CandidateType& I, const reco::Candidate::LorentzVector* V1, const reco::Candidate::LorentzVector* V2) {
  if(!V1 || !V2) return;
  reco::Candidate::LorentzVector V(*V1+*V2);
  I.pt          = V.pt();
  I.eta         = V.eta();
  I.phi         = V.phi();
  I.et          = V.Et();
  I.p           = V.P();
  I.energy      = V.energy();
  I.mass        = V.mass();
  I.dEta        = fabs( V1->eta() - V2->eta() );
  I.dPhi        = M_PI - fabs(fabs(V1->phi()-V2->phi()) - M_PI);
  I.dR          = sqrt( I.dEta*I.dEta + I.dPhi*I.dPhi );
  I.twist       = fabs(atan( I.dPhi/I.dEta ));
//  I.angle       = (R->daughter(0)->momentum().unit()).Dot(R->daughter(1)->momentum().unit()); // acos()
//  I.ptBalance   = ((R->daughter(0)->pt()-R->daughter(1)->pt()) / (R->daughter(0)->pt()+R->daughter(1)->pt())/2.);
//  I.centrality  = (R->daughter(0)->pt()+R->daughter(0)->pt()) / (R->daughter(0)->p()+R->daughter(0)->p());// /R->mass();
//  I.charge      = R->charge();
}
*/

void ObjectsFormat::FillLorentzType(LorentzType& I, const reco::Candidate::LorentzVector* V) {
  I.pt          = V->pt();
  I.eta         = V->eta();
  I.phi         = V->phi();
  I.energy      = V->energy();
  I.mass        = V->mass();
}

void ObjectsFormat::ResetLorentzType(LorentzType& I) {
  I.pt          = -1.;
  I.eta         = -9.;
  I.phi         = -9.;
  I.energy      = -1.;
  I.mass        = -1.;
}

std::string ObjectsFormat::ListLorentzType() {return "pt/F:eta/F:phi/F:energy/F:mass/F";}


//************************//
//       GenParticles     //
//***********************//

void ObjectsFormat::FillGenPType(GenPType& I, const reco::GenParticle* R) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.rapidity    = R->rapidity();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.status      = R->status();
    I.radius      = R->mother()? sqrt(pow(R->vx() - R->mother()->vx(),2) + pow(R->vy() - R->mother()->vy(),2) + pow(R->vz() - R->mother()->vz(),2)) : -1000.;
    I.motherid    = R->mother()? R->mother()->pdgId() : 0;
    I.vx          = R->vx();
    I.vy          = R->vy();
    I.vz          = R->vz();
}


void ObjectsFormat::ResetGenPType(GenPType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.rapidity    = -99999.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.charge      = 0;
    I.pdgId       = 0;
    I.status      = 0;
    I.radius      = -1.;
    I.motherid    = 0;
    I.vx          = -99;
    I.vy          = -99;
    I.vz          = -99;

}

std::string ObjectsFormat::ListGenPType() {return "pt/F:eta/F:rapidity/F:phi/F:mass/F:energy/F:charge/I:pdgId/I:status/I:radius/F:motherid/I:vx/F:vy/F:vz/F";}



//************************//
//   TriggerObjectType    //
//***********************//

void ObjectsFormat::FillTriggerObjectType(TriggerObjectType& I, const pat::TriggerObjectStandAlone* R) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.charge      = R->charge();
}


void ObjectsFormat::ResetTriggerObjectType(TriggerObjectType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.charge      = 0;
}

std::string ObjectsFormat::ListTriggerObjectType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:charge/I";}




//*******************//
//    Calo Jets      //
//*******************//

void ObjectsFormat::FillCaloJetType(CaloJetType& I, const reco::CaloJet* R, bool isMC, bool isMatched) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.emEnergyFraction = R->emEnergyFraction();
    I.emEnergyInEB     = R->emEnergyInEB();
    I.emEnergyInEE     = R->emEnergyInEE();
    I.emEnergyInHF     = R->emEnergyInHF();
    I.energyFractionHadronic = R->energyFractionHadronic();
    I.hadEnergyInHB    = R->hadEnergyInHB();
    I.hadEnergyInHE    = R->hadEnergyInHE();
    I.hadEnergyInHF    = R->hadEnergyInHF();
    I.hadEnergyInHO    = R->hadEnergyInHO();
    I.longLived        = R->longLived();
    I.maxEInEmTowers   = R->maxEInEmTowers();
    I.maxEInHadTowers  = R->maxEInHadTowers();
    I.n60              = R->n60();
    I.n90              = R->n90();
    I.nConstituents    = R->nConstituents();
    I.nPasses          = R->nPasses();
    I.isGenMatched     = isMatched;
}

void ObjectsFormat::ResetCaloJetType(CaloJetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.emEnergyFraction = -1.;
    I.emEnergyInEB     = -1.;
    I.emEnergyInEE     = -1.;
    I.emEnergyInHF     = -1.;
    I.energyFractionHadronic = -1.;
    I.hadEnergyInHB    = -1.;
    I.hadEnergyInHE    = -1.;
    I.hadEnergyInHF    = -1.;
    I.hadEnergyInHO    = -1.;
    I.longLived        = false;
    I.maxEInEmTowers   = -1.;
    I.maxEInHadTowers  = -1.;
    I.n60              = -1;
    I.n90              = -1;
    I.nConstituents    = -1;
    I.nPasses          = -1;
    I.isGenMatched        = false;
}

std::string ObjectsFormat::ListCaloJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:emEnergyFraction/F:emEnergyInEB/F:emEnergyInEE/F:emEnergyInHF/F:energyFractionHadronic/F:hadEnergyInHB/F:hadEnergyInHE/F:hadEnergyInHF/F:hadEnergyInHO/F:longLived/O:maxEInEmTowers/F:maxEInHadTowers/F:n60/I:n90/I:nConstituents/I:nPasses/I:isGenMatched/O";}


//**************************//
//      PFCandidates       //
//************************//

void ObjectsFormat::FillPFCandidateType(PFCandidateType& I, const pat::PackedCandidate* C, int jetIdx, int fatJetIdx, int vtxIdx) {
  I.pt              = C->pt();
  I.eta             = C->eta();
  I.phi             = C->phi();
  I.energy          = C->energy();
  I.mass            = C->mass();
  I.isTrack         = C->charge() ? true : false;
  I.jetIndex        = jetIdx;
  I.fatJetIndex     = fatJetIdx;
  I.pvIndex         = vtxIdx;
  I.px              = C->px();
  I.py              = C->py();
  I.pz              = C->pz();
  I.dxy             = C->dxy();
  I.dxyError        = C->dxyError();
  I.dz              = C->dzAssociatedPV();
  I.dzError         = C->dzError();
  I.nHits           = C->numberOfHits();
  I.nPixelHits      = C->numberOfPixelHits();
  I.lostInnerHits   = C->lostInnerHits();
  I.charge          = C->charge();
  //I.hcalFraction    = C->hcalFraction();
  I.POCA_x          = C->vx();
  I.POCA_y          = C->vy();
  I.POCA_z          = C->vz();
  I.POCA_phi        = C->phiAtVtx();
  I.pdgId           = C->pdgId();
  //I.longLived       = C->longLived();
  I.ptError         = C->pseudoTrack().ptError();
  I.etaError        = C->pseudoTrack().etaError();
  I.phiError        = C->pseudoTrack().phiError();
  I.theta           = C->pseudoTrack().theta();
  I.thetaError      = C->pseudoTrack().thetaError();
  //I.innerDetId      = C->pseudoTrack().innerDetId();
  //I.innerPosition_x = C->pseudoTrack().innerPosition().x();
  //I.innerPosition_y = C->pseudoTrack().innerPosition().y();
  //I.innerPosition_z = C->pseudoTrack().innerPosition().z();
  //I.innerMomentum_x = C->pseudoTrack().innerMomentum().x();
  //I.innerMomentum_y = C->pseudoTrack().innerMomentum().y();
  //I.innerMomentum_z = C->pseudoTrack().innerMomentum().z();
  I.chi2            = C->pseudoTrack().chi2();
  I.ndof            = C->pseudoTrack().ndof();
  I.normalizedChi2  = C->pseudoTrack().normalizedChi2();
}


void ObjectsFormat::ResetPFCandidateType(PFCandidateType& I) {
  I.pt              = -1.;
  I.eta             = -9.;
  I.phi             = -9.;
  I.energy          = -1.;
  I.mass            = -1.;
  I.isTrack         = false;
  I.jetIndex        = -1;
  I.pvIndex         = -1;
  I.px              = -99.;
  I.py              = -99.;
  I.pz              = -99.;
  I.dxy             = -99.;
  I.dxyError        = -99.;
  I.dz              = -99.;
  I.dzError         = -99.;
  I.nHits           = -1;
  I.nPixelHits      = -1;
  I.lostInnerHits   = -9;
  I.charge          = -99;
  //I.hcalFraction    = -1.;
  I.POCA_x          = -99.;
  I.POCA_y          = -99.;
  I.POCA_z          = -99.;
  I.POCA_phi        = -9.;
  I.pdgId           = 0;
  //I.longLived       = -1;
  I.ptError         = -1.;
  I.etaError        = -1.;
  I.phiError        = -1.;
  I.theta           = -9.;
  I.thetaError      = -1.;
  //I.innerDetId      = -1;
  //I.innerPosition_x = -9999.;
  //I.innerPosition_y = -9999.;
  //I.innerPosition_z = -9999.;
  //I.innerMomentum_x = -9999.;
  //I.innerMomentum_y = -9999.;
  //I.innerMomentum_z = -9999.;
  I.chi2            = -1.;
  I.ndof            = -1;
  I.normalizedChi2  = -1.;
}

std::string ObjectsFormat::ListPFCandidateType() {return"pt/F:eta/F:phi/F:energy/F:mass/F:isTrack/O:jetIndex/I:pvIndex/I:px/F:py/F:pz/F:dxy/F:dxyError/F:dz/F:dzError/F:nHits/I:nPixelHits/I:lostInnerHits/I:charge/I:POCA_x/F:POCA_y/F:POCA_z/F:POCA_phi/I:pdgId/I:ptError/F:etaError/F:phiError/F:theta/F:thetaError/F:chi2/F:ndof/F:normalizedChi2/F";}


//**************************//
//        Vertices         //
//************************//

void ObjectsFormat::FillPrimVertexType(VertexType& I, const reco::Vertex* R) {
  if(!R) return;
  I.chi2        = R->chi2();
  I.ndof        = R->ndof();
  I.x           = R->x();
  I.y           = R->y();
  I.z           = R->z();
  I.px          = R->p4().px();
  I.py          = R->p4().py();
  I.pz          = R->p4().pz();
  I.mass        = R->p4().mass();
  I.isValid     = R->isValid();
  I.isFake      = R->isFake();
  I.pt          = -1.;
}


void ObjectsFormat::FillSecVertexType(VertexType& I, const reco::VertexCompositePtrCandidate* R) {
  if(!R) return;
  I.chi2        = R->vertexChi2();
  I.ndof        = R->vertexNdof();
  I.x           = R->vx();
  I.y           = R->vy();
  I.z           = R->vz();
  I.px          = R->px();
  I.py          = R->py();
  I.pz          = R->pz();
  I.mass        = R->mass();
  I.isValid     = false;
  I.isFake      = false;
  I.pt          = R->pt();
}

void ObjectsFormat::FillBtagSecVertexType(VertexType& I, const reco::CandSecondaryVertexTagInfo * R, const reco::CandIPTagInfo * IP, unsigned int index, float dRSVJet, int jetIndex) {
  if(!R) return;
  if(R->nVertices() < 1) return;
  I.chi2        = R->secondaryVertex(index).vertexChi2();
  I.ndof        = R->secondaryVertex(index).vertexNdof();
  I.x           = R->secondaryVertex(index).vx();
  I.y           = R->secondaryVertex(index).vy();
  I.z           = R->secondaryVertex(index).vz();
  I.px          = R->secondaryVertex(index).px();
  I.py          = R->secondaryVertex(index).py();
  I.pz          = R->secondaryVertex(index).pz();
  I.mass        = R->secondaryVertex(index).mass();
  I.isValid     = false;
  I.isFake      = false;
  I.pt          = R->secondaryVertex(index).pt();
  I.flightDist2d = R->flightDistance(index, true).value();
  I.flightDist2dErr = R->flightDistance(index).error();
  I.flightDist2dSig = R->flightDistance(index, true).value()/R->flightDistance(index, true).error();
  I.flightDist3d = R->flightDistance(index).value();
  I.flightDist3dErr = R->flightDistance(index).error();
  I.flightDist3dSig = R->flightDistance(index).value()/R->flightDistance(index).error();
  I.nVertexTracks = R->nVertexTracks(index);
  I.dRSVJet = dRSVJet;
  I.jetIndex = jetIndex;
  if (R->taggingVariables().checkTag(reco::btau::trackJetPt)) I.trackJetPt= R->taggingVariables().get(reco::btau::trackJetPt);
  if (R->taggingVariables().checkTag(reco::btau::tau1_vertexMass)) I.vertexMass = R->taggingVariables().get(reco::btau::tau1_vertexMass);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip3dSig_3)) I.trackSip3dSig_3 = R->taggingVariables().get(reco::btau::trackSip3dSig_3);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip3dSig_2)) I.trackSip3dSig_2 = R->taggingVariables().get(reco::btau::trackSip3dSig_2);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip3dSig_1)) I.trackSip3dSig_1 = R->taggingVariables().get(reco::btau::trackSip3dSig_1);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip3dSig_0)) I.trackSip3dSig_0 = R->taggingVariables().get(reco::btau::trackSip3dSig_0);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip2dSigAboveCharm)) I.trackSip2dSigAboveCharm_0 = R->taggingVariables().get(reco::btau::trackSip2dSigAboveCharm);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip2dSigAboveBottom_0)) I.trackSip2dSigAboveBottom_0 = R->taggingVariables().get(reco::btau::trackSip2dSigAboveBottom_0);
  // if (R->taggingVariables().checkTag(reco::btau::trackSip2dSigAboveBottom_1)) I.trackSip2dSigAboveBottom_0 = R->taggingVariables().get(reco::btau::trackSip2dSigAboveBottom_1);

  if (IP->taggingVariables().checkTag(reco::btau::trackEta)) I.trackEta = IP->taggingVariables().get(reco::btau::trackEta,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackMomentum)) I.trackMomentum = IP->taggingVariables().get(reco::btau::trackMomentum,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackPhi)) I.trackPhi = IP->taggingVariables().get(reco::btau::trackPhi,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackPtRel)) I.trackPtRel = IP->taggingVariables().get(reco::btau::trackPtRel,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackPPar)) I.trackPPar = IP->taggingVariables().get(reco::btau::trackPPar,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackEtaRel)) I.trackEtaRel = IP->taggingVariables().get(reco::btau::trackEtaRel,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackDeltaR)) I.trackDeltaR = IP->taggingVariables().get(reco::btau::trackDeltaR,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackPtRatio)) I.trackPtRatio = IP->taggingVariables().get(reco::btau::trackPtRatio,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackPParRatio)) I.trackPParRatio = IP->taggingVariables().get(reco::btau::trackPParRatio,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackSip2dVal)) I.trackSip2dVal = IP->taggingVariables().get(reco::btau::trackSip2dVal,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackSip2dSig)) I.trackSip2dSig = IP->taggingVariables().get(reco::btau::trackSip2dSig,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackSip3dVal)) I.trackSip3dVal = IP->taggingVariables().get(reco::btau::trackSip3dVal,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackSip3dSig)) I.trackSip3dSig = IP->taggingVariables().get(reco::btau::trackSip3dSig,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackDecayLenVal)) I.trackDecayLenVal = IP->taggingVariables().get(reco::btau::trackDecayLenVal,false);
  if (IP->taggingVariables().checkTag(reco::btau::trackDecayLenSig)) I.trackDecayLenSig = IP->taggingVariables().get(reco::btau::trackDecayLenSig,false);
  if (R->taggingVariables().checkTag(reco::btau::trackJetPt)) I.trackJetPt = IP->taggingVariables().get(reco::btau::trackJetPt);
  // if (IP->taggingVariables().checkTag()) I. = IP->taggingVariables().get();

  // if (IP->taggingVariables().checkTag()) I. = IP->taggingVariables().get();

}



void ObjectsFormat::ResetVertexType(VertexType& I) {
  I.chi2        = -1.;
  I.ndof        = -1.;
  I.x           = -1.;
  I.y           = -1.;
  I.z           = -1.;
  I.px          = -1.;
  I.py          = -1.;
  I.pz          = -1.;
  I.mass        = -1.;
  I.isValid     = false;
  I.isFake      = false;
  I.pt          = -1.;
  I.flightDist2d = -100.;
  I.flightDist2dErr = -100.;
  I.flightDist2dSig = -100.;
  I.flightDist3d = -100.;
  I.flightDist3dErr = -100.;
  I.flightDist3dSig = -100.;
  I.nVertexTracks = -1;
  I.dRSVJet = -100.;
  I.jetIndex = -1;
  I.vertexMass = -1.;
  I.trackEta = -100.;
  I.trackMomentum = -100.;
  I.trackPhi = -100.;
  I.trackPtRel = -100.;
  I.trackPPar = -100.;
  I.trackEtaRel = -100.;
  I.trackDeltaR = -100.;
  I.trackPtRatio = -100.;
  I.trackPParRatio = -100.;
  I.trackSip2dVal = -100.;
  I.trackSip2dSig = -100.;
  I.trackSip3dVal = -100.;
  I.trackSip3dSig = -100.;
  I.trackDecayLenVal = -100.;
  I.trackDecayLenSig = -100.;
  I.trackJetDistVal = -100.;
  I.trackJetDistSig = -100.;
  I.trackChi2 = -100.;
  I.trackNTotalHits = -1;
  I.trackNPixelHits = -1;
  I.trackJetPt = -100.;
}

std::string ObjectsFormat::ListVertexType() {return "chi2/F:ndof/F:x/F:y/F:z/F:px/F:py/F:pz/F:mass/F:isValid/O:isFake/O:pt/F:flightDist2d/F:flightDist2dErr/F:flightDist2dSig/F:flightDist3d/F:flightDist3dErr/F:flightDist3dSig/F:nVertexTracks/I:dRSVJet/F:jetIndex/I:vertexMass/F:trackEta/F:trackMomentum/F:trackPhi/F:trackPtRel/F:trackPPar/F:trackEtaRel/F:trackDeltaR/F:trackPtRatio/F:trackPParRatio/F:trackSip2dVal/F:trackSip2dSig/F:trackSip3dVal/F:trackSip3dSig/F:trackDecayLenVal/F:trackDecayLenSig/F:trackJetDistVal/F:trackJetDistSig/F:trackChi2/F:trackNTotalHits/I:trackNPixelHits/I:trackJetPt/F";}

//*****************************//
//        SimplifiedJets       //
//*****************************//

void ObjectsFormat::FillSimplifiedJetType(SimplifiedJetType& I, const pat::Jet* R, bool isMC) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.cHadE       = R->chargedHadronEnergy();
    I.nHadE       = R->neutralHadronEnergy();
    I.cHadEFrac   = R->hasUserFloat("cHadEFrac") ? R->userFloat("cHadEFrac") : -1.;//R->chargedHadronEnergyFraction();
    I.nHadEFrac   = R->hasUserFloat("nHadEFrac") ? R->userFloat("nHadEFrac") : -1.;//R->neutralHadronEnergyFraction();
    I.nEmE        = R->neutralEmEnergy();
    I.nEmEFrac    = R->hasUserFloat("nEmEFrac") ? R->userFloat("nEmEFrac") : -1.;//R->neutralEmEnergyFraction();
    I.cEmE        = R->chargedEmEnergy();
    I.cEmEFrac    = R->hasUserFloat("cEmEFrac") ? R->userFloat("cEmEFrac") : -1.;//R->chargedEmEnergyFraction();
    I.cmuE        = R->chargedMuEnergy();
    I.cmuEFrac    = R->hasUserFloat("cmuEFrac") ? R->userFloat("cmuEFrac") : -1.;//R->chargedMuEnergyFraction();
    I.muE         = R->muonEnergy();
    I.muEFrac     = R->hasUserFloat("muEFrac") ? R->userFloat("muEFrac") : -1.;//R->muonEnergyFraction();
    I.eleE        = R->electronEnergy();
    I.eleEFrac    = R->hasUserFloat("eleEFrac") ? R->userFloat("eleEFrac") : -1.;//R->electronEnergyFraction();
    I.eleMulti    = R->electronMultiplicity();
    I.photonE     = R->photonEnergy();
    I.photonEFrac = R->hasUserFloat("photonEFrac") ? R->userFloat("photonEFrac") : -1.;//R->photonEnergyFraction();
    I.photonMulti = R->photonMultiplicity();
    I.cHadMulti   = R->chargedHadronMultiplicity();
    I.nHadMulti   = R->neutralHadronMultiplicity();
    I.npr         = R->chargedMultiplicity() + R->neutralMultiplicity();
    I.cMulti      = R->chargedMultiplicity();
    I.nMulti      = R->neutralMultiplicity();
    I.isLoose     = R->hasUserInt("isLoose") ? R->userInt("isLoose") : false;
    I.isMedium    = false;
    I.isTight     = R->hasUserInt("isTight") ? R->userInt("isTight") : false;
    I.isGenMatched = R->hasUserInt("isGenMatched") ? R->userInt("isGenMatched") : 0;
    if (R->tagInfoLabels().size() > 0){
      if(R->hasTagInfo("pfSecondaryVertex")){
	I.nSelectedTracks = R->tagInfoCandSecondaryVertex("pfSecondaryVertex")->nSelectedTracks();
      }
    }
    I.nConstituents = R->hasUserInt("nConstituents") ? R->userInt("nConstituents") : -1;
    I.TriggerMatched_VBFJet = R->hasUserInt("TriggerMatched_VBFJet") ? R->userInt("TriggerMatched_VBFJet") : 0;
    I.TriggerMatched_DisplacedJet = R->hasUserInt("TriggerMatched_DisplacedJet") ? R->userInt("TriggerMatched_DisplacedJet") : 0;
    I.TriggerMatched_TripleJet50 = R->hasUserInt("TriggerMatched_TripleJet50") ? R->userInt("TriggerMatched_TripleJet50") : 0;
}

void ObjectsFormat::ResetSimplifiedJetType(SimplifiedJetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.cHadE         = -1.;
    I.nHadE         = -1.;
    I.cHadEFrac     = -1.;
    I.nHadEFrac     = -1.;
    I.nEmE          = -1.;
    I.nEmEFrac      = -1.;
    I.cEmE          = -1.;
    I.cEmEFrac      = -1.;
    I.cmuE          = -1.;
    I.cmuEFrac      = -1.;
    I.muE           = -1.;
    I.muEFrac       = -1.;
    I.eleE          = -1.;
    I.eleEFrac      = -1.;
    I.eleMulti      = -1.;
    I.photonE       = -1.;
    I.photonEFrac   = -1.;
    I.photonMulti   = -1.;
    I.cHadMulti     = -1.;
    I.nHadMulti     = -1.;
    I.npr           = -1.;
    I.cMulti        = -1.;
    I.nMulti        = -1.;
    I.isLoose     = false;
    I.isMedium    = false;
    I.isTight     = false;
    I.isGenMatched = 0;
    I.nSelectedTracks = 0;
    I.nConstituents = -1;
    I.TriggerMatched_VBFJet = 0;
    I.TriggerMatched_DisplacedJet = 0;
    I.TriggerMatched_TripleJet50 = 0;
}

std::string ObjectsFormat::ListSimplifiedJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:cHadE/F:nHadE/F:cHadEFrac/F:nHadEFrac/F:nEmE/F:nEmEFrac/F:cEmE/F:cEmEFrac/F:cmuE/F:cmuEFrac/F:muE/F:muEFrac/F:eleE/F:eleEFrac/F:eleMulti/F:photonE/F:photonEFrac/F:photonMulti/F:cHadMulti/F:nHadMulti/F:npr/F:cMulti/F:nMulti/F:isLoose/O:isMedium/O:isTight/O:isGenMatched/I:nSelectedTracks/I:nConstituents/I:TriggerMatched_VBFJet/I:TriggerMatched_DisplacedJet/I:TriggerMatched_TripleJet50/I";}
