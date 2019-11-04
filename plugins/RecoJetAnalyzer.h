#ifndef RECOJETANALYZER_H
#define RECOJETANALYZER_H

#include <iostream>
#include <fstream>

#include <cmath>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "RecoilCorrector.h" // From: https://github.com/cms-met/MetTools/tree/master/RecoilCorrections
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//Reco classes
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "DataFormats/Common/interface/Ptr.h"

#include "FWCore/Utilities/interface/transform.h"

#include "TFile.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "FWCore/Utilities/interface/TypeID.h"




class RecoJetAnalyzer {
    public:
        RecoJetAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~RecoJetAnalyzer();
        virtual std::vector<reco::PFJet> FillJetVector(const edm::Event&);
        virtual void CorrectJet(reco::PFJet&, float, float, bool);
        virtual double CorrectMass(reco::PFJet&, float, float, bool);
        //virtual void CorrectPuppiMass(reco::PFJet&, bool);
        virtual void CleanJetsFromMuons(std::vector<reco::PFJet>&, std::vector<reco::Muon>&, float);
        virtual void CleanJetsFromElectrons(std::vector<reco::PFJet>&, std::vector<reco::GsfElectron>&, float);
        virtual void CleanFatJetsFromAK4(std::vector<reco::PFJet>&, std::vector<reco::PFJet>&, float);
        //virtual void AddVariables(std::vector<reco::PFJet>&, reco::PFMET&);
        //virtual void GenMatcher(std::vector<reco::PFJet>&, std::vector<reco::GenParticle>&, std::string);
        //virtual int GetNBJets(std::vector<reco::PFJet>&);
        virtual reco::PFMET FillMetVector(const edm::Event&);
	//virtual float GetMetTriggerEfficiency(reco::PFMET&);
        //virtual void ApplyRecoilCorrections(reco::PFMET&, const reco::Candidate::LorentzVector*, const reco::Candidate::LorentzVector*, int);
        virtual float CalculateHT(const edm::Event&, int, float, float);
        virtual bool isLooseJet(reco::PFJet&);
        virtual bool isTightJet(reco::PFJet&);
        virtual bool isTightLepVetoJet(reco::PFJet&);
        //virtual std::vector<float> ReshapeBtagDiscriminator(reco::PFJet&);
      
    private:

	edm::EDGetTokenT<std::vector<reco::PFJet> > JetToken;
        edm::EDGetTokenT<std::vector<reco::PFMET> > MetToken;
        //edm::EDGetTokenT<edm::ValueMap<float>> QGToken;
        int JetId;
        float Jet1Pt, Jet2Pt, JetEta;
        //bool AddQG;
        bool RecalibrateJets;
        bool RecalibrateMass;
        //bool RecalibratePuppiMass; 
	//std::string SoftdropPuppiMassString;
        bool SmearJets;
        std::string JECUncertaintyMC;
        std::string JECUncertaintyDATA;
        std::vector<std::string> JetCorrectorMC;
        std::vector<std::string> JetCorrectorDATA;
        std::vector<std::string> MassCorrectorMC;
        std::vector<std::string> MassCorrectorDATA;
        //std::string MassCorrectorPuppi;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;        
        edm::EDGetTokenT<double> RhoToken;
        //bool UseReshape;
        //std::string BTag;
        //int Jet1BTag, Jet2BTag;
        //std::string BTagDB;
        //bool UseRecoil;
        //std::string RecoilMCFile;
        //std::string RecoilDataFile;
        //std::string MetTriggerFileName;
        std::string JerName_res;
        std::string JerName_sf;
	//std::vector<std::string> BTagNames;
        float Rparameter;
        
        //TFile* PuppiCorrFile;
        //TF1* PuppiJECcorr_gen;
        //TF1* PuppiJECcorr_reco_0eta1v3;
        //TF1* PuppiJECcorr_reco_1v3eta2v5;

	//TFile* MetTriggerFile;
	//TH1F* MetTriggerHisto;
	//bool isMetTriggerFile;

        // JEC Uncertainty
        JetCorrectionUncertainty* jecUncMC;
        JetCorrectionUncertainty* jecUncDATA;
        
        boost::shared_ptr<FactorizedJetCorrector> jetCorrMC;
        boost::shared_ptr<FactorizedJetCorrector> jetCorrDATA;
        boost::shared_ptr<FactorizedJetCorrector> massCorrMC;
        boost::shared_ptr<FactorizedJetCorrector> massCorrDATA;
        
        // Btag calibrations
        //BTagCalibration       * calib;
	//std::map < int , BTagEntry::JetFlavor > flavour_map; 
	//std::map< BTagEntry::JetFlavor, std::vector<std::string>> syst_map; 
	//std::map<std::string, BTagCalibrationReader> cr_map;
	//std::string sf_mode;

        
        //JME
        JME::JetResolution              * resolution;
        JME::JetResolutionScaleFactor   * resolution_sf;        
        
        // Recoil corrections
        //RecoilCorrector* recoilCorr;
};

#endif
