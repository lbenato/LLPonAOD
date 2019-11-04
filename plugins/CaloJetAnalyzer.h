#ifndef CALOJETANALYZER_H
#define CALOJETANALYZER_H

#include <iostream>
#include <fstream>

#include <cmath>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
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


#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "DataFormats/Common/interface/Ptr.h"

#include "TFile.h"
#include "TH2.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "FWCore/Utilities/interface/TypeID.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"


class CaloJetAnalyzer {
    public:
        CaloJetAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~CaloJetAnalyzer();
        virtual std::vector<reco::CaloJet> FillJetVector(const edm::Event&);
        virtual void CorrectJet(reco::CaloJet&, float, float, bool);
        //virtual void CorrectMass(reco::CaloJet&, float, float, bool);
        //virtual void CorrectPuppiMass(reco::CaloJet&, bool);
        virtual void CleanJetsFromMuons(std::vector<reco::CaloJet>&, std::vector<pat::Muon>&, float);
        virtual void CleanJetsFromElectrons(std::vector<reco::CaloJet>&, std::vector<pat::Electron>&, float);
        //virtual void AddVariables(std::vector<reco::CaloJet>&, reco::PFMET&);
        virtual std::map<std::string,float> GenMatcherCalo(std::vector<reco::CaloJet>&, std::vector<reco::GenParticle>&, std::string);
        //virtual int GetNBJets(std::vector<reco::CaloJet>&);
        //virtual reco::PFMET FillMetVector(const edm::Event&);
	//virtual float GetMetTriggerEfficiency(reco::PFMET&);
        //virtual void ApplyRecoilCorrections(reco::PFMET&, const reco::Candidate::LorentzVector*, const reco::Candidate::LorentzVector*, int);
        //virtual bool isLooseJet(reco::CaloJet&);
        //virtual bool isTightJet(reco::CaloJet&);
        //virtual bool isTightLepVetoJet(reco::CaloJet&);
        //virtual std::vector<float> ReshapeBtagDiscriminator(reco::CaloJet&);
      
    private:

	edm::EDGetTokenT<std::vector<reco::CaloJet> > CaloJetToken;
        //edm::EDGetTokenT<std::vector<reco::PFMET> > CaloMetToken;
        //edm::EDGetTokenT<edm::ValueMap<float>> QGToken;
        //int JetId;
        float Jet1Pt, Jet2Pt, JetEta;
        //bool AddQG;
        bool RecalibrateJets, RecalibrateMass;
	//bool RecalibratePuppiMass;
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
