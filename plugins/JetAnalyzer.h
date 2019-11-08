#ifndef JETANALYZER_H
#define JETANALYZER_H

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
#include "RecoilCorrector.h" // From: https://github.com/cms-met/MetTools/tree/master/RecoilCorrections
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//Reco met
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


class JetAnalyzer {
    public:
        JetAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~JetAnalyzer();
        virtual std::vector<pat::Jet> FillJetVector(const edm::Event&);
        virtual void CorrectJet(pat::Jet&, float, float, bool);
        virtual void CorrectMass(pat::Jet&, float, float, bool);
        virtual void CorrectPuppiMass(pat::Jet&, bool);
        virtual void CleanJetsFromMuons(std::vector<pat::Jet>&, std::vector<pat::Muon>&, float);
        virtual void CleanJetsFromElectrons(std::vector<pat::Jet>&, std::vector<pat::Electron>&, float);
        virtual void CleanFatJetsFromAK4(std::vector<pat::Jet>&, std::vector<pat::Jet>&, float);
        virtual void AddVariables(std::vector<pat::Jet>&, pat::MET&);
        virtual void GenMatcher(std::vector<pat::Jet>&, std::vector<reco::GenParticle>&, std::string);
        virtual int GetNBJets(std::vector<pat::Jet>&);
        virtual pat::MET FillMetVector(const edm::Event&);
        virtual reco::PFMET FillRecoMetVector(const edm::Event&);
	virtual float GetMetTriggerEfficiency(pat::MET&);
        virtual void ApplyRecoilCorrections(pat::MET&, const reco::Candidate::LorentzVector*, const reco::Candidate::LorentzVector*, int);
        virtual float CalculateHT(const edm::Event&, int, float, float);
        virtual bool isLooseJet(pat::Jet&);
        virtual bool isTightJet(pat::Jet&);
        virtual bool isTightLepVetoJet(pat::Jet&);
        virtual std::vector<float> ReshapeBtagDiscriminator(pat::Jet&);
      
    private:

	edm::EDGetTokenT<std::vector<pat::Jet> > JetToken;
        edm::EDGetTokenT<std::vector<pat::MET> > MetToken;
        edm::EDGetTokenT<std::vector<reco::PFMET> > RecoMetToken;
        edm::EDGetTokenT<edm::ValueMap<float>> QGToken;
        int JetId;
        float Jet1Pt, Jet2Pt, JetEta;
        bool AddQG, RecalibrateJets, RecalibrateMass, RecalibratePuppiMass; 
	std::string SoftdropPuppiMassString;
        bool SmearJets;
        std::string JECUncertaintyMC;
        std::string JECUncertaintyDATA;
        std::vector<std::string> JetCorrectorMC;
        std::vector<std::string> JetCorrectorDATA;
        std::vector<std::string> MassCorrectorMC;
        std::vector<std::string> MassCorrectorDATA;
        std::string MassCorrectorPuppi;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;        
        edm::EDGetTokenT<double> RhoToken;
        bool UseReshape;
        std::string BTag;
        int Jet1BTag, Jet2BTag;
        std::string BTagDB;
        bool UseRecoil;
        std::string RecoilMCFile;
        std::string RecoilDataFile;
        std::string MetTriggerFileName;
        std::string JerName_res;
        std::string JerName_sf;
	std::vector<std::string> BTagNames;
        float Rparameter;
        
        TFile* PuppiCorrFile;
        TF1* PuppiJECcorr_gen;
        TF1* PuppiJECcorr_reco_0eta1v3;
        TF1* PuppiJECcorr_reco_1v3eta2v5;

	TFile* MetTriggerFile;
	TH1F* MetTriggerHisto;
	bool isMetTriggerFile;

        // JEC Uncertainty
        JetCorrectionUncertainty* jecUncMC;
        JetCorrectionUncertainty* jecUncDATA;
        
        boost::shared_ptr<FactorizedJetCorrector> jetCorrMC;
        boost::shared_ptr<FactorizedJetCorrector> jetCorrDATA;
        boost::shared_ptr<FactorizedJetCorrector> massCorrMC;
        boost::shared_ptr<FactorizedJetCorrector> massCorrDATA;
        
        // Btag calibrations
        BTagCalibration       * calib;
	std::map < int , BTagEntry::JetFlavor > flavour_map; 
	std::map< BTagEntry::JetFlavor, std::vector<std::string>> syst_map; 
	std::map<std::string, BTagCalibrationReader> cr_map;
	std::string sf_mode;

	//	std::vector<edm::EDGetTokenT<edm::View<reco::BaseTagInfo> > > BTagInfos_;
        
        //JME
        JME::JetResolution              * resolution;
        JME::JetResolutionScaleFactor   * resolution_sf;        
        
        // Recoil corrections
        RecoilCorrector* recoilCorr;
};

#endif
