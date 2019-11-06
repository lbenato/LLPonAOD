#ifndef MUONANALYZER_H
#define MUONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
//#include "rochcor2016.h"
//#include "RoccoR.h"

class MuonAnalyzer {
    public:
        MuonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~MuonAnalyzer();
        virtual std::vector<pat::Muon> FillMuonVector(const edm::Event&);
        virtual void AddVariables(std::vector<pat::Muon>&, pat::MET&);
        virtual bool IsTrackerHighPtMuon(pat::Muon&, const reco::Vertex*);
        virtual std::vector<float> FixTrackerIsolation(pat::Muon&, pat::Muon&);
        virtual std::string GetMuon1Id(pat::Muon&);
        virtual float GetMuonTrkSF(pat::Muon&);
        virtual float GetMuonTrkSFError(pat::Muon&);
        virtual float GetMuonIdSF(pat::Muon&, int);
        virtual float GetMuonIdSFError(pat::Muon&, int);
        virtual float GetMuonIsoSF(pat::Muon&, int);
        virtual float GetMuonIsoSFError(pat::Muon&, int);
        virtual float GetDoubleMuonTriggerSF(pat::Muon&, pat::Muon&);
        virtual float GetDoubleMuonTriggerSFError(pat::Muon&, pat::Muon&);
        virtual float GetMuonTriggerSFIsoMu24(pat::Muon&);
        virtual float GetMuonTriggerSFErrorIsoMu24(pat::Muon&);
        virtual float GetMuonTriggerSFMu50(pat::Muon&);
        virtual float GetMuonTriggerSFErrorMu50(pat::Muon&);
        virtual TH1F* ConvertTGraph(TGraphAsymmErrors*);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Muon> > MuonToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
        std::string MuonTrkFileName;
        std::string MuonIdFileName;
        std::string MuonIsoFileName;
        std::string MuonTrkHighptFileName;
        std::string MuonTriggerFileName;
        std::string DoubleMuonTriggerFileName;
        int Muon1Id, Muon2Id, Muon1Iso, Muon2Iso;
        float Muon1Pt, Muon2Pt;
        float MuonPtMax, MuonPtMin;
        bool UseTuneP, DoRochester;
        
        //rochcor2016 *rmcor;
        
        bool isMuonTriggerFile, isDoubleMuonTriggerFile, isMuonTrkFile, isMuonIdFile, isMuonIsoFile, isMuonTrkHighptFile;
        
        TFile* MuonTriggerFile;
        TFile* DoubleMuonTriggerFile;
        TFile* MuonTrkFile;
        TFile* MuonIdFile;
        TFile* MuonTrkHighptFile;
        TFile* MuonIsoFile;
        
        TH2F* MuonTriggerLt20;
        TH2F* MuonTriggerGt20;
        TH2F* MuonTriggerIsoMu24;
        TH2F* MuonTriggerMu50;
        
        TGraphAsymmErrors* MuonTrkGraph;
        TH1F* MuonTrk;
        TH2F* MuonIdLoose;
        TH2F* MuonIdMedium;
        TH2F* MuonIdTight;
        TH2F* MuonIdHighpt;
        TH2F* MuonIdTrkHighpt;
        TH2F* MuonIsoLoose;
        TH2F* MuonIsoTight;
        TH2F* MuonIsoHighpt;
}; 


#endif
