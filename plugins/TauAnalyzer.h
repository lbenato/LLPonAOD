#ifndef TAUANALYZER_H
#define TAUANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"


#include "TFile.h"
#include "TH2.h"

class TauAnalyzer {
    public:
        TauAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~TauAnalyzer();
        virtual std::vector<pat::Tau> FillTauVector(const edm::Event&);
        virtual void CleanTausFromMuons(std::vector<pat::Tau>&, std::vector<pat::Muon>&, float);
        virtual void CleanTausFromElectrons(std::vector<pat::Tau>&, std::vector<pat::Electron>&, float);
        virtual void CleanTausFromRecoElectrons(std::vector<pat::Tau>&, std::vector<reco::GsfElectron>&, float);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Tau> > TauToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
        float TauPt, TauEta;
	      int TauIdByDecayMode, TauIdByDeltaBetaIso, TauIdByMVAIso, TauIdByMuonRejection, TauIdByElectronRejection;
};


#endif
