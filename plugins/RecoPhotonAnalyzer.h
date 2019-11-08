#ifndef RECOPHOTONANALYZER_H
#define RECOPHOTONANALYZER_H

#include <iostream>
#include <fstream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "TFile.h"
#include "TH2.h"

class RecoPhotonAnalyzer {
    public:
        RecoPhotonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~RecoPhotonAnalyzer();
        virtual std::vector<reco::Photon> FillPhotonVector(const edm::Event&);
        virtual void PlotPhotons(std::vector<reco::Photon>&, std::map<std::string, TH1F*>&, float);
        virtual float GetPhotonIdSFLoose(reco::Photon&);
        virtual float GetPhotonIdSFLooseError(reco::Photon&);
        virtual float GetPhotonIdSFMedium(reco::Photon&);
        virtual float GetPhotonIdSFMediumError(reco::Photon&);
        virtual float GetPhotonIdSFTight(reco::Photon&);
        virtual float GetPhotonIdSFTightError(reco::Photon&);
        virtual float GetPhotonIdSFMVANonTrigMedium(reco::Photon&);
        virtual float GetPhotonIdSFMVANonTrigMediumError(reco::Photon&);
        //virtual bool isLoosePhoton(reco::Photon&, const reco::Vertex*);
        virtual void CleanPhotonsFromMuons(std::vector<reco::Photon>&, std::vector<pat::Muon>&, float);
        virtual void CleanPhotonsFromElectrons(std::vector<reco::Photon>&, std::vector<pat::Electron>&, float);
      
    private:
      
        edm::EDGetTokenT<std::vector<reco::Photon> > PhotonToken;
        edm::EDGetTokenT<reco::VertexCollection> VertexToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> PhoLooseIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> PhoMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> PhoTightIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool>> PhoMVANonTrigMediumIdMapToken;
        edm::EDGetTokenT<EcalRecHitCollection> PhoEcalRecHitCollectionToken;        
        std::string PhoLooseIdFileName;
        std::string PhoMediumIdFileName;
        std::string PhoTightIdFileName;
        std::string PhoMVANonTrigMediumIdFileName;
        int PhotonId;
        float PhotonPt;
        bool isPhoLooseIdFile, isPhoMediumIdFile, isPhoTightIdFile, isPhoMVANonTrigMediumIdFile;
        TFile* PhoLooseIdFile;
        TFile* PhoMediumIdFile;
        TFile* PhoTightIdFile;
        TFile* PhoMVANonTrigMediumIdFile;
        TH2F* PhotonIdLoose;
        TH2F* PhotonIdMedium;
        TH2F* PhotonIdTight;
        TH2F* PhotonIdMVANonTrigMedium;

};


#endif
