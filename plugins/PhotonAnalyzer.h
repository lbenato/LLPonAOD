#ifndef PHOTONANALYZER_H
#define PHOTONANALYZER_H

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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "TFile.h"
#include "TH2.h"

class PhotonAnalyzer {
    public:
        PhotonAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~PhotonAnalyzer();
        virtual std::vector<pat::Photon> FillPhotonVector(const edm::Event&);
        virtual void PlotPhotons(std::vector<pat::Photon>&, std::map<std::string, TH1F*>&, float);
        virtual float GetPhotonIdSFLoose(pat::Photon&);
        virtual float GetPhotonIdSFLooseError(pat::Photon&);
        virtual float GetPhotonIdSFMedium(pat::Photon&);
        virtual float GetPhotonIdSFMediumError(pat::Photon&);
        virtual float GetPhotonIdSFTight(pat::Photon&);
        virtual float GetPhotonIdSFTightError(pat::Photon&);
        virtual float GetPhotonIdSFMVANonTrigMedium(pat::Photon&);
        virtual float GetPhotonIdSFMVANonTrigMediumError(pat::Photon&);
        //virtual bool isLoosePhoton(pat::Photon&, const reco::Vertex*);
        virtual void CleanPhotonsFromMuons(std::vector<pat::Photon>&, std::vector<pat::Muon>&, float);
        virtual void CleanPhotonsFromElectrons(std::vector<pat::Photon>&, std::vector<pat::Electron>&, float);
      
    private:
      
        edm::EDGetTokenT<std::vector<pat::Photon> > PhotonToken;
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
