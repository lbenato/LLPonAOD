#ifndef STANDALONEMUONSANALYZER_H
#define STANDALONEMUONSANALYZER_H

#include <iostream>
#include <fstream>

#include <cmath>
#include <map>


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//Track header
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"

#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Objects.h"
#include "ObjectsFormat.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Riostream.h"





class StandAloneMuonsAnalyzer {
    public:
        StandAloneMuonsAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~StandAloneMuonsAnalyzer();
        virtual std::vector<reco::Track> FillStandAloneMuonsVector(const edm::Event&);
        virtual std::map<std::string,float> GenMatcherStandAloneMuons(std::vector<reco::Track>&,std::vector<reco::GenParticle>&, std::string);
        
    private:
        edm::EDGetTokenT<std::vector<reco::Track>> StandAloneMuonsToken;
        
        
        
        
        
        
};

#endif
        
   
