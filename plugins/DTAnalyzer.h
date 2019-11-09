#ifndef DTANALYZER_H
#define DTANALYZER_H

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

// DT Segment Collection
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecClusterCollection.h"

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





class DTAnalyzer {
    public:
        DTAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~DTAnalyzer();
        virtual std::vector<DTRecSegment4D> FillDTSegment4DVector(const edm::Event&);
        virtual std::map<std::string,float> GenMatcherDTSegments4D(std::vector<GlobalPoint>&,std::vector<reco::GenParticle>&, std::string);
        virtual std::vector<GlobalPoint> FillGlobalPointDT4DSegmentVector(const edm::Event&,const edm::EventSetup&,std::vector<DTRecSegment4D>&);
        
    private:
        edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken;
        
        
        
        
        
        
};

#endif
        
   