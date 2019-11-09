#ifndef CSCANALYZER_H
#define CSCANALYZER_H

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

// CSC Segment Collection
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

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





class CSCAnalyzer {
    public:
        CSCAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~CSCAnalyzer();
        virtual std::vector<CSCSegment> FillCSCSegmentVector(const edm::Event&);
        virtual std::map<std::string,float> GenMatcherCSCSegments(std::vector<GlobalPoint>&,std::vector<reco::GenParticle>&, std::string);
        virtual std::vector<GlobalPoint> FillGlobalPointCSCSegmentVector(const edm::Event&,const edm::EventSetup&,std::vector<CSCSegment>&);
        
    private:
        edm::EDGetTokenT<CSCSegmentCollection> CSCSegmentToken;
        
        
        
        
        
        
};

#endif
        
   
