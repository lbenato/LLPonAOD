#ifndef PILEUPANALYZER_H
#define PILEUPANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TFile.h"
#include "TH1.h"

class PileupAnalyzer {
    public:
        PileupAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
        ~PileupAnalyzer();
        virtual float GetPUWeight(const edm::Event&);
        virtual float GetPUWeightUp(const edm::Event&);
        virtual float GetPUWeightDown(const edm::Event&);
        virtual float GetPV(const edm::Event&);

      
    private:
        edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUToken;
        edm::EDGetTokenT<std::vector<reco::Vertex> > PVToken;
        edm::LumiReWeighting* LumiWeights;
        edm::LumiReWeighting* LumiWeightsUp;
        edm::LumiReWeighting* LumiWeightsDown;
        
        std::string DataFileName;
        std::string DataFileNameUp;
        std::string DataFileNameDown;
        std::string MCFileName;
        std::string DataName;
        std::string MCName;
};

#endif
