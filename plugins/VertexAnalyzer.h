#ifndef VERTEXANALYZER_H
#define VERTEXANALYZER_H

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
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "TFile.h"
#include "TH1.h"

class VertexAnalyzer {
 public:
  VertexAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
  ~VertexAnalyzer();
  virtual std::vector<reco::Vertex> FillPvVector(const edm::Event&);
  virtual std::vector<reco::VertexCompositePtrCandidate> FillSvVector(const edm::Event&);

      
 private:
  edm::EDGetTokenT<std::vector<reco::Vertex> > PVToken;
  edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > SVToken;

};

#endif
