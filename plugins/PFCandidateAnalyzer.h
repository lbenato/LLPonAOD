#ifndef PFCANDIDATEANALYZER_H
#define PFCANDIDATEANALYZER_H

#include <iostream>
#include <cmath>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TFile.h"
#include "TH1.h"

class PFCandidateAnalyzer {
 public:
  PFCandidateAnalyzer(edm::ParameterSet&, edm::ConsumesCollector&&);
  ~PFCandidateAnalyzer();
  virtual std::vector<pat::PackedCandidate> FillPFCandidateVector(const edm::Event&);
  virtual std::vector<pat::PackedCandidate> FillLostTrackVector(const edm::Event&);

      
 private:
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > PFCandidateToken;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > LostTrackToken;

};

#endif
