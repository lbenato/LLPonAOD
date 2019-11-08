#include "PFCandidateAnalyzer.h"


PFCandidateAnalyzer::PFCandidateAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
  PFCandidateToken(CColl.consumes<std::vector<pat::PackedCandidate> >(PSet.getParameter<edm::InputTag>("pfCandidates"))),
  LostTrackToken(CColl.consumes<std::vector<pat::PackedCandidate> >(PSet.getParameter<edm::InputTag>("lostTracks")))
{   
  std::cout << " --- PFCandidateAnalyzer initialization ---" << std::endl;
  std::cout << std::endl;
}

PFCandidateAnalyzer::~PFCandidateAnalyzer() {

}


// ---------- PFCandidates and LostTracks ----------

std::vector<pat::PackedCandidate> PFCandidateAnalyzer::FillPFCandidateVector(const edm::Event& iEvent) {
  edm::Handle<std::vector<pat::PackedCandidate>> PFCandidates;
  iEvent.getByToken(PFCandidateToken, PFCandidates);
  return *PFCandidates;
}


std::vector<pat::PackedCandidate> PFCandidateAnalyzer::FillLostTrackVector(const edm::Event& iEvent) {
  edm::Handle<std::vector<pat::PackedCandidate>> LostTracks;
  iEvent.getByToken(LostTrackToken,LostTracks);
  return *LostTracks;
}
