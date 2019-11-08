#include "VertexAnalyzer.h"


VertexAnalyzer::VertexAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
  PVToken(CColl.consumes<std::vector<reco::Vertex> >(PSet.getParameter<edm::InputTag>("primaryVertices"))),
  SVToken(CColl.consumes<std::vector<reco::VertexCompositePtrCandidate> >(PSet.getParameter<edm::InputTag>("secondaryVertices")))
{   
  std::cout << " --- VertexAnalyzer initialization ---" << std::endl;
  std::cout << std::endl;
}

VertexAnalyzer::~VertexAnalyzer() {

}


// ---------- Vertices ----------

std::vector<reco::Vertex> VertexAnalyzer::FillPvVector(const edm::Event& iEvent) {
  edm::Handle<reco::VertexCollection> PVCollection;
  iEvent.getByToken(PVToken, PVCollection);
  std::vector<reco::Vertex> Vect;
  for (std::vector<reco::Vertex>::const_iterator it = PVCollection->begin(); it != PVCollection->end(); ++it){
    Vect.push_back(*it);
  }
  return Vect;
}

std::vector<reco::VertexCompositePtrCandidate> VertexAnalyzer::FillSvVector(const edm::Event& iEvent) {
  edm::Handle<reco::VertexCompositePtrCandidateCollection> SVCollection;
  iEvent.getByToken(SVToken, SVCollection);
  return *SVCollection;
}
