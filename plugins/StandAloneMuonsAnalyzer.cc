#include "StandAloneMuonsAnalyzer.h"


StandAloneMuonsAnalyzer::StandAloneMuonsAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
  StandAloneMuonsToken(CColl.consumes<reco::TrackCollection>(PSet.getParameter<edm::InputTag>("standaloneMuons")))
{
    
    std::cout << " --- StandAloneMuonsAnalyzer initialization ---" << std::endl;

    // std::cout << "  sample            :\t" << Sample << std::endl;
    // if(ApplyEWK) std::cout << "  EWK file          :\t" << EWKFileName << std::endl;
    std::cout << std::endl;
}

StandAloneMuonsAnalyzer::~StandAloneMuonsAnalyzer() {

}



std::vector<reco::Track> StandAloneMuonsAnalyzer::FillStandAloneMuonsVector(const edm::Event& iEvent) {
    
    std::vector<reco::Track> Vect;
    
    // Declare and open collections
    edm::Handle<reco::TrackCollection> StandAloneMuonsCollection;
    iEvent.getByToken(StandAloneMuonsToken, StandAloneMuonsCollection); 
    
    // Iterate over StandAloneMuonsCollection and save DT segments in vect     
    for (reco::TrackCollection::const_iterator track_it = StandAloneMuonsCollection->begin(); track_it != StandAloneMuonsCollection->end();track_it++) {
        reco::Track standalonemuon = *track_it;
        Vect.push_back(standalonemuon);
        
    }
    
    
    return Vect;
}


std::map<std::string,float> StandAloneMuonsAnalyzer::GenMatcherStandAloneMuons(std::vector<reco::Track>& StandAloneMuon, std::vector<reco::GenParticle>& Quarks, std::string label) {
    std::map<std::string,float> match_map;
    
     for(unsigned int j = 0; j < StandAloneMuon.size(); j++) {
        for(unsigned int q = 0; q < Quarks.size(); q++) {
	  match_map.insert(std::make_pair(("dR_j"+std::to_string(j+1)+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(StandAloneMuon[j].eta(),StandAloneMuon[j].phi(),Quarks[q].eta(),Quarks[q].phi())) ));
        }
    }
    return match_map;
    
}



