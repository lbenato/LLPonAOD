#include "DTAnalyzer.h"


DTAnalyzer::DTAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
dtSegmentToken(CColl.consumes<DTRecSegment4DCollection>(PSet.getParameter<edm::InputTag>("dtsegments")))
{
    
    std::cout << " --- DTAnalyzer initialization ---" << std::endl;

    // std::cout << "  sample            :\t" << Sample << std::endl;
    // if(ApplyEWK) std::cout << "  EWK file          :\t" << EWKFileName << std::endl;
    std::cout << std::endl;
}

DTAnalyzer::~DTAnalyzer() {

}



std::vector<DTRecSegment4D> DTAnalyzer::FillDTSegment4DVector(const edm::Event& iEvent) {
    
    std::vector<DTRecSegment4D> Vect;
    
    // Declare and open collections
    edm::Handle<DTRecSegment4DCollection> DT4DCollection;
    iEvent.getByToken(dtSegmentToken, DT4DCollection); 
    
    // Iterate over DT4DCollection and save DT segments in vect     
    for (DTRecSegment4DCollection::const_iterator DTSegment = DT4DCollection->begin(); DTSegment != DT4DCollection->end();DTSegment++) {
        DTRecSegment4D segment = *DTSegment;
        Vect.push_back(segment);
        
    }
    
    
    return Vect;
}


std::vector<GlobalPoint> DTAnalyzer::FillGlobalPointDT4DSegmentVector(const edm::Event& iEvent,const edm::EventSetup& iSetup,std::vector<DTRecSegment4D>& Segment){

    std::vector<GlobalPoint> Vect;
    
    edm::Handle<DTRecSegment4DCollection> DT4DCollection;
    iEvent.getByToken(dtSegmentToken, DT4DCollection); 
    
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
    
    for (DTRecSegment4DCollection::const_iterator DTSegment = DT4DCollection->begin(); DTSegment != DT4DCollection->end();DTSegment++) {
        const GeomDet* geomDet = theTrackingGeometry->idToDet(DTSegment->geographicalId());
        GlobalPoint dtpoint = geomDet->toGlobal(DTSegment->localPosition());
        Vect.push_back(dtpoint);
        
    }
    return Vect;
    
}

std::map<std::string,float> DTAnalyzer::GenMatcherDTSegments4D(std::vector<GlobalPoint>& GlobalPoint, std::vector<reco::GenParticle>& Quarks, std::string label) {
    std::map<std::string,float> match_map;
    
     for(unsigned int j = 0; j < GlobalPoint.size(); j++) {
        for(unsigned int q = 0; q < Quarks.size(); q++) {
	  //std::cout << "jet: " << j << " quark: " << q << " delta R: " <<fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi()) ) << std::endl;
	  //std::cout << ("dR_q"+std::to_string(q)).c_str() << std::endl;
	  //std::cout <<  ("dR_"+label+std::to_string(q+1)).c_str() << std::endl;
	  match_map.insert(std::make_pair(("dR_j"+std::to_string(j+1)+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(GlobalPoint[j].eta(),GlobalPoint[j].phi(),Quarks[q].eta(),Quarks[q].phi())) ));
	  //Jets[j].addUserFloat(("dR_"+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi())) );
            //Jets[j].addUserFloat("quark_index", fabs(reco::deltaPhi(Jets[j].phi(), Jets[0].phi())));
        }
    }
    return match_map;
    
}



