#include "CSCAnalyzer.h"


CSCAnalyzer::CSCAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
CSCSegmentToken(CColl.consumes<CSCSegmentCollection>(PSet.getParameter<edm::InputTag>("cscsegments")))
{
    
    std::cout << " --- CSCAnalyzer initialization ---" << std::endl;
    // std::cout << "  sample            :\t" << Sample << std::endl;
    // if(ApplyEWK) std::cout << "  EWK file          :\t" << EWKFileName << std::endl;
    std::cout << std::endl;
}

CSCAnalyzer::~CSCAnalyzer() {

}



std::vector<CSCSegment> CSCAnalyzer::FillCSCSegmentVector(const edm::Event& iEvent) {
    
    std::vector<CSCSegment> Vect;
    
    // Declare and open collections
    edm::Handle<CSCSegmentCollection> CSCCollection;
    iEvent.getByToken(CSCSegmentToken, CSCCollection); 
    
    // Iterate over CSC4DCollection and save CSC segments in vect     
    for (CSCSegmentCollection::const_iterator CSCSegment_ = CSCCollection->begin(); CSCSegment_ != CSCCollection->end();CSCSegment_++) {
        CSCSegment segment = *CSCSegment_;
        Vect.push_back(segment);
        
    }
    
    
    return Vect;
}


std::vector<GlobalPoint> CSCAnalyzer::FillGlobalPointCSCSegmentVector(const edm::Event& iEvent,const edm::EventSetup& iSetup,std::vector<CSCSegment>& Segment){

    std::vector<GlobalPoint> Vect;
    
    edm::Handle<CSCSegmentCollection> CSCCollection;
    iEvent.getByToken(CSCSegmentToken, CSCCollection); 
    
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
    
    for (CSCSegmentCollection::const_iterator CSCSegment = CSCCollection->begin(); CSCSegment != CSCCollection->end();CSCSegment++) {
        const GeomDet* geomDet = theTrackingGeometry->idToDet(CSCSegment->geographicalId());
        GlobalPoint dtpoint = geomDet->toGlobal(CSCSegment->localPosition());
        Vect.push_back(dtpoint);
        
    }
    return Vect;
    
}

std::map<std::string,float> CSCAnalyzer::GenMatcherCSCSegments(std::vector<GlobalPoint>& GlobalPoint, std::vector<reco::GenParticle>& Quarks, std::string label) {
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



