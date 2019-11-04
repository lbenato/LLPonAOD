#include "CounterAnalyzer.h"



CounterAnalyzer::CounterAnalyzer(const edm::ParameterSet& iConfig):
    LheToken(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheProduct"))),
    PythiaLOSample(iConfig.getParameter<bool>("pythiaLOSample"))
{
    //now do what ever initialization is needed
    usesResource("TFileService");
    Weight = fs->make<TH1F>("c_nEvents", "Event Counter", 1, 0., 1.); Weight->Sumw2();

    //Additional plots, not needed yet

    //NPartons = fs->make<TH1F>("c_lhePartons", "Event Counter", 5, 0., 5.); NPartons->Sumw2();
    //NBPartons = fs->make<TH1F>("c_lheBPartons", "Event Counter", 3, 0., 3.); NBPartons->Sumw2();
    //LheHT = fs->make<TH1F>("c_lheHT", "Event Counter", 400, 0., 4000.); LheHT->Sumw2();
    //LhePtZ = fs->make<TH1F>("c_lhePtZ", "Event Counter", 100, 0., 1000.); LhePtZ->Sumw2();
    
    //float binNp[6] = {0., 1., 2., 3., 4., 5.};
    //float binNb[4] = {0., 1., 2., 3.};
    //float binHT[9] = {0., 100., 200., 400., 600., 800., 1200., 2500., 4000.};
    //Bin = fs->make<TH3F>("c_bin", "Event Counter", 8, binHT, 5, binNp, 3, binNb); Bin->Sumw2();
    
    std::cout << " --- CounterAnalyzer initialization ---" << std::endl;
    if(PythiaLOSample) std::cout << "  Pythia LO sample" << std::endl;
    std::cout << std::endl;
}


CounterAnalyzer::~CounterAnalyzer() {
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void CounterAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) { 
    float weight(1.);
    
    //int lhePartons(0), lheBPartons(0);
    //float lheHT(0.), lhePtZ(0.), pt(0.);
    
    if(!iEvent.isRealData() && !PythiaLOSample) {
        // Declare and open collection
        edm::Handle<LHEEventProduct> LheEventCollection;
        iEvent.getByToken(LheToken, LheEventCollection);
        weight = LheEventCollection.product()->originalXWGTUP();
        weight = weight > 0. ? 1. : -1.;
        
	//
	//Additional plots, not needed yet
	//
        //const lhef::HEPEUP hepeup = LheEventCollection->hepeup();
        
        //for(int i = 0; i < hepeup.NUP; ++i) {
	    //int id=abs(hepeup.IDUP[i]);
            //// Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries in this event
            ////reco::Candidate::LorentzVector P4(hepeup.PUP[i][0], hepeup.PUP[i][1], hepeup.PUP[i][2], hepeup.PUP[i][3]);
            //pt = sqrt(hepeup.PUP[i][0]*hepeup.PUP[i][0] + hepeup.PUP[i][1]*hepeup.PUP[i][1]);
            //if(hepeup.ISTUP[i]==1 && (id<6 || id==21)) {
	        //lheHT += pt; //P4.pt() 
                //lhePartons++;
                //if(id==5) lheBPartons++;
	    //}
            //if(hepeup.ISTUP[i]==2 && (abs(hepeup.IDUP[i])==23 || abs(hepeup.IDUP[i])==24)) lhePtZ = pt;
	//}
    }
    Weight->Fill(0., weight);
    
    //Additional plots, not needed yet
    //NPartons->Fill(lhePartons, weight);
    //NBPartons->Fill(lheBPartons, weight);
    //LheHT->Fill(lheHT, weight);
    //LhePtZ->Fill(lhePtZ, weight);
    //Bin->Fill(lheHT, lhePartons, lheBPartons, weight);
}


// ------------ method called once each job just before starting event loop  ------------
void CounterAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void CounterAnalyzer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CounterAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CounterAnalyzer);
