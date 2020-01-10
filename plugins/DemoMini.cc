// -*- C++ -*-
//
// Package:    Analyzer/DemoMini
// Class:      DemoMini
// 
/**\class DemoMini DemoMini.cc Analyzer/DemoMini/plugins/DemoMini.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lisa Benato
//         Created:  Mon, 26 Aug 2019 08:53:05 GMT
//
//


//Gen info stuff from: https://github.com/lviliani/LHEWeightsReader/blob/master/plugins/LHEWeightsDumper.cc

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Pat Jet classes
#include "DataFormats/PatCandidates/interface/Jet.h"
//Reco Jet classes
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TH1.h"
#include "TH3.h"
#include <string>


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DemoMini : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DemoMini(const edm::ParameterSet&);
      ~DemoMini();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

    edm::EDGetTokenT<pat::JetCollection> jetToken;
    edm::EDGetTokenT<LHEEventProduct> lheToken_;
    edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;

    long int nJets;
    float weight, pos_weight, neg_weight, weight_norm;
    float lhe_weight, lhe_weight_norm;
    float j0_pt;
    TH1F* Weight;
    TH1F* Weight_norm;
    TH1F* lhe_Weight;
    TH1F* lhe_Weight_norm;
    TH1F* jet_pt_unweighted;
    TH1F* jet_pt_weighted;
    TH1F* jet_pt_weighted_norm;
    TH1F* jet_pt_pos_weighted;
    TH1F* jet_pt_neg_weighted;

    //Initialize tree                                                                                                                     
    //edm::Service<TFileService> fs;
    TTree* tree;

    bool PythiaLOSample;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoMini::DemoMini(const edm::ParameterSet& iConfig):
    PythiaLOSample(iConfig.getParameter<bool>("pythiaLOSample"))
{
   if(PythiaLOSample) std::cout << "  Pythia LO sample" << std::endl;

   edm::InputTag IT_jets = edm::InputTag("slimmedJets");
   jetToken = consumes<pat::JetCollection>(IT_jets);

   if(!PythiaLOSample) lheToken_ = consumes <LHEEventProduct,edm::InEvent> (edm::InputTag("externalLHEProducer"));//(lheProduct_);
   genEventToken_ = consumes <GenEventInfoProduct,edm::InEvent> (edm::InputTag("generator"));//(genEventProduct_);

   //now do what ever initialization is needed
   usesResource("TFileService");

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "tree");
   tree -> Branch("nJets" , &nJets , "nJets/L");
   tree -> Branch("j0_pt" , &j0_pt , "j0_pt/F");
   tree -> Branch("weight" , &weight , "weight/F");
   tree -> Branch("weight_norm" , &weight_norm , "weight_norm/F");
   tree -> Branch("lhe_weight" , &lhe_weight , "lhe_weight/F");
   tree -> Branch("lhe_weight_norm" , &lhe_weight_norm , "lhe_weight_norm/F");
   tree -> Branch("pos_weight" , &pos_weight , "pos_weight/F");
   tree -> Branch("neg_weight" , &neg_weight , "neg_weight/F");

   Weight = fs->make<TH1F>("nEvents", "Event Counter", 1, 0., 1.);
   Weight->Sumw2();
   Weight_norm = fs->make<TH1F>("nEvents_norm", "Event Counter", 1, 0., 1.);
   Weight_norm->Sumw2();

   lhe_Weight = fs->make<TH1F>("lhe_nEvents", "Event Counter", 1, 0., 1.);
   lhe_Weight->Sumw2();
   lhe_Weight_norm = fs->make<TH1F>("lhe_nEvents_norm", "Event Counter", 1, 0., 1.);
   lhe_Weight_norm->Sumw2();

   jet_pt_unweighted = fs->make<TH1F>("jet_pt_unweighted", "jet_pt_unweighted", 50, 15., 215.);
   jet_pt_unweighted->Sumw2();
   jet_pt_weighted = fs->make<TH1F>("jet_pt_weighted", "jet_pt_weighted", 50, 15., 215.);
   jet_pt_weighted->Sumw2();

   jet_pt_weighted_norm = fs->make<TH1F>("jet_pt_weighted_norm", "jet_pt_weighted_norm", 50, 15., 215.);
   jet_pt_weighted_norm->Sumw2();

   jet_pt_pos_weighted = fs->make<TH1F>("jet_pt_pos_weighted", "jet_pt_pos_weighted", 50, 15., 215.);
   jet_pt_pos_weighted->Sumw2();
   jet_pt_neg_weighted = fs->make<TH1F>("jet_pt_neg_weighted", "jet_pt_neg_weighted", 50, 15., 215.);
   jet_pt_neg_weighted->Sumw2();

}


DemoMini::~DemoMini()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoMini::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Initialize 
   nJets = 0;
   weight = 0.;
   weight_norm = 0.;
   lhe_weight = 0.;
   lhe_weight_norm = 0.;
   pos_weight = 0.;
   neg_weight = 0.;
   j0_pt = -1;



   edm::Handle<GenEventInfoProduct> genEventInfoProduct;
   iEvent.getByToken(genEventToken_, genEventInfoProduct);
   const GenEventInfoProduct& genEventInfo = *(genEventInfoProduct.product());

   weight = genEventInfo.weights()[0];
   weight_norm = weight > 0. ? 1. : -1.;
   Weight->Fill(0., weight);
   Weight_norm->Fill(0., weight_norm);
   if (weight > 0.) {pos_weight = pos_weight + fabs(weight);}
   else {neg_weight = neg_weight + fabs(weight);}

   //std::cout << "originalXWGTUP: " <<  lhe_event->originalXWGTUP() << std::endl;
   //for (unsigned int i=0; i<lhe_event->weights().size(); i++) {
   //  std::cout << "Weight " << i << " ID: " << lhe_event->weights()[i].id << " " <<  lhe_event->weights()[i].wgt << std::endl;
   //}
   //std::cout << "as per twiki: " << lhe_event->weights() << std::endl;

   if(!PythiaLOSample)
   {
     edm::Handle<LHEEventProduct> lhe_event;  
     iEvent.getByToken(lheToken_,lhe_event);
     lhe_weight = lhe_event->originalXWGTUP();
   }

   //weight = weight/fabs(weight);
   lhe_weight_norm = lhe_weight > 0. ? 1. : -1.;
   lhe_Weight->Fill(0., lhe_weight);
   lhe_Weight_norm->Fill(0., lhe_weight_norm);

   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken( jetToken, jets );


   for(pat::JetCollection::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet) {
      if(jet->pt()>15 and abs(jet->eta())<2.4)
        {
          nJets++;
          if(nJets==1) j0_pt = jet->pt();
          jet_pt_unweighted->Fill(jet->pt());
          jet_pt_weighted->Fill(jet->pt(),weight);
          jet_pt_weighted_norm->Fill(jet->pt(),weight_norm);
          if(weight>0) jet_pt_pos_weighted->Fill(jet->pt(),fabs(weight));
          if(weight<0) jet_pt_neg_weighted->Fill(jet->pt(),fabs(weight));
        }
   }

   //std::cout << nJets << std::endl;

   tree -> Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoMini::beginJob()
{
  //Tree branches                                                                                                                       
  //tree = fs->make<TTree>("tree", "tree");
  //
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoMini::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoMini::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoMini);
