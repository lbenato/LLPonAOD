// -*- C++ -*-
//
// Package:    Analyzer/Demo
// Class:      Demo
// 
/**\class Demo Demo.cc Analyzer/Demo/plugins/Demo.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lisa Benato
//         Created:  Mon, 26 Aug 2019 08:53:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Reco Jet classes
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <string>


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Demo : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Demo(const edm::ParameterSet&);
      ~Demo();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

    edm::EDGetTokenT<reco::PFJetCollection> jetToken;
    long int nJets;

    //Initialize tree                                                                                                                     
    //edm::Service<TFileService> fs;
    TTree* tree;

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
Demo::Demo(const edm::ParameterSet& iConfig)

{

   edm::InputTag IT_jets = edm::InputTag("ak4PFJetsCHS");
   jetToken = consumes<reco::PFJetCollection>(IT_jets);

   //now do what ever initialization is needed
   usesResource("TFileService");

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "tree");
   tree -> Branch("nJets" , &nJets , "nJets/L");

}


Demo::~Demo()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Demo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // Initialize 
   nJets = 0;

   edm::Handle<reco::PFJetCollection> jets;
   iEvent.getByToken( jetToken, jets );


   for(reco::PFJetCollection::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet) {
      if(jet->pt()>15 and abs(jet->eta())<2.4) nJets++;
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
Demo::beginJob()
{
  //Tree branches                                                                                                                       
  //tree = fs->make<TTree>("tree", "tree");
  //
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Demo::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Demo::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Demo);
