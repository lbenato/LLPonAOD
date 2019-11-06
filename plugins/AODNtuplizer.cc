// -*- C++ -*-
//
// Package:    Analyzer/AODNtuplizer
// Class:      AODNtuplizer
// 
/**\class AODNtuplizer AODNtuplizer.cc Analyzer/LLPonAOD/plugins/AODNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lisa Benato
//         Created:  Fri, 1 Nov 2019 09:48:39 GMT
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
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//Reco Jet classes
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//Pat classes
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <string>

#include "RecoJetAnalyzer.h"
#include "CaloJetAnalyzer.h"
#include "GenAnalyzer.h"
#include "PileupAnalyzer.h"
#include "TriggerAnalyzer.h"
#include "MuonAnalyzer.h"
#include "Objects.h"
#include "ObjectsFormat.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class AODNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AODNtuplizer(const edm::ParameterSet&);
      ~AODNtuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
    edm::ParameterSet GenPSet;
    edm::ParameterSet PileupPSet;
    edm::ParameterSet TriggerPSet;
    edm::ParameterSet CHSJetPSet;
    edm::ParameterSet CaloJetPSet;
    edm::ParameterSet VBFJetPSet;
    edm::ParameterSet MuonPSet;
    edm::EDGetTokenT<reco::PFJetCollection> jetToken;


    RecoJetAnalyzer* theCHSJetAnalyzer;
    CaloJetAnalyzer* theCaloJetAnalyzer;
    RecoJetAnalyzer* theVBFJetAnalyzer;
    GenAnalyzer* theGenAnalyzer;
    PileupAnalyzer* thePileupAnalyzer;
    TriggerAnalyzer* theTriggerAnalyzer;
    MuonAnalyzer* theMuonAnalyzer;

    double MinGenBpt, MaxGenBeta;
    double InvmassVBF, DetaVBF;//VBF tagging
    bool WriteGenVBFquarks, WriteGenHiggs, WriteGenBquarks, WriteGenLLPs;
    bool WriteOnlyTriggerEvents, WriteOnlyisVBFEvents;
    bool PerformPreFiringStudies;

    std::vector<JetType> CHSJets;
    //std::vector<JetType> ManualJets;
    std::vector<CaloJetType> CaloJets;
    //std::vector<LeptonType> Muons; //maybe later!
    std::vector<GenPType> GenVBFquarks;
    std::vector<GenPType> GenBquarks;
    std::vector<GenPType> GenLLPs;
    GenPType GenHiggs;

    MEtType MEt;
    CandidateType VBF;//VBF tagging

    std::map<std::string, bool> TriggerMap;
    //std::map<std::string, bool> MetFiltersMap;


    bool isMC;
    bool isVBF;
    long int EventNumber, LumiNumber, RunNumber, nPV, nSV;
    bool AtLeastOneTrigger;
    float EventWeight;
    float PUWeight, PUWeightUp, PUWeightDown;
    long int nJets;
    long int nCaloJets;
    long int nElectrons, nMuons, nTaus, nPhotons;
    long int nTightMuons;
    AddFourMomenta addP4;
    float HT;
    float MinJetMetDPhi;
    float MinJetMetDPhiAllJets;
    float m_pi, gen_b_radius;
    //MET filters
    bool BadPFMuonFlag, BadChCandFlag;
    //Pre-firing
    bool Prefired;
    long int nGenBquarks, nGenLL;
    //Initialize tree                                                                                                                     
    edm::Service<TFileService> fs;
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
AODNtuplizer::AODNtuplizer(const edm::ParameterSet& iConfig):
   GenPSet(iConfig.getParameter<edm::ParameterSet>("genSet")),
   PileupPSet(iConfig.getParameter<edm::ParameterSet>("pileupSet")),
   TriggerPSet(iConfig.getParameter<edm::ParameterSet>("triggerSet")),
   CHSJetPSet(iConfig.getParameter<edm::ParameterSet>("chsJetSet")),
   CaloJetPSet(iConfig.getParameter<edm::ParameterSet>("caloJetSet")),
   VBFJetPSet(iConfig.getParameter<edm::ParameterSet>("vbfJetSet")),
   MuonPSet(iConfig.getParameter<edm::ParameterSet>("muonSet")),
   MinGenBpt(iConfig.getParameter<double>("minGenBpt")),
   MaxGenBeta(iConfig.getParameter<double>("maxGenBeta")),
   InvmassVBF(iConfig.getParameter<double>("invmassVBF")),
   DetaVBF(iConfig.getParameter<double>("detaVBF")),
   WriteGenVBFquarks(iConfig.getParameter<bool>("writeGenVBFquarks")),
   WriteGenHiggs(iConfig.getParameter<bool>("writeGenHiggs")),
   WriteGenBquarks(iConfig.getParameter<bool>("writeGenBquarks")),
   WriteGenLLPs(iConfig.getParameter<bool>("writeGenLLPs")),
   WriteOnlyTriggerEvents(iConfig.getParameter<bool>("writeOnlyTriggerEvents")),
   WriteOnlyisVBFEvents(iConfig.getParameter<bool>("writeOnlyisVBFEvents")),
   PerformPreFiringStudies(iConfig.getParameter<bool>("performPreFiringStudies"))

{

   theCHSJetAnalyzer      = new RecoJetAnalyzer(CHSJetPSet, consumesCollector());
   theCaloJetAnalyzer     = new CaloJetAnalyzer(CaloJetPSet, consumesCollector());
   theVBFJetAnalyzer      = new RecoJetAnalyzer(VBFJetPSet, consumesCollector());
   theGenAnalyzer         = new GenAnalyzer(GenPSet, consumesCollector());
   thePileupAnalyzer      = new PileupAnalyzer(PileupPSet, consumesCollector());
   theTriggerAnalyzer     = new TriggerAnalyzer(TriggerPSet, consumesCollector());
   theMuonAnalyzer        = new MuonAnalyzer(MuonPSet, consumesCollector());

   std::vector<std::string> TriggerList(TriggerPSet.getParameter<std::vector<std::string> >("paths"));
   for(unsigned int i = 0; i < TriggerList.size(); i++) TriggerMap[ TriggerList[i] ] = false;
   //for(unsigned int i = 0; i < TriggerList.size(); i++) PrescalesTriggerMap[ TriggerList[i] ] = -1;
   //std::vector<std::string> MetFiltersList(TriggerPSet.getParameter<std::vector<std::string> >("metpaths"));
   //for(unsigned int i = 0; i < MetFiltersList.size(); i++) MetFiltersMap[ MetFiltersList[i] ] = false;
   //std::vector<std::string> L1FiltersList(TriggerPSet.getParameter<std::vector<std::string> >("l1filters"));
   //for(unsigned int i = 0; i < L1FiltersList.size(); i++) L1FiltersMap[ L1FiltersList[i] ] = false;

   edm::InputTag IT_jets = edm::InputTag("ak4PFJetsCHS");
   jetToken = consumes<reco::PFJetCollection>(IT_jets);

   //now do what ever initialization is needed
   usesResource("TFileService");


}


AODNtuplizer::~AODNtuplizer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   delete theCHSJetAnalyzer;
   delete theCaloJetAnalyzer;
   delete theVBFJetAnalyzer;
   delete theGenAnalyzer;
   delete thePileupAnalyzer;
   delete theTriggerAnalyzer;
   delete theMuonAnalyzer;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AODNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   // Initialize types
   ObjectsFormat::ResetGenPType(GenHiggs);
   ObjectsFormat::ResetCandidateType(VBF);


   nJets = nCaloJets = 0;
   nElectrons = nMuons = nTaus = nPhotons = 999;
   nTightMuons = 0;
   isMC = false;
   isVBF = false;
   EventNumber = LumiNumber = RunNumber = nPV = 0;
   EventWeight = PUWeight = PUWeightDown = PUWeightUp = 1.;
   HT = 0.;
   nGenBquarks = nGenLL = 0;
   m_pi = 0.;
   gen_b_radius = -1.;
   Prefired = false;



   //Event info                                                                
   isMC = !iEvent.isRealData();
   EventNumber = iEvent.id().event();
   LumiNumber = iEvent.luminosityBlock();
   RunNumber = iEvent.id().run();

   //Not needed anymore
   edm::Handle<reco::PFJetCollection> JetColl;
   iEvent.getByToken( jetToken, JetColl );

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // Trigger and MET filters
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------

   //if(isVerbose) std::cout << "Trigger and met filters" << std::endl;
   theTriggerAnalyzer->FillTriggerMap(iEvent, TriggerMap);
   //theTriggerAnalyzer->FillMetFiltersMap(iEvent, MetFiltersMap);
   BadPFMuonFlag = theTriggerAnalyzer->GetBadPFMuonFlag(iEvent);
   BadChCandFlag = theTriggerAnalyzer->GetBadChCandFlag(iEvent);
   //theTriggerAnalyzer->FillL1FiltersMap(iEvent, L1FiltersMap);

   // 27 Sep 2018: saving only events that fired at least one trigger, to reduce output size
   for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++)
      {
	if(it->second)
	  {
	    ////std::cout << "Trigger fired!!! " << it->first << std::endl;
	    AtLeastOneTrigger = true;
	  }
	//if(AtLeastOneTrigger) break;// no need to go through everything; once one trigger is fired, event can be saved
      }

   ////if(!AtLeastOneTrigger && WriteOnlyTriggerEvents) std::cout << "This event can be rejected" << std::endl;
   if(!AtLeastOneTrigger && WriteOnlyTriggerEvents) return;

   //Pre-firing
   if(PerformPreFiringStudies)
     {
         Prefired = theTriggerAnalyzer->EvaluatePrefiring(iEvent);
     }

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // HT
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   HT = theCHSJetAnalyzer->CalculateHT(iEvent,3,15,3.);

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // Muons
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   //if(isVerbose) std::cout << "Muons" << std::endl;
   std::vector<pat::Muon> MuonVect = theMuonAnalyzer->FillMuonVector(iEvent);
   std::vector<pat::Muon> TightMuonVect;
   for(unsigned int a = 0; a<MuonVect.size(); a++)
      {
	if(MuonVect.at(a).hasUserInt("isTight") && MuonVect.at(a).userInt("isTight")>0 && MuonVect.at(a).hasUserFloat("pfIso04") && MuonVect.at(a).userFloat("pfIso04")<0.15)//tight iso for muons
	  {
	    TightMuonVect.push_back(MuonVect.at(a));
	    nTightMuons++;
	  }
      }
   nMuons = MuonVect.size();

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // Missing Energy
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   //if(isVerbose) std::cout << "MET" << std::endl;
   reco::PFMET MET = theCHSJetAnalyzer->FillMetVector(iEvent);


   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // Gen particles
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------

   GenVBFquarks.clear();
   GenLLPs.clear();
   GenBquarks.clear();

   std::vector<reco::GenParticle> GenVBFVect = theGenAnalyzer->FillVBFGenVector(iEvent);
   std::vector<reco::GenParticle> GenHiggsVect = theGenAnalyzer->FillGenVectorByIdAndStatus(iEvent,25,22);
   std::vector<reco::GenParticle> GenLongLivedVect = theGenAnalyzer->FillGenVectorByIdAndStatus(iEvent,9000006,22);
   std::vector<reco::GenParticle> GenBquarksVect;

   nGenLL = GenLongLivedVect.size();

   if(nGenLL>0)
      {
	GenBquarksVect = theGenAnalyzer->FillGenVectorByIdStatusAndMotherAndKin(iEvent,5,23,9000006,float(MinGenBpt),float(MaxGenBeta));
      }
   else
      {
	GenBquarksVect = theGenAnalyzer->FillGenVectorByIdAndStatusAndKin(iEvent,5,23,float(MinGenBpt),float(MaxGenBeta));
      }

   nGenBquarks = GenBquarksVect.size();

   for(unsigned int i = 0; i < GenVBFVect.size(); i++) GenVBFquarks.push_back( GenPType() );
   for(unsigned int i = 0; i < GenLongLivedVect.size(); i++) GenLLPs.push_back( GenPType() );
   for(unsigned int i = 0; i < GenBquarksVect.size(); i++) GenBquarks.push_back( GenPType() );

   if(nGenBquarks>0) gen_b_radius = GenBquarksVect.at(0).mother()? sqrt(pow(GenBquarksVect.at(0).vx() - GenBquarksVect.at(0).mother()->vx(),2) + pow(GenBquarksVect.at(0).vy() - GenBquarksVect.at(0).mother()->vy(),2) + pow(GenBquarksVect.at(0).vz() - GenBquarksVect.at(0).mother()->vz(),2)) : -1.;
   if(nGenLL>0) m_pi = GenLongLivedVect.at(0).mass();

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // Pu weight and number of vertices
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   //if(isVerbose) std::cout << "Pile-up" << std::endl;
   PUWeight     = thePileupAnalyzer->GetPUWeight(iEvent);//calculates pileup weights
   PUWeightUp   = thePileupAnalyzer->GetPUWeightUp(iEvent);//syst uncertainties due to pileup
   PUWeightDown = thePileupAnalyzer->GetPUWeightDown(iEvent);//syst uncertainties due to pileup
   nPV = thePileupAnalyzer->GetPV(iEvent);//calculates number of vertices

   EventWeight *= PUWeight;

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // VBF Jets
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------

   std::vector<reco::PFJet> VBFJetsVect = theVBFJetAnalyzer->FillJetVector(iEvent);
   reco::CompositeCandidate theVBF;
   std::vector<pat::Jet> VBFPairJetsVect;
   
   float delta_eta_reco (0.), curr_delta_eta_reco(0.) ;
   int j1(-1), j2(-1);
   float curr_mjj(0.);
   reco::CompositeCandidate current_VBF;

   if(VBFJetsVect.size()>=2) {
      for(unsigned int a = 0; a<VBFJetsVect.size(); a++) {
        //find the VBF pair
        for(unsigned int b = 1; b<VBFJetsVect.size(); b++) {
	  //if(b!=a and VBFJetsVect.at(a).pt()>=30 and VBFJetsVect.at(b).pt()>=30 and VBFJetsVect.at(a).userInt("isLoose")>0 and VBFJetsVect.at(b).userInt("isLoose")>0 and (VBFJetsVect.at(a).eta()*VBFJetsVect.at(b).eta())<0)//if looser thresholds applied; we will use that for JER-JEC effects
	  if(b!=a and (VBFJetsVect.at(a).eta()*VBFJetsVect.at(b).eta())<0)//currently we are fine with what is defined in vbfJetSet
            {
	      //calculate delta eta
              curr_delta_eta_reco = abs(VBFJetsVect.at(a).eta() - VBFJetsVect.at(b).eta());
	      current_VBF.clearDaughters();
	      current_VBF.addDaughter(VBFJetsVect.at(a));
	      current_VBF.addDaughter(VBFJetsVect.at(b));
	      addP4.set(current_VBF);
	      curr_mjj = current_VBF.mass();
              if(curr_delta_eta_reco>delta_eta_reco and curr_delta_eta_reco>DetaVBF and curr_mjj>InvmassVBF)
                {
                  delta_eta_reco = curr_delta_eta_reco;
                  j1=std::min(a,b);
                  j2=std::max(a,b);
                }
            }
        }
      }

      if(j1>=0 && j2>=0)//otherwise, if indices are -1, theVBF seg faults
	{
	  theVBF.addDaughter(VBFJetsVect.at(j1));
	  theVBF.addDaughter(VBFJetsVect.at(j2));
	  addP4.set(theVBF);
	  VBFPairJetsVect.push_back(VBFJetsVect.at(j1));
	  VBFPairJetsVect.push_back(VBFJetsVect.at(j2));
	}
    }

    if(theVBF.pt()>0 && theVBF.mass()>InvmassVBF && abs(theVBF.daughter(1)->eta() - theVBF.daughter(0)->eta())>DetaVBF) {isVBF = true;}
    if(WriteOnlyisVBFEvents)//set in cfg file
      {
	if(!isVBF) return;
      }

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // AK4 CHS jets
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   std::vector<reco::PFJet> CHSJetsVect = theCHSJetAnalyzer->FillJetVector(iEvent);

   //Filling Jet structure manually, without filling a vector first. Used as cross-check.
   //for(reco::PFJetCollection::const_iterator it=JetColl->begin(); it!=JetColl->end(); ++it) {
   //   if(it->pt()>15 and abs(it->eta())<2.4) 
   //	{
   //	  reco::Jet jet=*it;
   //     ManualJets.push_back( JetType() );
   //     ObjectsFormat::ResetJetType(ManualJets[nJets]);
   //     ObjectsFormat::FillJetType(ManualJets[nJets], &jet, isMC);
   //     nJets++;
   //   }
   //}

   //std::cout << " --------------- " << std::endl;
   //std::cout<<nJets<<std::endl;
   //std::cout << ManualJets.size() << std::endl;

   //Remove jets tagged as VBF from the list of potential signal
   for(unsigned int r = 0; r<CHSJetsVect.size(); r++)
      {
	for(unsigned int s = 0; s<VBFPairJetsVect.size(); s++)
	  {
	    if(VBFPairJetsVect[s].pt()==CHSJetsVect[r].pt() && isVBF) //if jets aren't tagged as VBF jets, don't remove them
	      {
		CHSJetsVect.erase(CHSJetsVect.begin()+r);
	      }
	  }//VBF jet pair removed
    }

   //Gen matching: to be performed!

   nJets = CHSJetsVect.size();

   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   // AK4 Calo Jets
   //------------------------------------------------------------------------------------------
   //------------------------------------------------------------------------------------------
   std::vector<reco::CaloJet> CaloJetsVect = theCaloJetAnalyzer->FillJetVector(iEvent);
   nCaloJets = CaloJetsVect.size();
   // for gen matching, to be filled later
   std::vector<bool> caloGenMatched;
   for(unsigned int i = 0; i < CaloJetsVect.size(); i++) caloGenMatched.push_back(false);//to be implemented later


   // Fill objects
   if (WriteGenVBFquarks) for(unsigned int i = 0; i < GenVBFVect.size(); i++) ObjectsFormat::FillGenPType(GenVBFquarks[i], &GenVBFVect[i]);
   if (WriteGenHiggs) for(unsigned int i = 0; i < GenHiggsVect.size(); i++) ObjectsFormat::FillGenPType(GenHiggs, &GenHiggsVect[i]);
   if (WriteGenLLPs) for(unsigned int i = 0; i < GenLongLivedVect.size(); i++) ObjectsFormat::FillGenPType(GenLLPs[i], &GenLongLivedVect[i]);
   if (WriteGenBquarks) for(unsigned int i = 0; i < GenBquarksVect.size(); i++) ObjectsFormat::FillGenPType(GenBquarks[i], &GenBquarksVect[i]);
   ObjectsFormat::FillMEtType(MEt, &MET, isMC);
   ObjectsFormat::FillCandidateType(VBF, &theVBF, isMC);

   for(unsigned int i = 0; i < CHSJetsVect.size(); i++) CHSJets.push_back( JetType() );
   for(unsigned int i = 0; i < CHSJetsVect.size(); i++){
     ObjectsFormat::FillJetType(CHSJets[i], &CHSJetsVect[i], isMC);
   }
   for(unsigned int i = 0; i < CaloJetsVect.size(); i++) CaloJets.push_back( CaloJetType() );
   for(unsigned int i = 0; i < CaloJetsVect.size(); i++){ ObjectsFormat::FillCaloJetType(CaloJets[i], &CaloJetsVect[i], isMC, caloGenMatched[i]);}
   

   //Fill tree
   tree -> Fill();
   //ManualJets.clear();
   CHSJets.clear();
   CaloJets.clear();


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
AODNtuplizer::beginJob()
{
  //Tree branches                                                                                                                       
  //tree = fs->make<TTree>("tree", "tree");
  //

   tree = fs->make<TTree>("tree", "tree");
   tree -> Branch("isMC" , &isMC, "isMC/O");
   tree -> Branch("EventNumber" , &EventNumber , "EventNumber/L");
   tree -> Branch("LumiNumber" , &LumiNumber , "LumiNumber/L");
   tree -> Branch("RunNumber" , &RunNumber , "RunNumber/L");
   tree -> Branch("EventWeight", &EventWeight, "EventWeight/F");
   tree -> Branch("PUWeight", &PUWeight, "PUWeight/F");
   tree -> Branch("PUWeightUp", &PUWeightUp, "PUWeightUp/F");
   tree -> Branch("PUWeightDown", &PUWeightDown, "PUWeightDown/F");
   tree -> Branch("AtLeastOneTrigger" , &AtLeastOneTrigger , "AtLeastOneTrigger/O");
   tree -> Branch("Prefired" , &Prefired , "Prefired/O");
   tree -> Branch("nPV" , &nPV , "nPV/L");
   tree -> Branch("isVBF" , &isVBF, "isVBF/O");
   tree -> Branch("HT" , &HT , "HT/F");
   tree -> Branch("nGenBquarks" , &nGenBquarks , "nGenBquarks/L");
   tree -> Branch("nGenLL" , &nGenLL , "nGenLL/L");
   tree -> Branch("gen_b_radius" , &gen_b_radius , "gen_b_radius/F");
   tree -> Branch("m_pi" , &m_pi , "m_pi/F");
   tree -> Branch("nMuons", &nMuons, "nMuons/I");
   tree -> Branch("nTightMuons", &nTightMuons, "nTightMuons/I");
   tree -> Branch("Flag_BadPFMuon", &BadPFMuonFlag, "Flag_BadPFMuon/O");
   tree -> Branch("Flag_BadChCand", &BadChCandFlag, "Flag_BadChCand/O");
   tree -> Branch("nJets" , &nJets , "nJets/L");
   tree -> Branch("nCaloJets" , &nCaloJets , "nCaloJets/L");
   tree -> Branch("Flag_BadPFMuon", &BadPFMuonFlag, "Flag_BadPFMuon/O");
   tree -> Branch("Flag_BadChCand", &BadChCandFlag, "Flag_BadChCand/O");
   // Set trigger branches
   for(auto it = TriggerMap.begin(); it != TriggerMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());
   //for(auto it = MetFiltersMap.begin(); it != MetFiltersMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());
   //for(auto it = L1FiltersMap.begin(); it != L1FiltersMap.end(); it++) tree->Branch(it->first.c_str(), &(it->second), (it->first+"/O").c_str());

   //tree -> Branch("ManualJets", &ManualJets);
   tree -> Branch("GenHiggs", &GenHiggs.pt, ObjectsFormat::ListGenPType().c_str());
   tree -> Branch("GenLLPs", &GenLLPs);
   tree -> Branch("GenBquarks", &GenBquarks);
   tree -> Branch("MEt", &MEt.pt, ObjectsFormat::ListMEtType().c_str());
   tree -> Branch("CHSJets", &CHSJets);
   tree -> Branch("CaloJets", &CaloJets);
   tree -> Branch("VBFPair", &VBF.pt, ObjectsFormat::ListCandidateType().c_str());
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AODNtuplizer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AODNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODNtuplizer);
