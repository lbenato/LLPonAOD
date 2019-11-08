#include "TriggerAnalyzer.h"


TriggerAnalyzer::TriggerAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TriggerToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("trigger"))),
    TriggerList(PSet.getParameter<std::vector<std::string> >("paths")),
    MetFiltersToken(CColl.consumes<edm::TriggerResults>(PSet.getParameter<edm::InputTag>("metfilters"))),
    MetFiltersList(PSet.getParameter<std::vector<std::string> >("metpaths")),
    PrescalesToken(CColl.consumes<pat::PackedTriggerPrescales>(PSet.getParameter<edm::InputTag>("prescales"))),
    L1MinPrescalesToken(CColl.consumes<pat::PackedTriggerPrescales>(PSet.getParameter<edm::InputTag>("l1Minprescales"))),
    L1MaxPrescalesToken(CColl.consumes<pat::PackedTriggerPrescales>(PSet.getParameter<edm::InputTag>("l1Maxprescales"))),
    TriggerObjectToken(CColl.consumes<std::vector<pat::TriggerObjectStandAlone> >(PSet.getParameter<edm::InputTag>("objects"))),
    BadPFMuonFilterToken(CColl.consumes<bool>(PSet.getParameter<edm::InputTag>("badPFMuonFilter"))),
    BadChCandFilterToken(CColl.consumes<bool>(PSet.getParameter<edm::InputTag>("badChCandFilter"))),
    L1GtToken(CColl.consumes<BXVector<GlobalAlgBlk>>(PSet.getParameter<edm::InputTag>("l1Gt"))),//Pre-Firing
    L1FiltersList(PSet.getParameter<std::vector<std::string> >("l1filters"))
{
    std::cout << " --- TriggerAnalyzer initialization ---" << std::endl;
    std::cout << "  HLT paths:" << std::endl;
    for(unsigned int i = 0; i < TriggerList.size(); i++) std::cout << "    " << TriggerList[i] << "*" << std::endl;
    std::cout << "  Met filters:" << std::endl;
    for(unsigned int i = 0; i < MetFiltersList.size(); i++) std::cout << "    " << MetFiltersList[i] << "*" << std::endl;
    std::cout << std::endl;
    std::cout << "  L1 filters:" << std::endl;
    for(unsigned int i = 0; i < L1FiltersList.size(); i++) std::cout << "    " << L1FiltersList[i] << "*" << std::endl;
    std::cout << std::endl;

    //l1GtUtils_ = new l1t::L1TGlobalUtil(PSet,CColl);


}

TriggerAnalyzer::~TriggerAnalyzer() {

}



// ---------- TRIGGER ----------

void TriggerAnalyzer::FillTriggerMap(const edm::Event& iEvent, std::map<std::string, bool>& Map, std::map<std::string, int>& PrescalesMap, bool& verbosity) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(PrescalesToken, triggerPrescales);

    edm::Handle<pat::PackedTriggerPrescales> l1MinTriggerPrescales;
    iEvent.getByToken(L1MinPrescalesToken, l1MinTriggerPrescales);

    edm::Handle<pat::PackedTriggerPrescales> l1MaxTriggerPrescales;
    iEvent.getByToken(L1MaxPrescalesToken, l1MaxTriggerPrescales);

    //edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjectCollection;
    //iEvent.getByToken(TriggerObjectToken, triggerObjectCollection);

    //for(unsigned int j=0, in=trigNames.size(); j < in; j++) std::cout << trigNames.triggerName(j) << std::endl;
    
    // Get Trigger index
    for(unsigned int i = 0; i < TriggerList.size(); i++) {
        Map[TriggerList[i]] = false;
        PrescalesMap[TriggerList[i]] = -1;
        for(unsigned int j=0, in=trigNames.size(); j < in; j++) {
            if(trigNames.triggerName(j).find(TriggerList[i]) != std::string::npos) {
                unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
                if(hltTriggerResults->accept(index)) Map[TriggerList[i]] = true;
                //if(hltTriggerResults->accept(index)) PrescalesMap[TriggerList[i]] = triggerPrescales->getPrescaleForIndex(index);
                PrescalesMap[TriggerList[i]] = triggerPrescales->getPrescaleForIndex(index); //they must be filled even if the trigger is not fired!!!!
		//if(verbosity) std::cout << "Trigger: " << TriggerList[i] << " , prescale: " << triggerPrescales->getPrescaleForIndex(index) << std::endl;
		//std::cout << "L1 Min prescale: " << l1MinTriggerPrescales->getPrescaleForIndex(index) << " , L1 Max prescale: " << l1MaxTriggerPrescales->getPrescaleForIndex(index) << std::endl;
            }
        }
    }

    /*
    std::cout << "\n TRIGGER OBJECTS " << std::endl;
    for(std::vector<pat::TriggerObjectStandAlone>::const_iterator it=triggerObjectCollection->begin(); it!=triggerObjectCollection->end(); ++it)
      {
	pat::TriggerObjectStandAlone obj=*it;
	obj.unpackPathNames(trigNames);

	std::vector< std::string > pathNamesAll = obj.pathNames(false);
	std::vector< std::string > pathNamesLast = obj.pathNames(true);
	if(pathNamesAll.size()>0)
	  {

	    for(unsigned int i = 0; i < TriggerList.size(); i++)
	      {
		for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
		  {
		    if(pathNamesAll[h].find(TriggerList[i]) != std::string::npos)
		      {
			std::cout << "\tTrigger accomplished: " << TriggerList[i]  << std::endl;
			std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
			//// Print trigger object collection and type
			std::cout << "\t   Collection: " << obj.collection() << std::endl;
			std::cout << "\t   Type IDs:   ";
			for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
			std::cout << std::endl;
			//// Print associated trigger filters
			std::cout << "\t   Filters:    ";
			for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
			std::cout << std::endl;

			// Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
			// definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
			// means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
			std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
			for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
			  {
			    bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
			    bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
			    bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
			    bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
			    std::cout << "   " << pathNamesAll[h];
			    if (isBoth) std::cout << "(L,3)";
			    if (isL3 && !isBoth) std::cout << "(*,3)";
			    if (isLF && !isBoth) std::cout << "(L,*)";
			    if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
			  }
			std::cout << std::endl;
			std::cout << std::endl;
		      }

		  }
	      }
	    
	  }


      }
    
    */

}



std::vector<pat::TriggerObjectStandAlone> TriggerAnalyzer::FillTriggerObjectVector(const edm::Event& iEvent, std::string& TrigName) {
//std::vector<pat::TriggerObjectStandAlone> TriggerAnalyzer::FillTriggerObjectVector(const edm::Event& iEvent) {


    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);

    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjectCollection;
    iEvent.getByToken(TriggerObjectToken, triggerObjectCollection);

    std::vector<pat::TriggerObjectStandAlone> Vect;

    //std::cout << "\n TRIGGER OBJECTS " << std::endl;
    for(std::vector<pat::TriggerObjectStandAlone>::const_iterator it=triggerObjectCollection->begin(); it!=triggerObjectCollection->end(); ++it)
      {
	pat::TriggerObjectStandAlone obj=*it;
	obj.unpackPathNames(trigNames);

	std::vector< std::string > pathNamesAll = obj.pathNames(false);
	std::vector< std::string > pathNamesLast = obj.pathNames(true);
	if(pathNamesAll.size()>0)
	  {

		for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
		  {
		    //if(pathNamesAll[h].find("HLT_VBF_DisplacedJet40_VTightID_Hadronic_v") != std::string::npos)
		    if(pathNamesAll[h].find(TrigName) != std::string::npos)
		      {
			/*
			std::cout << "\tTrigger accomplished: " << pathNamesAll[h] << std::endl;
			std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
			//// Print trigger object collection and type
			std::cout << "\t   Collection: " << obj.collection() << std::endl;
			std::cout << "\t   Type IDs:   ";
			for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
			std::cout << std::endl;
			//// Print associated trigger filters
			std::cout << "\t   Filters:    ";
			for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
			std::cout << std::endl;
			*/
			bool type = false;
			for (unsigned l = 0; l < obj.filterIds().size(); ++l)
			  {
			    if(obj.filterIds()[l]>0) type =true;
			  }
			if(type) Vect.push_back(obj);

			/*
			// Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
			// definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
			// means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
			std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
			for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
			  {
			    bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
			    bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
			    bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
			    bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
			    std::cout << "   " << pathNamesAll[h];
			    if (isBoth) std::cout << "(L,3)";
			    if (isL3 && !isBoth) std::cout << "(*,3)";
			    if (isLF && !isBoth) std::cout << "(L,*)";
			    if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
			  }
			std::cout << std::endl;
			std::cout << std::endl;
			*/
		      }

		  }
	      
	    
	  }
	

      }
    
    return Vect;

}


void TriggerAnalyzer::FillMetFiltersMap(const edm::Event& iEvent, std::map<std::string, bool>& Map) { //, edm::EDGetTokenT<edm::TriggerResults>& iToken, std::vector<std::string>& iList) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(MetFiltersToken, hltTriggerResults);//(iToken, hltTriggerResults);//
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    
    //for(unsigned int j=0, in=trigNames.size(); j < in; j++) std::cout << trigNames.triggerName(j) << std::endl;
    
    // Get Trigger index
    for(unsigned int i = 0; i < MetFiltersList.size(); i++) {
        Map[MetFiltersList[i]] = false;
        for(unsigned int j=0, in=trigNames.size(); j < in; j++) {
            if(trigNames.triggerName(j).find(MetFiltersList[i]) != std::string::npos) {
                unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
                if(hltTriggerResults->accept(index)) Map[MetFiltersList[i]] = true;
            }
        }
    }
}

////////////////////////////////////////////////
void TriggerAnalyzer::FillL1FiltersMap(const edm::Event& iEvent, std::map<std::string, bool>& Map) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    
    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjectCollection;
    iEvent.getByToken(TriggerObjectToken, triggerObjectCollection);


    //std::cout << "\n TRIGGER OBJECTS " << std::endl;

    for(unsigned int i = 0; i < L1FiltersList.size(); i++) {
      Map[L1FiltersList[i]] = false;

      for(std::vector<pat::TriggerObjectStandAlone>::const_iterator it=triggerObjectCollection->begin(); it!=triggerObjectCollection->end(); ++it)
	{
	  pat::TriggerObjectStandAlone obj=*it;
	  obj.unpackPathNames(trigNames);

	  std::vector< std::string > pathNamesAll = obj.pathNames(false);
	  std::vector< std::string > pathNamesLast = obj.pathNames(true);
	  if(pathNamesAll.size()>0)
	    {

	      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h)
		{





		  
		  //if(pathNamesAll[h].find("HLT_VBF_DisplacedJet40_VTightID_Hadronic_v") != std::string::npos)//new!NO!IN THIS WAY IT SAVES THE FILTERS ONLY IF THE MAIN PATH IS FIRED! BUUUG!

		    //{//new! and to be debugged....

		      //bool type = false;

		      for (unsigned h = 0; h < obj.filterLabels().size(); ++h)
			{
			  /////if(obj.filterLabels()[h]==) then.... fill the map;
			  if(obj.filterLabels()[h].find(L1FiltersList[i]) != std::string::npos) {
			    //std::cout << "Looping over L1 filters, this found: " << std::endl;
			    //std::cout << obj.filterLabels()[h] << std::endl;
			    Map[L1FiltersList[i]] = true;
			    //unsigned int index = trigNames.triggerIndex(trigNames.triggerName(j));
			    //if(hltTriggerResults->accept(index)) Map[MetFiltersList[i]] = true;
			
			  }
		      
			}



		 //}//new!















		}
	      
	    
	    }
	

	}

    }
    
}

bool TriggerAnalyzer::GetBadPFMuonFlag(const edm::Event& iEvent) { 

    edm::Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(BadPFMuonFilterToken, ifilterbadPFMuon);
    bool filterbadPFMuon = *ifilterbadPFMuon;

    return filterbadPFMuon;
}

bool TriggerAnalyzer::GetBadChCandFlag(const edm::Event& iEvent) { 

    edm::Handle<bool> ifilterbadChCand;
    iEvent.getByToken(BadChCandFilterToken, ifilterbadChCand);
    bool filterbadChCand = *ifilterbadChCand;

    return filterbadChCand;
}

//Pre-Firing
bool TriggerAnalyzer::EvaluatePrefiring(const edm::Event& iEvent) { 

    edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
    iEvent.getByToken(L1GtToken, l1GtHandle);
    bool prefire = l1GtHandle->begin(-1)->getFinalOR();
    return prefire;
}


void TriggerAnalyzer::L1Bits(const edm::Event& iEvent) {

  //(CColl.consumes<BXVector<GlobalAlgBlk>>(PSet.getParameter<edm::InputTag>("l1Gt"))),//Pre-Firing
    //cfg is the ParameterSet

  return;
}





void TriggerAnalyzer::Debug(const edm::Event& iEvent) {

    edm::Handle<edm::TriggerResults> hltTriggerResults;
    iEvent.getByToken(TriggerToken, hltTriggerResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltTriggerResults);
    //open trigger bits


    edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjectCollection;
    iEvent.getByToken(TriggerObjectToken, triggerObjectCollection);
    std::vector<pat::TriggerObjectStandAlone> Vect;
    //open trigger objects

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(PrescalesToken, triggerPrescales);
    //open trigger prescales


    std::cout << "\n == TRIGGER PATHS= " << std::endl;
    for (unsigned int i = 0, n = hltTriggerResults->size(); i < n; ++i) {
      std::cout << "Trigger " << trigNames.triggerName(i) <<
	", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
	": " << (hltTriggerResults->accept(i) ? "PASS" : "fail (or not run)")
                << std::endl;
    }



    std::cout << "\n TRIGGER OBJECTS " << std::endl;
    for (pat::TriggerObjectStandAlone obj : *triggerObjectCollection) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(trigNames);
      std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
      std::cout << "\t   Type IDs:   ";
      for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
      std::cout << std::endl;
      // Print associated trigger filters
      std::cout << "\t   Filters:    ";
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
      std::cout << std::endl;
      std::vector< std::string > pathNamesAll = obj.pathNames(false);
      std::vector< std::string > pathNamesLast = obj.pathNames(true);
      // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
      // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
      // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
      std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
	bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
	bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
	bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
	std::cout << "   " << pathNamesAll[h];
	if (isBoth) std::cout << "(L,3)";
	if (isL3 && !isBoth) std::cout << "(*,3)";
	if (isLF && !isBoth) std::cout << "(L,*)";
	if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;


}
