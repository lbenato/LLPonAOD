#include "TauAnalyzer.h"
  
TauAnalyzer::TauAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    TauToken(CColl.consumes<std::vector<pat::Tau> >(PSet.getParameter<edm::InputTag>("taus"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    TauPt(PSet.getParameter<double>("taupt")),
    TauEta(PSet.getParameter<double>("taueta")),
    TauIdByDecayMode(PSet.getParameter<int>("tauIdByDecayMode")),
    TauIdByDeltaBetaIso(PSet.getParameter<int>("tauIdByDeltaBetaIso")),
    TauIdByMVAIso(PSet.getParameter<int>("tauIdByMVAIso")),
    TauIdByMuonRejection(PSet.getParameter<int>("tauIdByMuonRejection")),
    TauIdByElectronRejection(PSet.getParameter<int>("tauIdByElectronRejection"))

{

    std::cout << " --- TauAnalyzer initialization ---" << std::endl;
    std::cout << "  tau pT            :\t" << TauPt << std::endl;
    std::cout << "  tau eta           :\t" << TauEta << std::endl;
    std::cout << std::endl;

}

TauAnalyzer::~TauAnalyzer() {
}



std::vector<pat::Tau> TauAnalyzer::FillTauVector(const edm::Event& iEvent) {
    //bool isMC(!iEvent.isRealData());
    float PtTh(TauPt), EtaTh(TauEta);
    std::vector<pat::Tau> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Tau> > TauCollection;
    iEvent.getByToken(TauToken, TauCollection);
        
    //edm::Handle<reco::VertexCollection> PVCollection;
    //iEvent.getByToken(VertexToken, PVCollection);
    //const reco::Vertex* vertex=&PVCollection->front();
    
    // Loop on Tau collection
    for(std::vector<pat::Tau>::const_iterator it=TauCollection->begin(); it!=TauCollection->end(); ++it) {
        pat::Tau tau=*it;
        // Pt and eta
        if(tau.pt()<PtTh || fabs(tau.eta())>EtaTh) continue;
        float pfIso = ( tau.chargedHadronIso() + std::max(tau.neutralHadronIso() + tau.photonIso() - 0.5*tau.puChargedHadronIso(), 0.) ) / tau.pt();

        //Tau id by Decay Mode
        if(TauIdByDecayMode==1 && !tau.tauID("decayModeFinding")) continue;
        if(TauIdByDecayMode==2 && !tau.tauID("decayModeFindingNewDMs")) continue;

        //Tau id by Delta Beta Iso
        if(TauIdByDeltaBetaIso==1 && !tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) continue;
        if(TauIdByDeltaBetaIso==2 && !tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) continue;
        if(TauIdByDeltaBetaIso==3 && !tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")) continue;

        //Tau id by MVA Iso New Decay Mode
        if(TauIdByMVAIso==1 && !tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT")) continue;
        if(TauIdByMVAIso==2 && !tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT")) continue;
        if(TauIdByMVAIso==3 && !tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT")) continue;
        if(TauIdByMVAIso==4 && !tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT")) continue;
        if(TauIdByMVAIso==5 && !tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT")) continue;

        //Tau id by Muon Rejection
        if(TauIdByMuonRejection==1 && !tau.tauID("againstMuonLoose3")) continue;
        if(TauIdByMuonRejection==2 && !tau.tauID("againstMuonTight3")) continue;

        //Tau id by Electron Rejection
        if(TauIdByElectronRejection==1 && !tau.tauID("againstElectronVLooseMVA6")) continue;
        if(TauIdByElectronRejection==2 && !tau.tauID("againstElectronLooseMVA6")) continue;
        if(TauIdByElectronRejection==3 && !tau.tauID("againstElectronMediumMVA6")) continue;
        if(TauIdByElectronRejection==4 && !tau.tauID("againstElectronTightMVA6")) continue;


        //if(IdTh==1 && !isPassLoose) continue;

        //Tau CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2  
        //bool isPassLoose = (*LooseIdDecisions)[tauRef];
        //bool isPassMedium = (*MediumIdDecisions)[tauRef];
        //bool isPassTight = (*TightIdDecisions)[tauRef];
        //bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[tauRef];

        //if(IdTh==1 && !isPassLoose) continue;
        //if(IdTh==2 && !isPassMedium) continue;
        //if(IdTh==3 && !isPassTight) continue;
        //if(IdTh==4 && !isPassMVANonTrigMedium) continue;

        //Fill user float
        tau.addUserFloat("pfIso", pfIso);
        //tau.addUserInt("isLoose", isPassLoose ? 1 : 0);
        //tau.addUserInt("isMedium", isPassMedium ? 1 : 0);
        //tau.addUserInt("isTight", isPassTight ? 1 : 0);        
      
        Vect.push_back(tau); // Fill vector
    }
    return Vect;
}


void TauAnalyzer::CleanTausFromMuons(std::vector<pat::Tau>& Taus, std::vector<pat::Muon>& Muons, float angle) {
    for(unsigned int m = 0; m < Muons.size(); m++) {
        for(unsigned int j = 0; j < Taus.size(); ) {
            if(deltaR(Taus[j], Muons[m]) < angle) Taus.erase(Taus.begin() + j);
            else j++;
        }
    }
}

void TauAnalyzer::CleanTausFromElectrons(std::vector<pat::Tau>& Taus, std::vector<pat::Electron>& Electrons, float angle) {
    for(unsigned int e = 0; e < Electrons.size(); e++) {
        for(unsigned int j = 0; j < Taus.size(); ) {
            if(deltaR(Taus[j], Electrons[e]) < angle) Taus.erase(Taus.begin() + j);
            else j++;
        }
    }
}

void TauAnalyzer::CleanTausFromRecoElectrons(std::vector<pat::Tau>& Taus, std::vector<reco::GsfElectron>& Electrons, float angle) {
    for(unsigned int e = 0; e < Electrons.size(); e++) {
        for(unsigned int j = 0; j < Taus.size(); ) {
            if(deltaR(Taus[j], Electrons[e]) < angle) Taus.erase(Taus.begin() + j);
            else j++;
        }
    }
}
