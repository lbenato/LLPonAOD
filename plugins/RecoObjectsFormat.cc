#include "RecoObjectsFormat.h"

//*******************//                                                         
//        RecoJets       //                                                         
//*******************//                                                         

void RecoObjectsFormat::FillRecoJetType(RecoJetType& I, const reco::Jet* R, bool isMC) {
  if(!R) return;
  I.pt          = R->pt();
  I.eta         = R->eta();
  I.phi         = R->phi();
  I.mass        = R->mass();
}

void RecoObjectsFormat::ResetRecoJetType(RecoJetType& I) {
  I.pt          = -1.;
  I.eta         = -9.;
  I.phi         = -9.;
  I.mass        = -1.;
}

std::string RecoObjectsFormat::ListRecoJetType() {return "pt/F:eta/F:phi/F:mass/F";}

//*******************//
//    Reco MEt       //
//*******************//

void RecoObjectsFormat::FillRecoMEtType(RecoMEtType& I, const reco::PFMET* R, bool isMC) {
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
}

void RecoObjectsFormat::ResetRecoMEtType(RecoMEtType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
}

std::string RecoObjectsFormat::ListRecoMEtType() {return "pt/F:eta/F:phi/F";}

//*******************//                                                         
//  RecoElectrons    //                                                         
//*******************//                                                         

void RecoObjectsFormat::FillRecoElectronType(RecoLeptonType& I, const reco::GsfElectron* R, bool isMC) {
  if(!R) return;
  I.pt          = R->pt();
  I.eta         = R->eta();
  I.phi         = R->phi();
  I.mass        = R->mass();
}

void RecoObjectsFormat::ResetRecoElectronType(RecoLeptonType& I) {
  I.pt          = -1.;
  I.eta         = -9.;
  I.phi         = -9.;
  I.mass        = -1.;
}

std::string RecoObjectsFormat::ListRecoElectronType() {return "pt/F:eta/F:phi/F:mass/F";}


//*******************//
//    Calo Jets      //
//*******************//

/*
void RecoObjectsFormat::FillCaloJetType(CaloJetType& I, const reco::CaloJet* R, bool isMC, bool isMatched) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.emEnergyFraction = R->emEnergyFraction();
    I.emEnergyInEB     = R->emEnergyInEB();
    I.emEnergyInEE     = R->emEnergyInEE();
    I.emEnergyInHF     = R->emEnergyInHF();
    I.energyFractionHadronic = R->energyFractionHadronic();
    I.hadEnergyInHB    = R->hadEnergyInHB();
    I.hadEnergyInHE    = R->hadEnergyInHE();
    I.hadEnergyInHF    = R->hadEnergyInHF();
    I.hadEnergyInHO    = R->hadEnergyInHO();
    I.longLived        = R->longLived();
    I.maxEInEmTowers   = R->maxEInEmTowers();
    I.maxEInHadTowers  = R->maxEInHadTowers();
    I.n60              = R->n60();
    I.n90              = R->n90();
    I.nConstituents    = R->nConstituents();
    I.nPasses          = R->nPasses();
    I.isGenMatched     = isMatched;
}

void RecoObjectsFormat::ResetCaloJetType(CaloJetType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.emEnergyFraction = -1.;
    I.emEnergyInEB     = -1.;
    I.emEnergyInEE     = -1.;
    I.emEnergyInHF     = -1.;
    I.energyFractionHadronic = -1.;
    I.hadEnergyInHB    = -1.;
    I.hadEnergyInHE    = -1.;
    I.hadEnergyInHF    = -1.;
    I.hadEnergyInHO    = -1.;
    I.longLived        = false;
    I.maxEInEmTowers   = -1.;
    I.maxEInHadTowers  = -1.;
    I.n60              = -1;
    I.n90              = -1;
    I.nConstituents    = -1;
    I.nPasses          = -1;
    I.isGenMatched        = false;
}

std::string RecoObjectsFormat::ListCaloJetType() {return "pt/F:eta/F:phi/F:mass/F:energy/F:emEnergyFraction/F:emEnergyInEB/F:emEnergyInEE/F:emEnergyInHF/F:energyFractionHadronic/F:hadEnergyInHB/F:hadEnergyInHE/F:hadEnergyInHF/F:hadEnergyInHO/F:longLived/O:maxEInEmTowers/F:maxEInHadTowers/F:n60/I:n90/I:nConstituents/I:nPasses/I:isGenMatched/O";}

void RecoObjectsFormat::FillMEtType(MEtType& I, const reco::PFMET* R, bool isMC) {
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.phi         = R->phi();
    I.sign        = -1.;//R->metSignificance();//not available for PFMET
    I.ptCalo      = -1.;//R->caloMETPt();//not available for PFMET
}

void RecoObjectsFormat::ResetMEtType(MEtType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.phi         = -9.;
    I.sign        = -1.;
    I.ptCalo      = -1.;
}


std::string RecoObjectsFormat::ListMEtType() {return "pt/F:eta/F:phi/F:sign/F:ptCalo/F";}
*/

//************************//
//       GenParticles     //
//***********************//
/*
void RecoObjectsFormat::FillGenPType(GenPType& I, const reco::GenParticle* R) {
    if(!R) return;
    I.pt          = R->pt();
    I.eta         = R->eta();
    I.rapidity    = R->rapidity();
    I.phi         = R->phi();
    I.mass        = R->mass();
    I.energy      = R->energy();
    I.charge      = R->charge();
    I.pdgId       = R->pdgId();
    I.status      = R->status();
    I.radius      = R->mother()? sqrt(pow(R->vx() - R->mother()->vx(),2) + pow(R->vy() - R->mother()->vy(),2) + pow(R->vz() - R->mother()->vz(),2)) : -1000.;
    I.motherid    = R->mother()? R->mother()->pdgId() : 0;
    I.vx          = R->vx();
    I.vy          = R->vy();
    I.vz          = R->vz();
}


void RecoObjectsFormat::ResetGenPType(GenPType& I) {
    I.pt          = -1.;
    I.eta         = -9.;
    I.rapidity    = -99999.;
    I.phi         = -9.;
    I.mass        = -1.;
    I.energy      = -1.;
    I.charge      = 0;
    I.pdgId       = 0;
    I.status      = 0;
    I.radius      = -1.;
    I.motherid    = 0;
    I.vx          = -99;
    I.vy          = -99;
    I.vz          = -99;

}

std::string RecoObjectsFormat::ListGenPType() {return "pt/F:eta/F:rapidity/F:phi/F:mass/F:energy/F:charge/I:pdgId/I:status/I:radius/F:motherid/I:vx/F:vy/F:vz/F";}




void RecoObjectsFormat::FillCandidateType(CandidateType& I, reco::CompositeCandidate* R, bool isMC) {
  if(!R) return;
  if(R->numberOfDaughters() == 0) return;
  I.pt          = R->pt();
  I.eta         = R->eta();
  I.phi         = R->phi();
  I.mass        = R->mass();
  I.tmass       = R->numberOfDaughters()>1 ? sqrt( 2.*R->daughter(0)->pt()*R->daughter(1)->pt()*(1.-cos(deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi())) ) ) : (R->numberOfDaughters()==1 ? R->mt() : -1.);
  I.dR          = R->numberOfDaughters()>1 ? deltaR(*R->daughter(0), *R->daughter(1)) : -1.;
  I.dEta        = R->numberOfDaughters()>1 ? fabs( R->daughter(0)->eta() - R->daughter(1)->eta() ) : -1.;
  I.dPhi        = R->numberOfDaughters()>1 ? fabs( deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi()) ) : -1.;
  I.twist       = R->numberOfDaughters()>1 ? fabs(atan( deltaPhi(R->daughter(0)->phi(), R->daughter(1)->phi())/fabs(R->daughter(0)->eta()-R->daughter(1)->eta()) )) : -1.;
}

void RecoObjectsFormat::ResetCandidateType(CandidateType& I) {
  I.pt          = -1.;
  I.eta         = -9.;
  I.phi         = -9.;
  I.mass        = -1.;
  I.tmass       = -1.;
  I.dR          = -1.;
  I.dEta        = -1.;
  I.dPhi        = -1.;
  I.twist       = -1.;
}

std::string RecoObjectsFormat::ListCandidateType() {return "pt/F:eta/F:phi/F:mass/F:tmass/F:dR/F:dEta/F:dPhi/F:twist/F";}

*/
