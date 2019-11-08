#include "PhotonAnalyzer.h"
  
PhotonAnalyzer::PhotonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    PhotonToken(CColl.consumes<std::vector<pat::Photon> >(PSet.getParameter<edm::InputTag>("photons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    PhoLooseIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoLooseIdMap"))),
    PhoMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoMediumIdMap"))),
    PhoTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoTightIdMap"))),
    PhoMVANonTrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("phoMVANonTrigMediumIdMap"))),
    PhoEcalRecHitCollectionToken(CColl.consumes<EcalRecHitCollection>(PSet.getParameter<edm::InputTag>("phoEcalRecHitCollection"))),
    PhoLooseIdFileName(PSet.getParameter<std::string>("phoLooseIdFileName")),
    PhoMediumIdFileName(PSet.getParameter<std::string>("phoMediumIdFileName")),
    PhoTightIdFileName(PSet.getParameter<std::string>("phoTightIdFileName")),
    PhoMVANonTrigMediumIdFileName(PSet.getParameter<std::string>("phoMVANonTrigMediumIdFileName")),
    PhotonId(PSet.getParameter<int>("photonid")),
    PhotonPt(PSet.getParameter<double>("photonpt"))

{

    isPhoLooseIdFile = isPhoMediumIdFile = isPhoTightIdFile = isPhoMVANonTrigMediumIdFile = false;

    PhoLooseIdFile=new TFile(PhoLooseIdFileName.c_str(), "READ");
    if(!PhoLooseIdFile->IsZombie()) {
      PhotonIdLoose=(TH2F*)PhoLooseIdFile->Get("EGamma_SF2D");
      isPhoLooseIdFile=true;
    }
    else {
      throw cms::Exception("PhotonAnalyzer", "No PhoLooseId Weight File");
      return;
    }

    PhoMediumIdFile=new TFile(PhoMediumIdFileName.c_str(), "READ");
    if(!PhoMediumIdFile->IsZombie()) {
      PhotonIdMedium=(TH2F*)PhoMediumIdFile->Get("EGamma_SF2D");
      isPhoMediumIdFile=true;
    }
    else {
      throw cms::Exception("PhotonAnalyzer", "No PhoMediumId Weight File");
      return;
    }

    PhoTightIdFile=new TFile(PhoTightIdFileName.c_str(), "READ");
    if(!PhoTightIdFile->IsZombie()) {
      PhotonIdTight=(TH2F*)PhoTightIdFile->Get("EGamma_SF2D");
      isPhoTightIdFile=true;
    }
    else {
      throw cms::Exception("PhotonAnalyzer", "No PhoTightId Weight File");
      return;
    }

    PhoMVANonTrigMediumIdFile=new TFile(PhoMVANonTrigMediumIdFileName.c_str(), "READ");
    if(!PhoMVANonTrigMediumIdFile->IsZombie()) {
      PhotonIdMVANonTrigMedium=(TH2F*)PhoMVANonTrigMediumIdFile->Get("EGamma_SF2D");
      isPhoMVANonTrigMediumIdFile=true;
    }
    else {
      throw cms::Exception("PhotonAnalyzer", "No PhoMVANonTrigMediumId Weight File");
      return;
    }
    
    
    std::cout << " --- PhotonAnalyzer initialization ---" << std::endl;
    std::cout << "  photon pT         :\t" << PhotonPt << std::endl;
    std::cout << "  photon id         :\t" << PhotonId << std::endl;
    std::cout << std::endl;
}

PhotonAnalyzer::~PhotonAnalyzer() {
  PhoLooseIdFile->Close();
  PhoMediumIdFile->Close();
  PhoTightIdFile->Close();
  PhoMVANonTrigMediumIdFile->Close();
}



std::vector<pat::Photon> PhotonAnalyzer::FillPhotonVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(PhotonId);
    float PtTh(PhotonPt);
    std::vector<pat::Photon> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Photon> > PhoCollection;
    iEvent.getByToken(PhotonToken, PhoCollection);
    
    //edm::Handle<std::vector<pat::Conversion> > EleConv;
    //iEvent.getByToken(edm::InputTag("patConversions"), EleConv);
    
    //edm::Handle<reco::VertexCollection> PVCollection;
    //iEvent.getByToken(VertexToken, PVCollection);
    //const reco::Vertex* vertex=&PVCollection->front();

    edm::Handle<EcalRecHitCollection> _ebrechits;
    iEvent.getByToken(PhoEcalRecHitCollectionToken, _ebrechits);    

    //value map for ID 2015-2016
    edm::Handle<edm::ValueMap<bool> > LooseIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > TightIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVANonTrigMediumIdDecisions;
    iEvent.getByToken(PhoLooseIdMapToken, LooseIdDecisions);
    iEvent.getByToken(PhoMediumIdMapToken, MediumIdDecisions);
    iEvent.getByToken(PhoTightIdMapToken, TightIdDecisions);
    iEvent.getByToken(PhoMVANonTrigMediumIdMapToken, MVANonTrigMediumIdDecisions);
    unsigned int phIdx = 0;

    // Loop on Photon collection
    for(std::vector<pat::Photon>::const_iterator it=PhoCollection->begin(); it!=PhoCollection->end(); ++it) {
        pat::Photon ph=*it;
  pat::PhotonRef phRef(PhoCollection,phIdx);
  
  
        // Corrections for Ele Smearing (on data only) -- MORIOND 2017
        double Ecorr=1;
        if(!isMC) {
            DetId detid = ph.superCluster()->seed()->seed();
            const EcalRecHit * rh = NULL;
            if (detid.subdetId() == EcalBarrel) {
                auto rh_i =  _ebrechits->find(detid);
                            if( rh_i != _ebrechits->end()) rh =  &(*rh_i);
                            else rh = NULL;
                    } 
            if(rh==NULL) Ecorr=1;
            else{
            if(rh->energy() > 200 && rh->energy()<300)  Ecorr=1.0199;
            else if(rh->energy()>300 && rh->energy()<400) Ecorr=  1.052;
            else if(rh->energy()>400 && rh->energy()<500) Ecorr = 1.015;
            }
        }
        ph.setP4(reco::Candidate::LorentzVector(ph.px(), ph.py(), ph.pz(), Ecorr*ph.energy() ));      
  
  
        // Pt and eta
        if(ph.pt()<PtTh || fabs(ph.eta())>2.5) continue;
        float pfIso = ( ph.chargedHadronIso() + std::max(ph.neutralHadronIso() + ph.photonIso() - 0.5*ph.puChargedHadronIso(), 0.) ) / ph.pt();

        //Photon CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2  
        bool isPassLoose = (*LooseIdDecisions)[phRef];
        bool isPassMedium = (*MediumIdDecisions)[phRef];
        bool isPassTight = (*TightIdDecisions)[phRef];
        bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[phRef];

        if(IdTh==1 && !isPassLoose) continue;
        if(IdTh==2 && !isPassMedium) continue;
        if(IdTh==3 && !isPassTight) continue;
        if(IdTh==4 && !isPassMVANonTrigMedium) continue;

        //Fill user float
        ph.addUserFloat("pfIso", pfIso);
        ph.addUserInt("isLoose", isPassLoose ? 1 : 0);
        ph.addUserInt("isMedium", isPassMedium ? 1 : 0);
        ph.addUserInt("isTight", isPassTight ? 1 : 0);
        ph.addUserInt("isMVANonTrigMedium", isPassMVANonTrigMedium ? 1 : 0);
        
        ++phIdx;
      
        Vect.push_back(ph); // Fill vector
    }
    return Vect;
}

void PhotonAnalyzer::PlotPhotons(std::vector<pat::Photon>& PhotonVect, std::map<std::string, TH1F*>& Hist, float EventWeight) {
    if(!(PhotonVect.size()>=2)) return;
    if(!(PhotonVect[0].pt()>80. && PhotonVect[1].pt()>80. && (PhotonVect[0].isEB() || PhotonVect[1].isEB()) )) return;
    float diphoton = (PhotonVect[0].p4()+PhotonVect[1].p4()).mass();
    Hist["p_mass"]->Fill(diphoton, EventWeight);
    if(PhotonVect[0].isEB() && PhotonVect[1].isEB()) Hist["p_massEBEB"]->Fill(diphoton, EventWeight);
    else Hist["p_massEBEE"]->Fill(diphoton, EventWeight);
    Hist["p_massLow"]->Fill(diphoton, EventWeight);
    Hist["p_pt1"]->Fill(PhotonVect[0].pt(), EventWeight);
    Hist["p_pt2"]->Fill(PhotonVect[1].pt(), EventWeight);
    Hist["p_eta1"]->Fill(PhotonVect[0].eta(), EventWeight);
    Hist["p_eta2"]->Fill(PhotonVect[1].eta(), EventWeight);
}

float PhotonAnalyzer::GetPhotonIdSFLoose(pat::Photon& el) {
  if(!isPhoLooseIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdLoose->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdLoose->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdLoose->GetBinContent(PhotonIdLoose->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFLooseError(pat::Photon& el) {
  if(!isPhoLooseIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdLoose->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdLoose->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdLoose->GetBinError(PhotonIdLoose->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFMedium(pat::Photon& el) {
  if(!isPhoMediumIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdMedium->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdMedium->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdMedium->GetBinContent(PhotonIdMedium->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFMediumError(pat::Photon& el) {
  if(!isPhoMediumIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdMedium->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdMedium->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdMedium->GetBinError(PhotonIdMedium->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFTight(pat::Photon& el) {
  if(!isPhoTightIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdTight->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdTight->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdTight->GetBinContent(PhotonIdTight->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFTightError(pat::Photon& el) {
  if(!isPhoTightIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdTight->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdTight->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdTight->GetBinError(PhotonIdTight->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFMVANonTrigMedium(pat::Photon& el) {
  if(!isPhoMVANonTrigMediumIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdMVANonTrigMedium->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdMVANonTrigMedium->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdMVANonTrigMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdMVANonTrigMedium->GetBinContent(PhotonIdMVANonTrigMedium->FindBin(abseta, pt));
}

float PhotonAnalyzer::GetPhotonIdSFMVANonTrigMediumError(pat::Photon& el) {
  if(!isPhoMVANonTrigMediumIdFile) return 1.;
  double pt = std::min( std::max( PhotonIdMVANonTrigMedium->GetYaxis()->GetXmin(), el.pt() ) , PhotonIdMVANonTrigMedium->GetYaxis()->GetXmax() - 0.000001 );
  double abseta = std::min( PhotonIdMVANonTrigMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
  return PhotonIdMVANonTrigMedium->GetBinError(PhotonIdMVANonTrigMedium->FindBin(abseta, pt));
}

void PhotonAnalyzer::CleanPhotonsFromMuons(std::vector<pat::Photon>& Photons, std::vector<pat::Muon>& Muons, float angle) {
    for(unsigned int m = 0; m < Muons.size(); m++) {
        for(unsigned int j = 0; j < Photons.size(); ) {
            if(deltaR(Photons[j], Muons[m]) < angle) Photons.erase(Photons.begin() + j);
            else j++;
        }
    }
}

void PhotonAnalyzer::CleanPhotonsFromElectrons(std::vector<pat::Photon>& Photons, std::vector<pat::Electron>& Electrons, float angle) {
    for(unsigned int e = 0; e < Electrons.size(); e++) {
        for(unsigned int j = 0; j < Photons.size(); ) {
            if(deltaR(Photons[j], Electrons[e]) < angle) Photons.erase(Photons.begin() + j);
            else j++;
        }
    }
}


/*bool PhotonAnalyzer::isLoosePhoton(pat::Photon& el, const reco::Vertex* vertex) {
    return true;
}*/
