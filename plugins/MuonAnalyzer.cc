#include "MuonAnalyzer.h"
  
MuonAnalyzer::MuonAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    MuonToken(CColl.consumes<std::vector<pat::Muon> >(PSet.getParameter<edm::InputTag>("muons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    MuonTrkFileName(PSet.getParameter<std::string>("muonTrkFileName")),
    MuonIdFileName(PSet.getParameter<std::string>("muonIdFileName")),
    MuonIsoFileName(PSet.getParameter<std::string>("muonIsoFileName")),
    MuonTrkHighptFileName(PSet.getParameter<std::string>("muonTrkHighptFileName")),
    MuonTriggerFileName(PSet.getParameter<std::string>("muonTriggerFileName")),
    DoubleMuonTriggerFileName(PSet.getParameter<std::string>("doubleMuonTriggerFileName")), //obsolete
    Muon1Id(PSet.getParameter<int>("muon1id")),
    Muon2Id(PSet.getParameter<int>("muon2id")),
    Muon1Iso(PSet.getParameter<int>("muon1iso")),
    Muon2Iso(PSet.getParameter<int>("muon2iso")),
    Muon1Pt(PSet.getParameter<double>("muon1pt")),
    Muon2Pt(PSet.getParameter<double>("muon2pt")),
    UseTuneP(PSet.getParameter<bool>("useTuneP")),
    DoRochester(PSet.getParameter<bool>("doRochester"))
{
    isMuonTriggerFile = isDoubleMuonTriggerFile = isMuonIdFile = isMuonTrkFile = isMuonTrkHighptFile = false;
    
    // Double Muon trigger: obsolete!
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
    DoubleMuonTriggerFile=new TFile(DoubleMuonTriggerFileName.c_str(), "READ");
    if(!DoubleMuonTriggerFile->IsZombie()) {
        MuonTriggerLt20=(TH2F*)DoubleMuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_10To20_&_Mu2_20ToInfty_with_SYST_uncrt");
        MuonTriggerGt20=(TH2F*)DoubleMuonTriggerFile->Get("DATA_over_MC_Mu17Mu8_Tight_Mu1_20ToInfty_&_Mu2_20ToInfty_with_SYST_uncrt");
        for(int i=1; i<=MuonTriggerGt20->GetNbinsX(); i++) {
            for(int j=1; j<=MuonTriggerGt20->GetNbinsY(); j++) {
                if(j>i) {
                    if(MuonTriggerGt20->GetBinContent(i, j)>0.) std::cout << " - MuonAnalyzer Warning: Trying to symmetrize diagonal matrix in bin " << i << ", " << j << std::endl;
                    MuonTriggerGt20->SetBinContent(i, j, MuonTriggerGt20->GetBinContent(j, i));
                    MuonTriggerGt20->SetBinError(i, j, MuonTriggerGt20->GetBinError(j, i));
                }
            }
        }
        isDoubleMuonTriggerFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No Double Muon Trigger Weight File");
        return;
    }

    //Single Muon Trigger, 2016
    //NOTE -> SF APPLIED AS PER-EVENT WEIGHTS
    MuonTriggerFile=new TFile(MuonTriggerFileName.c_str(), "READ");
    if(!MuonTriggerFile->IsZombie()) {
        MuonTriggerIsoMu24=(TH2F*)MuonTriggerFile->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio");
        MuonTriggerMu50   =(TH2F*)MuonTriggerFile->Get("Mu50_OR_TkMu50_PtEtaBins/pt_abseta_ratio");
        isMuonTriggerFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No Muon Trigger Weight File");
        return;
    }
    
    //Muon tracker eff
    // FIXME -> STILL ICHEP-2016 -> TO BE UPDATED ?
    MuonTrkFile=new TFile(MuonTrkFileName.c_str(), "READ");
    if(!MuonTrkFile->IsZombie()) {
        MuonTrkGraph=(TGraphAsymmErrors*)MuonTrkFile->Get("ratio_eff_eta3_dr030e030_corr");
        MuonTrk=(TH1F*)ConvertTGraph(MuonTrkGraph);
        isMuonTrkFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonTrk Weight File");
        return;
    }

    //Muon id, 2016
    //NOTE -> SF APPLIED AS PER-EVENT WEIGHTS
    MuonIdFile=new TFile(MuonIdFileName.c_str(), "READ");
    if(!MuonIdFile->IsZombie()) {
        MuonIdLoose =(TH2F*)MuonIdFile->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
        MuonIdMedium=(TH2F*)MuonIdFile->Get("MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
        MuonIdTight =(TH2F*)MuonIdFile->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
        MuonIdHighpt=(TH2F*)MuonIdFile->Get("MC_NUM_HighPtIDP_DEN_genTracks_PAR_newpt_eta/pair_ne_ratio"); // done wrt tune-p pt
        isMuonIdFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonId Weight File");
        return;
    }

    //Muon iso, 2016
    //NOTE -> SF APPLIED AS PER-EVENT WEIGHTS
    MuonIsoFile=new TFile(MuonIsoFileName.c_str(), "READ");
    if(!MuonIsoFile->IsZombie()) {
        MuonIsoHighpt=(TH2F*)MuonIsoFile->Get("tkLooseISO_highptID_newpt_eta/pair_ne_ratio");
        MuonIsoLoose=(TH2F*)MuonIsoFile->Get("LooseISO_LooseID_pt_eta/pt_abseta_ratio");
        MuonIsoTight=(TH2F*)MuonIsoFile->Get("TightISO_TightID_pt_eta/pt_abseta_ratio");
        isMuonIsoFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonIso Weight File");
        return;
    }

    //Muon custom TrackerHighPt id, 2016
    // FIXME -> STILL ICHEP-2016 -> TO BE UPDATED ?
    MuonTrkHighptFile=new TFile(MuonTrkHighptFileName.c_str(), "READ");
    if(!MuonTrkHighptFile->IsZombie()) {
        MuonIdTrkHighpt=(TH2F*)MuonTrkHighptFile->Get("scalefactor");
        isMuonTrkHighptFile=true;
    }
    else {
        throw cms::Exception("MuonAnalyzer", "No MuonTrkHighpt Weight File");
        return;
    }

    //rmcor = new rochcor2016();

    std::cout << " --- MuonAnalyzer initialization ---" << std::endl;
    std::cout << "  mu Id  [1, 2]     :\t" << Muon1Id << "\t" << Muon2Id << std::endl;
    std::cout << "  mu Iso [1, 2]     :\t" << Muon1Iso << "\t" << Muon2Iso << std::endl;
    std::cout << "  mu pT  [1, 2]     :\t" << Muon1Pt << "\t" << Muon2Pt << std::endl;
    //std::cout << "  DoRochester       :\t" << DoRochester << std::endl;
    std::cout << std::endl;

}

MuonAnalyzer::~MuonAnalyzer() {
    DoubleMuonTriggerFile->Close();
    MuonTriggerFile->Close();
    MuonTrkFile->Close();
    MuonIdFile->Close();
    MuonIsoFile->Close();
    MuonTrkHighptFile->Close();
}





std::vector<pat::Muon> MuonAnalyzer::FillMuonVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(Muon1Id), IsoTh(Muon1Iso);
    float PtTh(Muon1Pt);
    std::vector<pat::Muon> Vect;
    // Declare and open collections
    edm::Handle<std::vector<pat::Muon> > MuonCollection;
    iEvent.getByToken(MuonToken, MuonCollection);
    
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    
    
    // Loop on Muon collection
    for(std::vector<pat::Muon>::const_iterator it=MuonCollection->begin(); it!=MuonCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Muon2Id;
            IsoTh=Muon2Iso;
            PtTh=Muon2Pt;
        }
        pat::Muon mu=*it;
        // Pt and eta
        if (UseTuneP)
            mu.setP4(reco::Candidate::PolarLorentzVector(mu.tunePMuonBestTrack()->pt(), mu.eta(), mu.phi(), mu.mass()));
        // Apply Rochester corrections
	/*
        if (DoRochester){
            TLorentzVector * mup4 = new TLorentzVector ();        
            mup4->SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
            if(!isMC)
                rmcor->momcor_mc(*mup4, float(mu.charge()), mu.innerTrack()->hitPattern().trackerLayersWithMeasurement(), float(1.0));
            else
                rmcor->momcor_data(*mup4, float(mu.charge()), 0, float(1.0));      
            mu.setP4(reco::Candidate::PolarLorentzVector(mup4->Pt(), mu.eta(), mu.phi(), mu.mass()));
            delete mup4;
        }
	*/
        if(mu.pt()<PtTh || fabs(mu.eta())>2.4) continue;
        // Muon Quality ID 2015-2016: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        if(IdTh==0 && !IsTrackerHighPtMuon(mu, vertex)) continue;
        if(IdTh==1 && !mu.isLooseMuon()) continue;
        if(IdTh==2 && !mu.isMediumMuon()) continue;
        if(IdTh==3 && !mu.isTightMuon(*vertex)) continue;
        if(IdTh==4 && !mu.isHighPtMuon(*vertex)) continue;
        // Isolation 
        float pfIso03 = (mu.pfIsolationR03().sumChargedHadronPt + std::max(mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5*mu.pfIsolationR03().sumPUPt, 0.) ) / mu.pt(); // PF-based pt for PFIso
        float pfIso04 = (mu.pfIsolationR04().sumChargedHadronPt + std::max(mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5*mu.pfIsolationR04().sumPUPt, 0.) ) / mu.pt(); // PF-based pt for PFIso
	      // Tracker iso corrected with by-hand subtraction
        float trkIso = mu.trackIso();
        // Subtrack all muons from iso cone
        for(auto mit=MuonCollection->begin(); mit!=MuonCollection->end(); ++mit) if(mit!=it && deltaR(*mit, mu)<0.3 && IsTrackerHighPtMuon(mu, vertex)) trkIso -= mit->innerTrack()->pt();
        //if(Vect.size() == 0 && std::next(it, 1)!=MuonCollection->end() && deltaR(*std::next(it, 1), mu) < 0.3) trkIso -= std::next(it, 1)->pt();
        //if(Vect.size() == 1 && deltaR(Vect[0], mu) < 0.3) trkIso -= Vect[0].tunePMuonBestTrack()->pt();
        if(trkIso < 0.) trkIso = 0.;
        trkIso /= mu.pt();
        // Muon Isolation working point 2015-2016: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        if(IsoTh==0 && trkIso>0.1) continue;
        if(IsoTh==1 && pfIso04>0.25) continue;
        if(IsoTh==2 && pfIso04>0.15) continue;
        // Add userFloat
        mu.addUserFloat("inTrkPt", mu.innerTrack().isNonnull() ? mu.innerTrack()->pt() : -1);
        mu.addUserFloat("trkIso", trkIso);
        mu.addUserFloat("pfIso03", pfIso03);
        mu.addUserFloat("pfIso04", pfIso04);
        mu.addUserFloat("dxy", mu.muonBestTrack()->dxy(vertex->position()));
        mu.addUserFloat("dz", mu.muonBestTrack()->dz(vertex->position()));
        mu.addUserInt("isTrackerHighPt", IsTrackerHighPtMuon(mu, vertex) ? 1 : 0);
        mu.addUserInt("isLoose", mu.isLooseMuon() ? 1 : 0);
        mu.addUserInt("isMedium", mu.isMediumMuon() ? 1 : 0);
        mu.addUserInt("isTight", mu.isTightMuon(*vertex) ? 1 : 0);
        mu.addUserInt("isHighPt", mu.isHighPtMuon(*vertex) ? 1 : 0);
        // Fill vector
        Vect.push_back(mu);
    }
    return Vect;
}

void MuonAnalyzer::AddVariables(std::vector<pat::Muon>& Vect, pat::MET& MET) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
        Vect[i].addUserFloat("dPhi_met", fabs(reco::deltaPhi(Vect[i].phi(), MET.phi())));
    }
}


bool MuonAnalyzer::IsTrackerHighPtMuon(pat::Muon& mu, const reco::Vertex* vertex) {
    
    if (! (mu.isMuon()) ) return false;
    if (! (mu.isTrackerMuon()) ) return false;
    if (! (mu.tunePMuonBestTrack().isNonnull()) ) return false;
    if (! (mu.numberOfMatchedStations() > 1) ) return false;
    if (! (mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) ) return false;
    if (! (mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) return false;
    if (! (mu.tunePMuonBestTrack()->ptError()/mu.pt() < 0.3) ) return false;
    if (! (fabs(mu.tunePMuonBestTrack()->dxy(vertex->position()) ) < 0.2) ) return false;
    if (! (fabs(mu.tunePMuonBestTrack()->dz(vertex->position()) ) < 0.5) ) return false;

    return true;
    
}

std::vector<float> MuonAnalyzer::FixTrackerIsolation(pat::Muon& mu1, pat::Muon& mu2){
    std::vector<float> FixedTrackIso;
    if(mu1.innerTrack().isNonnull() && mu2.innerTrack().isNonnull()){
        if(deltaR(mu1,mu2)<0.3){
	    FixedTrackIso.push_back(std::max(mu1.trackIso() - mu2.innerTrack()->pt(), 0.) / mu1.tunePMuonBestTrack()->pt() );
	    FixedTrackIso.push_back(std::max(mu2.trackIso() - mu1.innerTrack()->pt(), 0.) / mu2.tunePMuonBestTrack()->pt() );
            return FixedTrackIso;
        }
        else{
            FixedTrackIso.push_back(mu1.trackIso());
            FixedTrackIso.push_back(mu2.trackIso());
            return FixedTrackIso;
	}
    }
    else{
        FixedTrackIso.push_back(-1.);
        FixedTrackIso.push_back(-1.);
        return FixedTrackIso;
    }
}


std::string MuonAnalyzer::GetMuon1Id(pat::Muon& mu){
    if(Muon1Id==0) return "isTrackerHighPt";
    if(Muon1Id==1) return "isLoose";
    if(Muon1Id==2) return "isMedium";
    if(Muon1Id==3) return "isTight";
    if(Muon1Id==4) return "isHighPt";
    else return "";
}

//ID
float MuonAnalyzer::GetMuonIdSF(pat::Muon& mu, int id) {
    if(id==0 && isMuonTrkHighptFile){
        double pt = std::min( std::max( MuonIdTrkHighpt->GetYaxis()->GetXmin(), mu.pt() ) , MuonIdTrkHighpt->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdTrkHighpt->GetBinContent( MuonIdTrkHighpt->FindBin(abseta,pt) );
    }
    if(id==1 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdLoose->GetBinContent( MuonIdLoose->FindBin(pt,abseta) );
    }
    if(id==2 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdMedium->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdMedium->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdMedium->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdMedium->GetBinContent( MuonIdMedium->FindBin(pt,abseta) );
    }
    if(id==3 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdTight->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdTight->GetBinContent( MuonIdTight->FindBin(pt,abseta) );
    }
    if(id==4 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdHighpt->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdHighpt->GetBinContent( MuonIdHighpt->FindBin(pt,abseta) );
    }
    else return 1.;
}

float MuonAnalyzer::GetMuonIdSFError(pat::Muon& mu, int id) {
    if(id==0 && isMuonTrkHighptFile){
        double pt = std::min( std::max( MuonIdTrkHighpt->GetYaxis()->GetXmin(), mu.pt() ) , MuonIdTrkHighpt->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdTrkHighpt->GetBinError( MuonIdTrkHighpt->FindBin(abseta,pt) );
    }
    if(id==1 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdLoose->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdLoose->GetBinError( MuonIdLoose->FindBin(pt,abseta) );
    }
    if(id==2 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdMedium->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdMedium->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdMedium->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdMedium->GetBinError( MuonIdMedium->FindBin(pt,abseta) );
    }
    if(id==3 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdTight->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdTight->GetBinError( MuonIdTight->FindBin(pt,abseta) );
    }
    if(id==4 && isMuonIdFile){
        double pt = std::min( std::max( MuonIdHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIdHighpt->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIdHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIdHighpt->GetBinError( MuonIdHighpt->FindBin(pt,abseta) );
    }
    else return 1.;
}

//TRK
float MuonAnalyzer::GetMuonTrkSF(pat::Muon& mu) {
    if(isMuonTrkFile){
        double eta = 0.;       
        if (mu.eta() > 0)
            eta = std::min( MuonTrk->GetXaxis()->GetXmax() - 0.000001 , mu.eta() );
        else
            eta = std::max( MuonTrk->GetXaxis()->GetXmin() + 0.000001 , mu.eta() );
        return MuonTrk->GetBinContent( MuonTrk->FindBin(eta) );
    }
    else return 1.;
}

float MuonAnalyzer::GetMuonTrkSFError(pat::Muon& mu) {
    if(isMuonTrkFile){
        double eta = 0.;
        if (mu.eta() > 0)
            eta = std::min( MuonTrk->GetXaxis()->GetXmax() - 0.000001 , mu.eta() );
        else
            eta = std::max( MuonTrk->GetXaxis()->GetXmin() + 0.000001 , mu.eta() );
        return MuonTrk->GetBinError( MuonTrk->FindBin(eta) );
    }
    else return 1.;
}

//ISO
float MuonAnalyzer::GetMuonIsoSF(pat::Muon& mu, int id) {
    if(id==0 && isMuonIsoFile){
        double pt = std::min( std::max( MuonIsoHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoHighpt->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoHighpt->GetBinContent( MuonIsoHighpt->FindBin(pt,abseta) );
    }
    if(id==1 && isMuonIsoFile){
        double pt = std::min( std::max( MuonIsoLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoLoose->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoLoose->GetBinContent( MuonIsoLoose->FindBin(pt,abseta) );
    }
    if(id==2 && isMuonIsoFile){
        double pt = std::min( std::max( MuonIsoTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoTight->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoTight->GetBinContent( MuonIsoTight->FindBin(pt,abseta) );
    }
    else return 1.;
}

float MuonAnalyzer::GetMuonIsoSFError(pat::Muon& mu, int id) {
    if(!isMuonIsoFile) return 1.;
    if(id==0){
        double pt = std::min( std::max( MuonIsoHighpt->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoHighpt->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoHighpt->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoHighpt->GetBinError( MuonIsoHighpt->FindBin(pt,abseta) );
    }
    if(id==1){
        double pt = std::min( std::max( MuonIsoLoose->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoLoose->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoLoose->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoLoose->GetBinError( MuonIsoLoose->FindBin(pt,abseta) );
    }
    if(id==2){
        double pt = std::min( std::max( MuonIsoTight->GetXaxis()->GetXmin(), mu.pt() ) , MuonIsoTight->GetXaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( MuonIsoTight->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
        return MuonIsoTight->GetBinError( MuonIsoTight->FindBin(pt,abseta) );
    }
    else return 1.;
}


//obsolete
float MuonAnalyzer::GetDoubleMuonTriggerSF(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isDoubleMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.tunePMuonBestTrack()->pt()<20.) return MuonTriggerLt20->GetBinContent(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinContent(MuonTriggerGt20->FindBin(eta2, eta1));
}

//obsolete
float MuonAnalyzer::GetDoubleMuonTriggerSFError(pat::Muon& mu1, pat::Muon& mu2) {
    if(!isDoubleMuonTriggerFile) return 1.;
    float eta1=fabs(mu1.eta());
    float eta2=fabs(mu2.eta());
    // Muon POG enumeration is inverted 1 <-> 2
    if(mu2.tunePMuonBestTrack()->pt()<20.) return MuonTriggerLt20->GetBinError(MuonTriggerLt20->FindBin(eta2, eta1));
    return MuonTriggerGt20->GetBinError(MuonTriggerGt20->FindBin(eta2, eta1));
}

float MuonAnalyzer::GetMuonTriggerSFIsoMu24(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerIsoMu24->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerIsoMu24->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerIsoMu24->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerIsoMu24->GetBinContent( MuonTriggerIsoMu24->FindBin(pt, abseta) );
}

float MuonAnalyzer::GetMuonTriggerSFErrorIsoMu24(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerIsoMu24->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerIsoMu24->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerIsoMu24->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerIsoMu24->GetBinError( MuonTriggerIsoMu24->FindBin(pt,abseta) );
}

float MuonAnalyzer::GetMuonTriggerSFMu50(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerMu50->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerMu50->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerMu50->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerMu50->GetBinContent( MuonTriggerMu50->FindBin(pt, abseta) );
}

float MuonAnalyzer::GetMuonTriggerSFErrorMu50(pat::Muon& mu) {
    if(!isMuonTriggerFile) return 1.;
    double pt = std::min( std::max( MuonTriggerMu50->GetXaxis()->GetXmin(), mu.pt() ) , MuonTriggerMu50->GetXaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( MuonTriggerMu50->GetYaxis()->GetXmax() - 0.000001 , fabs(mu.eta()) );
    return MuonTriggerMu50->GetBinError( MuonTriggerMu50->FindBin(pt,abseta) );
}




TH1F* MuonAnalyzer::ConvertTGraph(TGraphAsymmErrors* g) {
    int n=g->GetN();
    float x[n+1];
    for(int i=0; i<n; i++) x[i]=g->GetX()[i]-g->GetEXlow()[i];
    x[n]=g->GetX()[n-1]+g->GetEXhigh()[n-1];
    
    TH1F* h=new TH1F(g->GetName(), g->GetTitle(), n, x); h->Sumw2();
    for(int i=0; i<n; i++) {
      h->SetBinContent(i+1, g->GetY()[i]);
      h->SetBinError(i+1, g->GetEYhigh()[i]);
    }

    return h;
}
