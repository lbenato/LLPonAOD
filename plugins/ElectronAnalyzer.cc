#include "ElectronAnalyzer.h"
#include "EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class.cc"


ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    ElectronToken(CColl.consumes<std::vector<pat::Electron> >(PSet.getParameter<edm::InputTag>("electrons"))),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    EleVetoIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleVetoIdMap"))),
    EleLooseIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleLooseIdMap"))),
    EleMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMediumIdMap"))),
    EleTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleTightIdMap"))),
    EleHEEPIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleHEEPIdMap"))),
    EleMVANonTrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVANonTrigMediumIdMap"))),
    EleMVANonTrigTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVANonTrigTightIdMap"))),
    EleMVATrigMediumIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVATrigMediumIdMap"))),
    EleMVATrigTightIdMapToken(CColl.consumes<edm::ValueMap<bool>>(PSet.getParameter<edm::InputTag>("eleMVATrigTightIdMap"))),
    EleEcalRecHitCollectionToken(CColl.consumes<EcalRecHitCollection>(PSet.getParameter<edm::InputTag>("eleEcalRecHitCollection"))),
    EleSingleTriggerIsoFileName(PSet.getParameter<std::string>("eleSingleTriggerIsoFileName")),
    EleSingleTriggerFileName(PSet.getParameter<std::string>("eleSingleTriggerFileName")),
    EleVetoIdFileName(PSet.getParameter<std::string>("eleVetoIdFileName")),
    EleLooseIdFileName(PSet.getParameter<std::string>("eleLooseIdFileName")),
    EleMediumIdFileName(PSet.getParameter<std::string>("eleMediumIdFileName")),
    EleTightIdFileName(PSet.getParameter<std::string>("eleTightIdFileName")),
    EleMVATrigMediumIdFileName(PSet.getParameter<std::string>("eleMVATrigMediumIdFileName")),
    EleMVATrigTightIdFileName(PSet.getParameter<std::string>("eleMVATrigTightIdFileName")),
    EleRecoEffFileName(PSet.getParameter<std::string>("eleRecoEffFileName")),
    EleScaleSmearCorrectionName(PSet.getParameter<std::string>("eleScaleSmearCorrectionName")),
    Electron1Id(PSet.getParameter<int>("electron1id")),
    Electron2Id(PSet.getParameter<int>("electron2id")),
    //Electron1Iso(PSet.getParameter<int>("electron1iso")),
    //Electron2Iso(PSet.getParameter<int>("electron2iso")),
    Electron1Pt(PSet.getParameter<double>("electron1pt")),
    Electron2Pt(PSet.getParameter<double>("electron2pt"))
{
    isEleVetoIdFile = isEleLooseIdFile = isEleMediumIdFile = isEleTightIdFile = isEleRecoEffFile = isEleMVATrigMediumIdFile = isEleMVATrigTightIdFile = isEleTriggerFile = isEleSingleTriggerFile = isEleSingleTriggerIsoFile = false;    

    // Electron SingleTriggerIso
    EleSingleTriggerIsoFile=new TFile(EleSingleTriggerIsoFileName.c_str(), "READ");
    if(!EleSingleTriggerIsoFile->IsZombie()) {
      ElectronTriggerEle27Tight=(TH2F*)EleSingleTriggerIsoFile->Get("Ele27_WPTight/eleTrigEff_Ele27Tight");//X:eta;Y:pt
      isEleSingleTriggerIsoFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No SingleTriggerIso File");
      return;
    }

    // Electron SingleTrigger
    EleSingleTriggerFile=new TFile(EleSingleTriggerFileName.c_str(), "READ");
    if(!EleSingleTriggerFile->IsZombie()) {
      ElectronTriggerEle105=(TH2F*)EleSingleTriggerFile->Get("sf");//X:pt;Y:eta
      isEleSingleTriggerFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No SingleTrigger File");
      return;
    }

    // Electron reco SF 2015-2016
    EleRecoEffFile=new TFile(EleRecoEffFileName.c_str(), "READ");
    if(!EleRecoEffFile->IsZombie()) {
      ElectronRecoEff=(TH2F*)EleRecoEffFile->Get("EGamma_SF2D");
      isEleRecoEffFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleRecoEff Weight File");
      return;
    }

    // Electron id SF, 2015-2016
    EleVetoIdFile=new TFile(EleVetoIdFileName.c_str(), "READ");
    if(!EleVetoIdFile->IsZombie()) {
      ElectronIdVeto=(TH2F*)EleVetoIdFile->Get("EGamma_SF2D");
      isEleVetoIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleVetoId Weight File");
      return;
    }

    EleLooseIdFile=new TFile(EleLooseIdFileName.c_str(), "READ");
    if(!EleLooseIdFile->IsZombie()) {
      ElectronIdLoose=(TH2F*)EleLooseIdFile->Get("EGamma_SF2D");
      isEleLooseIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleLooseId Weight File");
      return;
    }

    EleMediumIdFile=new TFile(EleMediumIdFileName.c_str(), "READ");
    if(!EleMediumIdFile->IsZombie()) {
      ElectronIdMedium=(TH2F*)EleMediumIdFile->Get("EGamma_SF2D");
      isEleMediumIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleMediumId Weight File");
      return;
    }

    EleTightIdFile=new TFile(EleTightIdFileName.c_str(), "READ");
    if(!EleTightIdFile->IsZombie()) {
      ElectronIdTight=(TH2F*)EleTightIdFile->Get("EGamma_SF2D");
      isEleTightIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleTightId Weight File");
      return;
    }

    EleMVATrigMediumIdFile=new TFile(EleMVATrigMediumIdFileName.c_str(), "READ");
    if(!EleMVATrigMediumIdFile->IsZombie()) {
      ElectronIdMVATrigMedium=(TH2F*)EleMVATrigMediumIdFile->Get("EGamma_SF2D");
      isEleMVATrigMediumIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleMVATrigMediumId Weight File");
      return;
    }

    EleMVATrigTightIdFile=new TFile(EleMVATrigTightIdFileName.c_str(), "READ");
    if(!EleMVATrigTightIdFile->IsZombie()) {
      ElectronIdMVATrigTight=(TH2F*)EleMVATrigTightIdFile->Get("EGamma_SF2D");
      isEleMVATrigTightIdFile=true;
    }
    else {
      throw cms::Exception("ElectronAnalyzer", "No EleMVATrigTightId Weight File");
      return;
    }
    
    eScaleSmearer = new EnergyScaleCorrection_class(EleScaleSmearCorrectionName);

    std::cout << " --- ElectronAnalyzer initialization ---" << std::endl;
    std::cout << "  jet Id [1, 2]     :\t" << Electron1Id << "\t" << Electron2Id << std::endl;
    std::cout << "  jet pT [1, 2]     :\t" << Electron1Pt << "\t" << Electron2Pt << std::endl;
    std::cout << "  veto Id file      :\t" << EleVetoIdFileName << std::endl;
    std::cout << "  loose Id file     :\t" << EleLooseIdFileName << std::endl;
    std::cout << "  medium Id file    :\t" << EleMediumIdFileName << std::endl;
    std::cout << "  tight Id file     :\t" << EleTightIdFileName << std::endl;
    std::cout << "  mva medium Id file:\t" << EleMVATrigMediumIdFileName << std::endl;
    std::cout << "  mva tight Id file :\t" << EleMVATrigTightIdFileName << std::endl;
    std::cout << "  reco eff file     :\t" << EleRecoEffFileName << std::endl;
    std::cout << std::endl;
}

ElectronAnalyzer::~ElectronAnalyzer() {
    EleSingleTriggerFile->Close();
    EleSingleTriggerIsoFile->Close();
    EleVetoIdFile->Close();
    EleLooseIdFile->Close();
    EleMediumIdFile->Close();
    EleTightIdFile->Close();
    EleMVATrigMediumIdFile->Close();
    EleMVATrigTightIdFile->Close();
    EleRecoEffFile->Close();
}





std::vector<pat::Electron> ElectronAnalyzer::FillElectronVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    int IdTh(Electron1Id);//, IsoTh(Electron1Iso);
    float PtTh(Electron1Pt);
    std::vector<pat::Electron> Vect;
    // Declare and open collection
    edm::Handle<std::vector<pat::Electron> > EleCollection;
    iEvent.getByToken(ElectronToken, EleCollection);
    
    //edm::Handle<std::vector<pat::Conversion> > EleConv;
    //iEvent.getByToken(edm::InputTag("patConversions"), EleConv);
    
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    const reco::Vertex* vertex=&PVCollection->front();
    
    edm::Handle<EcalRecHitCollection> _ebrechits;
    iEvent.getByToken(EleEcalRecHitCollectionToken, _ebrechits);

    //value map for ID 2015-2016
    edm::Handle<edm::ValueMap<bool> > VetoIdDecisions;
    edm::Handle<edm::ValueMap<bool> > LooseIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > TightIdDecisions;
    edm::Handle<edm::ValueMap<bool> > HEEPIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVANonTrigMediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVANonTrigTightIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVATrigMediumIdDecisions;
    edm::Handle<edm::ValueMap<bool> > MVATrigTightIdDecisions;
    iEvent.getByToken(EleVetoIdMapToken, VetoIdDecisions);
    iEvent.getByToken(EleLooseIdMapToken, LooseIdDecisions);
    iEvent.getByToken(EleMediumIdMapToken, MediumIdDecisions);
    iEvent.getByToken(EleTightIdMapToken, TightIdDecisions);
    iEvent.getByToken(EleHEEPIdMapToken, HEEPIdDecisions);
    iEvent.getByToken(EleMVANonTrigMediumIdMapToken, MVANonTrigMediumIdDecisions);
    iEvent.getByToken(EleMVANonTrigTightIdMapToken, MVANonTrigTightIdDecisions);
    iEvent.getByToken(EleMVATrigMediumIdMapToken, MVATrigMediumIdDecisions);
    iEvent.getByToken(EleMVATrigTightIdMapToken, MVATrigTightIdDecisions);
    
    unsigned int elIdx = 0;

    if (isMC){
        eScaleSmearer->doSmearings = true;
    }
    else{
        eScaleSmearer->doScale = true;
    }

    // Loop on Electron collection
    for(std::vector<pat::Electron>::const_iterator it=EleCollection->begin(); it!=EleCollection->end(); ++it) {
        if(Vect.size()>0) {
            IdTh=Electron2Id;
            //IsoTh=Electron2Iso;
            PtTh=Electron2Pt;
        }
        pat::Electron el=*it;
        pat::ElectronRef elRef(EleCollection, elIdx);
                        
        // Pt and eta
        if(el.pt()<PtTh || fabs(el.eta())>2.5) continue;
        // PF (?) Isolation R=0.4 https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPFBasedIsolation#for_PAT_electron_users_using_sta
        float pfIso04 = ( el.chargedHadronIso() + std::max(el.neutralHadronIso() + el.photonIso() - 0.5*el.puChargedHadronIso(), 0.) ) / el.pt();
        float pfIso03 = ( el.pfIsolationVariables().sumChargedHadronPt + std::max(el.pfIsolationVariables().sumNeutralHadronEt + el.pfIsolationVariables().sumPhotonEt - 0.5*el.pfIsolationVariables().sumPUPt, 0.) ) / el.pt();
        //if(IsoTh==1 && pfIso03>0.15) continue;
        //Electron CutBased and HEEP ID 2015-2016, https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
        bool isPassVeto = (*VetoIdDecisions)[elRef];
        bool isPassLoose = (*LooseIdDecisions)[elRef];
        bool isPassMedium = (*MediumIdDecisions)[elRef];
        bool isPassTight = (*TightIdDecisions)[elRef];
        bool isPassHEEP = (*HEEPIdDecisions)[elRef];
        bool isPassMVANonTrigMedium = (*MVANonTrigMediumIdDecisions)[elRef];
        bool isPassMVANonTrigTight = (*MVANonTrigTightIdDecisions)[elRef];
        bool isPassMVATrigMedium = (*MVATrigMediumIdDecisions)[elRef];
        bool isPassMVATrigTight = (*MVATrigTightIdDecisions)[elRef];

        if(IdTh==0 && !isPassVeto) continue;
        if(IdTh==1 && !isPassLoose) continue;
        if(IdTh==2 && !isPassMedium) continue;
        if(IdTh==3 && !isPassTight) continue;
        if(IdTh==4 && !isPassHEEP) continue;
        if(IdTh==5 && !isPassMVANonTrigMedium) continue;
        if(IdTh==6 && !isPassMVANonTrigTight) continue;
        if(IdTh==7 && !isPassMVATrigMedium) continue;
        if(IdTh==8 && !isPassMVATrigTight) continue;
        // Fill userFloat
        el.addUserFloat("pfIso03", pfIso03);
        el.addUserFloat("pfIso04", pfIso04);
        el.addUserFloat("trkIso", el.pfIsolationVariables().sumChargedHadronPt);
        el.addUserFloat("dxy", el.gsfTrack()->dxy( vertex->position() ));
        el.addUserFloat("dz", el.gsfTrack()->dz( vertex->position() ));
        el.addUserInt("isVeto", isPassVeto ? 1 : 0);
        el.addUserInt("isLoose", isPassLoose ? 1 : 0);
        el.addUserInt("isMedium", isPassMedium ? 1 : 0);
        el.addUserInt("isTight", isPassTight ? 1 : 0);
        el.addUserInt("isHEEP", isPassHEEP ? 1 : 0);
        el.addUserInt("isMVANonTrigMedium", isPassMVANonTrigMedium ? 1 : 0);
        el.addUserInt("isMVANonTrigTight", isPassMVANonTrigTight ? 1 : 0);
        el.addUserInt("isMVATrigMedium", isPassMVATrigMedium ? 1 : 0);
        el.addUserInt("isMVATrigTight", isPassMVATrigTight ? 1 : 0);

        float scale         = eScaleSmearer->ScaleCorrection           (iEvent.id().run(), el.isEB(), el.r9(), el.superCluster()->eta(), el.et(), 0);
        float sigma         = eScaleSmearer->getSmearingSigma          (iEvent.id().run(), el.isEB(), el.r9(), el.superCluster()->eta(), el.et(), 0, 0, 0);
        float error_scale   = eScaleSmearer->ScaleCorrectionUncertainty(iEvent.id().run(), el.isEB(), el.r9(), el.superCluster()->eta(), el.et(), 0);
        float error_sigma_u = eScaleSmearer->getSmearingSigma          (iEvent.id().run(), el.isEB(), el.r9(), el.superCluster()->eta(), el.et(), 0, 1,  0);
        float error_sigma_d = eScaleSmearer->getSmearingSigma          (iEvent.id().run(), el.isEB(), el.r9(), el.superCluster()->eta(), el.et(), 0, -1, 0);

        el.addUserFloat("SSscale",          scale);
        el.addUserFloat("SSsigma",          sigma);
        el.addUserFloat("SSscaleUnc",       error_scale);
        el.addUserFloat("SSsigmaUncUp",     error_sigma_u);
        el.addUserFloat("SSsigmaUncDown",   error_sigma_d);

        TRandom3 * rgen_ = new TRandom3(0);

        if(isMC){
            float smear = rgen_->Gaus(1, sigma);
            float smearUp = rgen_->Gaus(1, error_sigma_u);
            float smearDown = rgen_->Gaus(1, error_sigma_d);
            el.addUserFloat("SScorr",               smear);
            el.addUserFloat("ptSScorr",             el.pt()*smear);
            el.addUserFloat("ptSScorrUncUp",        el.pt()*smearUp);
            el.addUserFloat("ptSScorrUncDown",      el.pt()*smearDown);
            el.addUserFloat("energySScorr",         el.energy()*smear);
            el.addUserFloat("energySScorrUncUp",    el.energy()*smearUp);
            el.addUserFloat("energySScorrUncDown",  el.energy()*smearDown);
        }
        else{
            el.addUserFloat("SScorr",               scale);
            el.addUserFloat("ptSScorr",             el.pt()*scale);
            el.addUserFloat("ptSScorrUncUp",        el.pt()*(1+error_scale));
            el.addUserFloat("ptSScorrUncDown",      el.pt()*(1-error_scale));
            el.addUserFloat("energySScorr",         el.energy()*scale);
            el.addUserFloat("energySScorrUncUp",    el.energy()*(1+error_scale));
            el.addUserFloat("energySScorrUncDown",  el.energy()*(1-error_scale));            
        }

        ++elIdx;
        // Fill vector
        Vect.push_back(el);

    }
    return Vect;
}


void ElectronAnalyzer::AddVariables(std::vector<pat::Electron>& Vect, pat::MET& MET) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
        Vect[i].addUserFloat("dPhi_met", fabs(reco::deltaPhi(Vect[i].phi(), MET.phi())));
    }
}


float ElectronAnalyzer::GetElectronIdSF(pat::Electron& el, int id) {
    if(id==0 && isEleVetoIdFile){
        double pt = std::min( std::max( ElectronIdVeto->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdVeto->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdVeto->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdVeto->GetBinContent(ElectronIdVeto->FindBin(abseta, pt));
    }
    if(id==1 && isEleLooseIdFile){
        double pt = std::min( std::max( ElectronIdLoose->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdLoose->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdLoose->GetBinContent(ElectronIdLoose->FindBin(abseta, pt));
    }
    if(id==2 && isEleMediumIdFile){
        double pt = std::min( std::max( ElectronIdMedium->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMedium->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMedium->GetBinError(ElectronIdMedium->FindBin(abseta, pt));
    }
    if(id==3 && isEleTightIdFile){
        double pt = std::min( std::max( ElectronIdTight->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdTight->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdTight->GetBinContent(ElectronIdTight->FindBin(abseta, pt));
    }
    if(id==7 && isEleMVATrigMediumIdFile){
        double pt = std::min( std::max( ElectronIdMVATrigMedium->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMVATrigMedium->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMVATrigMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMVATrigMedium->GetBinContent(ElectronIdMVATrigMedium->FindBin(abseta, pt));
    }
    if(id==8 && isEleMVATrigTightIdFile){
        double pt = std::min( std::max( ElectronIdMVATrigTight->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMVATrigTight->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMVATrigTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMVATrigTight->GetBinContent(ElectronIdMVATrigTight->FindBin(abseta, pt));
    }
    else{
        return 1.;
    }
}

float ElectronAnalyzer::GetElectronIdSFError(pat::Electron& el, int id) {
    if(id==0 && isEleVetoIdFile){
        double pt = std::min( std::max( ElectronIdVeto->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdVeto->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdVeto->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdVeto->GetBinError(ElectronIdVeto->FindBin(abseta, pt));
    }
    if(id==1 && isEleLooseIdFile){
        double pt = std::min( std::max( ElectronIdLoose->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdLoose->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdLoose->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdLoose->GetBinError(ElectronIdLoose->FindBin(abseta, pt));
    }
    if(id==2 && isEleMediumIdFile){
        double pt = std::min( std::max( ElectronIdMedium->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMedium->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMedium->GetBinError(ElectronIdMedium->FindBin(abseta, pt));
    }
    if(id==3 && isEleTightIdFile){
        double pt = std::min( std::max( ElectronIdTight->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdTight->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdTight->GetBinError(ElectronIdTight->FindBin(abseta, pt));
    }
    if(id==7 && isEleMVATrigMediumIdFile){
        double pt = std::min( std::max( ElectronIdMVATrigMedium->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMVATrigMedium->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMVATrigMedium->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMVATrigMedium->GetBinError(ElectronIdMVATrigMedium->FindBin(abseta, pt));
    }
    if(id==8 && isEleMVATrigTightIdFile){
        double pt = std::min( std::max( ElectronIdMVATrigTight->GetYaxis()->GetXmin(), el.pt() ) , ElectronIdMVATrigTight->GetYaxis()->GetXmax() - 0.000001 );
        double abseta = std::min( ElectronIdMVATrigTight->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
        return ElectronIdMVATrigTight->GetBinError(ElectronIdMVATrigTight->FindBin(abseta, pt));
    }
    else{
        return 1.;
    }
}

float ElectronAnalyzer::GetDoubleElectronTriggerSF(pat::Electron& el1, pat::Electron& el2) { //obsolete!
    if(!isEleTriggerFile) return 1.;
    float pt1=el1.pt(), pt2=el2.pt();
    float eta1(fabs(el1.eta())), eta2(fabs(el2.eta()));
    if(pt1>=EleTriggerPtMax) pt1=EleTriggerPtMax-1.;
    if(pt2>=EleTriggerPtMax) pt2=EleTriggerPtMax-1.;
    
    float effDATAHighEle1 = EleTriggerDATAHighLeg->GetBinContent(EleTriggerDATAHighLeg->FindBin(pt1, eta1));
    float effDATALowEle2 = EleTriggerDATALowLeg->GetBinContent(EleTriggerDATALowLeg->FindBin(pt2, eta2));
    float effDATAHighEle2 = EleTriggerDATAHighLeg->GetBinContent(EleTriggerDATAHighLeg->FindBin(pt2, eta2));
    float effDATALowEle1 = EleTriggerDATALowLeg->GetBinContent(EleTriggerDATALowLeg->FindBin(pt1, eta1));
    
    float effMCHighEle1 = EleTriggerMCHighLeg->GetBinContent(EleTriggerMCHighLeg->FindBin(pt1, eta1));
    float effMCLowEle2 = EleTriggerMCLowLeg->GetBinContent(EleTriggerMCLowLeg->FindBin(pt2, eta2));
    float effMCHighEle2 = EleTriggerMCHighLeg->GetBinContent(EleTriggerMCHighLeg->FindBin(pt2, eta2));
    float effMCLowEle1 = EleTriggerMCLowLeg->GetBinContent(EleTriggerMCLowLeg->FindBin(pt1, eta1));
    
    float effDATA = effDATAHighEle1*effDATALowEle2 + effDATALowEle1*effDATAHighEle2 - effDATAHighEle1*effDATAHighEle2;
    float effMC = effMCHighEle1*effMCLowEle2 + effMCLowEle1*effMCHighEle2 - effMCHighEle1*effMCHighEle2;
    return effDATA/effMC;
}

float ElectronAnalyzer::GetElectronRecoEffSF(pat::Electron& el) {
    if(!isEleRecoEffFile) return 1.;
    double pt = std::min( std::max( ElectronRecoEff->GetYaxis()->GetXmin(), el.pt() ) , ElectronRecoEff->GetYaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( ElectronRecoEff->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
    return ElectronRecoEff->GetBinContent(ElectronRecoEff->FindBin(abseta, pt));
}

float ElectronAnalyzer::GetElectronRecoEffSFError(pat::Electron& el) {
    if(!isEleRecoEffFile) return 1.;
    double pt = std::min( std::max( ElectronRecoEff->GetYaxis()->GetXmin(), el.pt() ) , ElectronRecoEff->GetYaxis()->GetXmax() - 0.000001 );
    double abseta = std::min( ElectronRecoEff->GetXaxis()->GetXmax() - 0.000001 , fabs(el.eta()) );
    return ElectronRecoEff->GetBinError(ElectronRecoEff->FindBin(abseta, pt));
}

float ElectronAnalyzer::GetElectronTriggerSFEle105(pat::Electron& ele) {
    if(!isEleSingleTriggerFile) return 1.;
    double pt = std::min( std::max( ElectronTriggerEle105->GetYaxis()->GetXmin(), ele.pt() ) , ElectronTriggerEle105->GetYaxis()->GetXmax() - 0.000001 );
    double eta = 0.;
    if (ele.eta() > 0){
        eta = std::min( ElectronTriggerEle105->GetXaxis()->GetXmax() - 0.000001 , ele.eta() );
    }
    else{
        eta = std::max( ElectronTriggerEle105->GetXaxis()->GetXmin() + 0.000001 , ele.eta() );
    }
    return ElectronTriggerEle105->GetBinContent( ElectronTriggerEle105->FindBin(eta, pt) );
}

float ElectronAnalyzer::GetElectronTriggerSFErrorEle105(pat::Electron& ele) {
    if(!isEleSingleTriggerFile) return 1.;
    double pt = std::min( std::max( ElectronTriggerEle105->GetYaxis()->GetXmin(), ele.pt() ) , ElectronTriggerEle105->GetYaxis()->GetXmax() - 0.000001 );
    double eta = 0.;
    if (ele.eta() > 0)
        eta = std::min( ElectronTriggerEle105->GetXaxis()->GetXmax() - 0.000001 , ele.eta() );
    else
        eta = std::max( ElectronTriggerEle105->GetXaxis()->GetXmin() + 0.000001 , ele.eta() );
    return ElectronTriggerEle105->GetBinError( ElectronTriggerEle105->FindBin(eta, pt) );
}

float ElectronAnalyzer::GetElectronTriggerSFEle27Tight(pat::Electron& ele) {
    if(!isEleSingleTriggerIsoFile) return 1.;
    double pt = std::min( std::max( ElectronTriggerEle27Tight->GetYaxis()->GetXmin(), ele.pt() ) , ElectronTriggerEle27Tight->GetYaxis()->GetXmax() - 0.000001 );
    double eta = 0.;
    if (ele.eta() > 0)
        eta = std::min( ElectronTriggerEle27Tight->GetXaxis()->GetXmax() - 0.000001 , ele.eta() );
    else
        eta = std::max( ElectronTriggerEle27Tight->GetXaxis()->GetXmin() + 0.000001 , ele.eta() );
        
    return ElectronTriggerEle27Tight->GetBinContent( ElectronTriggerEle27Tight->FindBin(eta, pt) );
}

float ElectronAnalyzer::GetElectronTriggerSFErrorEle27Tight(pat::Electron& ele) {
  if(!isEleSingleTriggerIsoFile) return 1.;
  double pt = std::min( std::max( ElectronTriggerEle27Tight->GetYaxis()->GetXmin(), ele.pt() ) , ElectronTriggerEle27Tight->GetYaxis()->GetXmax() - 0.000001 );
  double eta = 0.;
  if (ele.eta() > 0)
    eta = std::min( ElectronTriggerEle27Tight->GetXaxis()->GetXmax() - 0.000001 , ele.eta() );
  else
    eta = std::max( ElectronTriggerEle27Tight->GetXaxis()->GetXmin() + 0.000001 , ele.eta() );

  return ElectronTriggerEle27Tight->GetBinError( ElectronTriggerEle27Tight->FindBin(eta, pt) );
}
