#include "GenAnalyzer.h"


GenAnalyzer::GenAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    GenToken(CColl.consumes<GenEventInfoProduct>(PSet.getParameter<edm::InputTag>("genProduct"))),
    LheToken(CColl.consumes<LHEEventProduct>(PSet.getParameter<edm::InputTag>("lheProduct"))),
    GenParticlesToken(CColl.consumes<std::vector<reco::GenParticle> >(PSet.getParameter<edm::InputTag>("genParticles"))),
    ParticleList(PSet.getParameter<std::vector<int> >("pdgId")),
    ParticleStatus(PSet.getParameter<std::vector<int> >("status")),
    SampleDYJetsToLL(PSet.getParameter<std::vector<std::string> >("samplesDYJetsToLL")),
    SampleZJetsToNuNu(PSet.getParameter<std::vector<std::string> >("samplesZJetsToNuNu")),
    SampleWJetsToLNu(PSet.getParameter<std::vector<std::string> >("samplesWJetsToLNu")),
    SampleDir(PSet.getParameter<std::string>("samplesDir")),
    Sample(PSet.getParameter<std::string>("sample")),
    EWKFileName(PSet.getParameter<std::string>("ewkFile")),
    ApplyEWK(PSet.getParameter<bool>("applyEWK")),
    ApplyTopPtReweigth(PSet.getParameter<bool>("applyTopPtReweigth")),
    PythiaLOSample(PSet.getParameter<bool>("pythiaLOSample"))
{
    for(unsigned int i = 0; i < SampleDYJetsToLL.size(); i++) {
        Files[SampleDYJetsToLL[i]] = new TFile((SampleDir+SampleDYJetsToLL[i]+".root").c_str(), "READ");
        hPartons[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lhePartons");
        hBPartons[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lheBPartons");
        hHT[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lheHT");
        hPtV[SampleDYJetsToLL[i]] = (TH1F*)Files[SampleDYJetsToLL[i]]->Get("counter/c_lhePtZ");
    }
    for(unsigned int i = 0; i < SampleZJetsToNuNu.size(); i++) {
        Files[SampleZJetsToNuNu[i]] = new TFile((SampleDir+SampleZJetsToNuNu[i]+".root").c_str(), "READ");
        hPartons[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lhePartons");
        hBPartons[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lheBPartons");
        hHT[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lheHT");
        hPtV[SampleZJetsToNuNu[i]] = (TH1F*)Files[SampleZJetsToNuNu[i]]->Get("counter/c_lhePtZ");
    }
    for(unsigned int i = 0; i < SampleWJetsToLNu.size(); i++) {
        Files[SampleWJetsToLNu[i]] = new TFile((SampleDir+SampleWJetsToLNu[i]+".root").c_str(), "READ");
        hPartons[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lhePartons");
        hBPartons[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lheBPartons");
        hHT[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lheHT");
        hPtV[SampleWJetsToLNu[i]] = (TH1F*)Files[SampleWJetsToLNu[i]]->Get("counter/c_lhePtZ");
    }
    
    EWKFile = new TFile(EWKFileName.c_str(), "READ");
    
    fZEWK = (TF1*)EWKFile->Get("z_ewkcorr/z_ewkcorr_func");
    fWEWK = (TF1*)EWKFile->Get("w_ewkcorr/w_ewkcorr_func");
    
    /*
    Sample=sample;
    isDYFile=false;
    Sum=Num=NULL;
    // Read root file
    DYFile=new TFile("data/DYWeight.root", "READ");
    DYFile->cd();
    if(!DYFile->IsZombie()) {
      Num=(TH3F*)DYFile->Get(Sample.c_str());
      Sum=(TH3F*)DYFile->Get("Sum");
      if(Sum && Num && Sum->GetEntries()>0 && Num->GetEntries()>0) {
        npMax=Sum->GetXaxis()->GetBinCenter(Sum->GetNbinsX());
        ptMax=Sum->GetYaxis()->GetBinCenter(Sum->GetNbinsY());
        htMax=Sum->GetZaxis()->GetBinCenter(Sum->GetNbinsZ());
        isDYFile=true;
      }
      else std::cout << " - GenAnalyzer Warning: Drell-Yan initialization failed, check rootfile" << std::endl;
      
      
    }
    else std::cout << " - GenAnalyzer Warning: No Drell-Yan File" << std::endl;
    */
    
    // PU reweighting
//    LumiWeights=new edm::LumiReWeighting("data/MC_True.root", "data/Prod6.root", "S10", "pileup");
    
    std::cout << " --- GenAnalyzer initialization ---" << std::endl;
    std::cout << "  sample            :\t" << Sample << std::endl;
    if(ApplyEWK) std::cout << "  EWK file          :\t" << EWKFileName << std::endl;
    std::cout << std::endl;
}

GenAnalyzer::~GenAnalyzer() {
//    delete LumiWeights;
    for(auto const &it : Files) it.second->Close();
    EWKFile->Close();
}

// ---------- GEN WEIGHTS ----------

std::map<int, float> GenAnalyzer::FillWeightsMap(const edm::Event& iEvent) {
    std::map<int, float> Weights;
    Weights[-1] = 1.; // EventWeight
    if(iEvent.isRealData() or PythiaLOSample) return Weights;
    // Declare and open collection
    edm::Handle<GenEventInfoProduct> GenEventCollection;
    iEvent.getByToken(GenToken, GenEventCollection);
    // Declare and open collection
    edm::Handle<LHEEventProduct> LheEventCollection;
    iEvent.getByToken(LheToken, LheEventCollection);
    const LHEEventProduct* Product = LheEventCollection.product();
    
    Weights[-1] = GenEventCollection->weight()/fabs(GenEventCollection->weight());
     
    for(unsigned int i = 0; i < Product->weights().size(); i++) {
        Weights[ i ] = Product->weights()[i].wgt / Product->originalXWGTUP();
    }

    return Weights;
}


// ---------- GEN PARTICLES ----------

std::vector<reco::GenParticle> GenAnalyzer::FillGenVector(const edm::Event& iEvent) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
//        std::cout << it->pdgId() << "  " << it->status() << "  " << it->pt() << "  " << it->eta() << "  " << it->phi() << "  " << it->mass() << std::endl;
//        if(it->numberOfDaughters()>0) std::cout << "  " << it->daughter(0)->pdgId() << "  " << it->daughter(0)->status() << "  " << it->daughter(0)->pt() << std::endl;
//        if(it->numberOfDaughters()>1) std::cout << "  " << it->daughter(1)->pdgId() << "  " << it->daughter(1)->status() << "  " << it->daughter(1)->pt() << std::endl;
        for(unsigned int i = 0; i < ParticleList.size(); i++) {
            if(abs(it->pdgId()) == ParticleList[i]) Vect.push_back(*it); // Fill vector
        }
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdAndStatus(const edm::Event& iEvent, int partid, int partstatus) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      //for(unsigned int i = 0; i < ParticleList.size(); i++) {

      //for(unsigned int s = 0; s < ParticleStatus.size(); s++) {
            if(abs(it->pdgId()) == partid && (it->status()) == partstatus) Vect.push_back(*it); // Fill vector
	    //  }
	    //}
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdListAndStatusAndMotherList(const edm::Event& iEvent, std::vector<int> partidlist, int partstatus, std::vector<int> motheridlist) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {

      for(unsigned int s = 0; s < partidlist.size(); s++) {

	for(unsigned int m = 0; m < motheridlist.size(); m++) {

	    if(abs(it->pdgId()) == s && (it->status()) == partstatus && abs(it->mother()->pdgId()) == m) Vect.push_back(*it); // Fill vector

	}

      }
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}



std::vector<reco::GenParticle> GenAnalyzer::FillVBFGenVector(const edm::Event& iEvent){//, std::vector<int> partidlist, int partstatus, std::vector<int> motheridlist, std::vector<int> motherindexlist) {

    std::vector<reco::GenParticle> Vect;
    std::vector<reco::GenParticle> VBFmothers;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    // stop at fourth particle
    int iter = 0;
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      //std::cout << "Particle n. " << iter << std::endl;
      if(iter==2 or iter==3) 
	{
	  //std::cout << "VBF mother: pdgId " << it->pdgId() << " ; status "  << it->status() << " ; pz "<< it->pz() << std::endl;
	  VBFmothers.push_back(*it);
	}

      if(iter>3)
	{
	  break;
	}

      iter++;

    }

    
    std::vector<int> partidlist = {1,2,3,4,5,6,21};
    int partstatus = 23;
    //std::vector<int> motheridlist;

    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {


      for(unsigned int m = 0; m < VBFmothers.size(); m++) {

	for(unsigned int s = 0; s < partidlist.size(); s++) {
	  if(abs(it->pdgId()) == partidlist.at(s) && it->status() == partstatus && it->mother()->pdgId() == VBFmothers.at(m).pdgId() && it->mother()->p() == VBFmothers.at(m).p()) 
	    {
	      //std::cout << "Scattered quarks/gluons due to VBF: pdgId " << it->pdgId() << " ; status " << it->status() << " ; pt " << it->pt()  << " ; eta " << it->eta()  << " ; mother pz " << it->mother()->pz() << std::endl;
	      Vect.push_back(*it); // Fill vector
	    }
	}

      }


    }
    //    std::cout << "\n\n\n" << std::endl;
    return Vect;
}



std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdStatusAndMother(const edm::Event& iEvent, int partid, int partstatus, int motherid) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      if(abs(it->pdgId()) == partid && (it->status()) == partstatus && fabs(it->mother()->pdgId()) == motherid) Vect.push_back(*it); // Fill vector
    }
    return Vect;
}

//////////////////////////////////
std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdAndStatusAndKin(const edm::Event& iEvent, int partid, int partstatus, float pt, float eta) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
        if(abs(it->pdgId()) == partid && (it->status()) == partstatus && (it->pt())>pt && fabs(it->eta())<fabs(eta) ) Vect.push_back(*it); // Fill vector
    }
//    std::cout << "\n\n\n" << std::endl;
    return Vect;
}

std::vector<reco::GenParticle> GenAnalyzer::FillGenVectorByIdStatusAndMotherAndKin(const edm::Event& iEvent, int partid, int partstatus, int motherid, float pt, float eta) {

    std::vector<reco::GenParticle> Vect;

    // check if is real data
    isRealData = iEvent.isRealData();
    if(isRealData or PythiaLOSample) return Vect;
    // fill collection for this event 
    iEvent.getByToken(GenParticlesToken, GenCollection);
    // Loop on Gen Particles collection
    for(std::vector<reco::GenParticle>::const_iterator it = GenCollection->begin(); it != GenCollection->end(); ++it) {
      if(abs(it->pdgId()) == partid && (it->status()) == partstatus && fabs(it->mother()->pdgId()) == motherid && (it->pt())>pt && fabs(it->eta())<fabs(eta)) Vect.push_back(*it); // Fill vector
    }
    return Vect;
}

//////////////////////////////////

std::map<std::string, float> GenAnalyzer::FillLheMap(const edm::Event& iEvent) {
    std::map<std::string, float> Var;
    if(iEvent.isRealData() or PythiaLOSample) return Var;
    
    int lhePartons(0), lheBPartons(0);
    float lheHT(0.), lhePtZ(0.), lhePtW(0.), pt(0.);
    
    // Declare and open collection
    edm::Handle<LHEEventProduct> LheEventCollection;
    iEvent.getByToken(LheToken, LheEventCollection);
    const lhef::HEPEUP hepeup = LheEventCollection->hepeup();
    
    for(int i = 0; i < hepeup.NUP; ++i) {
        int id=abs(hepeup.IDUP[i]);
        // Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle entries in this event
        //reco::Candidate::LorentzVector P4(hepeup.PUP[i][0], hepeup.PUP[i][1], hepeup.PUP[i][2], hepeup.PUP[i][3]);
        pt = sqrt(hepeup.PUP[i][0]*hepeup.PUP[i][0] + hepeup.PUP[i][1]*hepeup.PUP[i][1]);
        if(hepeup.ISTUP[i]==1 && (id<6 || id==21)) {
            lheHT += pt; //P4.pt() 
            lhePartons++;
            if(id==5) lheBPartons++;
        }
        if(hepeup.ISTUP[i]==2 && abs(hepeup.IDUP[i])==23) lhePtZ = pt;
        if(hepeup.ISTUP[i]==2 && abs(hepeup.IDUP[i])==24) lhePtW = pt;
    }
    
    Var["lhePartons"] = lhePartons;
    Var["lheBPartons"] = lheBPartons;
    Var["lheHT"] = lheHT;
    Var["lhePtZ"] = lhePtZ;
    Var["lhePtW"] = lhePtW;
    return Var;
}

reco::Candidate* GenAnalyzer::FindGenParticle(std::vector<reco::GenParticle>& Vect, int pdg) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
        if(Vect[i].pdgId() == pdg) return FindLastDaughter(dynamic_cast<reco::Candidate*>(&Vect[i]));
    }
    return NULL;
}

// Recursive function to find the last particle in the chain before decay: e.g. 23 -> 23 -> *23* -> 13 -13. Returns a reco Candidate
reco::Candidate* GenAnalyzer::FindLastDaughter(reco::Candidate* p) {
    if(p->numberOfDaughters() <= 0 || !p->daughter(0)) return p;
    if(p->daughter(0)->pdgId() != p->pdgId()) return p;
    return FindLastDaughter(p->daughter(0));
}

reco::GenParticle* GenAnalyzer::FindGenParticleGenByIds(std::vector<reco::GenParticle>& Vect, std::vector<int> pdgs, int motherId) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
        for(unsigned int j=0; j<pdgs.size();j++) {
            if( Vect[i].pdgId() == pdgs[j] ) {
                if(motherId<=0) return FindLastDaughterGen(&Vect[i]);
                else {
                    if(Vect[i].mother() && Vect[i].mother()->pdgId() == motherId) return &Vect[i];
                }
            }
        }
    }
    return NULL;
}

// Recursive function to find the last particle in the chain before decay: e.g. 23 -> 23 -> *23* -> 13 -13.  Returns a GenParticle
reco::GenParticle* GenAnalyzer::FindLastDaughterGen(reco::GenParticle* p) {
    if(p->numberOfDaughters() <= 0 || !p->daughter(0)) return p;
    if(p->daughter(0)->pdgId() != p->pdgId()) return p;
    return FindLastDaughterGen( dynamic_cast<reco::GenParticle*>( p->daughter(0) ));
}

// Recursive function to find the last particle in the chain before decay: e.g. 23 -> 23 -> *23* -> 13 -13.  Returns a GenParticle
const reco::GenParticle* GenAnalyzer::FindLastDaughterGen(const reco::GenParticle* p) {
    if(p->numberOfDaughters() <= 0 || !p->daughter(0)) return p;
    if(p->daughter(0)->pdgId() != p->pdgId()) return p;
    return FindLastDaughterGen( dynamic_cast<const reco::GenParticle*>( p->daughter(0) ));
}

// Some particles radiate other particles with the same pdgId but a different status. Method to find a mother with different pdgId
const reco::Candidate* GenAnalyzer::FindMother(reco::GenParticle* p) {
    int pId = p->pdgId();
    const reco::Candidate* mom = p->mother();
    while (mom != 0 && mom->pdgId() == pId)
        mom = mom->mother();
    return mom;
}

// In Pythia 8, final state heavy bosons with kinematical info enclosed have status 62.
reco::Candidate* GenAnalyzer::FindGenParticleByIdAndStatus(std::vector<reco::GenParticle>& Vect, int pdg, int stat) {
    for(unsigned int i = 0; i < Vect.size(); i++) {
      if( (fabs(Vect[i].pdgId()) == pdg) && (Vect[i].status() == stat) ) return dynamic_cast<reco::Candidate*>(&Vect[i]);
    }
    return NULL;
}


// ---------- Pileup ----------

float GenAnalyzer::GetPUWeight(const edm::Event& iEvent) {
  //  int nPT(0);
  //  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchool2012PileupReweighting
  //  if(!iEvent.isRealData()) {
  //    edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
  //    iEvent.getByLabel(edm::InputTag("PileupSummaryInfo"), PUInfo);
  //    for(std::vector<PileupSummaryInfo>::const_iterator pvi=PUInfo->begin(), pvn=PUInfo->end(); pvi!=pvn; ++pvi) {
  //      if(pvi->getBunchCrossing()==0) nPT=pvi->getTrueNumInteractions(); // getPU_NumInteractions();
  //    }
  //    return LumiWeights->weight( nPT );
  //  }
    return 1.;
}


// ---------- PDF ----------

std::pair<float, float> GenAnalyzer::GetQ2Weight(const edm::Event& iEvent) {
  //  float Q, id1, id2, x1, x2;
  //  // Open LHE event
  //  edm::Handle<LHEEventProduct> lheProduct;
  //  iEvent.getByLabel(edm::InputTag("source"), lheProduct);
  //  // Access event info from LHE
  //  if(lheProduct.isValid()) {
  //    const lhef::HEPEUP hepeup=lheProduct->hepeup();
  //    // PDF
  //    Q   = hepeup.SCALUP;
  //    // id of the particle: 0 is for gluons
  //    id1 = hepeup.IDUP[0]==21 ? 0 : hepeup.IDUP[0];
  //    id2 = hepeup.IDUP[1]==21 ? 0 : hepeup.IDUP[1];
  //    x1  = fabs(hepeup.PUP[0][2]/6500.);
  //    x2  = fabs(hepeup.PUP[1][2]/6500.);
  //    //xfx1 = hepeup.XPDWUP.first;
  //    //xfx2 = hepeup.XPDWUP.second;
  //  }
  //  // Access event info from Gen if LHE info are not available
  //  else {
  //    edm::Handle<GenEventInfoProduct> genProduct;
  //    iEvent.getByLabel(edm::InputTag("generator"), genProduct);
  //    Q   = genProduct->pdf()->scalePDF;
  //    id1 = genProduct->pdf()->id.first;
  //    id2 = genProduct->pdf()->id.second;
  //    x1  = genProduct->pdf()->x.first;
  //    x2  = genProduct->pdf()->x.second;
  //    //pdf1 = genProduct->pdf()->xPDF.first;
  //    //pdf2 = genProduct->pdf()->xPDF.second;
  //  }
  //  //LHAPDF::usePDFMember(1, 0);
  //  float pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
  //  float pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
  //  // Q2 scale
  //  float newqcd1_up = LHAPDF::xfx(1, x1, Q*2, id1)/x1;
  //  float newqcd2_up = LHAPDF::xfx(1, x2, Q*2, id2)/x2;
  //  float newqcd1_down = LHAPDF::xfx(1, x1, Q/2, id1)/x1;
  //  float newqcd2_down = LHAPDF::xfx(1, x2, Q/2, id2)/x2;
  //  float Q2ScaleWeightUp = (newqcd1_up/pdf1)*(newqcd2_up/pdf2);
  //  float Q2ScaleWeightDown = (newqcd1_down/pdf1)*(newqcd2_down/pdf2);
    float Q2ScaleWeightUp = 1.;
    float Q2ScaleWeightDown = 1.;
    return std::pair<float, float>(Q2ScaleWeightUp, Q2ScaleWeightDown);
}

//float GenAnalyzer::GetPDFWeight(const edm::Event& iEvent) {
//  float Q, id1, id2, x1, x2;
//  // Open LHE event
//  edm::Handle<LHEEventProduct> lheProduct;
//  iEvent.getByLabel(edm::InputTag("source"), lheProduct);
//  // Access event info from LHE
//  if(lheProduct.isValid()) {
//    const lhef::HEPEUP hepeup=lheProduct->hepeup();
//    // PDF
//    Q   = hepeup.SCALUP;
//    // id of the particle: 0 is for gluons
//    id1 = hepeup.IDUP[0]==21 ? 0 : hepeup.IDUP[0];
//    id2 = hepeup.IDUP[1]==21 ? 0 : hepeup.IDUP[1];
//    x1  = fabs(hepeup.PUP[0][2]/6500.);
//    x2  = fabs(hepeup.PUP[1][2]/6500.);
//    //xfx1 = hepeup.XPDWUP.first;
//    //xfx2 = hepeup.XPDWUP.second;
//  }
//  // Access event info from Gen if LHE info are not available
//  else {
//    edm::Handle<GenEventInfoProduct> genProduct;
//    iEvent.getByLabel(edm::InputTag("generator"), genProduct);
//    Q   = genProduct->pdf()->scalePDF;
//    id1 = genProduct->pdf()->id.first;
//    id2 = genProduct->pdf()->id.second;
//    x1  = genProduct->pdf()->x.first;
//    x2  = genProduct->pdf()->x.second;
//    //pdf1 = genProduct->pdf()->xPDF.first;
//    //pdf2 = genProduct->pdf()->xPDF.second;
//  }
//  //LHAPDF::usePDFMember(1, 0);
//  float pdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
//  float pdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
//  // New PDF, if not working type <scramv1 setup lhapdffull> to enable more than one LHAPDF set
//  //LHAPDF::usePDFMember(2, 0);
//  float newpdf1 = LHAPDF::xfx(2, x1, Q, id1)/x1;
//  float newpdf2 = LHAPDF::xfx(2, x2, Q, id2)/x2;
//  return (newpdf1/pdf1)*(newpdf2/pdf2);
//}

// ---------- LHE Event ----------

float GenAnalyzer::GetStitchWeight(std::map<std::string, float> Map) {

    if(Map.size() <= 1) return 1.;
    if(Sample=="" || (SampleDYJetsToLL.size()==0 && SampleZJetsToNuNu.size()==0 && SampleWJetsToLNu.size()==0)) return 1.;
    if(Sample.find("JetsToLL") == std::string::npos && Sample.find("JetsToNuNu") == std::string::npos && Sample.find("JetsToLNu") == std::string::npos) return 1.;
    if(hPartons.find(Sample) == hPartons.end()) return 1.;
    
    // Find bins
    float StitchWeight(1.);
    int binPartons = hPartons[Sample]->FindBin(Map["lhePartons"]);
    int binBPartons = hBPartons[Sample]->FindBin(Map["lheBPartons"]);
    int binHT = hHT[Sample]->FindBin(Map["lheHT"]);
    int binPtZ = hPtV[Sample]->FindBin(Map["lhePtZ"] > 0. ? Map["lhePtZ"] : Map["lhePtW"]);
    
    // Calculate numerator and denominator
    if(Sample.find("JetsToLL") != std::string::npos) {
        float num(0.), den(0.);
        for(unsigned int i = 0; i < SampleDYJetsToLL.size(); i++) {
            den += hPartons[SampleDYJetsToLL[i]]->GetBinContent(binPartons);
            den += hBPartons[SampleDYJetsToLL[i]]->GetBinContent(binBPartons);
            den += hHT[SampleDYJetsToLL[i]]->GetBinContent(binHT);
            den += hPtV[SampleDYJetsToLL[i]]->GetBinContent(binPtZ);
        }
        num += hPartons[Sample]->GetBinContent(binPartons);
        num += hBPartons[Sample]->GetBinContent(binBPartons);
        num += hHT[Sample]->GetBinContent(binHT);
        num += hPtV[Sample]->GetBinContent(binPtZ);
        if(!(den != den || num != num) && den > 0.) StitchWeight = num / den;
    }
    return StitchWeight;
}

/// returning the topPt weight for each top/antitop
/// the total ttbar weight is implemented as sqrt(w_t1*w_t2)
float GenAnalyzer::GetTopPtWeight(float toppt) {
    if(!ApplyTopPtReweigth) return 1.;
    if(toppt <= 0) return 1.;
    return (toppt<400) ? sqrt(exp(0.0615-0.0005*toppt)) : sqrt(0.87066325126911626); // 0.87066325126911626 = exp(0.0615-0.0005*400)
}

float GenAnalyzer::GetZewkWeight(float zpt) {
    if(!ApplyEWK) return 1.;
    if(zpt <= 0) return 1.;
    return fZEWK->Eval(zpt);
}

float GenAnalyzer::GetWewkWeight(float wpt) {
    if(!ApplyEWK) return 1.;
    if(wpt <= 0) return 1.;
    return fWEWK->Eval(wpt);
}


// returns a vector with the gen particles from the decays of the particles in the list
std::vector<reco::GenParticle> GenAnalyzer::PartonsFromDecays(const std::vector<int> & pdgIds)  {
  
  // vector to save the partons
  std::vector<reco::GenParticle> partons;

  if(isRealData or PythiaLOSample) return partons;

  for (auto & genParticle : *GenCollection ) {
    // check if pdgid in vector
    if (std::find(pdgIds.begin(), pdgIds.end(), genParticle.pdgId()) != pdgIds.end()) {
      if (genParticle.numberOfDaughters() > 1) {
        if (genParticle.daughter(0) && genParticle.daughter(1)) {
          // emplace genParticles in output vector
          partons.emplace_back(*FindLastDaughterGen(dynamic_cast<const reco::GenParticle*>(genParticle.daughter(0))));
          partons.emplace_back(*FindLastDaughterGen(dynamic_cast<const reco::GenParticle*>(genParticle.daughter(1))));
        }
      }
    }
  }

  return partons;
}

// returns a vector with the gen particles from the decays of the particles in the list
// saves decaying particles in vector passed a reference
std::vector<reco::GenParticle> GenAnalyzer::PartonsFromDecays(const std::vector<int> & pdgIds, std::vector<reco::GenParticle> & genDecay)  {
  
  // vector to save the partons
  std::vector<reco::GenParticle> partons;

  if(isRealData or PythiaLOSample) return partons;

  for (auto & genParticle : *GenCollection ) {
    // check if pdgid in vector
    if (std::find(pdgIds.begin(), pdgIds.end(), genParticle.pdgId()) != pdgIds.end()) {
      if (genParticle.numberOfDaughters() > 1) {
        if (genParticle.daughter(0) && genParticle.daughter(1)) {
          genDecay.emplace_back(genParticle);
          // emplace genParticles in output vector
          partons.emplace_back(*FindLastDaughterGen(dynamic_cast<const reco::GenParticle*>(genParticle.daughter(0))));
          partons.emplace_back(*FindLastDaughterGen(dynamic_cast<const reco::GenParticle*>(genParticle.daughter(1))));
        }
      }
    }
  }

  return partons;
}

// returns a vector with the first n particles which have pdgIDs
std::vector<reco::GenParticle> GenAnalyzer::FirstNGenParticles(const std::vector<int> & pdgIds, std::size_t n)  {
  
  // vector to save the rtons
  std::vector<reco::GenParticle> particles;

  if(isRealData or PythiaLOSample) return particles;

  for (auto & genParticle : *GenCollection ) {
    // check if pdgid in vector
    if (std::find(pdgIds.begin(), pdgIds.end(), genParticle.pdgId()) != pdgIds.end()) {
          particles.emplace_back(genParticle);
    }
    if (particles.size() == n) return particles;
  }

  return particles;
}
