//#include "RecoilCorrector.hh" // From: https://github.com/cms-met/MetTools/tree/master/RecoilCorrections

//Important: for correcting calo jets: https://github.com/cms-jet/JMEDAS/blob/master/plugins/JetCorrectionsOnTheFly.cc


#include "CaloJetAnalyzer.h"


CaloJetAnalyzer::CaloJetAnalyzer(edm::ParameterSet& PSet, edm::ConsumesCollector&& CColl):
    CaloJetToken(CColl.consumes<std::vector<reco::CaloJet> >(PSet.getParameter<edm::InputTag>("jets"))),
    //CaloMetToken(CColl.consumes<std::vector<reco::PFMET> >(PSet.getParameter<edm::InputTag>("met"))),
    //QGToken(CColl.consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    //JetId(PSet.getParameter<int>("jetid")),
    Jet1Pt(PSet.getParameter<double>("jet1pt")),
    Jet2Pt(PSet.getParameter<double>("jet2pt")),
    JetEta(PSet.getParameter<double>("jeteta")),
    //AddQG(PSet.getParameter<bool>("addQGdiscriminator")),
    RecalibrateJets(PSet.getParameter<bool>("recalibrateJets")),
    RecalibrateMass(PSet.getParameter<bool>("recalibrateMass")),
    //RecalibratePuppiMass(PSet.getParameter<bool>("recalibratePuppiMass")),
    SmearJets(PSet.getParameter<bool>("smearJets")),
    JECUncertaintyMC(PSet.getParameter<std::string>("jecUncertaintyMC")),
    JECUncertaintyDATA(PSet.getParameter<std::string>("jecUncertaintyDATA")),
    JetCorrectorMC(PSet.getParameter<std::vector<std::string> >("jecCorrectorMC")),
    JetCorrectorDATA(PSet.getParameter<std::vector<std::string> >("jecCorrectorDATA")),
    MassCorrectorMC(PSet.getParameter<std::vector<std::string> >("massCorrectorMC")),
    MassCorrectorDATA(PSet.getParameter<std::vector<std::string> >("massCorrectorDATA")),
    //MassCorrectorPuppi(PSet.getParameter<std::string>("massCorrectorPuppi")),
    VertexToken(CColl.consumes<reco::VertexCollection>(PSet.getParameter<edm::InputTag>("vertices"))),
    RhoToken(CColl.consumes<double>(PSet.getParameter<edm::InputTag>("rho"))),
    //UseReshape(PSet.getParameter<bool>("reshapeBTag")),
    //BTag(PSet.getParameter<std::string>("btag")),
    //Jet1BTag(PSet.getParameter<int>("jet1btag")),
    //Jet2BTag(PSet.getParameter<int>("jet2btag")),
    //BTagDB(PSet.getParameter<std::string>("btagDB")),
    //UseRecoil(PSet.getParameter<bool>("metRecoil")),
    //RecoilMCFile(PSet.getParameter<std::string>("metRecoilMC")),
    //RecoilDataFile(PSet.getParameter<std::string>("metRecoilData")),
    //MetTriggerFileName(PSet.getParameter<std::string>("metTriggerFileName")),
    JerName_res(PSet.getParameter<std::string>("jerNameRes")),
    JerName_sf(PSet.getParameter<std::string>("jerNameSf"))
{
  
    jecUncMC = new JetCorrectionUncertainty(JECUncertaintyMC);
    jecUncDATA = new JetCorrectionUncertainty(JECUncertaintyDATA);

    //isMetTriggerFile = false;
    
    if(RecalibrateJets) {
        std::vector<JetCorrectorParameters> jetParMC;
        for ( std::vector<std::string>::const_iterator payloadBegin = JetCorrectorMC.begin(), payloadEnd = JetCorrectorMC.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
            //std::cout << *ipayload << "\n";
            jetParMC.push_back(JetCorrectorParameters(*ipayload));
        }    
        std::vector<JetCorrectorParameters> jetParDATA;
        for ( std::vector<std::string>::const_iterator payloadBegin = JetCorrectorDATA.begin(), payloadEnd = JetCorrectorDATA.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
            //std::cout << *ipayload << "\n";
            jetParDATA.push_back(JetCorrectorParameters(*ipayload));
        }
        // Make the FactorizedJetCorrector
        jetCorrMC = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(jetParMC) );
        jetCorrDATA = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(jetParDATA) );
    }
    
    if(RecalibrateMass) {
        std::vector<JetCorrectorParameters> massParMC;
        for ( std::vector<std::string>::const_iterator payloadBegin = MassCorrectorMC.begin(), payloadEnd = MassCorrectorMC.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
            massParMC.push_back(JetCorrectorParameters(*ipayload));
        }    
        std::vector<JetCorrectorParameters> massParDATA;
        for ( std::vector<std::string>::const_iterator payloadBegin = MassCorrectorDATA.begin(), payloadEnd = MassCorrectorDATA.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
            massParDATA.push_back(JetCorrectorParameters(*ipayload));
        }
        // Make the FactorizedJetCorrector
        massCorrMC = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(massParMC) );
        massCorrDATA = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(massParDATA) );
    }
    
    if(SmearJets) {
        resolution    = new JME::JetResolution(JerName_res);
        resolution_sf = new JME::JetResolutionScaleFactor(JerName_sf);
        if (JerName_res.find("AK8") != std::string::npos)
            Rparameter = 0.8;
        else 
            Rparameter = 0.4;
    }
    
    /*
    if(RecalibratePuppiMass) {
        PuppiCorrFile = new TFile(MassCorrectorPuppi.c_str(), "READ");
        PuppiJECcorr_gen = (TF1*)PuppiCorrFile->Get("puppiJECcorr_gen");
        PuppiJECcorr_reco_0eta1v3 = (TF1*)PuppiCorrFile->Get("puppiJECcorr_reco_0eta1v3");
        PuppiJECcorr_reco_1v3eta2v5 = (TF1*)PuppiCorrFile->Get("puppiJECcorr_reco_1v3eta2v5");
    }
    */

    // BTag calibrator
    /*
    if(UseReshape) {
        calib           = new BTagCalibration("CSVv2", BTagDB);
	
	// Set up readers for systematics. This code is largely thanks to Martino & Pablo in
	// https://github.com/cms-hh-pd/alp_analysis/blob/master/interface/BTagFilterOperator.h
	
	// Map of flavor type
	flavour_map = {{5, BTagEntry::FLAV_B},
		       {4, BTagEntry::FLAV_C},
		       {0, BTagEntry::FLAV_UDSG}};
	// Systematics to use for each flavor type
	syst_map = {{BTagEntry::FLAV_B, {"up_jes","down_jes",
					 "up_lf","down_lf",
					 "up_hfstats1", "down_hfstats1",
					 "up_hfstats2", "down_hfstats2"}},
                    {BTagEntry::FLAV_C, {"up_cferr1","down_cferr1",
                                         "up_cferr2", "down_cferr2"}},
                    {BTagEntry::FLAV_UDSG, {"up_jes","down_jes",
                                            "up_hf","down_hf",
                                            "up_lfstats1", "down_lfstats1",
					    "up_lfstats2", "down_lfstats2"}}};
	
	sf_mode = "iterativefit";
	
	// Load the reader with each systematic type.
	cr_map.emplace("central",
		       BTagCalibrationReader{BTagEntry::OP_RESHAPING,
			       "central", {}});
	for (const auto & kv : flavour_map) 
	    cr_map.at("central").load(*calib, kv.second, sf_mode);
	// for every flavour
	for (const auto & kv : syst_map) {
	    auto & syst_vector = kv.second;
	    // for every systematic relevant per flavour
	    for (const auto & syst : syst_vector) {
		auto it = cr_map.find(syst);
		if (it ==cr_map.end()) {
		    // return iterator as first pair element
		    it = cr_map.emplace(syst,
					BTagCalibrationReader{BTagEntry::OP_RESHAPING,
						syst, {}})
			.first;
		}
		// load calibration for this flavour and reader
		it->second.load(*calib, kv.first, sf_mode);
	    }
	}
    }
    */

    // Recoil Corrector
    /*
    if(UseRecoil) {
        recoilCorr = new RecoilCorrector(RecoilMCFile);
        recoilCorr->addDataFile(RecoilDataFile);
        recoilCorr->addMCFile(RecoilMCFile);
    }
    */

    /*
    MetTriggerFile=new TFile(MetTriggerFileName.c_str(), "READ");
    if(!MetTriggerFile->IsZombie()) {
        MetTriggerHisto=(TH1F*)MetTriggerFile->Get("SingleMuAll_numOR");
        isMetTriggerFile=true;
    }
    else {
        throw cms::Exception("CaloJetAnalyzer", "No Met Trigger File");
        return;
    }
    */

    std::cout << " --- CaloJetAnalyzer initialization ---" << std::endl;
    //std::cout << "  jet collection    :\t" << JetToken << std::endl;
    //std::cout << "  jet Id            :\t" << JetId << std::endl;
    std::cout << "  jet pT [1, 2]     :\t" << Jet1Pt << "\t" << Jet2Pt << std::endl;
    std::cout << "  jet eta           :\t" << JetEta << std::endl;
    //std::cout << "  b-tagging algo    :\t" << BTag << std::endl;
    //std::cout << "  b-tag cut [1, 2]  :\t" << Jet1BTag << "\t" << Jet2BTag << std::endl;
    //std::cout << "  apply recoil corr :\t" << (UseRecoil ? "YES" : "NO") << std::endl;
    //std::cout << "  recoil file MC    :\t" << RecoilMCFile << std::endl;
    //std::cout << "  recoil file Data  :\t" << RecoilDataFile << std::endl;
    std::cout << std::endl;
}

CaloJetAnalyzer::~CaloJetAnalyzer() {
  //if(RecalibratePuppiMass) PuppiCorrFile->Close();

  //if(UseReshape) {
  //delete calib;
  //}
    delete jecUncMC;
    delete jecUncDATA;
    //if(UseRecoil) delete recoilCorr;
    //MetTriggerFile->Close();
}





std::vector<reco::CaloJet> CaloJetAnalyzer::FillJetVector(const edm::Event& iEvent) {
    bool isMC(!iEvent.isRealData());
    //int BTagTh(Jet1BTag);
    float PtTh(Jet1Pt), EtaTh(JetEta);
    std::vector<reco::CaloJet> Vect;
    // Declare and open collection
    edm::Handle<std::vector<reco::CaloJet> > PFJetsCollection;
    iEvent.getByToken(CaloJetToken, PFJetsCollection);
    
    //// Open QG value maps
    //edm::Handle<edm::ValueMap<float>> QGHandle;
    //if(AddQG) iEvent.getByToken(QGToken, QGHandle);

    // Vertex collection
    edm::Handle<reco::VertexCollection> PVCollection;
    iEvent.getByToken(VertexToken, PVCollection);
    
    // Rho handle
    edm::Handle<double> rho_handle;
    iEvent.getByToken(RhoToken, rho_handle);
 
    // Loop on Jet collection
    for(std::vector<reco::CaloJet>::const_iterator it=PFJetsCollection->begin(); it!=PFJetsCollection->end(); ++it) {

        if(Vect.size()>0) {
            PtTh=Jet2Pt;
            //BTagTh=Jet2BTag;
        }
        reco::CaloJet jet=*it;
        int idx=it-PFJetsCollection->begin();
        //jet.addUserInt("Index", idx);
        reco::CaloJetRef jetRef(PFJetsCollection, idx);

	//First of all, jet id selections
        //// Quality cut
        //if(JetId==1 && !isLooseJet(jet)) continue;
        //if(JetId==2 && !isTightJet(jet)) continue;
        //if(JetId==3 && !isTightLepVetoJet(jet)) continue;
        //// b-tagging
        //if(BTagTh==1 && jet.bDiscriminator(BTag)<BTagTh) continue;
        //// Fill jet scale uncertainty
        //jet.addUserInt("isLoose", isLooseJet(jet) ? 1 : 0);
        //jet.addUserInt("isTight", isTightJet(jet) ? 1 : 0);
        //jet.addUserInt("isTightLepVeto", isTightLepVetoJet(jet) ? 1 : 0);

	
        if(RecalibrateJets) CorrectJet(jet, *rho_handle, PVCollection->size(), isMC);

        // JEC Uncertainty
        if (!isMC){
            jecUncDATA->setJetEta(jet.eta());
            jecUncDATA->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
            //jet.addUserFloat("JESUncertainty", jecUncDATA->getUncertainty(true));
        } else {
            jecUncMC->setJetEta(jet.eta());
            jecUncMC->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
            //jet.addUserFloat("JESUncertainty", jecUncMC->getUncertainty(true));
        }

	/*
	//std::cout << "JES uncertainty: " << jet.userFloat("JESUncertainty") <<std::endl;
        // PUPPI soft drop mass for AK8 jets
        if(jet.hasSubjets("SoftDropPuppi")) {
//            TLorentzVector puppiSoftdrop, puppiSoftdropSubjet;
//            auto const & sdSubjetsPuppi = jet.subjets("SoftDropPuppi");
//            for (auto const & it : sdSubjetsPuppi) {
//                puppiSoftdropSubjet.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(), it->mass());
//                puppiSoftdrop += puppiSoftdropSubjet;
//            }
            reco::Particle::LorentzVector puppiSoftdrop;
            for (auto const & it : jet.subjets("SoftDropPuppi")) puppiSoftdrop += it->correctedP4(0);
            jet.addUserFloat("ak8PFJetsPuppiSoftDropPt", puppiSoftdrop.pt());
            jet.addUserFloat("ak8PFJetsPuppiSoftDropEta", puppiSoftdrop.eta());
            jet.addUserFloat("ak8PFJetsPuppiSoftDropPhi", puppiSoftdrop.phi());
            jet.addUserFloat("ak8PFJetsPuppiSoftDropEnergy", puppiSoftdrop.energy());
            jet.addUserFloat("ak8PFJetsPuppiSoftDropMass", puppiSoftdrop.mass());
            
            float tau21 = jet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2")/jet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
            float ddt = tau21 + 0.063 * log( jet.userFloat("ak8PFJetsPuppiSoftDropMass")*jet.userFloat("ak8PFJetsPuppiSoftDropMass")/jet.userFloat("ak8PFJetsPuppiSoftDropPt") );
            jet.addUserFloat("ddtTau21", ddt);
        }
        
        if(RecalibrateMass) CorrectMass(jet, *rho_handle, PVCollection->size(), isMC);
        if(RecalibratePuppiMass) CorrectPuppiMass(jet, isMC);
        
	*/

	/*  
        // JER NEW IMPLEMENTATION
	
        if(SmearJets) {
            JME::JetParameters TheJetParameters;
            TheJetParameters.setJetPt(jet.pt());
            TheJetParameters.setJetEta(jet.eta());
            TheJetParameters.setRho(*rho_handle);

            float smearFactor = 1.;
            float smearFactorUp = 1.;
            float smearFactorDown = 1.;

            if(isMC) {
                float JERresolution = resolution->getResolution(TheJetParameters);
                float JERsf         = resolution_sf->getScaleFactor(TheJetParameters);
                float JERsfUp       = resolution_sf->getScaleFactor(TheJetParameters, Variation::UP);
                float JERsfDown     = resolution_sf->getScaleFactor(TheJetParameters, Variation::DOWN);
                //std::cout << "JERresolution " << JERresolution << "\n";
                //std::cout << "JERsf         " << JERsf << "\n";
                //std::cout << "JERsfUp       " << JERsfUp << "\n";
                //std::cout << "JERsfDown     " << JERsfDown << "\n";
                const reco::GenJet* genJet=jet.genJet();
                if(genJet) {
                    if ( ( sqrt( pow(jet.eta() - genJet->eta(),2) + pow(jet.phi() - genJet->phi(),2) ) < 0.5*Rparameter )  &&
                         fabs( jet.pt() - genJet->pt()) < 3.*JERresolution*jet.pt() ) { // (DeltaR < R/2) AND (DeltaPt < 3*PtRes)
                        smearFactor = max(0.,genJet->pt()+JERsf*(jet.pt() - genJet->pt()))/jet.pt();
                        smearFactorUp = max(0.,genJet->pt()+JERsfUp*(jet.pt() - genJet->pt()))/jet.pt();
                        smearFactorDown = max(0.,genJet->pt()+JERsfDown*(jet.pt() - genJet->pt()))/jet.pt();
                    }  
                    else {
                        TRandom3 rnd(0);
                        smearFactor = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsf*JERsf-1.)));
                        smearFactorUp = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsfUp*JERsfUp-1.)));
                        smearFactorDown = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsfDown*JERsfDown-1.)));
                    }
                }
		else {
		    TRandom3 rnd(0);
		    smearFactor = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsf*JERsf-1.)));
		    smearFactorUp = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsfUp*JERsfUp-1.)));
		    smearFactorDown = 1. + rnd.Gaus(0.,JERresolution*sqrt(max(0.,JERsfDown*JERsfDown-1.)));
                }
            }        
            //std::cout << "Rparameter      " << Rparameter << "\n";
            //std::cout << "smearFactor     " << smearFactor << "\n";
            //std::cout << "smearFactorUp   " << smearFactorUp << "\n";
            //std::cout << "smearFactorDown " << smearFactorDown << "\n";
	    reco::CaloJet jetJERUp = jet;
	    reco::CaloJet jetJERDown = jet;
            jet.setP4(jet.p4() * smearFactor);
            jetJERUp.setP4(jet.p4() * smearFactorUp);
            jetJERDown.setP4(jet.p4() * smearFactorDown);
            jet.addUserFloat("ptJERUp", jetJERUp.pt());
            jet.addUserFloat("etaJERUp", jetJERUp.eta());
            jet.addUserFloat("phiJERUp", jetJERUp.phi());
            jet.addUserFloat("energyJERUp", jetJERUp.energy());
            jet.addUserFloat("ptJERDown", jetJERDown.pt());
            jet.addUserFloat("etaJERDown", jetJERDown.eta());
            jet.addUserFloat("phiJERDown", jetJERDown.phi());
            jet.addUserFloat("energyJERDown", jetJERDown.energy());

            //jet.addUserFloat("smearFactor", smearFactor);
            //jet.addUserFloat("smearFactorUp", smearFactorUp);
            //jet.addUserFloat("smearFactorDown", smearFactorDown);           
        }        
        // JER NEW IMPLEMENTATION        

*/        

        // Pt and eta cut
        if(jet.pt()<PtTh || fabs(jet.eta())>EtaTh) continue;
	/*
	std::vector<float> reshapedDiscriminator = ReshapeBtagDiscriminator(jet);
        jet.addUserFloat("ReshapedDiscriminator", reshapedDiscriminator[0]);
        jet.addUserFloat("ReshapedDiscriminatorUp", reshapedDiscriminator[1]);
        jet.addUserFloat("ReshapedDiscriminatorDown", reshapedDiscriminator[2]);
	*/

        // CSV reshaping for soft drop subjets
	/*
        if(jet.hasSubjets("SoftDrop")) {
            auto const & sdSubjets = jet.subjets("SoftDrop");
            short nsj = 1;
            for (auto const & it : sdSubjets) {
                reco::CaloJet subjet = it;
		std::vector<float> reshapedDiscriminatorSubjet = ReshapeBtagDiscriminator(subjet);
                jet.addUserFloat(Form("ReshapedDiscriminator%d",nsj), reshapedDiscriminatorSubjet[0]);
                jet.addUserFloat(Form("ReshapedDiscriminatorUp%d",nsj), reshapedDiscriminatorSubjet[1]);
                jet.addUserFloat(Form("ReshapedDiscriminatorDown%d",nsj), reshapedDiscriminatorSubjet[2]);
                ++nsj;
            }
        }
        
        //QG tagger for AK4 jets
        if(AddQG && jet.nSubjetCollections()<=0) {
            jet.addUserFloat("QGLikelihood", (*QGHandle)[jetRef]);
        }
        */

        Vect.push_back(jet); // Fill vector
    }
    return Vect;
}


/////////////////////////////////////////////
void CaloJetAnalyzer::CorrectJet(reco::CaloJet& jet, float rho, float nPV, bool isMC) {
    double corr(1.);
    reco::Candidate::LorentzVector uncorrJet = jet.p4();//pat::Jet has correctedP4(0); calo do not; see: https://github.com/cms-jet/JMEDAS/blob/master/plugins/JetCorrectionsOnTheFly.cc
    
    if(!isMC) {
        jetCorrDATA->setJetEta( uncorrJet.Eta() );
        jetCorrDATA->setJetPt ( uncorrJet.Pt() );
        jetCorrDATA->setJetE  ( uncorrJet.E() );
        jetCorrDATA->setJetA  ( jet.jetArea() );
        jetCorrDATA->setRho   ( rho );
        jetCorrDATA->setNPV   ( nPV );
        corr = jetCorrDATA->getCorrection();
    }
    else {
        jetCorrMC->setJetEta( uncorrJet.Eta() );
        jetCorrMC->setJetPt ( uncorrJet.Pt() );
        jetCorrMC->setJetE  ( uncorrJet.E() );
        jetCorrMC->setJetA  ( jet.jetArea() );
        jetCorrMC->setRho   ( rho );
        jetCorrMC->setNPV   ( nPV );
        corr = jetCorrMC->getCorrection();
    }
    
    reco::Candidate::LorentzVector corrJet(uncorrJet);
    jet.setP4(corrJet * corr);
}

/*
void CaloJetAnalyzer::CorrectMass(reco::CaloJet& jet, float rho, float nPV, bool isMC) {
    double corr(1.);
    reco::Candidate::LorentzVector uncorrJet = jet.correctedP4(0);
    
    if(!isMC) {
        massCorrDATA->setJetEta( uncorrJet.Eta() );
        massCorrDATA->setJetPt ( uncorrJet.Pt() );
        massCorrDATA->setJetE  ( uncorrJet.E() );
        massCorrDATA->setJetA  ( jet.jetArea() );
        massCorrDATA->setRho   ( rho );
        massCorrDATA->setNPV   ( nPV );
        corr = massCorrDATA->getCorrection();
    }
    else {
        massCorrMC->setJetEta( uncorrJet.Eta() );
        massCorrMC->setJetPt ( uncorrJet.Pt() );
        massCorrMC->setJetE  ( uncorrJet.E() );
        massCorrMC->setJetA  ( jet.jetArea() );
        massCorrMC->setRho   ( rho );
        massCorrMC->setNPV   ( nPV );
        corr = massCorrMC->getCorrection();
    }
    if(jet.hasUserFloat("ak8PFJetsCHSPrunedMass")) jet.addUserFloat("ak8PFJetsCHSPrunedMassCorr", jet.userFloat("ak8PFJetsCHSPrunedMass") * corr);
    if(jet.hasUserFloat("ak8PFJetsCHSSoftDropMass")) jet.addUserFloat("ak8PFJetsCHSSoftDropMassCorr", jet.userFloat("ak8PFJetsCHSSoftDropMass") * corr);
    if(jet.hasUserFloat("ak8PFJetsPrunedMass")) jet.addUserFloat("ak8PFJetsPrunedMassCorr", jet.userFloat("ak8PFJetsPrunedMass") * corr);
    if(jet.hasUserFloat("ak8PFJetsSoftDropMass")) jet.addUserFloat("ak8PFJetsSoftDropMassCorr", jet.userFloat("ak8PFJetsSoftDropMass") * corr);
    //if(jet.hasUserFloat("ak8PFJetsPuppiSoftDropMass")) jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorr", jet.userFloat("ak8PFJetsPuppiSoftDropMass") * corr);
}
*/


/*
void CaloJetAnalyzer::CorrectPuppiMass(reco::CaloJet& jet, bool isMC) {
    bool hasInfo( jet.hasUserFloat("ak8PFJetsPuppiSoftDropMass") && jet.hasUserFloat("ak8PFJetsPuppiSoftDropPt") && jet.hasUserFloat("ak8PFJetsPuppiSoftDropEta") );
    float corr(1.), genCorr(1.), recoCorr(1.);
    if(hasInfo && jet.userFloat("ak8PFJetsPuppiSoftDropMass") > 0.) {
        genCorr = PuppiJECcorr_gen->Eval( jet.userFloat("ak8PFJetsPuppiSoftDropPt") );
        if(fabs(jet.userFloat("ak8PFJetsPuppiSoftDropEta")) <= 1.3) 
            recoCorr = PuppiJECcorr_reco_0eta1v3->Eval( jet.userFloat("ak8PFJetsPuppiSoftDropPt") );
        else if(fabs(jet.userFloat("ak8PFJetsPuppiSoftDropEta")) > 1.3 ) 
            recoCorr = PuppiJECcorr_reco_1v3eta2v5->Eval( jet.userFloat("ak8PFJetsPuppiSoftDropPt") );
        corr = genCorr * recoCorr;
    }
    if(corr < 0.) corr = 0.;
    jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorr", jet.userFloat("ak8PFJetsPuppiSoftDropMass") * corr);
    jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrNotSmeared", jet.userFloat("ak8PFJetsPuppiSoftDropMass") * corr);

    if(isMC){
      float JMSSf  = 1.;//Moriond17
        float JMSUnc = 0.0094;//Moriond17
        float JESUnc = jet.userFloat("JESUncertainty");    
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMS", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")       * JMSSf);
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSUp", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")     * (JMSSf + sqrt(JMSUnc*JMSUnc + JESUnc*JESUnc) ) );
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMSDown", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")   * (JMSSf - sqrt(JMSUnc*JMSUnc + JESUnc*JESUnc) ) );
        
        float JMRSf   = 1.;//Moriond17
        float JMRUnc  = 0.20;//Moriond17
        TRandom3 rnd(0);
        float smearJMR    = rnd.Gaus(1.,JMRSf-1.);
        float smearJMRUp    = rnd.Gaus(1.,(JMRSf-1.)*(1. + JMRUnc/JMRSf));
	float smearJMRDown    = rnd.Gaus(1.,(JMRSf -1.)*(1. - JMRUnc/JMRSf));
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMR", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")       * smearJMR);
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRUp", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")     * smearJMRUp);
        jet.addUserFloat("ak8PFJetsPuppiSoftDropMassCorrJMRDown", jet.userFloat("ak8PFJetsPuppiSoftDropMassCorr")   * smearJMRDown);
    }
}
*/

void CaloJetAnalyzer::CleanJetsFromMuons(std::vector<reco::CaloJet>& Jets, std::vector<pat::Muon>& Muons, float angle) {
    for(unsigned int m = 0; m < Muons.size(); m++) {
        for(unsigned int j = 0; j < Jets.size(); ) {
            if(deltaR(Jets[j], Muons[m]) < angle) Jets.erase(Jets.begin() + j);
            else j++;
        }
    }
}

void CaloJetAnalyzer::CleanJetsFromElectrons(std::vector<reco::CaloJet>& Jets, std::vector<pat::Electron>& Electrons, float angle) {
    for(unsigned int e = 0; e < Electrons.size(); e++) {
        for(unsigned int j = 0; j < Jets.size(); ) {
            if(deltaR(Jets[j], Electrons[e]) < angle) Jets.erase(Jets.begin() + j);
            else j++;
        }
    }
}

/*
void CaloJetAnalyzer::AddVariables(std::vector<reco::CaloJet>& Jets, reco::PFMET& MET) {
    for(unsigned int j = 0; j < Jets.size(); j++) {
        Jets[j].addUserFloat("dPhi_met", fabs(reco::deltaPhi(Jets[j].phi(), MET.phi())));
        Jets[j].addUserFloat("dPhi_Jet1", fabs(reco::deltaPhi(Jets[j].phi(), Jets[0].phi())));
    }
}
*/

/*
void CaloJetAnalyzer::GenMatcherCalo(std::vector<reco::CaloJet>& Jets, std::vector<reco::GenParticle>& Quarks, std::string label) {
    for(unsigned int j = 0; j < Jets.size(); j++){
        for(unsigned int q = 0; q < Quarks.size(); q++) {
	  //std::cout << "jet: " << j << " quark: " << q << " delta R: " <<fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi()) ) << std::endl;
	  //std::cout << ("dR_q"+std::to_string(q)).c_str() << std::endl;
	  //std::cout <<  ("dR_"+label+std::to_string(q+1)).c_str() << std::endl;
	  Jets[j].addUserFloat(("dR_"+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi())) );

	  if(Jets[j].hasSubjets("SoftDrop"))
	     {
	       if(Jets[j].subjets("SoftDrop").size() > 0) Jets[j].addUserFloat(("dR_"+label+std::to_string(q+1)+"_sj1").c_str(), fabs(reco::deltaR(Jets[j].subjets("SoftDrop")[0]->eta(),Jets[j].subjets("SoftDrop")[0]->phi(),Quarks[q].eta(),Quarks[q].phi())) );
	       if(Jets[j].subjets("SoftDrop").size() > 1) Jets[j].addUserFloat(("dR_"+label+std::to_string(q+1)+"_sj2").c_str(), fabs(reco::deltaR(Jets[j].subjets("SoftDrop")[1]->eta(),Jets[j].subjets("SoftDrop")[1]->phi(),Quarks[q].eta(),Quarks[q].phi())) );
	     }

            //Jets[j].addUserFloat("quark_index", fabs(reco::deltaPhi(Jets[j].phi(), Jets[0].phi())));
	}
    }
}
*/


/* nice prototype, do not delete!*/
std::map<std::string,float> CaloJetAnalyzer::GenMatcherCalo(std::vector<reco::CaloJet>& Jets, std::vector<reco::GenParticle>& Quarks, std::string label) {
    std::map<std::string,float> match_map;
    for(unsigned int j = 0; j < Jets.size(); j++) {
        for(unsigned int q = 0; q < Quarks.size(); q++) {
	  //std::cout << "jet: " << j << " quark: " << q << " delta R: " <<fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi()) ) << std::endl;
	  //std::cout << ("dR_q"+std::to_string(q)).c_str() << std::endl;
	  //std::cout <<  ("dR_"+label+std::to_string(q+1)).c_str() << std::endl;
	  match_map.insert(std::make_pair(("dR_j"+std::to_string(j+1)+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi())) ));
	  //Jets[j].addUserFloat(("dR_"+label+std::to_string(q+1)).c_str(), fabs(reco::deltaR(Jets[j].eta(),Jets[j].phi(),Quarks[q].eta(),Quarks[q].phi())) );
            //Jets[j].addUserFloat("quark_index", fabs(reco::deltaPhi(Jets[j].phi(), Jets[0].phi())));
        }
    }
    return match_map;
}


/*
int CaloJetAnalyzer::GetNBJets(std::vector<reco::CaloJet>& Jets) {
    int n(0);
    for(unsigned int i = 0; i < Jets.size(); i++) if(abs(Jets[i].hadronFlavour()) == 5) n++;
    return n;
}
*/

/*
float CaloJetAnalyzer::GetMetTriggerEfficiency(reco::PFMET& MET) {
    if(!isMetTriggerFile) return 1.;
    double pt = std::min( std::max( MetTriggerHisto->GetXaxis()->GetXmin(), MET.pt() ) , MetTriggerHisto->GetXaxis()->GetXmax() - 0.000001 );
    return(MetTriggerHisto->Interpolate(pt));
}
*/

/*
reco::PFMET CaloJetAnalyzer::FillMetVector(const edm::Event& iEvent) {
    
    edm::Handle<std::vector<reco::PFMET> > MetCollection;
    iEvent.getByToken(CaloMetToken, MetCollection);
    reco::PFMET MEt = MetCollection->front();
    //MEt.addUserFloat("ptShiftJetResUp", MEt.shiftedPt(reco::PFMET::METUncertainty::JetResUp));
    //MEt.addUserFloat("ptShiftJetResDown", MEt.shiftedPt(reco::PFMET::METUncertainty::JetResDown));
    //MEt.addUserFloat("ptShiftJetEnUp", MEt.shiftedPt(reco::PFMET::METUncertainty::JetEnUp));
    //MEt.addUserFloat("ptShiftJetEnDown", MEt.shiftedPt(reco::PFMET::METUncertainty::JetEnDown));
    //MEt.addUserFloat("ptShiftUnclusteredEnUp", MEt.shiftedPt(reco::PFMET::METUncertainty::UnclusteredEnUp));
    //MEt.addUserFloat("ptShiftUnclusteredEnDown", MEt.shiftedPt(reco::PFMET::METUncertainty::UnclusteredEnDown));
    //MEt.addUserFloat("ptRaw", MEt.uncorPt());
    //MEt.addUserFloat("phiRaw", MEt.uncorPhi());
    return MEt;
}
*/

/*
void CaloJetAnalyzer::ApplyRecoilCorrections(reco::PFMET& MET, const reco::Candidate::LorentzVector* GenV, const reco::Candidate::LorentzVector* RecoV, int nJets) {
    double MetPt(MET.pt()), MetPhi(MET.phi()), MetPtScaleUp(MET.pt()), MetPhiScaleUp(MET.phi()), MetPtScaleDown(MET.pt()), MetPhiScaleDown(MET.phi()), MetPtResUp(MET.pt()), MetPhiResUp(MET.phi()), MetPtResDown(MET.pt()), MetPhiResDown(MET.phi());
    double GenPt(0.), GenPhi(0.), LepPt(0.), LepPhi(0.), LepPx(0.), LepPy(0.);
    double RecoilX(0.), RecoilY(0.), Upara(0.), Uperp(0.);
    
    if(GenV) {
        GenPt = GenV->pt();
        GenPhi = GenV->phi();
    }
    else {
        throw cms::Exception("CaloJetAnalyzer", "GenV boson is null. No Recoil Correction can be derived");
        return;
    }
    
    if(RecoV) {
        LepPt = RecoV->pt();
        LepPhi = RecoV->phi();
        LepPx = RecoV->px();
        LepPy = RecoV->py();
        RecoilX = - MET.px() - LepPx;
        RecoilY = - MET.py() - LepPy;
        Upara = (RecoilX*LepPx + RecoilY*LepPy) / LepPt;
        Uperp = (RecoilX*LepPy - RecoilY*LepPx) / LepPt;
    }
    
    // Apply Recoil Corrections
    if(UseRecoil) {
        recoilCorr->CorrectType2(MetPt,          MetPhi,          GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0,  0, nJets);
        recoilCorr->CorrectType2(MetPtScaleUp,   MetPhiScaleUp,   GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  3,  0, nJets);
        recoilCorr->CorrectType2(MetPtScaleDown, MetPhiScaleDown, GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp, -3,  0, nJets);
        recoilCorr->CorrectType2(MetPtResUp,     MetPhiResUp,     GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0,  3, nJets);
        recoilCorr->CorrectType2(MetPtResDown,   MetPhiResDown,   GenPt, GenPhi, LepPt, LepPhi, Upara, Uperp,  0, -3, nJets);
    }
    
    // Set userFloats for systematics
    MET.addUserFloat("ptScaleUp", MetPtScaleUp);
    MET.addUserFloat("ptScaleDown", MetPtScaleDown);
    MET.addUserFloat("ptResUp", MetPtResUp);
    MET.addUserFloat("ptResDown", MetPtResDown);
    
    // Set new P4
    MET.setP4(reco::Candidate::PolarLorentzVector(MetPt, MET.eta(), MetPhi, MET.mass()));
}
*/

// // https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
// float CaloJetAnalyzer::GetResolutionRatio(float eta) {
//     eta=fabs(eta);
//     if(eta>=0.0 && eta<0.5) return 1.122; 
//     if(eta>=0.5 && eta<0.8) return 1.167;
//     if(eta>=0.8 && eta<1.1) return 1.168;
//     if(eta>=1.1 && eta<1.3) return 1.029;
//     if(eta>=1.3 && eta<1.7) return 1.115;
//     if(eta>=1.7 && eta<1.9) return 1.041;
//     if(eta>=1.9 && eta<2.1) return 1.167;
//     if(eta>=2.1 && eta<2.3) return 1.094;
//     if(eta>=2.3 && eta<2.5) return 1.168;
//     if(eta>=2.5 && eta<2.8) return 1.266;
//     if(eta>=2.8 && eta<3.0) return 1.595;
//     if(eta>=3.0 && eta<3.2) return 0.998;
//     if(eta>=3.2 && eta<5.0) return 1.226;
//     return -1.;
// }
// float CaloJetAnalyzer::GetResolutionErrorUp(float eta) {
//     eta=fabs(eta);
//     if(eta>=0.0 && eta<0.5) return 1.122 + 0.026; 
//     if(eta>=0.5 && eta<0.8) return 1.167 + 0.048;
//     if(eta>=0.8 && eta<1.1) return 1.168 + 0.046;
//     if(eta>=1.1 && eta<1.3) return 1.029 + 0.066;
//     if(eta>=1.3 && eta<1.7) return 1.115 + 0.030;
//     if(eta>=1.7 && eta<1.9) return 1.041 + 0.062;
//     if(eta>=1.9 && eta<2.1) return 1.167 + 0.086;
//     if(eta>=2.1 && eta<2.3) return 1.094 + 0.093;
//     if(eta>=2.3 && eta<2.5) return 1.168 + 0.120;
//     if(eta>=2.5 && eta<2.8) return 1.266 + 0.132;
//     if(eta>=2.8 && eta<3.0) return 1.595 + 0.175;
//     if(eta>=3.0 && eta<3.2) return 0.998 + 0.066;
//     if(eta>=3.2 && eta<5.0) return 1.226 + 0.145;
//     return -1.;
// }
// float CaloJetAnalyzer::GetResolutionErrorDown(float eta) {
//     eta=fabs(eta);
//     if(eta>=0.0 && eta<0.5) return 1.122 - 0.026; 
//     if(eta>=0.5 && eta<0.8) return 1.167 - 0.048;
//     if(eta>=0.8 && eta<1.1) return 1.168 - 0.046;
//     if(eta>=1.1 && eta<1.3) return 1.029 - 0.066;
//     if(eta>=1.3 && eta<1.7) return 1.115 - 0.030;
//     if(eta>=1.7 && eta<1.9) return 1.041 - 0.062;
//     if(eta>=1.9 && eta<2.1) return 1.167 - 0.086;
//     if(eta>=2.1 && eta<2.3) return 1.094 - 0.093;
//     if(eta>=2.3 && eta<2.5) return 1.168 - 0.120;
//     if(eta>=2.5 && eta<2.8) return 1.266 - 0.132;
//     if(eta>=2.8 && eta<3.0) return 1.595 - 0.175;
//     if(eta>=3.0 && eta<3.2) return 0.998 - 0.066;
//     if(eta>=3.2 && eta<5.0) return 1.226 - 0.145;
//     return -1.;
// }

/*
// PFJet Quality ID 2015-2016: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
bool CaloJetAnalyzer::isLooseJet(reco::CaloJet& jet) {
    if(fabs(jet.eta())<=2.7){ /// |eta| < 2.7
        if(jet.neutralHadronEnergyFraction()>=0.99) return false;
        if(jet.neutralEmEnergyFraction()>=0.99) return false;
        if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
        if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedMultiplicity()<=0) return false;
            if(jet.chargedEmEnergyFraction()>=0.99) return false;
        }
    }
    else{ /// |eta| > 2.7
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
            if(jet.neutralMultiplicity()<=2) return false;
        }
        else{ /// |eta| > 3.0
            if(jet.neutralMultiplicity()<=10) return false;
        }
    }
    return true;
}

bool CaloJetAnalyzer::isTightJet(reco::CaloJet& jet) {
    if(fabs(jet.eta())<=2.7){ /// |eta| < 2.7
        if(jet.neutralHadronEnergyFraction()>=0.90) return false;
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
        if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedMultiplicity()<=0) return false;
            if(jet.chargedEmEnergyFraction()>=0.99) return false;
        }
    }
    else{ /// |eta| > 2.7
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
            if(jet.neutralMultiplicity()<=2) return false;
        }
        else{ /// |eta| > 3.0
            if(jet.neutralMultiplicity()<=10) return false;
        }
    }
    return true;
}

bool CaloJetAnalyzer::isTightLepVetoJet(reco::CaloJet& jet) {
    if(fabs(jet.eta())<=2.7){ /// |eta| < 2.7
        if(jet.neutralHadronEnergyFraction()>=0.90) return false;
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
        if(jet.muonEnergyFraction()>=0.80) return false;
        if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
            if(jet.chargedHadronEnergyFraction()<=0.) return false;
            if(jet.chargedMultiplicity()<=0) return false;
            if(jet.chargedEmEnergyFraction()>=0.99) return false;
        }
    }
    else{ /// |eta| > 2.7
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
            if(jet.neutralMultiplicity()<=2) return false;
        }
        else{ /// |eta| > 3.0
            if(jet.neutralMultiplicity()<=10) return false;
        }
    }
    return true;
}


std::vector<float> CaloJetAnalyzer::ReshapeBtagDiscriminator(reco::CaloJet& jet) {
    float pt(jet.pt()), eta(fabs(jet.eta())), discr(jet.bDiscriminator(BTag));
    int hadronFlavour_ = std::abs(jet.hadronFlavour());
    std::vector<float> reshapedDiscr(3, discr);
    
    if(UseReshape) {
        BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
        if (hadronFlavour_ == 5) jf = BTagEntry::FLAV_B;
	else if (hadronFlavour_ == 4) jf = BTagEntry::FLAV_C;
	else if (hadronFlavour_ == 0) jf = BTagEntry::FLAV_UDSG;

	auto central_sf = cr_map.at("central").eval_auto_bounds("central", jf, eta, pt, discr);
	// default to 1 rather than 0 if out of bounds
	if (central_sf == 0) central_sf = 1.0;
	
	// Get the systematic shifts. For the time being just add up the differences from the central
	// value in quadrature and take that as the overall systematic.
	float up_total2 = 0;
	float down_total2 = 0;
	for (const auto & syst : syst_map.at(jf)) {
	    auto syst_sf = cr_map.at(syst).eval_auto_bounds(syst, jf, eta, pt, discr);
	    // default to 1, as above
	    if (syst_sf == 0) syst_sf = 1.0;

	    if (syst.find("up") != std::string::npos) {
		up_total2 += (syst_sf-central_sf)*(syst_sf-central_sf);
	    } else if (syst.find("down") != std::string::npos) {
		down_total2 += (syst_sf-central_sf)*(syst_sf-central_sf);
	    } else {
		std::cerr << "Unknown systematic " << syst << " -- don't know if this is up or down!" << std::endl;
	    }
	}
	float up_sf = central_sf - sqrt(up_total2);
	float down_sf = central_sf + sqrt(down_total2);

	reshapedDiscr[0] = discr*central_sf;
	reshapedDiscr[1] = discr*up_sf;
	reshapedDiscr[2] = discr*down_sf;
	  
	//std::cout << Form("pt, eta, b-tag, flav, reshapedDiscr : %f, %f, %f, %d, %f, %f, %f\n",
	// 		  pt, eta, discr, jf, reshapedDiscr[0], reshapedDiscr[1], reshapedDiscr[2]);
    }
    return reshapedDiscr;
}
*/
