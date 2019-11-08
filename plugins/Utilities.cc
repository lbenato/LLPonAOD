
#include "Utilities.h"

// --------------------   ANGULAR  -------------------

float Utilities::ReturnCosThetaStar(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theZ) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pZ(theZ.px(), theZ.py(), theZ.pz(), theZ.energy());
  // Boost the Z to the A rest frame
  pZ.Boost( -pA.BoostVector() );
  float value=pZ.CosTheta();
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnCosTheta1(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Reconstruct H in A rest frame
  TLorentzVector pHr = pB1 + pB2;
  // cos theta = H dot L1 / (|H|*|L1|)
  float value=pHr.Vect().Dot( pL1.Vect() ) / ( pHr.Vect().Mag()*pL1.Vect().Mag() );
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnCosTheta2(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Reconstruct Z in A rest frame
  TLorentzVector pZr = pL1 + pL2;
  // cos theta = Z dot B1 / (|Z|*|B1|)
  float value=pZr.Vect().Dot( pB1.Vect() ) / ( pZr.Vect().Mag()*pB1.Vect().Mag() );
  if(value!=value || isinf(value)) return 0.;
  return value;
}

float Utilities::ReturnPhi(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2, const reco::Candidate::LorentzVector& theB1, const reco::Candidate::LorentzVector& theB2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TLorentzVector pB1(theB1.px(), theB1.py(), theB1.pz(), theB1.energy());
  TLorentzVector pB2(theB2.px(), theB2.py(), theB2.pz(), theB2.energy());
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  pB1.Boost( -pA.BoostVector() );
  pB2.Boost( -pA.BoostVector() );
  // Build unit vectors orthogonal to the decay planes
  TVector3 Zplane=pL1.Vect().Cross( pL2.Vect() ); // L1 x L2
  TVector3 Hplane=pB1.Vect().Cross( pB2.Vect() ); // B1 x B2
  Zplane.SetMag(1.);
  Hplane.SetMag(1.);
  // Sign of Phi
  TLorentzVector pZr = pL1 + pL2;
  float sgn = pZr.Vect().Dot( Zplane.Cross(Hplane) );
  sgn/=fabs(sgn);
  
  float value=sgn * acos( Zplane.Dot(Hplane) );
  if(value!=value || isinf(value)) return 0.;
  return value;
}


float Utilities::ReturnPhi1(const reco::Candidate::LorentzVector& theA, const reco::Candidate::LorentzVector& theL1, const reco::Candidate::LorentzVector& theL2) {
  TLorentzVector pA(theA.px(), theA.py(), theA.pz(), theA.energy());
  TLorentzVector pL1(theL1.px(), theL1.py(), theL1.pz(), theL1.energy());
  TLorentzVector pL2(theL2.px(), theL2.py(), theL2.pz(), theL2.energy());
  TVector3 beamAxis(0., 0., 1.);
  // Boost objects to the A rest frame
  pL1.Boost( -pA.BoostVector() );
  pL2.Boost( -pA.BoostVector() );
  // Reconstruct Z in A rest frame
  TLorentzVector pZr = pL1 + pL2;
  // Build unit vectors orthogonal to the decay planes
  TVector3 Zplane=pL1.Vect().Cross( pL2.Vect() ); // L1 x L2
  TVector3 Bplane=beamAxis.Cross( pZr.Vect() ); // Beam x Z, beam/Z plane
  Zplane.SetMag(1.);
  Bplane.SetMag(1.);
  // Sign of Phi1
  float sgn = pZr.Vect().Dot( Zplane.Cross(Bplane) );
  sgn/=fabs(sgn);
  
  float value=sgn * acos( Zplane.Dot(Bplane) );
  if(value!=value || isinf(value)) return 0.;
  return value;
}


// Kinematics and reconstruction

/*

// -----------------------------------
// ---------- KINEMATIC FIT ----------
// -----------------------------------

float Utilities::PerformKinematicFit(pat::Jet* tJet1, pat::Jet* tJet2, reco::Candidate::LorentzVector* fJet1, reco::Candidate::LorentzVector* fJet2, float mass) {
  //  TLorentzVector b1, b2;
  //  b1.SetPtEtaPhiE(tJet1->pt(), tJet1->eta(), tJet1->phi(), tJet1->energy());
  //  b2.SetPtEtaPhiE(tJet2->pt(), tJet2->eta(), tJet2->phi(), tJet2->energy());
    
    TMatrixD m1(3,3);
    TMatrixD m2(3,3);
    m1.Zero();
    m2.Zero();

    //In this example the covariant matrix depends on the transverse energy and eta of the jets
    m1(0,0) = GetErrEt (tJet1->et(), tJet1->eta()); // et
    m1(1,1) = GetErrEta(tJet1->et(), tJet1->eta()); // eta
    m1(2,2) = GetErrPhi(tJet1->et(), tJet1->eta()); // phi
    m2(0,0) = GetErrEt (tJet2->et(), tJet2->eta()); // et
    m2(1,1) = GetErrEta(tJet2->et(), tJet2->eta()); // eta
    m2(2,2) = GetErrPhi(tJet2->et(), tJet2->eta()); // phi

  //  TFitParticleEtEtaPhi jet1("Jet1", "Jet1", &b1, &m1);
  //  TFitParticleEtEtaPhi jet2("Jet2", "Jet2", &b2, &m2);
    
  //  TVector3 b1_3=b1.Vect();
  //  TVector3 b2_3=b2.Vect();
    TVector3 b1(tJet1->px(), tJet1->py(), tJet1->pz());
    TVector3 b2(tJet2->px(), tJet2->py(), tJet2->pz());
    TFitParticlePtEtaPhi jet1("Jet1", "Jet1", &b1, tJet1->mass(), &m1 );
    TFitParticlePtEtaPhi jet2("Jet2", "Jet2", &b2, tJet2->mass(), &m2 );

  //  TFitParticleEScaledMomDev jet1("Jet1", "Jet1", &b1, &m1);
  //  TFitParticleEScaledMomDev jet2("Jet2", "Jet2", &b2, &m2);

    //vec1 and vec2 must make a W boson
    TFitConstraintM mCons1("hMassConstraint", "hMass-Constraint", 0, 0, mass);
    mCons1.addParticles1( &jet1, &jet2 );

    //Definition of the fitter
    //Add two constraints
    TKinFitter fitter("fitter", "fitter");
    fitter.addMeasParticle( &jet1 );
    fitter.addMeasParticle( &jet2 );

    fitter.addConstraint( &mCons1 );
    
    //Set convergence criteria
    fitter.setMaxNbIter( 30 );
    fitter.setMaxDeltaS( 1e-2 );
    fitter.setMaxF( 1e-1 );
    fitter.setVerbosity(1);

    // Perform the fit
    fitter.fit();
  //  fitter.print();
    
    
    float dPt1  = jet1.getCurr4Vec()->Pt()  - jet1.getIni4Vec()->Pt();
    float dEta1 = jet1.getCurr4Vec()->Eta() - jet1.getIni4Vec()->Eta();
    float dPhi1 = jet1.getCurr4Vec()->Phi() - jet1.getIni4Vec()->Phi();
    float dPt2  = jet2.getCurr4Vec()->Pt()  - jet2.getIni4Vec()->Pt();
    float dEta2 = jet2.getCurr4Vec()->Eta() - jet2.getIni4Vec()->Eta();
    float dPhi2 = jet2.getCurr4Vec()->Phi() - jet2.getIni4Vec()->Phi();

    float chi2( dPt1*dPt1/m1(0,0) + dEta1*dEta1/m1(1,1) + dPhi1*dPhi1/m1(2,2) + dPt2*dPt2/m2(0,0) + dEta2*dEta2/m2(1,1) + dPhi2*dPhi2/m2(2,2) ); //=fitter.getS();
    
    //float pchi2=TMath::Prob(chi2, fitter.getNDF());
    
    //Hist["k_chi2"]->Fill(chi2, EventWeight);
    //Hist["k_chi2Prob"]->Fill(pchi2, EventWeight);
    //Hist["k_deltaPt1" ]->Fill(dPt1, EventWeight);
    //Hist["k_deltaEta1"]->Fill(dEta1, EventWeight);
    //Hist["k_deltaPhi1"]->Fill(dPhi1, EventWeight);
    //Hist["k_deltaPt2" ]->Fill(dPt2, EventWeight);
    //Hist["k_deltaEta2"]->Fill(dEta2, EventWeight);
    //Hist["k_deltaPhi2"]->Fill(dPhi2, EventWeight);
    //if(tJet1->genParton()) {
      //Hist["k_pullPt1" ]->Fill((jet1.getCurr4Vec()->Pt()-tJet1->genParton()->pt())/sqrt(m1(0,0)), EventWeight);
      //Hist["k_pullEta1"]->Fill((jet1.getCurr4Vec()->Eta()-tJet1->genParton()->eta())/sqrt(m1(1,1)), EventWeight);
      //Hist["k_pullPhi1"]->Fill((jet1.getCurr4Vec()->Phi()-tJet1->genParton()->phi())/sqrt(m1(2,2)), EventWeight);
    //}
    //if(tJet2->genParton()) {
      //Hist["k_pullPt2" ]->Fill((jet2.getCurr4Vec()->Pt()-tJet2->genParton()->pt())/sqrt(m2(0,0)), EventWeight);
      //Hist["k_pullEta2"]->Fill((jet2.getCurr4Vec()->Eta()-tJet2->genParton()->eta())/sqrt(m2(1,1)), EventWeight);
      //Hist["k_pullPhi2"]->Fill((jet2.getCurr4Vec()->Phi()-tJet2->genParton()->phi())/sqrt(m2(2,2)), EventWeight);
    //}
    
    
    // Update objects
    fJet1->SetPxPyPzE(jet1.getCurr4Vec()->Px(), jet1.getCurr4Vec()->Py(), jet1.getCurr4Vec()->Pz(), jet1.getCurr4Vec()->Energy());
    fJet2->SetPxPyPzE(jet2.getCurr4Vec()->Px(), jet2.getCurr4Vec()->Py(), jet2.getCurr4Vec()->Pz(), jet2.getCurr4Vec()->Energy());
    
    return chi2;
}

*/

// Neutrino pz recovering by trading W mass

float Utilities::RecoverNeutrinoPz(const reco::Particle::LorentzVector* lep, const reco::Particle::LorentzVector* met) {
    // W kinematical reconstruction
    float pz = 0.;
    float a = pow(80.4,2) - pow(lep->mass(),2) + 2.*lep->px()*met->px() + 2.*lep->py()*met->py();
    float A = 4*( pow(lep->energy(),2) - pow(lep->pz(),2) );
    float B = -4*a*lep->pz();
    float C = 4*pow(lep->energy(),2) * (pow(met->px(),2)  + pow(met->py(),2)) - pow(a,2);
    float D = pow(B,2) - 4*A*C;
    // If there are real solutions, use the one with lowest pz                                            
    if (D>=0) {
        float s1 = (-B+sqrt(D))/(2*A);
        float s2 = (-B-sqrt(D))/(2*A);
        if(fabs(s1)<fabs(s2)) pz=s1;
        else pz=s2;
    }
    // Otherwise, use real part                                                                           
    else {
        pz = -B/(2*A);
    }
    return pz;
}






pat::CompositeCandidate Utilities::RecoilMassFormula(pat::CompositeCandidate& H, pat::MET& met){
    pat::CompositeCandidate X;
    AddFourMomenta addP4;
    X.addDaughter(H);
    X.addDaughter(met);
    addP4.set(X);
    reco::Particle::LorentzVector metp4 = met.p4();
    reco::Particle::LorentzVector Xp4;
    metp4.SetPz(-H.pz());
    Xp4 += metp4;
    Xp4.SetPz(0);
    float B = -2.*H.energy();
    float C = pow(H.mass(),2) - pow(90.18,2);
    float D = pow(B,2) - 4*1*C;
    float mX;
    if(D>0){
        float s1 = (-B+sqrt(D))/2.;
        float s2 = (-B-sqrt(D))/2.;
        if(fabs(s1)>fabs(s2)) mX = s1;
        else mX = s2;
    }
    else{
        mX = -B/2.;
    }
    Xp4.SetE(sqrt(pow(mX,2) + pow(Xp4.Px(),2) + pow(Xp4.Py(),2) + pow(Xp4.Pz(),2)));
    X.setP4(Xp4);
    X.setCharge(0);
    return X;
}

//Find mother id of a const reco::GenParticle, method used for genPartons MC truth
int Utilities::FindMotherId(const reco::GenParticle* p) {
  int pId = p->pdgId();
  const reco::Candidate* mom = p->mother();
  while (mom != 0 && mom->pdgId() == pId)
    mom = mom->mother();
  return mom->pdgId();
}

