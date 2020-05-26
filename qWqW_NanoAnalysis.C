#define qWqW_NanoAnalysis_cxx
#include "qWqW_NanoAnalysis.h"
#include <fstream>
#include <iostream>     // std::cout
#include <algorithm>    // std::sort
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>  // for high_resolution_clock

//Includes for HITFIT
#include "../TopQuarkAnalysis/TopHitFit/interface/fourvec.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/RunHitFit.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/Top_Decaykin.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/ElectronTranslatorBase.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/MuonTranslatorBase.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/JetTranslatorBase.h"
#include "../TopQuarkAnalysis/TopHitFit/interface/METTranslatorBase.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace std;
using namespace hitfit;

//Functions to help sorting
bool sortinrev(const pair<double, TLorentzVector> &a, const pair<double,TLorentzVector> &b)
{
  return (a.second.Pt() > b.second.Pt());
}

bool sortinrev2(const pair<double, int> &a, const pair<double,int> &b)
{
  return (a.first > b.first);
}

// ----------------------------------------------------------------------------
// MAIN EVENT LOOP
// ----------------------------------------------------------------------------
void qWqW_NanoAnalysis::Loop(bool isSignal, bool isTTbar, Long64_t skipevents, Long64_t maxevents)
{
  
  // ----------------------------------------------------------------------------
  // Variables for HITFIT
  // ____________________________________________________________________________
  
  // The following three initializers instantiate the translator to the HitFit objects
  ElectronTranslatorBase electronTranslator;
  MuonTranslatorBase muonTranslator;
  JetTranslatorBase jetTranslator;
  METTranslatorBase metTranslator;
  
  
  //TIMER->  Record start time
  auto start = std::chrono::high_resolution_clock::now();
  auto newStart = start;

  // ---------------------------------------------------------------------------
  // PUPPI W Corrections File
  // ---------------------------------------------------------------------------
  TFile* file = TFile::Open("puppiCorr.root","READ");
  TF1 *puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  TF1 *puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  TF1 *puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");
  file->Close();  

  // ---------------------------------------------------------------------------
  // OUTPUT FILE
  // ---------------------------------------------------------------------------
  TFile *outputFile;
  outputFile = new TFile(outputFileName,"RECREATE");
  outputFile->cd();

  TH1D *nEventsWeighted = new TH1D("nEventsWeighted","nEventsWeighted",3,-1,2);

  TTree *outputTree = new TTree("Events","Events");

  // ---------------------------------------------------------------------------
  // Define Branches
  // ---------------------------------------------------------------------------
  // Basics
  int isElectron;
  int isMuon;
  outputTree->Branch("event",&event,"event/L");
  outputTree->Branch("lumi",&luminosityBlock,"lumi/I");
  outputTree->Branch("run",&run,"run/I");
  outputTree->Branch("isElectron",&isElectron,"isElectron/I");
  outputTree->Branch("isMuon",&isMuon,"isMuon/I");

  // weights
  float pileupWeight, pileupWeightUp, pileupWeightDn, topPtWeight;
  float HTSF_Pol, HTSF_PolUp, HTSF_PolDn;
  float genTTbarMass;
  outputTree->Branch("genWeight",&genWeight,"genWeight/F");
  outputTree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
  outputTree->Branch("pileupWeightUp",&pileupWeightUp,"pileupWeightUp/F");
  outputTree->Branch("pileupWeightDn",&pileupWeightDn,"pileupWeightDn/F");
  outputTree->Branch("LHEScaleWeight",LHEScaleWeight,"LHEScaleWeight[9]/F");
  outputTree->Branch("LHEPdfWeight",LHEPdfWeight,"LHEPdfWeight[103]/F");
  outputTree->Branch("btagWeight_DeepCSVB",&btagWeight_DeepCSVB,"btagWeight_DeepCSVB/F");
  outputTree->Branch("HTSF_Pol",&HTSF_Pol,"HTSF_Pol/F");
  outputTree->Branch("HTSF_PolUp",&HTSF_PolUp,"HTSF_PolUp/F");
  outputTree->Branch("HTSF_PolDn",&HTSF_PolDn,"HTSF_PolDn/F");
  outputTree->Branch("topPtWeight",&topPtWeight,"topPtWeight/F");
  outputTree->Branch("genTTbarMass",&genTTbarMass,"genTTbarMass/F");

  //Lepton Information 
  float leppt, lepeta, lepphi, lepM;
  int lepId;
  outputTree->Branch("Lepton_pt",&leppt,"Lepton_pt/F");
  outputTree->Branch("Lepton_eta",&lepeta,"Lepton_eta/F");
  outputTree->Branch("Lepton_phi",&lepphi,"Lepton_phi/F");
  outputTree->Branch("Lepton_mass", &lepM,"Lepton_mass/F");
  outputTree->Branch("Lepton_Id",&lepId,"Lepton_Id/I");

  //AK4 Jets
  int NJets;
  int NJetsDeepCSVmed;
  int NJetsDeepCSVloose;
  float AK4HT;
  vector<float> Jet_pt_ptordered;
  vector<float> Jet_eta_ptordered;
  vector<float> Jet_phi_ptordered;
  vector<float> Jet_mass_ptordered;
  vector<float> Jet_bdisc_ptordered;  
  outputTree->Branch("NJets",&NJets, "NJets/I");
  outputTree->Branch("NJetsDeepCSVmed",&NJetsDeepCSVmed, "NJetsDeepCSVmed/I");
  outputTree->Branch("NJetsDeepCSVloose",&NJetsDeepCSVloose, "NJetsDeepCSVloose/I");
  outputTree->Branch("AK4HT",&AK4HT, "AK4HT/F");
  outputTree->Branch("Jet_pt",&Jet_pt_ptordered);
  outputTree->Branch("Jet_eta",&Jet_eta_ptordered);
  outputTree->Branch("Jet_phi",&Jet_phi_ptordered);
  outputTree->Branch("Jet_mass",&Jet_mass_ptordered);
  outputTree->Branch("Jet_bdisc",&Jet_bdisc_ptordered);

  //AK8 Jets
  int NJetsAK8;
  bool usesWtag;
  vector<float> FatJet_pt_ptordered;
  vector<float> FatJet_eta_ptordered;
  vector<float> FatJet_phi_ptordered;
  vector<float> FatJet_tau21_ptordered;
  vector<float> FatJet_mass_ptordered;
  vector<float> FatJet_sdmass_ptordered;
  vector<float> FatJet_sdmassraw_ptordered;
  vector<float> FatJet_sdmasscorr_ptordered;
  vector<int> FatJet_index1_ptordered;
  vector<int> FatJet_index2_ptordered;
  outputTree->Branch("usesWtag",&usesWtag,"usesWtag/O");
  outputTree->Branch("NFatJets",&NJetsAK8,"NFatJets/I");
  outputTree->Branch("FatJet_pt",&FatJet_pt_ptordered);
  outputTree->Branch("FatJet_eta",&FatJet_eta_ptordered);
  outputTree->Branch("FatJet_phi",&FatJet_phi_ptordered);
  outputTree->Branch("FatJet_mass",&FatJet_mass_ptordered);
  outputTree->Branch("FatJet_sdmass",&FatJet_sdmass_ptordered);
  outputTree->Branch("FatJet_sdmassraw",&FatJet_sdmassraw_ptordered);
  outputTree->Branch("FatJet_sdmasscorr",&FatJet_sdmasscorr_ptordered);
  outputTree->Branch("FatJet_index1",&FatJet_index1_ptordered);
  outputTree->Branch("FatJet_index2",&FatJet_index2_ptordered);
  outputTree->Branch("FatJet_tau21",&FatJet_tau21_ptordered);

  //Subjets 
  vector<float> SubJet_pt_fjordered;
  vector<float> SubJet_eta_fjordered;
  vector<float> SubJet_phi_fjordered;
  vector<float> SubJet_mass_fjordered;
  vector<float> SubJet_bdisc_fjordered;
  outputTree->Branch("SubJet_pt",&SubJet_pt_fjordered);
  outputTree->Branch("SubJet_eta",&SubJet_eta_fjordered);
  outputTree->Branch("SubJet_phi",&SubJet_phi_fjordered);
  outputTree->Branch("SubJet_mass",&SubJet_mass_fjordered);
  outputTree->Branch("Subjet_bdisc",&SubJet_bdisc_fjordered);
  
  // Hybrid List of Jets
  vector<float> HybridJet_pt;
  vector<float> HybridJet_eta;
  vector<float> HybridJet_phi;
  vector<float> HybridJet_mass;
  outputTree->Branch("HybridJet_pt",&HybridJet_pt);
  outputTree->Branch("HybridJet_eta",&HybridJet_eta);
  outputTree->Branch("HybridJet_phi",&HybridJet_phi);
  outputTree->Branch("HybridJet_mass",&HybridJet_mass);

  //MET -- defined in input tree
  outputTree->Branch("MET_pt",&MET_pt, "MET_pt/F");
  outputTree->Branch("MET_phi",&MET_phi, "MET_phi/F");
  
  //Masses
  float totalMass;
  float threeJetMass;
  outputTree->Branch("totalMass", &totalMass, "totalMass/F");
  outputTree->Branch("threeJetMass", &threeJetMass, "threeJetMass/F");  
  
  //HITFIT
  float tMass;
  float chi2;
  float probchi2;
  float sigmaM;
  float nJetsHitFit;
  float SLDivideST;
  outputTree->Branch("tMass",&tMass, "tMass/F");
  outputTree->Branch("chi2", &chi2, "chi2/F");
  outputTree->Branch("probchi2", &probchi2, "probchi2/F");
  outputTree->Branch("sigmaM", &sigmaM, "sigmaM/F");
  outputTree->Branch("SLDivideST",&SLDivideST,"SLDivideST/F");
  outputTree->Branch("nJetsHitFit",&nJetsHitFit, "nJetsHitFit/F");

  //Fitted Content
  float fittedLepPt;
  float fittedLepEta;
  float fittedLepPhi;
  float fittedLepMass;
  float fittedMETPt;
  float fittedMETEta;
  float fittedMETPhi;
  float fittedMETMass;
  vector<float> fittedJetPt;
  vector<float> fittedJetEta;
  vector<float> fittedJetPhi;
  vector<float> fittedJetMass;
  outputTree->Branch("fittedLepPt",&fittedLepPt, "fittedLepPt/F");
  outputTree->Branch("fittedLepEta",&fittedLepEta, "fittedLepEta/F");
  outputTree->Branch("fittedLepPhi",&fittedLepPhi, "fittedLepPhi/F");
  outputTree->Branch("fittedLepMass",&fittedLepMass, "fittedLepMass/F");
  outputTree->Branch("fittedMETPt",&fittedMETPt, "fittedMETPt/F");
  outputTree->Branch("fittedMETEta",&fittedMETEta, "fittedMETEta/F");
  outputTree->Branch("fittedMETPhi",&fittedMETPhi, "fittedMETPhi/F");
  outputTree->Branch("fittedMETMass",&fittedMETMass, "fittedMETMass/F");
  outputTree->Branch("fittedJetPt",&fittedJetPt);
  outputTree->Branch("fittedJetEta",&fittedJetEta);
  outputTree->Branch("fittedJetPhi",&fittedJetPhi);
  outputTree->Branch("fittedJetMass",&fittedJetMass);

  // ----------------------------------------------------------------------------
  // Define and initialize objects / cuts / efficiencies
  // ----------------------------------------------------------------------------
  
  float lepPtCut=55;
  float elEtaCut=2.5;
  float muEtaCut=2.4;
  float jetEtaCut=2.4;
  float ak8EtaCut=2.4;
  float ak8SJEtaCut=2.5;
  float jetPtCut=30;
  float ak8PtCut=200;
  float metCut=50;
  int   nAK8jetsCut=0;
  int overlap1=0;
  int overlap2=0;
  int highestOverlapping=0;
  
  
  //counters
  int passProbChisq=0;
  int numEventsElectron= 0;
  int numEventsMuon = 0;
  int numEvents=0;
  int numEvents2=0;
  int noMatch=0; 
  int noGoodMatch=0;
  int n_jetstotal        = 0;
  int n_jetsnearlep      = 0;
  int n_jetspassed       = 0;
  int npass_ThreeJets    = 0;
  int npass_trigger      = 0;
  int npass_mu500        = 0;
  int npass_met          = 0;
  int npass_ht           = 0;
  int npass_nAK8jets     = 0;
  int npass_nHjets       = 0;
  int npass_lepPt        = 0;
  int npass_ElEta        = 0;
  int npass_MuEta        = 0;
  int npass_all          = 0;
  int Nelectrons         = 0;
  int Nmuons             = 0;
  int passbwbw 	  = 0;
  int passflag  	  = 0;
  int passHLT 		  = 0;
  int passLepton 	  = 0;
  int passJetCut	  = 0;
  int passMET            = 0;
  int passIso   	  = 0;
  int passEl 		  = 0;
  int passMu 		  = 0;
  int passLooseEl	  = 0;
  int passTightEl	  = 0;
  int passElAll 	  = 0;
  int passHybridJetCut   = 0;
  int passHitFitJetCut   = 0;
  int passMuAll	  = 0;
  int passElAll2	  = 0;
  int passMuAll2 	  = 0;
  int passElAll3         = 0;
  int passMuAll3         = 0;
  int numPassed	  = 0;
  int passAnyChi2 = 0;
  
  
  //Counts for HITFIT errors
  int num999=0;
  int num998=0;
  int num997=0;
  int num996=0;
  int num995=0;
  int ntotal=0;
  int nfailed=0;
  int ngood=0;
  int good2M=0;
  int goodML=0;
  int good1M=0;
  int numLowLepB=0;
  int numLowHadB=0;
  int numLowHadW=0;
  int numBadW=0;  
  int numBadW2=0;
  int numLowNB=0;
  int numLowST=0;
  int numBadSLST=0;

  // pileup
  vector<float> pileupweight = { 0.000e+00,  1.250e+01,  5.025e+01,  1.900e+01,  1.197e+01,  9.118e+00,  6.647e+00,  4.832e+00,  3.633e+00,  2.747e+00,  
		   2.238e+00,  1.904e+00,  1.694e+00,  1.572e+00,  1.502e+00,  1.468e+00,  1.451e+00,  1.452e+00,  1.461e+00,  1.467e+00,  
		   1.458e+00,  1.436e+00,  1.401e+00,  1.352e+00,  1.302e+00,  1.250e+00,  1.205e+00,  1.167e+00,  1.138e+00,  1.116e+00,  
		   1.100e+00,  1.089e+00,  1.081e+00,  1.077e+00,  1.074e+00,  1.069e+00,  1.063e+00,  1.052e+00,  1.035e+00,  1.012e+00,  
		   9.794e-01,  9.367e-01,  8.867e-01,  8.272e-01,  7.609e-01,  6.925e-01,  6.209e-01,  5.505e-01,  4.824e-01,  4.193e-01,  
		   3.621e-01,  3.103e-01,  2.658e-01,  2.276e-01,  1.948e-01,  1.679e-01,  1.455e-01,  1.262e-01,  1.098e-01,  9.584e-02,  
		   8.398e-02,  7.342e-02,  6.469e-02,  5.654e-02,  4.900e-02,  4.233e-02,  3.634e-02,  3.081e-02,  2.661e-02,  2.266e-02,  
		   1.889e-02,  1.652e-02,  1.344e-02,  1.173e-02,  9.423e-03,  8.155e-03,  6.671e-03,  5.637e-03,  4.504e-03,  3.920e-03,  
		   3.007e-03,  2.567e-03,  1.909e-03,  1.425e-03,  1.084e-03,  8.974e-04,  8.267e-04,  4.805e-04,  2.814e-04,  1.588e-04,  
		   1.206e-04,  6.816e-05,  2.870e-05,  1.436e-05,  7.136e-06,  4.702e-06,  1.107e-06,  1.723e-06,  1.079e-06,  1.269e-07,  0.000e+00,  }; // TTToSemiLeptonic
  vector<float> pileupweightDn = { 0.000e+00,  1.447e+01,  5.931e+01,  2.226e+01,  1.397e+01,  1.063e+01,  7.865e+00,  5.773e+00,  4.367e+00,  3.362e+00,  
		     2.793e+00,  2.408e+00,  2.153e+00,  1.999e+00,  1.907e+00,  1.857e+00,  1.824e+00,  1.805e+00,  1.789e+00,  1.761e+00,  
		     1.710e+00,  1.645e+00,  1.572e+00,  1.493e+00,  1.423e+00,  1.358e+00,  1.306e+00,  1.264e+00,  1.232e+00,  1.206e+00,  
		     1.185e+00,  1.166e+00,  1.148e+00,  1.131e+00,  1.111e+00,  1.085e+00,  1.054e+00,  1.014e+00,  9.643e-01,  9.065e-01,  
		     8.402e-01,  7.663e-01,  6.892e-01,  6.091e-01,  5.296e-01,  4.550e-01,  3.850e-01,  3.223e-01,  2.671e-01,  2.202e-01,  
		     1.810e-01,  1.482e-01,  1.217e-01,  1.003e-01,  8.289e-02,  6.902e-02,  5.780e-02,  4.834e-02,  4.046e-02,  3.384e-02,  
		     2.830e-02,  2.351e-02,  1.962e-02,  1.619e-02,  1.322e-02,  1.075e-02,  8.680e-03,  6.919e-03,  5.624e-03,  4.508e-03,  
		     3.538e-03,  2.915e-03,  2.233e-03,  1.832e-03,  1.382e-03,  1.121e-03,  8.562e-04,  6.734e-04,  4.991e-04,  4.012e-04,  
		     2.833e-04,  2.218e-04,  1.506e-04,  1.023e-04,  7.060e-05,  5.282e-05,  4.385e-05,  2.290e-05,  1.202e-05,  6.059e-06,  
		     4.098e-06,  2.061e-06,  7.699e-07,  3.418e-07,  1.507e-07,  8.865e-08,  1.895e-08,  2.792e-08,  1.816e-08,  2.615e-09,  0.000e+00,  }; // TTToSemiLeptonic
  vector<float> pileupweightUp = { 0.000e+00,  1.091e+01,  4.288e+01,  1.632e+01,  1.032e+01,  7.892e+00,  5.683e+00,  4.086e+00,  3.057e+00,  2.285e+00,  
		     1.826e+00,  1.531e+00,  1.349e+00,  1.249e+00,  1.194e+00,  1.169e+00,  1.161e+00,  1.169e+00,  1.188e+00,  1.211e+00,  
		     1.226e+00,  1.233e+00,  1.230e+00,  1.210e+00,  1.183e+00,  1.147e+00,  1.112e+00,  1.080e+00,  1.053e+00,  1.034e+00,  
		     1.020e+00,  1.013e+00,  1.010e+00,  1.013e+00,  1.018e+00,  1.026e+00,  1.035e+00,  1.044e+00,  1.051e+00,  1.055e+00,  
		     1.054e+00,  1.045e+00,  1.030e+00,  1.004e+00,  9.672e-01,  9.247e-01,  8.727e-01,  8.155e-01,  7.534e-01,  6.903e-01,  
		     6.275e-01,  5.648e-01,  5.067e-01,  4.529e-01,  4.034e-01,  3.604e-01,  3.231e-01,  2.891e-01,  2.595e-01,  2.338e-01,  
		     2.119e-01,  1.921e-01,  1.761e-01,  1.607e-01,  1.459e-01,  1.325e-01,  1.199e-01,  1.073e-01,  9.796e-02,  8.821e-02,  
		     7.775e-02,  7.190e-02,  6.182e-02,  5.697e-02,  4.835e-02,  4.420e-02,  3.822e-02,  3.417e-02,  2.895e-02,  2.677e-02,  
		     2.188e-02,  1.996e-02,  1.591e-02,  1.278e-02,  1.049e-02,  9.408e-03,  9.417e-03,  5.965e-03,  3.820e-03,  2.364e-03,  
		     1.972e-03,  1.229e-03,  5.718e-04,  3.170e-04,  1.748e-04,  1.282e-04,  3.363e-05,  5.837e-05,  4.079e-05,  5.342e-06,  0.000e+00,  }; // TTToSemiLeptonic


  // WJets and DYJets reweighting for split HT samples
   TF1 *poly2 = new TF1("poly2","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2->SetParameter(0,    0.998174);  
   poly2->SetParameter(1, 8.40861e-05); 
   poly2->SetParameter(2,-6.63274e-07);
   poly2->SetParameter(3, 4.09272e-10); 
   poly2->SetParameter(4,-9.50233e-14); 
   poly2->SetParameter(5, 7.59648e-18); 
   poly2->SetParameter(6,0.402806);
   
   TF1 *poly2U = new TF1("poly2U","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2U->SetParameter(0,    0.998164);  
   poly2U->SetParameter(1, 8.51652e-05); 
   poly2U->SetParameter(2,-6.62919e-07);
   poly2U->SetParameter(3,  4.1205e-10); 
   poly2U->SetParameter(4,-9.57795e-14); 
   poly2U->SetParameter(5, 7.67371e-18); 
   poly2U->SetParameter(6,0.454456);
   
   TF1 *poly2D = new TF1("poly2D","max([6],[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x)",100,5000);
   poly2D->SetParameter(0,    0.998183);  
   poly2D->SetParameter(1, 8.30069e-05); 
   poly2D->SetParameter(2,-6.63629e-07);
   poly2D->SetParameter(3, 4.06495e-10); 
   poly2D->SetParameter(4,-9.42671e-14); 
   poly2D->SetParameter(5, 7.51924e-18); 
   poly2D->SetParameter(6,0.351156);


  //Lorentz vectors
  TLorentzVector jet_lv;
  TLorentzVector lepton_lv;
  TLorentzVector ak8_lv;
  vector<pair<int, TLorentzVector>> listOfJets;
  vector<float> listOfJetBDisc;
  vector<pair<float, int>> ptOrderedJets;
  vector<int> FatJet_sj2_ptordered;
  vector<int> FatJet_sj1_ptordered;
  vector<int> tPrimeID;
  vector<int> bPrimeID;
  vector<int> listofQuarkIDs;
  vector<int> listofBosonIDs;
  vector<unsigned int> quarks;
  vector<unsigned int> bosons;

  RunHitFit HitFit("../TopQuarkAnalysis/TopHitFit/data/setting/RunHitFitConfiguration.txt",80.4,80.4,0);//,hitfitLepWMass_,hitfitHadWMass_,hitfitTopMass_);
  
  // ----------------------------------------------------------------------------
  // RUN THE EVENT LOOP
  // ----------------------------------------------------------------------------
  
  cout << "RUNNING LOOP" << endl;
  Long64_t nentries = fChain->GetEntriesFast();
  
  if (fChain == 0) return;
  Long64_t nbytes = 0, nb = 0;
  //Run only specific events
  Long64_t runthrough = skipevents+maxevents;
  if (runthrough > nentries) runthrough = nentries;
  for (Long64_t jentry= skipevents; jentry< runthrough; jentry++) {

    //if (jentry > skipevents+10000) continue;

    if(jentry%1000 == 0){
      cout << "On entry " << jentry << " of " << runthrough << endl; 
    }
     
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;    

    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (Cut(ientry) != 1) continue;

    // --------------------------------------------------------------------------
    // Fill the counter histogram for total events processed
    // ----------------------------------------------------------------------
    nEventsWeighted->Fill(genWeight/fabs(genWeight)); // + or - 1
    
    // --------------------------------------------------------------------------
    // Find generated tops and save masses, get reweighting
    // ----------------------------------------------------------------------

    genTTbarMass = -999;
    topPtWeight = 1.0;    
    TLorentzVector top, antitop;
    bool gottop = false;
    bool gotantitop = false;	
    bool gottoppt = false;
    bool gotantitoppt = false;
    float toppt, antitoppt;
    // std::bitset<16> statbits; // saving this for info

    if(isTTbar){
      for(unsigned int p = 0; p < nGenPart; p++) {
	int id = GenPart_pdgId[p];
	if (abs(id) != 6) continue;
	if (GenPart_mass[p] < 10) continue; 

	// statbits = std::bitset<16>(GenPart_statusFlags[p]);	// saving this for info
	int motherid = GenPart_pdgId[GenPart_genPartIdxMother[p]];
	if(abs(motherid) != 6){
	  if (!gottop && id == 6){ 
	    top.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
	    gottop = true;
	  }
	  if (!gotantitop && id == -6){ 
	    antitop.SetPtEtaPhiM(GenPart_pt[p], GenPart_eta[p], GenPart_phi[p], GenPart_mass[p]);
	    gotantitop = true;
	  }
	}

	if(GenPart_status[p] == 62){
	  if (!gottoppt && id == 6){ 
	    toppt = GenPart_pt[p];
	    gottoppt = true;
	  }
	  if (!gotantitoppt && id == -6){ 
	    antitoppt = GenPart_pt[p];
	    gotantitoppt = true;
	  }
	}
      }
      if(gottop && gotantitop) genTTbarMass = (top+antitop).M();

      if(gottoppt && gotantitoppt){
	float SFtop = TMath::Exp(0.0615-0.0005*toppt);
	float SFantitop = TMath::Exp(0.0615-0.0005*antitoppt);
	topPtWeight = TMath::Sqrt(SFtop*SFantitop);
      }
    }	

    // FILTER events based on output file name and calculated tt mass:
    if(outputFileName.Contains("Mtt0") && genTTbarMass > 700) continue;
    if(outputFileName.Contains("Mtt700") && (genTTbarMass <= 700 || genTTbarMass > 1000)) continue;
    if(outputFileName.Contains("Mtt1000") && genTTbarMass <= 1000) continue; 

    // --------------------------------------------------------------------------
    // Check to make sure event isBWBW
    // ----------------------------------------------------------------------
    bool isBWBW=false;
    tPrimeID.clear();
    bPrimeID.clear();
    listofQuarkIDs.clear();
    listofBosonIDs.clear();
    quarks.clear();
    bosons.clear();
    
    if(isSignal){
      for(unsigned int  p = 0; p < nGenPart; p++) {
	int id=GenPart_pdgId[p];
	// find T' and B' particles 
	if(abs(id) != 8000001 && abs(id) != 8000002) {continue;}
	bool hasTdaughter = false;
	vector<unsigned int> daughters;
	daughters.clear();
	for(unsigned int  dau = 0; dau < nGenPart; dau++){
	  if(GenPart_genPartIdxMother[dau]!=p) continue;
	  daughters.push_back(dau);
	  if(abs(id) == 8000001 && abs(GenPart_pdgId[dau]) == 8000001) hasTdaughter=true;
	  if(abs(id) == 8000002 && abs(GenPart_pdgId[dau]) == 8000002) hasTdaughter=true;
	}
	if(hasTdaughter){continue;}
	int mother=GenPart_genPartIdxMother[p];
	int mother_id=GenPart_pdgId[mother];
	if(abs(id) == 8000001){
	  if(abs(mother_id)== 8000001){
	    tPrimeID.push_back( GenPart_pdgId[mother] );
	  }else{
	    tPrimeID.push_back( GenPart_pdgId[p] );
	  }
	}
	if(abs(id) == 8000002){
	  if(abs(mother_id)== 8000002){
	    bPrimeID.push_back( GenPart_pdgId[mother] );
	  }else{
	    bPrimeID.push_back( GenPart_pdgId[p] );
	  }
	}
	for(unsigned int j=0; j<daughters.size(); j++){
	  unsigned int d=daughters.at(j);
	  int dauId=GenPart_pdgId[d];
	  if(abs(dauId) == 5 || abs(dauId) == 6){
	    quarks.push_back(d);
	    listofQuarkIDs.push_back(dauId);
	  }
	  else if(abs(dauId) > 22 && abs(dauId) < 26){
	    bosons.push_back(d);
	    listofBosonIDs.push_back(dauId);
	  }
	  else continue;
	}
      }
      
      if(tPrimeID.size() > 0 && bPrimeID.size() > 0) cout << "Found both T' and B' " << endl;
      
      if(listofQuarkIDs.size()!=0 && listofQuarkIDs.size()!=2 ){
	std::cout << "More/less than 2 quarks stored: " << listofQuarkIDs.size() << std::endl;
	for(unsigned int i = 0; i < listofQuarkIDs.size(); i++){std::cout << "quark " << i << " = " << listofQuarkIDs.at(i) << std::endl;}
	int test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
	int sign=-1;
	if(test>0) sign=1;
	if(sign>0){
	  if(listofQuarkIDs.size()==4){
	    std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
	    std::swap(quarks.at(2),quarks.at(3));
	  }
	  std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
	  std::swap(quarks.at(1),quarks.at(2));
	  test = listofQuarkIDs.at(0)*listofQuarkIDs.at(1);
	  sign=-1;
	  if(test>0) sign=1;
	  if(sign < 0) std::cout << "Signs are fixed!" << std::endl;
	}
	if(listofQuarkIDs.size() > 3 && abs(listofQuarkIDs.at(3)) == 6){
	  std::swap(listofQuarkIDs.at(2),listofQuarkIDs.at(3));
	  std::swap(quarks.at(2),quarks.at(3));
	}
	if(listofQuarkIDs.size() > 2 && abs(listofQuarkIDs.at(2)) == 6){
	  std::swap(listofQuarkIDs.at(1),listofQuarkIDs.at(2));
	  std::swap(quarks.at(1),quarks.at(2));
	}
      }
      if(listofBosonIDs.size() != 0 && listofBosonIDs.size() != 2) std::cout << "More/less than 2 bosons stored: " << listofBosonIDs.size() << std::endl;
      
      //Check if isBWBW
      if(tPrimeID.size() > 1 && bPrimeID.size() == 0){
	if(abs(listofQuarkIDs.at(0)) == 5 && abs(listofQuarkIDs.at(1)) == 5){
	  if(abs(listofBosonIDs.at(0)) == 24 && abs(listofBosonIDs.at(1)) == 24)
	    isBWBW = true;
	}
      }
      
      if(!isBWBW){continue;}
    }
    passbwbw++;

    // --------------------------------------------------------------------------
    // Other Cleanup selections 
    // --------------------------------------------------------------------------
    int highPurity=0;
    if (nIsoTrack>=10){
      for(int i=0;i<nIsoTrack;i++){
	if(IsoTrack_isHighPurityTrack[i])
	  highPurity++;
      }
      if((double(highPurity)/nIsoTrack)<0.25){continue;}
    }
    passIso++;
    
     
    // ---------------------------------------------------------------------------
    // Flags
    // ---------------------------------------------------------------------------   
    if (!(Flag_HBHENoiseFilter==1&&Flag_HBHENoiseIsoFilter==1&&Flag_EcalDeadCellTriggerPrimitiveFilter==1&&Flag_goodVertices==1&&Flag_eeBadScFilter==1&&Flag_globalTightHalo2016Filter==1&&Flag_BadPFMuonFilter==1)){continue;}
    
    passflag++;
     
      
    // ----------------------------------------------------------------------------
    // Assign as electron or muon event
    // ----------------------------------------------------------------------------     
    isElectron =0;
    isMuon=0;
    lepId=-1;
    //Check Electron Events
    int numTight=0;
    int numLoose=0;
    int indexOfTight=0;
    for(int index=0;index<nElectron;index++){
      if( Electron_pt[index] > lepPtCut && 
	  fabs(Electron_eta[index]) < elEtaCut && 
	  Electron_cutBased[index] == 4 &&
	  ((Electron_eta[index] <= 1.479 && Electron_dxy[index] < 0.05 && Electron_dz[index] < 0.10) ||
	   (Electron_eta[index] > 1.479 && Electron_dxy[index] < 0.1 && Electron_dz[index] < 0.2)) && 
	  ((Electron_eta[index] < -2.5 || Electron_eta[index] > -1.479) || (Electron_phi[index] < -1.55 || Electron_phi[index] >- 0.9))){ // hole cut
	numTight++; //count number of tight electrons
	indexOfTight=index;
      }
      else if( Electron_pt[index] > 20 && 
	       fabs(Electron_eta[index]) < 2.5 && 
	       Electron_cutBased[index] >= 2 && 
	       ((Electron_eta[index] <= 1.479 && Electron_dxy[index] < 0.05 && Electron_dz[index] < 0.10) ||
		(Electron_eta[index] > 1.479 && Electron_dxy[index] < 0.1 && Electron_dz[index]<0.2)) && 
	       ((Electron_eta[index] < -2.5 || Electron_eta[index] > -1.479) || (Electron_phi[index] < -1.55 || Electron_phi[index] >- 0.9))){
	numLoose++; //count number of loose electrons`
      }
    } 
    //Check Muon Events
    int numTightMuon=0;
    int numLooseMuon=0;
    int indexMuon=0;
    for(int index=0;index<nMuon;index++){
      if((Muon_pt[index]>lepPtCut && fabs(Muon_eta[index])<muEtaCut && Muon_miniPFRelIso_all[index]<0.15 && Muon_tightId[index]==1)){
	numTightMuon++;
	indexMuon=index;
      }
      else if ((Muon_pt[index]>20 && fabs(Muon_eta[index])<2.5 && Muon_miniPFRelIso_all[index]<0.25 && Muon_looseId[index]==1))
	numLooseMuon++;
    }
        
    if (numTight==1 && numLoose==0 && numTightMuon==0 && numLooseMuon==0){
      isElectron=1;
      passEl++;
    } else if (numTight==0 && numLoose==0 && numTightMuon==1 && numLooseMuon==0){
      isMuon=1;
      passMu++;
    } else {continue;}
    passLepton++;
    //Set variables for pt and e of lepton
    leppt = 0;
    lepeta = 0;
    if(isElectron){leppt = Electron_pt[indexOfTight]; lepeta = Electron_eta[indexOfTight];}
    if(isMuon){leppt = Muon_pt[indexMuon]; lepeta = Muon_eta[indexMuon];}
    
    // ----------------------------------------------------------------------------
    // Triggers
    // ----------------------------------------------------------------------------
    if (isElectron && HLT_Ele38_WPTight_Gsf!=1) continue;
    if (isMuon && HLT_Mu50!=1){continue;}
    passHLT++;

    // ----------------------------------------------------------------------------
    // Pileup Weight
    // ----------------------------------------------------------------------------
    
    if(Pileup_nTrueInt > 99) Pileup_nTrueInt = 99;
    pileupWeight = pileupweight[Pileup_nTrueInt];
    pileupWeightUp = pileupweightUp[Pileup_nTrueInt];
    pileupWeightDn = pileupweightDn[Pileup_nTrueInt];

    // ----------------------------------------------------------------------------
    // Generator-level HT correction -- only apply to WJets / DYJets / QCD in plotting!
    // ----------------------------------------------------------------------------      
    
    HTSF_Pol = 1;
    HTSF_PolUp = 1;
    HTSF_PolDn = 1;

    // Piece-wise splice with a flat line. Uncertainty from upper/lower error bar fits
    HTSF_Pol = poly2->Eval(LHE_HT);
    HTSF_PolUp = poly2U->Eval(LHE_HT);
    HTSF_PolDn = poly2D->Eval(LHE_HT);
    
    // --------------------------------------------------------------------------
    // Lepton 4-vectors, calculate MT 
    // --------------------------------------------------------------------------
    
    //Set lepton 4-vectors
    lepM=0;
    lepphi=0;
    if (isMuon){
      lepM = 0.105658367;
      lepphi = Muon_phi[indexMuon];
      lepton_lv.SetPtEtaPhiM(Muon_pt[indexMuon],Muon_eta[indexMuon],Muon_phi[indexMuon],lepM);
      lepId = Muon_tightId[indexMuon];
    }
    else{
      lepM = 0.00051099891;
      lepphi = Electron_phi[indexOfTight];
      lepton_lv.SetPtEtaPhiM(Electron_pt[indexOfTight],Electron_eta[indexOfTight],Electron_phi[indexOfTight],lepM);
      lepId = Electron_cutBased[indexOfTight];
    }
    float MT_lepMet = sqrt(2*leppt*MET_pt*(1 - cos(lepphi - MET_phi)));
    
    
    //Cut out events without enough missing energy 
    if(MET_pt < metCut){continue;}
    passMET++;
    
    // ----------------------------------------------------------------------------
    // Loop over AK4 jets for calculations
    // ----------------------------------------------------------------------------
    //cout << "---------" << endl;
    NJets = 0;
    NJetsDeepCSVmed = 0;
    NJetsDeepCSVloose = 0;
    AK4HT = 0;
    int nskipped = 0;
    float deltaR1 = -1; float deltaR2 = -1;
    float ptRel1 = -1; float ptRel2 = -1;
    int ind1 = -1; int ind2 = -1;
    vector<pair<double,int>> jetptindpair;
    jetptindpair.clear();
    listOfJets.clear();
    for(int ijet=0; ijet < nJet; ijet++){
      
      // ----------------------------------------------------------------------------
      // Basic cuts
      // ----------------------------------------------------------------------------
      
      if(!(Jet_jetId[ijet] > 1 && Jet_pt[ijet]>jetPtCut && fabs(Jet_eta[ijet])<jetEtaCut)){continue;}

      if (Jet_btagDeepB[ijet] > 0.4184) NJetsDeepCSVmed += 1;
      if (Jet_btagDeepB[ijet] > 0.1241) NJetsDeepCSVloose += 1;
       
      //Checking overlapping with lepton -- ASKED ON B2G SELECTIONS HN about what's best
      jet_lv.SetPtEtaPhiM(Jet_pt[ijet],Jet_eta[ijet],Jet_phi[ijet],Jet_mass[ijet]);
      float deltaR = lepton_lv.DeltaR(jet_lv);
      if(deltaR < 0.4) continue;

      listOfJets.push_back(std::make_pair(ijet, jet_lv));
      jetptindpair.push_back(std::make_pair(Jet_pt[ijet],ijet));
      NJets+=1;
      AK4HT+=Jet_pt[ijet];   

    }
    if(nskipped > 1){
      std::cout << "------------------------------------------------------------" << std::endl;
      std::cout << "Weird thing! " << nskipped << " AK4s had DR < 0.4" << std::endl;
      std::cout << "dr1 = " << deltaR1 << ", dr2 = " << deltaR2 << std::endl;
      std::cout << "pt1 = " << ptRel1 << ", pt2 = " << ptRel2 << std::endl;
      std::cout << "------------------------------------------------------------" << std::endl;
    }

    // ----------------------------------------------------------------------------
    // Loop over AK8 jets for calculations and pt ordering pair
    // ----------------------------------------------------------------------------
    NJetsAK8 = 0;  //Part of output tree later
    vector<pair<double,int>> jetak8ptindpair;    
    jetak8ptindpair.clear();  
    nskipped = 0;      
    deltaR1 = -1; deltaR2 = -1;
    ptRel1 = -1; ptRel2 = -1;
    for(unsigned int ijet=0; ijet < nFatJet; ijet++){  
      // --------------------------------------------------------------------------
      // Basic cuts
      // --------------------------------------------------------------------------
      
      if(!(FatJet_jetId[ijet] > 1 && FatJet_pt[ijet]>ak8PtCut && fabs(FatJet_eta[ijet])<ak8EtaCut && 
	   fabs(SubJet_eta[FatJet_subJetIdx1[ijet]])<ak8SJEtaCut && fabs(SubJet_eta[FatJet_subJetIdx2[ijet]])<ak8SJEtaCut)) {continue;}
      
      //Checking if overlapping with lepton
      ak8_lv.SetPtEtaPhiM(FatJet_pt[ijet],FatJet_eta[ijet],FatJet_phi[ijet],FatJet_mass[ijet]);
      float deltaR = lepton_lv.DeltaR(ak8_lv);
      if(deltaR < 0.8) continue;

      // ----------------------------------------------------------------------------      
      // Counters and pt ordering pair                                                     
      // ----------------------------------------------------------------------------   
      
      NJetsAK8 += 1;
      jetak8ptindpair.push_back(std::make_pair(FatJet_pt[ijet],ijet));
    }
    if(nskipped > 1){
      std::cout << "------------------------------------------------------------" << std::endl;
      std::cout << "Weird thing! " << nskipped << " AK8s had DR < 0.8" << std::endl;
      std::cout << "dr1 = " << deltaR1 << ", dr2 = " << deltaR2 << std::endl;
      std::cout << "pt1 = " << ptRel1 << ", pt2 = " << ptRel2 << std::endl; 
      std::cout << "------------------------------------------------------------" << std::endl;
    }     
    
    if(NJets < 3 || (NJets == 3 && NJetsAK8 < 1)) continue;
    passJetCut++;
    
    // ----------------------------------------------------------------------------
    // PT Order the jet list
    // ----------------------------------------------------------------------------
    std::sort(jetptindpair.begin(), jetptindpair.end(), sortinrev2);
    Jet_mass_ptordered.clear();
    Jet_pt_ptordered.clear();
    Jet_eta_ptordered.clear();
    Jet_phi_ptordered.clear();
    Jet_bdisc_ptordered.clear();
    listOfJetBDisc.clear();
    for(unsigned int ijet=0; ijet < jetptindpair.size(); ijet++){
      Jet_mass_ptordered.push_back(Jet_mass[jetptindpair[ijet].second]);
      Jet_pt_ptordered.push_back(Jet_pt[jetptindpair[ijet].second]);
      Jet_eta_ptordered.push_back(Jet_eta[jetptindpair[ijet].second]);
      Jet_phi_ptordered.push_back(Jet_phi[jetptindpair[ijet].second]);
      Jet_bdisc_ptordered.push_back(Jet_btagDeepB[jetptindpair[ijet].second]);
      listOfJetBDisc.push_back(Jet_btagDeepB[jetptindpair[ijet].second]);
    }
    std::sort(listOfJets.begin(),listOfJets.end(),sortinrev);
      
    std::sort(jetak8ptindpair.begin(), jetak8ptindpair.end(), sortinrev2);
    FatJet_sdmass_ptordered.clear();
    FatJet_mass_ptordered.clear();
    FatJet_pt_ptordered.clear();
    FatJet_eta_ptordered.clear();
    FatJet_phi_ptordered.clear();
    FatJet_sj1_ptordered.clear();
    FatJet_sj2_ptordered.clear();
    FatJet_tau21_ptordered.clear();
    for(unsigned int ijet=0; ijet < jetak8ptindpair.size(); ijet++){
      FatJet_sdmass_ptordered.push_back(FatJet_msoftdrop[jetak8ptindpair[ijet].second]);
      FatJet_mass_ptordered.push_back(FatJet_mass[jetak8ptindpair[ijet].second]);
      FatJet_pt_ptordered.push_back(FatJet_pt[jetak8ptindpair[ijet].second]);
      FatJet_eta_ptordered.push_back(FatJet_eta[jetak8ptindpair[ijet].second]);
      FatJet_phi_ptordered.push_back(FatJet_phi[jetak8ptindpair[ijet].second]);
      FatJet_tau21_ptordered.push_back(FatJet_tau2[jetak8ptindpair[ijet].second]/FatJet_tau1[jetak8ptindpair[ijet].second]);
      FatJet_sj1_ptordered.push_back(FatJet_subJetIdx1[jetak8ptindpair[ijet].second]);
      FatJet_sj2_ptordered.push_back(FatJet_subJetIdx2[jetak8ptindpair[ijet].second]);
    }

    FatJet_index1_ptordered.clear();
    FatJet_index2_ptordered.clear();
    SubJet_bdisc_fjordered.clear();
    SubJet_mass_fjordered.clear();
    SubJet_pt_fjordered.clear();
    SubJet_eta_fjordered.clear();
    SubJet_phi_fjordered.clear();
    int nsubstotal = 0;
    for(unsigned int ijet = 0; ijet < FatJet_pt_ptordered.size(); ijet++){
      int sj1 = FatJet_sj1_ptordered[ijet];
      int sj2 = FatJet_sj2_ptordered[ijet];
      if(sj1 >= 0){
	FatJet_index1_ptordered.push_back(nsubstotal); // 1 per fatjet, marks location in longer vecs of sj2 
	SubJet_pt_fjordered.push_back(SubJet_pt[sj1]*(1 - SubJet_rawFactor[sj1]));
	SubJet_mass_fjordered.push_back(SubJet_mass[sj1]*(1 - SubJet_rawFactor[sj1]));
	SubJet_eta_fjordered.push_back(SubJet_eta[sj1]);
	SubJet_phi_fjordered.push_back(SubJet_phi[sj1]);
	SubJet_bdisc_fjordered.push_back(SubJet_btagDeepB[sj1]);
	nsubstotal++;
      }else{FatJet_index1_ptordered.push_back(-1);}
      if(sj2 >= 0){
	FatJet_index2_ptordered.push_back(nsubstotal); // 1 per fatjet, marks location in longer vecs of sj2 
	SubJet_pt_fjordered.push_back(SubJet_pt[sj2]*(1 - SubJet_rawFactor[sj2]));
	SubJet_mass_fjordered.push_back(SubJet_mass[sj2]*(1 - SubJet_rawFactor[sj2]));
	SubJet_eta_fjordered.push_back(SubJet_eta[sj2]);
	SubJet_phi_fjordered.push_back(SubJet_phi[sj2]);
	SubJet_bdisc_fjordered.push_back(SubJet_btagDeepB[sj2]);
	nsubstotal++;
      }else{FatJet_index2_ptordered.push_back(-1);}
    }
    

    // ----------------------------------------------------------------------------
    // Replace AK4s with AK8 subjets from W jets
    // ----------------------------------------------------------------------------

    usesWtag = false;
    FatJet_sdmassraw_ptordered.clear();
    FatJet_sdmasscorr_ptordered.clear();
    float mintau21 = 99;
    int mintau21jet = -1;
    TLorentzVector sj1_tlv;
    TLorentzVector sj2_tlv;
    for(unsigned int ijet=0; ijet < FatJet_pt_ptordered.size(); ijet++){

      float genCorr  = 1.;
      float recoCorr = 1.;
      float totalWeight = 1.;      
      genCorr =  puppisd_corrGEN->Eval(FatJet_pt_ptordered[ijet]);
      if(fabs(FatJet_eta_ptordered[ijet]) <= 1.3 ) recoCorr = puppisd_corrRECO_cen->Eval(FatJet_pt_ptordered[ijet]);
      else recoCorr = puppisd_corrRECO_for->Eval(FatJet_pt_ptordered[ijet]);      
      totalWeight = genCorr * recoCorr;

      int sj1 = FatJet_index1_ptordered[ijet];
      int sj2 = FatJet_index2_ptordered[ijet];
      if(sj1 < 0 || sj2 < 0){
	FatJet_sdmassraw_ptordered.push_back(-1);
	FatJet_sdmasscorr_ptordered.push_back(-1);
	continue;
      }

      sj1_tlv.SetPtEtaPhiM(SubJet_pt_fjordered[sj1],SubJet_eta_fjordered[sj1],SubJet_phi_fjordered[sj1],SubJet_mass_fjordered[sj1]);
      sj2_tlv.SetPtEtaPhiM(SubJet_pt_fjordered[sj2],SubJet_eta_fjordered[sj2],SubJet_phi_fjordered[sj2],SubJet_mass_fjordered[sj2]);

      float sdmasscorr = (sj1_tlv+sj2_tlv).M()*totalWeight;
      FatJet_sdmassraw_ptordered.push_back((sj1_tlv+sj2_tlv).M());
      FatJet_sdmasscorr_ptordered.push_back(sdmasscorr);

      if(!(sdmasscorr < 100 && sdmasscorr> 60)) continue;

      if(FatJet_tau21_ptordered[ijet] < mintau21){
	mintau21 = FatJet_tau21_ptordered[ijet];
	mintau21jet = ijet;
      }
    }

    TLorentzVector subjet1; 
    TLorentzVector subjet2; 
    if(mintau21jet > -1){

      int ijet = mintau21jet;
      ak8_lv.SetPtEtaPhiM(FatJet_pt_ptordered[ijet],FatJet_eta_ptordered[ijet],FatJet_phi_ptordered[ijet],FatJet_mass_ptordered[ijet]);
      subjet1.SetPtEtaPhiM(SubJet_pt[FatJet_sj1_ptordered[ijet]], SubJet_eta[FatJet_sj1_ptordered[ijet]], SubJet_phi[FatJet_sj1_ptordered[ijet]], SubJet_mass[FatJet_sj1_ptordered[ijet]]);
      subjet2.SetPtEtaPhiM(SubJet_pt[FatJet_sj2_ptordered[ijet]], SubJet_eta[FatJet_sj2_ptordered[ijet]], SubJet_phi[FatJet_sj2_ptordered[ijet]], SubJet_mass[FatJet_sj2_ptordered[ijet]]);
    
      //--------Count overlapping AK4--------------
      //Find indexes of overlapping AK4s
      TLorentzVector jet; // Jet currently checking for single
      TLorentzVector jet1; //Best of all the singles
      TLorentzVector jetPair1; //Best of all the pairs (first jet)
      TLorentzVector jetPair2; //Best of all the pairs (second jet)
      TLorentzVector testPair;  //Jet currently checking for pair
      TLorentzVector testPair2;
      int bestPair [2] = {-1,-1}; //Indexes of best pair
      int bestSingle = -1; //Index of best Single

      // Find the minimum DR for single jets
      float minDRsingle = 99999;
      for (int i = 0; i < listOfJets.size(); i++){
        jet = listOfJets[i].second;        
	if(jet.DeltaR(ak8_lv) < minDRsingle){
	  minDRsingle = jet.DeltaR(ak8_lv);
	  jet1 = jet;
	  bestSingle = i;
	}
      }

      // Find the minimum DR for a pair of jets
      float minDRpair = 99999;
      for(int i = 0; i < listOfJets.size(); i++){
	jet = listOfJets[i].second;
	for(int j = i+1; j < listOfJets.size(); j++){
	  testPair = listOfJets[j].second;
	  
	  if((jet+testPair).DeltaR(ak8_lv) < minDRpair){
	    minDRpair = (jet+testPair).DeltaR(ak8_lv);
	    jetPair1 = jet;
	    jetPair2 = testPair;
	    bestPair[0] = i;
	    bestPair[1] = j;
	  }
	}
      }
    
      if (minDRsingle <= minDRpair && minDRsingle < 0.05){
	usesWtag = true;
        overlap1++;
	
        listOfJetBDisc.erase(listOfJetBDisc.begin()+bestSingle);
        listOfJets.erase(listOfJets.begin()+bestSingle);
        listOfJetBDisc.push_back(SubJet_btagDeepB[FatJet_sj1_ptordered[ijet]]);
        listOfJetBDisc.push_back(SubJet_btagDeepB[FatJet_sj2_ptordered[ijet]]);
        listOfJets.push_back(std::make_pair(nJet+1, subjet1));
        listOfJets.push_back(std::make_pair(nJet+2, subjet2));

      }
      else if (minDRsingle > minDRpair && minDRpair < 0.05){
	usesWtag = true;
        overlap2++;

        listOfJetBDisc.erase(listOfJetBDisc.begin()+bestPair[1]);
        listOfJetBDisc.erase(listOfJetBDisc.begin()+bestPair[0]);
        listOfJets.erase(listOfJets.begin()+bestPair[1]);
        listOfJets.erase(listOfJets.begin()+bestPair[0]);
        listOfJetBDisc.push_back(SubJet_btagDeepB[FatJet_sj1_ptordered[ijet]]);
        listOfJetBDisc.push_back(SubJet_btagDeepB[FatJet_sj2_ptordered[ijet]]);
        listOfJets.push_back(std::make_pair(nJet+1, subjet1));
        listOfJets.push_back(std::make_pair(nJet+2, subjet2));

      }
    }   

    // ------------------------------
    // Find pt ordering of new list
    // ------------------------------

    ptOrderedJets.clear();
    for (unsigned int ijet = 0; ijet < listOfJets.size(); ijet++){
      ptOrderedJets.push_back(std::make_pair(listOfJets[ijet].second.Pt(),ijet));
    }
    std::sort(ptOrderedJets.begin(), ptOrderedJets.end(), sortinrev2);
 

    // ----------------------------------------------------------------------------
    // Check the hybrid jet list
    // ----------------------------------------------------------------------------
    if(ptOrderedJets[0].first>100){
      if(isElectron){ passElAll++;}
      if(isMuon){ passMuAll++;}
    }
    
    if((ptOrderedJets[1].first>70)){
      if(isMuon) passMuAll2++;
      if(isElectron) passElAll2++;
    }
    
    if(ptOrderedJets.size()>=4){
      if(isElectron){ passElAll3++;}
      if(isMuon){ passMuAll3++;}
    }
    
    if(!(ptOrderedJets[0].first>100 && ptOrderedJets[1].first>70 && ptOrderedJets.size()>=4)){continue;}
    
    passHybridJetCut++; // this is ALL PRESELECTION
    
    
    //////////////////////////////////////////////////////
    // HITFIT FURTHER JET SELECTION
    //////////////////////////////////////////////////////

    bool passHitFitJetPt = false;
    float stmax = ptOrderedJets[0].first + ptOrderedJets[1].first + ptOrderedJets[2].first + ptOrderedJets[3].first + leppt + MET_pt; 
    if(ptOrderedJets[0].first > 200 && ptOrderedJets[1].first > 100 && ptOrderedJets[2].first > 100 && stmax > 1000) passHitFitJetPt = true;
    
    float probChisq = -1;
    tMass = -1;
    chi2 = 999;
    sigmaM = -1;
    HybridJet_pt.clear();
    HybridJet_eta.clear();
    HybridJet_phi.clear();
    HybridJet_mass.clear();

    fittedJetPt.clear();
    fittedJetEta.clear();
    fittedJetPhi.clear();
    fittedJetMass.clear();
    fittedLepPt =-9;
    fittedLepEta=-9;
    fittedLepPhi=-9;
    fittedLepMass=-9;
    fittedMETPt =-9;
    fittedMETEta=-9;
    fittedMETPhi=-9;
    fittedMETMass=-9;
    tMass = -9;
    chi2 = 999;
    probchi2 = -9;
    sigmaM = 0;

    
    if(passHitFitJetPt) passHitFitJetCut += 1;
    
    // ------------------------------------------------------------------
    // HITFIT
    // ------------------------------------------------------------------
    
    //Clear internal state
    HitFit.clear();
    
    //Add lepton, jets, met
    HitFit.AddLepton(lepton_lv,isElectron);
    nJetsHitFit = ptOrderedJets.size();
    for(int i=0;(i < ptOrderedJets.size() && i < 5);i++){
      int isSubJet = 0;
      if(listOfJets[ptOrderedJets[i].second].first>nJet) isSubJet=1;
      HitFit.AddJet(listOfJets[ptOrderedJets[i].second].second,isSubJet);
      
      HybridJet_pt.push_back(listOfJets[ptOrderedJets[i].second].second.Pt());
      HybridJet_eta.push_back(listOfJets[ptOrderedJets[i].second].second.Eta());
      HybridJet_phi.push_back(listOfJets[ptOrderedJets[i].second].second.Phi());
      HybridJet_mass.push_back(listOfJets[ptOrderedJets[i].second].second.M());
    }
    double px = MET_pt * cos(MET_phi);
    double py = MET_pt * sin(MET_phi);
    HitFit.SetMet(px,py);
    
    //Number of all permutation of the event
    size_t nHitFit = 0;
    nHitFit = HitFit.FitAllPermutation();
    
    std::vector<hitfit::Fit_Result> hitfitResult = HitFit.GetFitAllPermutation();
    
    // mass of the top quark after fit
    std::vector<double> fittedTopMass;
    
    // the chi-square of the fit
    // a negative value means the fit of a particular permutation did not converge
    std::vector<double> fitChi2;
    
    float minChisq2M = 1e9;
    float minChisqML = 1e9;
    float minChisq1M = 1e9;
    float minChisqBadFit = 1e9;
    int indexMinChisq2M = -1;
    int indexMinChisqML = -1;
    int indexMinChisq1M = -1;
    int indexMinChisqBadFit=-1;
    double SLST [nHitFit];
    
    //cout << "HITFIT returned " << nHitFit << " fits" << endl;
    for (size_t fit = 0 ; fit != nHitFit ; ++fit) {
      ntotal++;
      hitfit::Lepjets_Event fittedEvent = hitfitResult[fit].ev();
      double fittedST =0;
      double fittedSL = 0;
      bool badFit =false;
      bool badFitST = false;
      bool badFitMW = false;
      bool badFitMW2 = false;
      bool badFitB = false;
      bool badFitSL = false;
      bool badFitJet1 = false;
      bool badFitJet2 = false;
      bool badFitJet3 = false;
      double WMass=0;
      double fittedWMass = 0;
      double st=0;
      TLorentzVector Wjet;
      Fourvec fittedWjet;
      int firstWindex = -1;
      int secondWindex = -1;
      float highWpt = 0;
      int bTagM=0;
      int bTagL=0;
      
      vector<int> jetTypes = hitfitResult[fit].jet_types();
      //cout << "--------------" << endl;
      for(int i=0;i<jetTypes.size();i++){	  
	if (jetTypes[i]==11){  // leptonic b
	  if (listOfJetBDisc[ptOrderedJets[i].second] > 0.4184) bTagM++;
	  else if (listOfJetBDisc[ptOrderedJets[i].second] > 0.1241) bTagL++;
	  st = st+listOfJets[ptOrderedJets[i].second].second.Pt();
	  fittedST = fittedST + fittedEvent.jet(i).p().perp();
	  fittedSL = fittedSL + fittedEvent.jet(i).p().pz();
	}
	else if (jetTypes[i]==12){ // hadronic b
	  if (listOfJetBDisc[ptOrderedJets[i].second] > 0.4184) bTagM++;
	  else if (listOfJetBDisc[ptOrderedJets[i].second] > 0.1241) bTagL++;
	  st = st+listOfJets[ptOrderedJets[i].second].second.Pt();
	  fittedST = fittedST + fittedEvent.jet(i).p().perp();
	  fittedSL = fittedSL + fittedEvent.jet(i).p().pz();
	}
	else if(jetTypes[i]==13 || jetTypes[i]==14){	    // hadronic W jets.
	  if (firstWindex < 0){
	    fittedWjet = fittedEvent.jet(i).p();
	    Wjet = listOfJets[ptOrderedJets[i].second].second;
	    firstWindex = i;
	  }
	  else{
	    secondWindex = i;
	    fittedWMass = (fittedWjet + fittedEvent.jet(i).p()).m();
	    WMass = (Wjet+listOfJets[ptOrderedJets[i].second].second).M();
	  }
	}
      }
      st = st+listOfJets[ptOrderedJets[firstWindex].second].second.Pt()+listOfJets[ptOrderedJets[secondWindex].second].second.Pt();
      fittedST = fittedST + fittedEvent.jet(firstWindex).p().perp() + fittedEvent.jet(secondWindex).p().perp();
      fittedSL = fittedSL + fittedEvent.jet(firstWindex).p().pz() + fittedEvent.jet(secondWindex).p().pz();
      
      jetTypes.clear();
      st = st + MET_pt + lepton_lv.Pt(); 
      fittedST = fittedST + fittedEvent.lep(0).p().perp() + fittedEvent.met().perp();
      fittedSL = fittedSL + fittedEvent.lep(0).p().pz() + fittedEvent.met().pz();
      SLST[fit]=fittedSL/fittedST;

      if (bTagM < 1){badFitB=true; badFit = true;}       
      if (fittedSL/fittedST >= 1.1){badFitSL=true; badFit = true;}     
      if(hitfitResult[fit].chisq() < 0) badFit = true;
      
      if (badFit==true){
	nfailed++;
	if(hitfitResult[fit].chisq()==-999){num999++;}
	else if(hitfitResult[fit].chisq()==-998){num998++;}
	else if(hitfitResult[fit].chisq()==-997){num997++;}
	else if(hitfitResult[fit].chisq()==-996){num996++;}
	else if(hitfitResult[fit].chisq()==-995){num995++;}
	if(badFitJet1) numLowHadB++;
	if(badFitJet2) numLowLepB++;
	if(badFitJet3) numLowHadW++;
	if(badFitST) numLowST++;
	if(badFitMW) numBadW++;
	if(badFitMW2) numBadW2++;
	if(badFitB) numLowNB++;
	if(badFitSL) numBadSLST++;
      }else ngood++;
      
      if (hitfitResult[fit].chisq() >= 0){
	if(!badFit){
	  if(bTagM == 2 && bTagL == 0){
	    good2M++;
	    if( hitfitResult[fit].chisq() < minChisq2M){
	      minChisq2M = hitfitResult[fit].chisq();
	      indexMinChisq2M = fit;
	    }
	  }else if(bTagM == 1 && bTagL == 1){
	    goodML++;
	    if(hitfitResult[fit].chisq() < minChisqML){
	      minChisqML = hitfitResult[fit].chisq();
	      indexMinChisqML = fit;
	    }
	  }else if(bTagM == 1 && bTagL == 0){
	    good1M++;
	    if(hitfitResult[fit].chisq() < minChisqML){
	      minChisq1M = hitfitResult[fit].chisq();
	      indexMinChisq1M = fit;
	    }
	  }else cout << "GOT WEIRD TAGS: btagM = " << bTagM << ", btagL = " << bTagL << endl;
	}else{
	  if(hitfitResult[fit].chisq() < minChisqBadFit){
	    minChisqBadFit = hitfitResult[fit].chisq();
	    indexMinChisqBadFit = fit;
	  }
	}
      }
    }
    
    int indexMinChisq = -1;
    if(indexMinChisq2M > -1) indexMinChisq = indexMinChisq2M;
    else if(indexMinChisqML > -1) indexMinChisq = indexMinChisqML;
    else if(indexMinChisq1M > -1) indexMinChisq = indexMinChisq1M;
    
    if(indexMinChisq > 0){
      passAnyChi2 += 1;
      probchi2 = exp(-0.5 * hitfitResult[indexMinChisq].chisq());
      if(probchi2 > 0.001) passProbChisq++;
      
      chi2 = hitfitResult[indexMinChisq].chisq();
      tMass = hitfitResult[indexMinChisq].mt();
      sigmaM = hitfitResult[indexMinChisq].sigmt();
      SLDivideST = SLST[indexMinChisq];
      fittedLepPt = hitfitResult[indexMinChisq].ev().lep(0).p().perp();
      fittedLepEta = hitfitResult[indexMinChisq].ev().lep(0).p().eta();
      fittedLepPhi = hitfitResult[indexMinChisq].ev().lep(0).p().phi();
      fittedLepMass = hitfitResult[indexMinChisq].ev().lep(0).p().m();
      fittedMETPt =  hitfitResult[indexMinChisq].ev().met().perp();
      fittedMETEta = hitfitResult[indexMinChisq].ev().met().eta();
      fittedMETPhi = hitfitResult[indexMinChisq].ev().met().phi();
      fittedMETMass = hitfitResult[indexMinChisq].ev().met().m();
      vector<int> jetTypes = hitfitResult[indexMinChisq].jet_types();
      if(jetTypes.size() < 4) cout << "WEIRD FITTED JET TYPES: " << jetTypes.size() << endl;
      for (int i=0;i<jetTypes.size();i++){	 
	if(jetTypes[i]==11 || jetTypes[i]== 12 || jetTypes[i]== 13 || jetTypes[i]== 14 ){	    	    
	  fittedJetPt.push_back(hitfitResult[indexMinChisq].ev().jet(i).p().perp());
	  fittedJetEta.push_back(hitfitResult[indexMinChisq].ev().jet(i).p().eta());
	  fittedJetPhi.push_back(hitfitResult[indexMinChisq].ev().jet(i).p().phi());
	  fittedJetMass.push_back(hitfitResult[indexMinChisq].ev().jet(i).p().m());
	}
      }
    }   

    // -------------------------------------------------------------------
    // Adding Masses
    // -------------------------------------------------------------------
    
    //Lepton + bJet + MET Masses
    TLorentzVector together;
    TLorentzVector MET_lv;
    MET_lv.SetPtEtaPhiM(MET_pt,0,MET_phi,0);
    double minMass = 999999999;
    totalMass=0;
    for (int i=0;i<listOfJets.size();i++){
      if(listOfJets[i].first>nJet){continue;}
      if(Jet_btagDeepB[listOfJets[i].first] <= 0.4184){continue;}
      if(listOfJets[i].second.Pt() <= 100) {continue;}
      if((listOfJets[i].second+lepton_lv).M() < minMass ){
	minMass = (listOfJets[i].second+lepton_lv).M();
	together = (listOfJets[i].second+lepton_lv);
      }
    }
    if (minMass == 999999999) {totalMass=-1;}
    else {totalMass = (together + MET_lv).M();}
    
    // 3 Jets
    TLorentzVector all3Jets;
    int num3Jets=0;
    double sum3JetMass=0;
    for(int i=0;i<listOfJets.size();i++){
      if (listOfJets[i].second.Pt() <= 30) {continue;}
      for(int j=(i+1);j<listOfJets.size();j++){
	if( listOfJets[j].second.Pt()<=30 ||(listOfJets[i].second.Pt() <= 100 && listOfJets[j].second.Pt() <= 100)) {continue;}
	for(int k=(j+1);k<listOfJets.size();k++){
	  if(listOfJets[k].second.Pt() <= 30) {continue;} //If not greater than 30 pt skip
	  if(listOfJets[i].second.Pt() <= 200 && listOfJets[j].second.Pt() <= 200 && listOfJets[k].second.Pt() <= 200){continue;} //If none greater then 200
	  if((listOfJets[i].second.Pt() <= 100 && listOfJets[k].second.Pt() <= 100)||(listOfJets[i].second.Pt() <= 100 && listOfJets[j].second.Pt() <= 100)||(listOfJets[k].second.Pt() <= 100 && listOfJets[j].second.Pt() <= 100)) {continue;} // If not two greater than 100 pt skip
	  all3Jets = listOfJets[i].second+listOfJets[j].second+listOfJets[k].second;
	  sum3JetMass = sum3JetMass + all3Jets.M();
	  num3Jets++;
	}
      }
    }
    threeJetMass =0;
    if(num3Jets>0){threeJetMass = sum3JetMass/num3Jets;}
    
    //Counting events
    if (isElectron){numEventsElectron++;}
    if (isMuon){numEventsMuon++;}
    
    //Filling output Tree
    outputTree->Fill();
    // ---------------------------------------------------------------------------
  }

  outputTree->Write();
  nEventsWeighted->Write();

  // Record end time
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  cout << "Elapsed time: " << elapsed.count() << endl;
  cout << "Pass bWbW Cut: "<<passbwbw<<endl;
  cout << "Pass isoTrack Cut: "<<passIso<<endl;
  cout << "Pass flag Cut: " <<passflag<<endl;
  cout << "Pass Lepton: "<<passLepton<<endl;
  cout << "\t Pass Electron: "<<passEl << endl;
  cout << "\t Pass Muon: "<<passMu<<endl;
  cout << "Pass HLT: "<<passHLT<<endl;
  cout << "Pass MET Cut: "<<passMET<<endl;
  cout << "Pass Jet Count Cut: "<<passJetCut<<endl;

  cout << "Pass Hybrid Jet List Cuts (AND of previous 3): "<<passHybridJetCut<<endl;
  cout << "\t Number of passed Muon events: " << numEventsMuon << endl;
  cout << "\t Number of passed Electron events: " << numEventsElectron << endl;
  cout << "Pass HITFIT Jet List Cuts (200,100,100): " << passHitFitJetCut<<endl;
  cout << "Number of events with Prob(chi squared) > 0.1%: "<<passProbChisq<<endl; 
  cout << "Number of events with any chi2 > 0: " <<passAnyChi2<<endl;
  cout << "Number of fits failed HITFIT out of: "<<ntotal << " = " << nfailed << endl;

  cout << "\t -998: "<<num998<< "/" << ntotal << endl;
  cout << "\t -997: "<<num997<< "/" << ntotal << endl;
  cout << "\t -996: "<<num996<< "/" << ntotal << endl;
  cout << "\t -995: "<<num995<< "/" << ntotal << endl;
  cout << "\t low lep B pT: "<<numLowLepB<< "/" << ntotal << endl;
  cout << "\t low had B pT: "<<numLowHadB<< "/" << ntotal << endl;
  cout << "\t low had W pT: "<<numLowHadW<< "/" << ntotal << endl;
  cout << "\t low ST: : "<<numLowST<< "/" << ntotal << endl;
  cout << "\t bad MW: : "<<numBadW<< "/" << ntotal << endl;
  cout << "\t bad umwhad from HF: : "<<numBadW2<< "/" << ntotal << endl;
  cout << "\t low btags: : "<<numLowNB<< "/" << ntotal << endl;
  cout << "\t bad SL/ST: : "<<numBadSLST<< "/" << ntotal << endl;
  cout << "Number of fits passed HITFIT out of: " <<ntotal << " = " << ngood << endl;
  cout << "\t good 2M: "<<good2M<<"/"<<ntotal << endl;
  cout << "\t good ML: "<<goodML<<"/"<<ntotal << endl;
  cout << "\t good 1M: "<<good1M<<"/"<<ntotal << endl;
   
}
