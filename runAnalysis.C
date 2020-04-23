#include "qWqW_NanoAnalysis.C"
#include "../TopQuarkAnalysis/TopHitFit/src/RunHitFit.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/PatElectronHitFitTranslator.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/PatMuonHitFitTranslator.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/PatJetHitFitTranslator.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/PatMETHitFitTranslator.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Base_Constrainer.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Constraint.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Fit_Result.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Lepjets_Event.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Resolution.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/fourvec.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Constraint_Intermed.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Fit_Result_Vec.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Lepjets_Event_Jet.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/gentop.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Chisq_Constrainer.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Defaults_Text.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Fit_Results.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Lepjets_Event_Lep.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Top_Decaykin.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/matutil.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Constrained_Top.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/EtaDepResElement.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Fourvec_Constrainer.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Objpair.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Top_Fit.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/EtaDepResolution.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Fourvec_Event.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Pair_Table.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Refcount.cc"
#include "../TopQuarkAnalysis/TopHitFit/src/Vector_Resolution.cc"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// PUT ALL YOUR INCLUDES

using namespace std;

void runAnalysis(TString inputFile, TString outputFile, Long64_t skipevents, Long64_t maxevents){

  bool isSignal = inputFile.Contains("TprimeTprime");
  bool isTTbar = (inputFile.Contains("TTTo") || inputFile.Contains("TT_Tune") || inputFile.Contains("_Mtt"));  
  
  qWqW_NanoAnalysis m = qWqW_NanoAnalysis(inputFile,outputFile); // instantiate your class with whatever it needs
  
  m.Loop(isSignal,isTTbar,skipevents,maxevents); // pass whatever arguments your loop function needs

}
