// modify line 3 for using MakeWario.C and CombineWarios.C with the Neutron or DP data, or both
// lines 5, 6, 7, 8 are obsolete
// contains the Hebing structure

std::string Source = "none";     // "neutron" or "none" or "*"
std::string recoil_type = "NR"; // "ER", "NR"
std::string cut = ""; // "", "_++", "_nruns"
std::string bump = "both"; // "both", "upper", "lower"
bool multi = false;   // multiple or single scattering

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/HybridResult.h"
#include "RooStats/HybridPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/UniformProposal.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/NumberCountingPdfFactory.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalFunction.h"
#include "RooStats/ProposalHelper.h"
#include "RooFitResult.h"
#include "TGraph2D.h"
#include "TTree.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

struct HEBING_STRUCT {
  // S1 variables
  int S1unixTime;
  unsigned int S1TTT; // TTT value
  float S1TTR;
  float S1LY;
  float S1F90;
  float S1F90Top;
  float S1F90Bottom;
  float S1COGX;
  float S1COGY;
  float S1maxFracPMT;
  int S1Nx;
  int nS1perS2;
  int nS1perS2max;

  // S2 variables
  int nS2;
  int matchedS2;
  // float S2LY[MAXMATCHEDS2];
  float S2LY;
  float S2maxLY;
  float S2maxDT;
  float S2maxrecoLSX;
  float S2maxrecoLSY;
  float S2maxrecoLSCHI2;
  float S2maxmaxFracPMT;
  float S2maxLoS1L;     // log10(S2maxL/S1L)
  float S2maxLoS1Lcorr; // log10(S2maxL/S1L) DRIFT TIME CORRECTED
  float S2maxQMatch;    // linearize TTR vs Drifttime relation
  float S2maxFracS2L;   // fraction of S2 total light of max S2
}; 