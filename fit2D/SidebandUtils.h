//
// Programmer: Created by Alex Stauffer on 07.09.21.
// File: SidebandUtils.h
// Purpose: Prototypes for SidebandUtils.cpp
//

#ifndef SIDEBANDUTILS_H
#define SIDEBANDUTILS_H

// Root header files
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooProdPdf.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TText.h"
#include "TPaveText.h"
#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/HybridResult.h"
#include "RooStats/HybridPlot.h"
#include <RooStats/LikelihoodInterval.h>
// General header files
#include <sys/stat.h>

// Namespaces
using namespace RooFit;
using namespace RooStats;
using namespace TMath;
using namespace std;

// Variables
extern std::string InputPath;
extern std::string OutputPath;

// Functions
std::vector<std::string> GetStdoutFromCommand(std::string cmd);
std::unique_ptr<RooDataSet> GetDataFromWorkspaces(std::string Path, std::string DataType, int MinRunNr, int MaxRunNr, int nr_S1bin);


#endif // SIDEBANDUTILS_H
