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

#include <sys/stat.h>

#define IS_ON_ESSOS 1
// 0 (no) or 1 (yes)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Similar to MakeWario.C, these files ("ReadData.C" and "ReadDataSimple.C") takes as input the Satoshi runs and perform
// cuts on the events and store histograms in Wario output files with the same naming scheme as the input Satoshi file
// However, these files are primarily though to prepare the data for ER rejection and NR event-counting analysis
// in S1 bins...

// Usage: "ReadData.C" will process all Satoshi runs one by one, takes several hours
// Usage: "ReadDataSimple.C" processes only one Satoshi run and is supposed to be processed in parallel
// in conjunction with "RunReadDataInParallel.sh", "RunReadDataInParallel.submit"
// and "RunReadDataInParallel_arguments.txt"

/////////////////////////////////////////////   USER INPUT   ///////////////////////////////////////////////////////////
// Run in batch-mode (no graphics):
// $ root -l -b ReadData.C

// Filepath In (where the Satoshi runs are located)
#if IS_ON_ESSOS == 1
  bool is_on_essos = true;
  // std::string InputPath = "/mnt/raid/users/duylai/Satoshi_DP_runs/";  // on essos
  std::string InputPath = "/mnt/raid/users/duylai/Satoshi_neutron_runs/";  // on essos
#elif IS_ON_ESSOS == 0
  bool is_on_essos = false;
  std::string InputPath = "/media/alex/Seagate_HDD/satoshiV2/";  // local
#endif

// Input details
std::string DataType = "*";  // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
// std::string DataType = "DP[12]";  // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
std::string Source = "*";  // "none" or "neutron" or "*"
int MinRunNr = 0;  // Minimum run number of runs to be processed
int MaxRunNr = 169;  // Maximum run number of runs to be processed

// Output path (where the Wario runs should be stored to)
#if IS_ON_ESSOS == 1
  std::string OutputPath("/mnt/raid/users/duylai/warios_bkgs_template/");  // on essos
#elif IS_ON_ESSOS == 0
  // std::string OutputPath("/media/alex/Seagate_HDD/WarioV2_Sideband/");  // local
  std::string OutputPath("/home/alex/Desktop/PhD/ArDM_Neutron_source/WarioV3_Sideband/");  // local
#endif

// Cut values
// S1 bin edges -> [) S1 bin edges; last one is open
// std::vector<float> S1L_bin_edges{0.0, 100.0};
std::vector<float> S1L_bin_edges{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
// Event quality cut
float S1maxFracPMT_cut = 0.3;  // S1maxFracPMT <= value
// Match quality cut
float QMatch_cut_low = 0.75;  // QMatch >= value
float QMatch_cut_high = 1.25;  // QMatch <= value
// Single-scatterer cut
float S2maxFracS2L_cut = 1.0;  // S2maxFracS2L == value (for single-scatterers)
// Sniper cut
int nS1perS2max_cut = 1;
// Fiducial volume cut
float radial_cut = 300.0;  // [mm] S2maxLSR <= value
float S2maxDT_cut_low = 0.2;  // [ms] S2maxDT >= value
float S2maxDT_cut_high = 1.0;  // [ms] S2maxDT <= value

// Options
bool SuppressInfo = true;  // whether or not to suppress all the annoying info messages
bool SuppressWarning = true;  // whether or not to suppress all the annoying warning messages
bool SuppressError = false;  // whether or not to suppress all the annoying error messages

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace RooFit;
using namespace RooStats;

std::vector<std::string> GetStdoutFromCommand(std::string cmd);
double calculate_n_disc(double F90_temp, double logS2LoS1L_temp);

void ReadData() {
  // Suppress text output
  if (SuppressInfo) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);  // Everything below specified level is killed
    // enum MsgLevel { DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5 };
    gErrorIgnoreLevel = kWarning;
    // > const Int_t kPrint = 0;
    // > const Int_t kInfo = 1000;
    // > const Int_t kWarning = 2000;
    // > const Int_t kError = 3000;
    // > const Int_t kBreak = 4000;
    // > const Int_t kSysError = 5000;
    // > const Int_t kFatal = 6000;
  }
  if (SuppressWarning) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);  // Everything below specified level is killed
    // enum MsgLevel { DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5 };
    gErrorIgnoreLevel = kError;
  }
  if (SuppressError) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);  // Everything below specified level is killed
    // enum MsgLevel { DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5 };
    gErrorIgnoreLevel = kFatal;
  }

  // Putting together cuts as strings
  std::string Event_quality_cut = "S1maxFracPMT <= " + std::to_string(S1maxFracPMT_cut);
  std::string Match_quality_cut = "S2maxQMatch >= " + std::to_string(QMatch_cut_low) + " && S2maxQMatch <= " + std::to_string(QMatch_cut_high);
  std::string Single_scatterer_cut = "S2maxFracS2L == " + std::to_string(S2maxFracS2L_cut);
  std::string Sniper_cut = "nS1perS2max == " + std::to_string(nS1perS2max_cut);
  std::string Fiducial_volume_cut = "sqrt(S2maxLSX*S2maxLSX + S2maxLSY*S2maxLSY) <= " + std::to_string(radial_cut) +
                                    " && S2maxDT >= " + std::to_string(S2maxDT_cut_low) + " && S2maxDT <= " + std::to_string(S2maxDT_cut_high);
  // S1L Bins
  std::map<int, std::string> map_S1_bin_cuts;
  for (std::vector<float>::size_type i = 0; i != S1L_bin_edges.size(); i++) {
    if (i != (S1L_bin_edges.size() - 1)) {
      map_S1_bin_cuts[i] = "S1LY >= " + std::to_string(S1L_bin_edges[i]) + " && S1LY < " + std::to_string(S1L_bin_edges[i + 1]);
    }
    // else {
      // map_S1_bin_cuts[i] = "S1LY >= " + std::to_string(S1L_bin_edges[i]);
    // }
  }

  // Potentially make target folder
  int status = mkdir(OutputPath.c_str(), 0777);

  // Initialize variables (must match names from TTree) - Variables of interest
  RooRealVar S1LY("S1LY", "S1LY", 0.0);  // want this to be unbounded -> only initial value
  RooRealVar S1F90("S1F90", "S1F90", 0, 1);  // can be bounded... -> initial value is center (0.5)
  RooRealVar S2maxLoS1Lcorr("S2maxLoS1Lcorr", "S2maxLoS1Lcorr", 0.0);  // initial value doesn't matter
  // Initialize variables (must match names from TTree) - Variables to cut
  RooRealVar S1maxFracPMT("S1maxFracPMT", "S1maxFracPMT", 0, 1);
  RooRealVar S1TTR("S1TTR", "S1TTR", 0, 1);
  RooRealVar S2maxQMatch("S2maxQMatch", "S2maxQMatch", 0, 4);
  RooRealVar S2maxFracS2L("S2maxFracS2L", "S2maxFracS2L", 0, 1);
  RooRealVar S2maxLSX("S2maxLSX", "S2maxLSX", -500, 500);
  RooRealVar S2maxLSY("S2maxLSY", "S2maxLSY", -500, 500);
  RooRealVar S2maxDT("S2maxDT", "S2maxDT", 0, 1.3);
  RooRealVar nS1perS2max("nS1perS2max", "nS1perS2max", 0, 10);
  RooRealVar S2maxmaxFracPMT("S2maxmaxFracPMT", "S2maxmaxFracPMT", 0, 1);
  // RooRealVar S1TTR("S1TTR", "S1TTR", 0, 1);
  // RooRealVar S2maxLY("S2maxLY", "S2maxLY", 0, 5000);
  // Data selection criteria
  std::string preselection = "S1F90 >= 0 && S1F90 <= 1 "
                             "&& S1maxFracPMT >= 0 && S1maxFracPMT <= 1 "
                            //  "&& S2maxQMatch >= 0 && S2maxQMatch <= 4 "
                             "&& S2maxFracS2L >= 0 && S2maxFracS2L <= 1 "
                             "&& S2maxLSX >= -600 && S2maxLSY >= -600 "
                             "&& S2maxDT >= 0 && S1TTR >= 0.15 "
                             "&& S2maxmaxFracPMT >= 0 && S2maxLoS1Lcorr >= -6";
  // Combining general cuts and preselection
  std::string data_of_interest = preselection + " && " + Event_quality_cut + " && " + Match_quality_cut
                                + " && " + Single_scatterer_cut + " && " + Sniper_cut
                                //  + " && " + Fiducial_volume_cut
                                 ;
  // Container of all variables; can only add up to 9 upon construction
  auto variables_all = new RooArgSet(S1LY, S1F90, S2maxLoS1Lcorr, S1maxFracPMT, S2maxQMatch,
          S2maxFracS2L, S2maxLSX, S2maxLSY, S2maxDT);
  variables_all->add(S1TTR);
  variables_all->add(S2maxmaxFracPMT);
  variables_all->add(nS1perS2max);
  auto variables_reduced = new RooArgSet(S1LY, S1F90, S2maxLoS1Lcorr, S2maxFracS2L, S2maxQMatch, nS1perS2max);
  auto variables_of_interest = new RooArgSet(S1LY, S1F90, S2maxLoS1Lcorr);  // NDisc will be added later

  // Extract list of filenames with/out source
  std::string command = "ls " + InputPath + "SRun*-" + DataType + "-" + Source + ".root";
  std::vector <std::string> pre_file_list = GetStdoutFromCommand(command);
  int run_ID_string_length = 4;
  int start = InputPath.length() + 4;  // getting the starting index of run number 00xx

  // Looping through the files and keeping the ones with appropriate run number
  // mainly to get the number of files to be processed for progress information...
  int files_to_be_processed = 0;
  std::vector <std::string> file_list;
  for (std::vector<std::string>::iterator t = pre_file_list.begin(); t != pre_file_list.end(); ++t) {
    const int run_nr_temp = std::stoi(t->substr(start, run_ID_string_length));
    if (run_nr_temp >= MinRunNr && run_nr_temp <= MaxRunNr) {
      file_list.push_back(t->c_str());
      files_to_be_processed++;
    }
  }

  cout << pre_file_list[pre_file_list.size()-1] << endl;

  // Actually processing the files
  std::cout << "=======================================================================" << std::endl;
  std::cout << "                           Processing:" << std::endl;
  std::cout << "=======================================================================" << std::endl;
  int processed_files = 0;
  for (std::vector<std::string>::iterator t = file_list.begin(); t != file_list.end(); ++t) {
    const int run_nr_temp = std::stoi(t->substr(start, run_ID_string_length));
    // const std::string ending_temp = t->substr(start, 18);  // extracting e. g. "0089-DP2-none.root"
    const std::string ending_temp = t->substr(start, 21);  // extracting e. g. "0089-DP2-neutron.root"
    const std::string name_temp = t->substr(0, start + ending_temp.length());  // removing additional characters
    std::cout << "------- File: " << name_temp << " (" << processed_files + 1 << "/" << files_to_be_processed
              << ") -------" << std::endl;
    std::string OutputFileName_temp = OutputPath + "WarioRun" + ending_temp;

    // Actual processing

    // Open the file and get the TTree
    TFile * DataFile = new TFile(name_temp.c_str());
    TTree * DataTree = (TTree*)DataFile->Get("HebingTree");
    // Preselected Data
    // First step: import all variables (unfortunately not possible otherwise)
    RooDataSet * DataPreselected_all = new RooDataSet("data_preselected", "data_preselected",
            *variables_all, Import(*DataTree), Cut(data_of_interest.c_str()));
    RooDataSet * DataPreselected = (RooDataSet *)DataPreselected_all->reduce(SelectVars(*variables_reduced));
    // Memory
    delete DataPreselected_all;
    delete DataTree;
    delete DataFile;
    /*
    // Initialise NDisc variable (and add it to variables of interest)
    RooRealVar NDisc("NDisc", "NDisc", 0 , 1);
    variables_of_interest->add(NDisc);
    // Empty dataset
    RooDataSet * DataNDisc = new RooDataSet("DataNDisc", "DataNDisc", NDisc);
    // Filling dataset by hand - keeps index
    RooArgSet * observables_temp = (RooArgSet *)DataPreselected->get();
    RooRealVar * F90_temp = (RooRealVar *)observables_temp->find(S1F90.GetName());
    RooRealVar * logS2_o_S1_temp = (RooRealVar *)observables_temp->find(S2maxLoS1Lcorr.GetName());
    for (int i = 0; i < DataPreselected->numEntries(); i++) {
      DataPreselected->get(i);
      NDisc = calculate_n_disc(F90_temp->getVal(), logS2_o_S1_temp->getVal());
      DataNDisc->add(NDisc);
    }
    // Adding NDisc to dataset
    DataPreselected->merge(DataNDisc);
    // Memory
    delete DataNDisc;
    */

    // Create a new empty workspace
    RooWorkspace * w = new RooWorkspace("w","workspace");
    // Produce datasets for all S1 bins - single- and multiple-scatterers separated
    RooDataSet * S1bin_slice;
    // and saving them to the workspace
    for (auto const& S1bin : map_S1_bin_cuts) {
      S1bin_slice = (RooDataSet *)DataPreselected->reduce(SelectVars(*variables_of_interest),
              Cut((S1bin.second).c_str()));
      S1bin_slice->SetName(("S1bin_" + to_string(S1bin.first) + "_single").c_str());
      if (S1bin.first != S1L_bin_edges.size()) {
        S1bin_slice->SetTitle(("S1 in [" + to_string(S1L_bin_edges[S1bin.first]) + ", "
        + to_string(S1L_bin_edges[S1bin.first + 1]) + "] p.e.").c_str());
      } else {
        S1bin_slice->SetTitle(("S1 >= " + to_string(S1L_bin_edges[S1bin.first]) + " p.e.").c_str());
      }
      // Import data into the workspace
    }
    std::cout << S1bin_slice->sumEntries() << std::endl;
    // Print workspace contents
    w->Print();
    // Memory
    delete S1bin_slice;
    delete DataPreselected;
    // Save the workspace into a ROOT file - unfortunately necessary to cut away weird characters from the "ls" output
    if (is_on_essos) {w->writeToFile(OutputFileName_temp.c_str());}
    else {w->writeToFile(OutputFileName_temp.substr(0, OutputFileName_temp.length() - 1).c_str());}

    // Processing done
    processed_files++;
  }
  std::cout << "\nDONE" << std::endl;
}


// This function takes as input a terminal command and return the standard output
std::vector<std::string> GetStdoutFromCommand(std::string cmd) {
  std::vector<std::string> data;
  FILE * stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL) data.push_back(buffer);
    pclose(stream);
  }
  return data;
}


// N_discriminator
double calculate_n_disc(double F90_temp, double logS2LoS1L_temp) {
  ////// LS Fit results
  const double center_F90 = 0.737841;
  const double distance_unit_F90 = 0.0652304;
  const double center_logS2LoS1L = 0.223837;
  const double distance_unit_logS2LoS1L = 0.469427;

  ////// LogLikelihood Fit results
  // const double center_F90 = 0.737728;
  // const double distance_unit_F90 = 0.0663849;
  // const double center_logS2LoS1L = 0.223794;
  // const double distance_unit_logS2LoS1L = 0.48069;

  // Calculation of N_discriminator | x and y are i.i.d. random variables each following a standard normal distribution
  double distance_in_x = (F90_temp - center_F90) / distance_unit_F90;
  double distance_in_y = (logS2LoS1L_temp - center_logS2LoS1L) / distance_unit_logS2LoS1L;
  double distance_in_r = TMath::Sqrt(distance_in_x*distance_in_x + distance_in_y*distance_in_y);
  // double pdf_probability = distance_in_r * TMath::Exp( - distance_in_r*distance_in_r / 2);
  double cdf_probability = 1. - TMath::Exp( - distance_in_r*distance_in_r / 2.);

  return cdf_probability;
}
