// Run the preselection of Satoshi runs into Wario files
//
// root -L MakeWario.C
//
// AR, June 2020
//
// adapted by Alex Stauffer
// USE ReadData instead, more uptodate and uses ROODatasets
// 
// Adapted by Duy Lai, 09.2022, requires param.h

#include "param.h"

/////////////////////////////////////////////   USER INPUT   ///////////////////////////////////////////////////////////

// Input details
string DataType = "*"; // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
// string DataType = "DP[12]"; // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
// int MinRunNr = 68;  // Minimum run number of runs to be "wariolised"
// int MaxRunNr = 75;  // Maximum run number of runs to be "wariolised"
int MinRunNr = 1;
int MaxRunNr = 169;
bool is_Hebing_output = false; // Otherwise it is assumed to be a Satoshi-Run

string InputPath, OutputPath;

void output() {
  if (Source == "neutron") {
    InputPath = "/mnt/raid/users/duylai/Satoshi_neutron_runs/";
    OutputPath = "/mnt/raid/users/duylai/matched_neutron/";
  } else if (Source == "none") {
    InputPath = "/mnt/raid/users/duylai/Satoshi_DP_runs/";
    OutputPath = "/mnt/raid/users/duylai/matched_DP/";
  }
  cout << "Producing warios at " << OutputPath << " ..." << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ProcessFile(const char* current_filename, const int runID, string final_path, string ending);
vector<string> GetStdoutFromCommand(string cmd);

// Main function which extracts filenames and calls processing
void MakeWario() {
  // Potentially make target folder
  output();
  int status = mkdir(OutputPath.c_str(), 0777);

  // Extract list of filenames
  string command = "ls " + InputPath + "SRun*-" + DataType + "-" + Source + ".root";
  if (is_Hebing_output) {
    command = "ls " + InputPath + "Run*.root";
  }
  vector<string> pre_file_list = GetStdoutFromCommand(command);
  int pre_string_length = 4;
  int run_ID_string_length = 4;
  if (is_Hebing_output) {
    pre_string_length = 3;
  }
  if (is_Hebing_output) {
    run_ID_string_length = 6;
  }
  int start = InputPath.length() + pre_string_length; // getting the starting index of run number 00xx

  // Looping through the files and keeping the ones with appropriate run number
  // mainly to get the number of files to be processed for progress information...
  int files_to_be_processed = 0;
  vector<string> file_list;
  for (vector<string>::iterator t = pre_file_list.begin(); t != pre_file_list.end(); ++t) {
    const int run_nr_temp = stoi(t->substr(start, run_ID_string_length));
    if (run_nr_temp >= MinRunNr && run_nr_temp <= MaxRunNr) {
      file_list.push_back(t->c_str());
      files_to_be_processed++;
    }
  }

  // Actually processing the files
  int processed_files = 0;
  for (vector<string>::iterator t = file_list.begin(); t != file_list.end(); ++t) {
    cout << "Bump ....." << bump << endl;
    const int run_nr_temp = stoi(t->substr(start, run_ID_string_length));
    const string ending_temp = t->substr(start);                             // extracting e. g. "0089-DP2-neutron.root"
    const string name_temp = t->substr(0, start + ending_temp.length() - 1); // removing additional characters
    cout << "=======================================================================" << endl;
    cout << "Processing file:" << endl;
    cout << name_temp << endl;
    cout << "File: " << processed_files + 1 << " out of " << files_to_be_processed << endl;
    cout << "=======================================================================" << endl;
    ProcessFile(name_temp.c_str(), run_nr_temp, OutputPath, ending_temp);
    processed_files++;
  }
}

// This function takes as input a terminal command and return the standard output
vector<string> GetStdoutFromCommand(string cmd) {
  vector<string> data;
  FILE* stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL)
        data.push_back(buffer);
    pclose(stream);
  }
  return data;
}

struct HEBING_STRUCT Hebing;

// Linearisation of TTR to Drift time for low drift fields
double linearised_DT(double x) {
  // QMatch - low drift field ~230 V/cm
  double QMatch_par_a = 0.5163192586635797;
  double QMatch_par_b = -0.04259468609403853; // 1/us
  double QMatch_par_c = -1.1997290781301084;  // 1/us**2
  double QMatch_par_d = 1.4317140748085602;   // 1/us**3
  double QMatch_par_e = -0.5083798612364056;  // 1/us**4

  return QMatch_par_a + QMatch_par_b * x + QMatch_par_c * x * x + QMatch_par_d * x * x * x + QMatch_par_e * x * x * x * x;
}

// Linearisation of TTR to Drift time for low drift fields
double linearised_DT_high(double x) {
  // QMatch - high drift field ~340 V/cm
  double QMatch_par_a = 0.5158214282713242;
  double QMatch_par_b = -0.0974761816076573; // 1/us
  double QMatch_par_c = -1.6179453473493883; // 1/us**2
  double QMatch_par_d = 2.3558627870750697;  // 1/us**3
  double QMatch_par_e = -1.0147690094106852; // 1/us**4

  return QMatch_par_a + QMatch_par_b * x + QMatch_par_c * x * x + QMatch_par_d * x * x * x + QMatch_par_e * x * x * x * x;
}

// Funtion of the ER region
double ER_F90_lower_limit(double x) {
  int stds = 5;
  double ER_par_a = 0.2882788009259751;
  double ER_par_b = 0.3028346934757867;
  double ER_par_c = -0.07320972612366412;
  double ER_par_d = 0.03054922411320943;
  double ER_par_e = 1.687045528047335;
  double ER_par_f = 3.805435737601408;

  return ER_par_a + ER_par_b * exp(ER_par_c * x) + stds * (ER_par_d + ER_par_e / (x + ER_par_f));
}

// Function of the NR region
double NR_F90_upper_limit(double x) {
  int stds = 3;
  double NR_par_a = 0.7367600167931209;
  double NR_par_b = -0.10129012255422735;
  double NR_par_c = -0.026119459318325253;
  double NR_par_d = 0.23776101022499113;
  double NR_par_e = 345.8223734905984;
  double NR_par_f = -2139.785720352299;

  return NR_par_a - NR_par_b * exp(NR_par_c * x) + stds * (NR_par_d + NR_par_e / (x + NR_par_f));
}

// Function of the NR region
double NR_F90_lower_limit(double x) {
  int stds = -3;
  double NR_par_a = 0.7367600167931209;
  double NR_par_b = -0.10129012255422735;
  double NR_par_c = -0.026119459318325253;
  double NR_par_d = 0.23776101022499113;
  double NR_par_e = 345.8223734905984;
  double NR_par_f = -2139.785720352299;

  return NR_par_a - NR_par_b * exp(NR_par_c * x) + stds * (NR_par_d + NR_par_e / (x + NR_par_f));
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
  double distance_in_r = TMath::Sqrt(distance_in_x * distance_in_x + distance_in_y * distance_in_y);
  // double pdf_probability = distance_in_r * TMath::Exp( - distance_in_r*distance_in_r / 2);
  double cdf_probability = 1. - TMath::Exp(-distance_in_r * distance_in_r / 2.);
  
  return cdf_probability;
}

// Actual processing of a given file(name), output in WarioRun file with same naming scheme
void ProcessFile(const char* current_filename, const int runID, string final_path, string ending) {
  // Open the right TTree in the given file(name)
  TFile oldfile(current_filename);
  TTree* oldtree;
  oldfile.GetObject("HebingTree", oldtree);

  // Get the total number of entries
  const auto nentries = oldtree->GetEntries();
  cout << "Total number of entries " << nentries << endl;

  //  Event *event = nullptr;
  //  oldtree->SetBranchAddress("event", &event);

  // Activate all branches
  oldtree->SetBranchStatus("*", 1);
  // Activate only needed branches
  //  for (auto activeBranchName : {"S1unixTime","nS2"})
  //    oldtree->SetBranchStatus(activeBranchName, 1);
  oldtree->SetBranchAddress("S1unixTime", &Hebing.S1unixTime);
  oldtree->SetBranchAddress("S1TTT", &Hebing.S1TTT);
  oldtree->SetBranchAddress("S1LY", &Hebing.S1LY);
  oldtree->SetBranchAddress("S1F90", &Hebing.S1F90);
  oldtree->SetBranchAddress("S1TTR", &Hebing.S1TTR);
  //  oldtree->SetBranchAddress("S1F90Top",&Hebing.S1F90Top,"S1F90Top/F");
  //  oldtree->SetBranchAddress("S1F90Bottom",&Hebing.S1F90Bottom,"S1F90Bottom/F");
  oldtree->SetBranchAddress("S1COGX", &Hebing.S1COGX);
  oldtree->SetBranchAddress("S1COGY", &Hebing.S1COGY);
  oldtree->SetBranchAddress("S1maxFracPMT", &Hebing.S1maxFracPMT);

  oldtree->SetBranchAddress("nS2", &Hebing.matchedS2);
  oldtree->SetBranchAddress("S2LY", &Hebing.S2LY);
  oldtree->SetBranchAddress("S2maxLY", &Hebing.S2maxLY);
  oldtree->SetBranchAddress("S2maxDT", &Hebing.S2maxDT);
  oldtree->SetBranchAddress("S2maxLSX", &Hebing.S2maxrecoLSX);
  oldtree->SetBranchAddress("S2maxLSY", &Hebing.S2maxrecoLSY);
  oldtree->SetBranchAddress("S2maxLSCHI2", &Hebing.S2maxrecoLSCHI2);
  oldtree->SetBranchAddress("S2maxmaxFracPMT", &Hebing.S2maxmaxFracPMT);
  oldtree->SetBranchAddress("S2maxLoS1L", &Hebing.S2maxLoS1L);
  oldtree->SetBranchAddress("S2maxLoS1Lcorr", &Hebing.S2maxLoS1Lcorr);
  oldtree->SetBranchAddress("S2maxQMatch", &Hebing.S2maxQMatch);
  oldtree->SetBranchAddress("S2maxFracS2L", &Hebing.S2maxFracS2L);
  oldtree->SetBranchAddress("S1Nx", &Hebing.S1Nx);
  oldtree->SetBranchAddress("nS2", &Hebing.nS2);
  oldtree->SetBranchAddress("nS1perS2", &Hebing.nS1perS2);
  oldtree->SetBranchAddress("nS1perS2max", &Hebing.nS1perS2max);

  // Create a new file and a clone of the structure of the old TTree for the output file
  string output_file = final_path + "WarioRun" + ending;
  TFile newfile(output_file.substr(0, output_file.length() - 1).c_str(), "recreate");
  auto newtree = oldtree->CloneTree(0);

  // Histos
  TH1F* myallS1LHisto = new TH1F("aS1L", "all S1L", 2000, 0., 2000.);
  myallS1LHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  myallS1LHisto->GetYaxis()->SetTitle("Events");

  TH2F* mymatchedlogS2maxLoS1vsF90Histo = new TH2F("matched_logS2maxLoS1L_vs_F90", "log10S2maxLoS1L vs F90",
          200, 0., 1., 200, -3., 5.);
  mymatchedlogS2maxLoS1vsF90Histo->GetXaxis()->SetTitle("S1 F90");
  mymatchedlogS2maxLoS1vsF90Histo->GetYaxis()->SetTitle("log10(S2maxLoS1L)");

  TH2F* myallF90vsS1LHisto = new TH2F("all_F90_vs_S1L", "All F90 vs S1L", 200, 0., 800.0, 200, 0., 1.);
  myallF90vsS1LHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  myallF90vsS1LHisto->GetYaxis()->SetTitle("S1 F90");
  
  TH1F* mypreS1LHisto = new TH1F("pre_S1L","Preselected S1 L", 800, 0., 800.);
  mypreS1LHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  mypreS1LHisto->GetYaxis()->SetTitle("Events");

  TH2F* mypreF90vsS1LHisto = new TH2F("pre_F90_vs_S1L", "Preselected F90 vs S1L", 200, 0., 100.0, 200, 0., 1.);
  mypreF90vsS1LHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  mypreF90vsS1LHisto->GetYaxis()->SetTitle("S1 F90");

  TH2F* myprelogS2maxLoS1LvsF90LHisto = new TH2F("pre_logS2maxLoS1L_vs_F90", "log10S2maxLoS1L vs F90", 200, 0., 1., 200, -3., 5.);
  myprelogS2maxLoS1LvsF90LHisto->GetXaxis()->SetTitle("S1 F90");
  myprelogS2maxLoS1LvsF90LHisto->GetYaxis()->SetTitle("log10(corrS2maxLoS1L)");

  TH2F* myRoI_F90vsS1LHisto = new TH2F("RoI_F90_vs_S1L", "F90 vs S1L", 200, 0., 100.0, 200, 0., 1.);
  myRoI_F90vsS1LHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  myRoI_F90vsS1LHisto->GetYaxis()->SetTitle("S1 F90");

  TH2F* myRoI_logS2maxLoS1LvsF90LHisto = new TH2F("RoI_logS2maxLoS1L_vs_F90", "log10S2maxLoS1L vs F90", 200, 0, 1., 200, -3., 5.);
  myRoI_logS2maxLoS1LvsF90LHisto->GetXaxis()->SetTitle("F90");
  myRoI_logS2maxLoS1LvsF90LHisto->GetYaxis()->SetTitle("log10(corrS2maxLoS1L)");

  TH1F* myQMatchHisto = new TH1F("QMatch", "QMatch", 100, 0., 4.0);
  myQMatchHisto->GetXaxis()->SetTitle("QMatch");
  myQMatchHisto->GetYaxis()->SetTitle("Events");

  TH1F* myS1maxFracPMTHisto = new TH1F("S1maxFracPMT", "S1maxFracPMT", 100, 0., 1.);
  myS1maxFracPMTHisto->GetXaxis()->SetTitle("S1maxFracPMT");
  myS1maxFracPMTHisto->GetYaxis()->SetTitle("Events");

  TH1F* myS2maxFracS2LHisto = new TH1F("S2maxFracS2L", "S2maxFracS2L", 100, 0.90, 1.);
  myS2maxFracS2LHisto->GetXaxis()->SetTitle("S2maxFracS2L");
  myS2maxFracS2LHisto->GetYaxis()->SetTitle("Events");

  TH1F* mycQMatchHisto = new TH1F("cQMatch", "Current QMatch", 100, 0., 4.0);
  mycQMatchHisto->GetXaxis()->SetTitle("QMatch");
  mycQMatchHisto->GetYaxis()->SetTitle("Events");

  TH2F* myQMatchvsDTHisto = new TH2F("QMatch_vs_DT", "QMatch vs Drift time", 200, 0., 1.4, 200, 0., 4.);
  myQMatchvsDTHisto->GetXaxis()->SetTitle("S2max DT [ms]");
  myQMatchvsDTHisto->GetYaxis()->SetTitle("QMatch");

  TH2F* mymatchedS2maxLSYvsS2maxLSXHisto = new TH2F("S2max_LS_position", "S2max LS position", 100, -500., 500., 100, -500., 500.);
  mymatchedS2maxLSYvsS2maxLSXHisto->GetXaxis()->SetTitle("S2maxLSX [mm]");
  mymatchedS2maxLSYvsS2maxLSXHisto->GetYaxis()->SetTitle("S2maxLSY [mm]");

  TH2F* myS2maxDTvsS2maxLSRadiusHisto = new TH2F("Radius_vs_DT", "Radius vs Drift time", 100, 0., 500., 100, 0., 1.4);
  myS2maxDTvsS2maxLSRadiusHisto->GetXaxis()->SetTitle("S2maxLSRadius [mm]");
  myS2maxDTvsS2maxLSRadiusHisto->GetYaxis()->SetTitle("S2maxDT [mm]");

  TH2F* myDTvsS1TTRHisto = new TH2F("DT_vs_TTR", "Drift time vs TTR", 200, 0, 1.4, 200, 0., 1.);
  myDTvsS1TTRHisto->GetXaxis()->SetTitle("S2max DT [ms]");
  myDTvsS1TTRHisto->GetYaxis()->SetTitle("S1TTR");

  TH1F* myS2maxDTHisto = new TH1F("S2maxDT", "Drift time", 100, 0., 1.4);
  myS2maxDTHisto->GetXaxis()->SetTitle("S2max DT [ms]");
  myS2maxDTHisto->GetYaxis()->SetTitle("Events");

  TH1F* myS2maxLYHisto = new TH1F("S2maxLY", "S2maxLY", 200, 0., 2000.);
  myS2maxLYHisto->GetXaxis()->SetTitle("S2maxLY [p.e]");
  myS2maxLYHisto->GetYaxis()->SetTitle("Events");

  TH1F* myS2LYHisto = new TH1F("S2LY", "S2LY", 200, 0., 2000.);
  myS2LYHisto->GetXaxis()->SetTitle("S2LY [p.e]");
  myS2LYHisto->GetYaxis()->SetTitle("Events");

  TH1I* myS1NxHisto = new TH1I("S1Nx", "S1Nx", 11, -0.5, 10.5);
  myS1NxHisto->GetXaxis()->SetTitle("S1Nx");
  myS1NxHisto->GetYaxis()->SetTitle("Events");

  TH1I* mynS2Histo = new TH1I("nS2", "nS2", 11, -0.5, 10.5);
  mynS2Histo->GetXaxis()->SetTitle("nS2");
  mynS2Histo->GetYaxis()->SetTitle("Events");

  TH1F* myS1LYHisto = new TH1F("S1LY", "S1LY", 100, 0., 100.);
  myS1LYHisto->GetXaxis()->SetTitle("S1 Total light [p.e.]");
  myS1LYHisto->GetYaxis()->SetTitle("Events");

  // Average S2maxFracS2L vs position plots
  int nr_bins = 100;
  int nr_global_bins = nr_bins*nr_bins;
  TH2F* myAverageS2maxFracS2LvsPositionHisto = new TH2F("Average_S2maxFracS2L_vs_position",
            "Average S2maxFracS2L vs position", nr_bins, 0., 500., nr_bins, 0., 1.4);
  myAverageS2maxFracS2LvsPositionHisto->GetXaxis()->SetTitle("S2maxLSRadius [mm]");
  myAverageS2maxFracS2LvsPositionHisto->GetYaxis()->SetTitle("S2maxDT [mm]");

  TH1F* mySumS2maxFracS2L = new TH1F("S2maxFracS2L_sums", "S2maxFracS2L sums", nr_global_bins, 0., nr_global_bins);
  mySumS2maxFracS2L->GetXaxis()->SetTitle("Global Bin number");
  mySumS2maxFracS2L->GetYaxis()->SetTitle("Sum");

  TH1F* myEntriesS2maxFracS2L = new TH1F("S2maxFracS2L_entries", "S2maxFracS2L entries", nr_global_bins, 0., nr_global_bins);
  myEntriesS2maxFracS2L->GetXaxis()->SetTitle("Global Bin number");
  myEntriesS2maxFracS2L->GetYaxis()->SetTitle("Entries");

  // Initialising cut statistics (counts)
  int n_mS1events = 0;
  int n_preselected = 0;
  int n_S2maxevents = 0;
  // int n_pmts2distrib = 0;
  // int n_qnmatchls2chi2 = 0;
  int n_S1PMT = 0;
  int n_overall_S1PMT = 0;
  int n_S2FracL = 0;
  int n_overall_S2FracL = 0;
  int n_QMatch = 0;
  int n_overall_QMatch = 0;
  int n_RoI_F90 = 0;
  int n_overall_RoI_F90 = 0;
  int n_n_disc = 0;
  int n_overall_n_disc = 0;
  int n_fiducial = 0;
  int n_overall_fiducial = 0;
  int nWarioCuts = 0;
  int count_hey = 0;
  int current_position_bin = 0;
  for (auto i : ROOT::TSeqI(nentries)) {
    if (i % 1000000 == 0)
      cout << 100.0 * (i) / (nentries) << "\%" << endl;
    oldtree->GetEntry(i);
    myallS1LHisto->Fill(Hebing.S1LY);
    // myallS1F90Histo->Fill(Hebing.S1F90);
    // myallS1TTRHisto->Fill(Hebing.S1TTR);
    // myallS1maxFracPMTHisto->Fill(Hebing.S1maxFracPMT);
    // if (Hebing.S1maxFracPMT <= 0.3) {n_overall_S1PMT++;}
    // mymatchednS2sHisto->Fill(Hebing.matchedS2);
    // if (Hebing.S1TTR < 0.15) {mymatchednS2slowTTRHisto->Fill(Hebing.matchedS2);}
    // if (Hebing.S1TTR >= 0.15) {mymatchednS2shighTTRHisto->Fill(Hebing.matchedS2);}

    // Calculating the RoI F90 bounds
    myallF90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
    double current_min_F90 = ER_F90_lower_limit(Hebing.S1LY);
    double current_max_F90 = NR_F90_upper_limit(Hebing.S1LY);
    if (Hebing.S1F90 >= current_min_F90 && Hebing.S1F90 <= current_max_F90) {n_overall_RoI_F90++;}

    // Matched S2
    if (Hebing.S1F90 >= 0.0 && Hebing.S1F90 <= 1.0 && Hebing.S1maxFracPMT >= 0.0 && Hebing.S1maxFracPMT <= 1.0
              && Hebing.S2maxFracS2L >= 0.0 && Hebing.S2maxFracS2L <= 1.0 && Hebing.S2maxDT >= 0.0
              && Hebing.S2maxrecoLSX >= -600.0 && Hebing.S2maxrecoLSY >= -600.0&& Hebing.S2maxmaxFracPMT >= 0.0
              && Hebing.S2maxLoS1Lcorr >= -6.0) {
      // if (Hebing.matchedS2 > 0){
      n_mS1events++;
      n_S2maxevents++;
      // double current_LSRadius = TMath::Power(Hebing.S2maxrecoLSX*Hebing.S2maxrecoLSX + Hebing.S2maxrecoLSY*Hebing.S2maxrecoLSY, 0.5);
      // Calculating the linearised QMatch variable
      // double current_QMatch;
      // High drift field runs have IDs 100-102
      // if (runID == 100 | runID == 101 | runID == 102) {
        // current_QMatch = Hebing.S1TTR / linearised_DT_high(Hebing.S2maxDT);
      // } else {
        // current_QMatch = Hebing.S1TTR / linearised_DT(Hebing.S2maxDT);
      // }
      // if (current_QMatch >= 0.75 && current_QMatch <= 1.3) {n_overall_QMatch++;}
      // mymatchedlogS2maxLoS1vsF90Histo->Fill(Hebing.S1F90, Hebing.S2maxLoS1L);
      // mycQMatchHisto->Fill(current_QMatch);

      // Preselection
      double current_LSRadius = TMath::Power(Hebing.S2maxrecoLSX * Hebing.S2maxrecoLSX
                + Hebing.S2maxrecoLSY * Hebing.S2maxrecoLSY, 0.5);
      if (Hebing.S1LY <= 100.0 && Hebing.S1maxFracPMT <= 0.3 && Hebing.S1TTR >= 0.15
      ) {
                // && current_QMatch >= 0.75 && current_QMatch <= 1.3
                // && current_LSRadius <= 300.0 && Hebing.S2maxDT >= 0.2 && Hebing.S2maxDT <= 1.0) {
        n_preselected++;
        mypreS1LHisto->Fill(Hebing.S1LY);
        mypreF90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
        myprelogS2maxLoS1LvsF90LHisto->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
        // RoI - NR and ER regions - F90
        // if ((recoil_type == "NR" && Hebing.S1LY <= 100.0 && Hebing.S1F90 >= current_min_F90 && Hebing.S1F90 <= current_max_F90) ||
                // (recoil_type == "ER" && Hebing.S1LY <= 100.0 && Hebing.S1F90 < current_min_F90)) {
        // if (recoil_type == "NR" && Hebing.S1LY <= 100.0 && Hebing.S1F90 >= current_min_F90 && Hebing.S1F90 <= current_max_F90) {
        // if (Hebing.S1LY <= 100.0 && Hebing.S1F90 < current_min_F90) {//} && Hebing.S1F90 <= current_max_F90) {
          n_RoI_F90++;
          if (Source == "neutron") {
            if ((recoil_type == "NR" && bump == "upper" && Hebing.S2maxLoS1Lcorr >= 0.2) || 
                  (recoil_type == "ER" && bump == "upper" && Hebing.S2maxLoS1Lcorr >= 0.5) || bump == "both") {
              if (Hebing.S2maxFracS2L == 1.0 && Hebing.S2maxQMatch >= 0.75 && Hebing.S2maxQMatch <= 1.25
                    && Hebing.nS1perS2max == 1) {
                myRoI_F90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
                myRoI_logS2maxLoS1LvsF90LHisto->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
                myS1maxFracPMTHisto->Fill(Hebing.S1maxFracPMT);
                myS2maxFracS2LHisto->Fill(Hebing.S2maxFracS2L);
                myQMatchvsDTHisto->Fill(Hebing.S2maxDT, Hebing.S2maxQMatch);
                mymatchedS2maxLSYvsS2maxLSXHisto->Fill(Hebing.S2maxrecoLSX, Hebing.S2maxrecoLSY);
                myS2maxDTvsS2maxLSRadiusHisto->Fill(current_LSRadius, Hebing.S2maxDT);
                myDTvsS1TTRHisto->Fill(Hebing.S2maxDT, Hebing.S1TTR);
                // Average histograms - adding values to the sum-histo and incrementing corresponding bin in the entries-histo
                current_position_bin = myAverageS2maxFracS2LvsPositionHisto->FindBin(current_LSRadius, Hebing.S2maxDT);
                mySumS2maxFracS2L->Fill(current_position_bin, Hebing.S2maxFracS2L);
                myEntriesS2maxFracS2L->Fill(current_position_bin);
                myS2maxDTHisto->Fill(Hebing.S2maxDT);
                myQMatchHisto->Fill(Hebing.S2maxQMatch);
                myS1NxHisto->Fill(Hebing.S1Nx);
                mynS2Histo->Fill(Hebing.nS2);
                myS2maxLYHisto->Fill(Hebing.S2maxLY);
                myS2LYHisto->Fill(Hebing.S2LY);
                myS1LYHisto->Fill(Hebing.S1LY);
                newtree->Fill();
              // if (Hebing.S2maxQMatch >= 0.8 && Hebing.S2maxQMatch <= 1.2
                  // && ((Hebing.nS2 == 1 && Hebing.S1Nx <= 3) || (Hebing.S1Nx == 0 && Hebing.nS2 <= 3))
              // ) {
                  // newtree->Fill();
              }
            }
            if ((recoil_type == "NR" && bump == "lower" && Hebing.S2maxLoS1Lcorr < 0.2) || 
                  (recoil_type == "ER" && bump == "lower" && Hebing.S2maxLoS1Lcorr < 0.5)) {
              if (Hebing.S2maxFracS2L == 1.0 &&  Hebing.S2maxQMatch >= 0.75 && Hebing.S2maxQMatch <= 1.25
                    && Hebing.nS1perS2max == 1) {
                myRoI_F90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
                myRoI_logS2maxLoS1LvsF90LHisto->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
                myS1maxFracPMTHisto->Fill(Hebing.S1maxFracPMT);
                myS2maxFracS2LHisto->Fill(Hebing.S2maxFracS2L);
                myQMatchvsDTHisto->Fill(Hebing.S2maxDT, Hebing.S2maxQMatch);
                mymatchedS2maxLSYvsS2maxLSXHisto->Fill(Hebing.S2maxrecoLSX, Hebing.S2maxrecoLSY);
                myS2maxDTvsS2maxLSRadiusHisto->Fill(current_LSRadius, Hebing.S2maxDT);
                myDTvsS1TTRHisto->Fill(Hebing.S2maxDT, Hebing.S1TTR);
                // Average histograms - adding values to the sum-histo and incrementing corresponding bin in the entries-histo
                current_position_bin = myS2maxDTvsS2maxLSRadiusHisto->FindBin(current_LSRadius, Hebing.S2maxDT);
                mySumS2maxFracS2L->Fill(current_position_bin, Hebing.S2maxFracS2L);
                myEntriesS2maxFracS2L->Fill(current_position_bin);
                myS2maxDTHisto->Fill(Hebing.S2maxDT);
                myQMatchHisto->Fill(Hebing.S2maxQMatch);
                myS1NxHisto->Fill(Hebing.S1Nx);
                mynS2Histo->Fill(Hebing.nS2);
                myS2maxLYHisto->Fill(Hebing.S2maxLY);
                myS2LYHisto->Fill(Hebing.S2LY);
                myS1LYHisto->Fill(Hebing.S1LY);
                newtree->Fill();
              // if (Hebing.S2maxQMatch >= 0.8 && Hebing.S2maxQMatch <= 1.2
                  // && ((Hebing.nS2 == 1 && Hebing.S1Nx <= 3) || (Hebing.S1Nx == 0 && Hebing.nS2 <= 3))
              // ) {
                  // newtree->Fill();
              }
            }
          } else {
            if ((recoil_type == "NR" && bump == "upper" && Hebing.S2maxLoS1Lcorr >= 0.2) || 
                  (recoil_type == "ER" && bump == "upper" && Hebing.S2maxLoS1Lcorr >= 0.5) || bump == "both") {
              if (Hebing.S2maxFracS2L == 1.0 &&  Hebing.S2maxQMatch >= 0.75 && Hebing.S2maxQMatch <= 1.25
                    && Hebing.nS1perS2max == 1) {
                myRoI_F90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
                myRoI_logS2maxLoS1LvsF90LHisto->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
                myS1maxFracPMTHisto->Fill(Hebing.S1maxFracPMT);
                myS2maxFracS2LHisto->Fill(Hebing.S2maxFracS2L);
                myQMatchvsDTHisto->Fill(Hebing.S2maxDT, Hebing.S2maxQMatch);
                mymatchedS2maxLSYvsS2maxLSXHisto->Fill(Hebing.S2maxrecoLSX, Hebing.S2maxrecoLSY);
                myS2maxDTvsS2maxLSRadiusHisto->Fill(current_LSRadius, Hebing.S2maxDT);
                myDTvsS1TTRHisto->Fill(Hebing.S2maxDT, Hebing.S1TTR);
                // Average histograms - adding values to the sum-histo and incrementing corresponding bin in the entries-histo
                current_position_bin = myS2maxDTvsS2maxLSRadiusHisto->FindBin(current_LSRadius, Hebing.S2maxDT);
                mySumS2maxFracS2L->Fill(current_position_bin, Hebing.S2maxFracS2L);
                myEntriesS2maxFracS2L->Fill(current_position_bin);
                myS2maxDTHisto->Fill(Hebing.S2maxDT);
                myQMatchHisto->Fill(Hebing.S2maxQMatch);
                myS1NxHisto->Fill(Hebing.S1Nx);
                mynS2Histo->Fill(Hebing.nS2);
                myS2maxLYHisto->Fill(Hebing.S2maxLY);
                myS2LYHisto->Fill(Hebing.S2LY);
                myS1LYHisto->Fill(Hebing.S1LY);
                newtree->Fill();
              }
            }
            if ((recoil_type == "NR" && bump == "lower" && Hebing.S2maxLoS1Lcorr < 0.2) || 
                  (recoil_type == "ER" && bump == "lower" && Hebing.S2maxLoS1Lcorr < 0.5)) {
              if (Hebing.S2maxFracS2L == 1.0 &&  Hebing.S2maxQMatch >= 0.75 && Hebing.S2maxQMatch <= 1.25
                    && Hebing.nS1perS2max == 1) {
                myRoI_F90vsS1LHisto->Fill(Hebing.S1LY, Hebing.S1F90);
                myRoI_logS2maxLoS1LvsF90LHisto->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
                myS1maxFracPMTHisto->Fill(Hebing.S1maxFracPMT);
                myS2maxFracS2LHisto->Fill(Hebing.S2maxFracS2L);
                myQMatchvsDTHisto->Fill(Hebing.S2maxDT, Hebing.S2maxQMatch);
                mymatchedS2maxLSYvsS2maxLSXHisto->Fill(Hebing.S2maxrecoLSX, Hebing.S2maxrecoLSY);
                myS2maxDTvsS2maxLSRadiusHisto->Fill(current_LSRadius, Hebing.S2maxDT);
                myDTvsS1TTRHisto->Fill(Hebing.S2maxDT, Hebing.S1TTR);
                // Average histograms - adding values to the sum-histo and incrementing corresponding bin in the entries-histo
                current_position_bin = myS2maxDTvsS2maxLSRadiusHisto->FindBin(current_LSRadius, Hebing.S2maxDT);
                mySumS2maxFracS2L->Fill(current_position_bin, Hebing.S2maxFracS2L);
                myEntriesS2maxFracS2L->Fill(current_position_bin);
                myS2maxDTHisto->Fill(Hebing.S2maxDT);
                myQMatchHisto->Fill(Hebing.S2maxQMatch);
                myS1NxHisto->Fill(Hebing.S1Nx);
                mynS2Histo->Fill(Hebing.nS2);
                myS2maxLYHisto->Fill(Hebing.S2maxLY);
                myS2LYHisto->Fill(Hebing.S2LY);
                myS1LYHisto->Fill(Hebing.S1LY);
                newtree->Fill();
              }
            }
          }
        // }
      }
    }
  }
  
  newfile.Write();

  // Absolute cut statistics
  cout << "======================================" << endl;
  cout << "Absolute cut statistics" << endl;
  cout << "======================================" << endl;
  cout << "Total number of S1 events " << nentries << endl;
  cout << "Total number of S1PMTFrac events " << n_overall_S1PMT << " (" << 100.0 * (n_overall_S1PMT) / (nentries) << " \%)" << endl;
  cout << "Total number of RoI - I (F90) events " << n_overall_RoI_F90 << " (" << 100.0 * (n_overall_RoI_F90) / (nentries) << " \%)" << endl;
  cout << "Total number of S2max events " << n_S2maxevents << endl;
  cout << "Total number of S2maxFracL events " << n_overall_S2FracL << " (" << 100.0 * (n_overall_S2FracL) / (n_S2maxevents) << " \%)" << endl;
  cout << "Total number of S2maxQmatch events " << n_overall_QMatch << " (" << 100.0 * (n_overall_QMatch) / (n_S2maxevents) << " \%)" << endl;
  cout << "Total number of RoI - II (n_disc) events " << n_overall_n_disc << " (" << 100.0 * (n_overall_n_disc) / (n_S2maxevents) << " \%)" << endl;
  cout << "Total number of fiducial volume events " << n_overall_fiducial << " (" << 100.0 * (n_overall_fiducial) / (n_S2maxevents) << " \%)" << endl;

  // Summary statistics
  cout << "======================================" << endl;
  cout << "Consecutive cut statistics" << endl;
  cout << "======================================" << endl;
  cout << "Total number of S1 events " << nentries << endl;
  cout << "Total number of matched S1 events " << n_mS1events << " (" << 100.0 * (n_mS1events) / (nentries)
            << " \%)" << endl;
  cout << "Total number of preselected events " << n_preselected << " (" << 100.0 * (n_preselected) / (nentries)
            << " \%)" << endl;
  cout << "Percentage reset to preselected events..." << endl;
  cout << "Total number of preselected events " << n_preselected << " (" << 100.0
            << " \%)" << endl;
  cout << "Total number of S1PMTFrac events " << n_S1PMT << " (" << 100.0 * (n_S1PMT) / (n_preselected) << " \%)"
            << endl;
  cout << "Total number of S2FracL events " << n_S2FracL << " (" << 100.0 * (n_S2FracL) / (n_preselected) << " \%)"
            << endl;
  cout << "Total number of QMatch events " << n_QMatch << " (" << 100.0 * (n_QMatch) / (n_preselected) << " \%)"
            << endl;
  cout << "Total number of RoI - I (F90) events " << n_RoI_F90 << " (" << 100.0 * (n_RoI_F90) / (n_preselected)
            << " \%)" << endl;
  cout << "Total number of RoI - II (n_disc) events " << n_n_disc << " (" << 100.0 * (n_n_disc) / (n_preselected)
            << " \%)" << endl;
  cout << "Total number of fiducial volume events " << n_fiducial << " (" << 100.0 * (n_fiducial) / (n_preselected)
            << " \%)" << endl;
}