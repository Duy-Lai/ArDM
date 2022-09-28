#include "param.h"

string InputPath;
string OutputPath;
string bump_name = "";

void input_output() {
  if (Source == "neutron") {
    if (bump == "both") {
      InputPath = "/mnt/raid/users/duylai/warios_sniper_" + recoil_type + cut + "/";
      OutputPath = "/home/duylai/ardm/figures_sniper_" + recoil_type + cut + "/";
      // InputPath = "/mnt/raid/users/duylai/warios_neutron_" + recoil_type + cut + "/";
      // OutputPath = "/home/duylai/ardm/figures_neutron_" + recoil_type + cut + "/";
    } else {
      if (bump == "upper") {
        InputPath = "/mnt/raid/users/duylai/warios_neutron_" + recoil_type + cut + "/1bump/";
      } else if (bump == "lower") {
        InputPath = "/mnt/raid/users/duylai/warios_neutron_" + recoil_type + cut + "/2bump/";
      }
      OutputPath = "/home/duylai/ardm/figures_neutron_" + recoil_type + cut + "/2bumps/";
    }
  } else if (Source == "none") {
    if (bump == "both") {
      InputPath = "/mnt/raid/users/duylai/warios_dp_" + recoil_type + cut + "/";
      OutputPath = "/home/duylai/ardm/figures_dp_" + recoil_type + cut + "/";
    } else {
      if (bump == "upper") {
        InputPath = "/mnt/raid/users/duylai/warios_dp_" + recoil_type + cut + "/1bump/";
      } else if (bump == "lower") {
        InputPath = "/mnt/raid/users/duylai/warios_dp_" + recoil_type + cut + "/2bump/";
      }
      OutputPath = "/home/duylai/ardm/figures_dp_" + recoil_type + cut + "/2bumps/";
    }
  } else {
    InputPath = "/home/duylai/ardm/warios_all/";
    OutputPath = "/home/duylai/ardm/figures_all/";
  }
  OutputPath += "runs/";
  cout << InputPath << endl;
}

// Input details
string DataType = "DP[12]";  // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
int MinRunNr = 0;  // 63 // Minimum run number of runs to be processed; Neutrons - 13 and partially 62 is weak source
int MaxRunNr = 124;  // Maximum run number of runs to be processed

string PlotFileType = "png";

// X and Y offset determines starting position of first plot (default x=0, y=0)
int XOffset = 10;  // Alex's Home Zürich: 10, Alex's Home Basel: 0
int YOffset = 35;  // Alex's Home Zürich: 35, Alex's Home Basel: 880

vector<string> GetStdoutFromCommand(string cmd);
double gauss_2D(double* x, double* par);

void old_calibrate_runs() {
  RooWorkspace* wks = new RooWorkspace("myWS");
  input_output();
  // Extract list of filenames
  string command = "ls " + InputPath + "WarioRun*-" + DataType + "-" + Source + ".root";
  vector<string> pre_file_list = GetStdoutFromCommand(command);
  // cout << pre_file_list[0] << endl;
  int start = InputPath.length() + 8; // getting the starting index of run number 00xx

  // Looping through the files and keeping the ones with appropriate run number
  // mainly to get the number of files to be processed for progress information...
  int files_to_be_processed = 0;
  string first_file;
  vector<string> consecutive_files;
  for (vector<string>::iterator t = pre_file_list.begin(); t != pre_file_list.end(); ++t) {
    const int run_nr_temp = stoi(t->substr(start, 4));
    if (run_nr_temp >= MinRunNr && run_nr_temp <= MaxRunNr) {
      consecutive_files.push_back(t->substr(0, t->length() - 1));
    }
  }
  
  // Plot positions
  int x_start = 80 + XOffset;
  int y_start = 0 + YOffset;
  float x_size = 1440;
  float y_size = 1200;
  float x_step = 0;
  float y_step = 0;

  ofstream data_file;
  // data_file.open(OutputPath + "fit_data.csv");
  // data_file << "name,amplitude,mean_F90,sigma_F90,mean_logS2maxLoS1L,sigma_logS2maxLoS1L,correlation_coeff,theta" << endl;
  data_file.open(OutputPath + "entries_bumps.csv");
  data_file << "run,upper,lower" << endl;

  int processed_files = 1;
  for (vector<string>::iterator t = consecutive_files.begin(); t != consecutive_files.end(); ++t) {
    const int run_nr_temp = stoi(t->substr(start, 4));
    const string ending_temp = t->substr(start);                             // extracting e. g. "0089-DP2-neutron.root"
    const string name_temp = t->substr(0, start + ending_temp.length()); // removing additional characters
    processed_files++;
    TFile* ToBeAdded = new TFile(name_temp.c_str());
    // auto RoI_logS2maxLoS1L_vs_F90_histo = (TH2I*)ToBeAdded->Get("RoI_logS2maxLoS1L_vs_F90");
    auto tree = (TTree*)ToBeAdded->Get("HebingTree");
    struct HEBING_STRUCT Hebing;
    tree->SetBranchAddress("S1F90", &Hebing.S1F90);
    tree->SetBranchAddress("S2maxLoS1Lcorr", &Hebing.S2maxLoS1Lcorr);
    // RoI: logS2maxLoS1L vs S1F90
    TH2F* RoI_logS2maxLoS1L_vs_F90_histo = new TH2F("RoI_logS2maxLoS1L_vs_F90", "log10S2maxLoS1L vs F90", 200, 0, 1., 200, -3., 5.);
    RoI_logS2maxLoS1L_vs_F90_histo->GetXaxis()->SetTitle("F90");
    RoI_logS2maxLoS1L_vs_F90_histo->GetYaxis()->SetTitle("log10(corrS2maxLoS1L)");
    int n_upper = 0;
    int n_lower = 0;
    for (auto i : ROOT::TSeqI(tree->GetEntries())) {
      tree->GetEntry(i);
      RoI_logS2maxLoS1L_vs_F90_histo->Fill(Hebing.S1F90, Hebing.S2maxLoS1Lcorr);
      if ((Hebing.S2maxLoS1Lcorr >= 0.5 && Source == "none") || (Hebing.S2maxLoS1Lcorr >= 0.2 && Source == "neutron"))
              ++n_upper;
      else if ((Hebing.S2maxLoS1Lcorr < 0.5 && Source == "none") || (Hebing.S2maxLoS1Lcorr < 0.2 && Source == "neutron"))
              ++n_lower;
    }
    data_file << run_nr_temp << "," << n_upper << "," << n_lower << endl;
    auto c11 = new TCanvas("c11", "c11", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    c11->SetLogz();
    c11->SetGridx();
    c11->SetGridy();
    RoI_logS2maxLoS1L_vs_F90_histo->DrawClone("colsz");
    c11->SaveAs((OutputPath + "RoI_logS2maxLoS1L_vs_F90_" + to_string(run_nr_temp) + "." + PlotFileType).c_str());
    delete c11;
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
      if (fgets(buffer, max_buffer, stream) != NULL) data.push_back(buffer);
    pclose(stream);
  }
  return data;
}

// 2D gaussian
double gauss_2D(double* x, double* par) {
  if (par[2] > 0 && par[4] > 0 && par[5] > 0) {
    double rxx = TMath::Power((x[0] - par[1]) / par[2], 2.0);
    double ryy = TMath::Power((x[1] - par[3]) / par[4], 2.0);
    double rxy = 2.0 * par[5] * (x[0] - par[1]) * (x[1] - par[3]) / par[2] / par[4];
    return TMath::InvPi() * par[0] / (2.0 * par[2] * par[4] * TMath::Sqrt(1 - par[5] * par[5])) * TMath::Exp(-(rxx + ryy - rxy) / 2.0 / (1 - par[5] * par[5]));
  }
  else {return 0.0;}
}