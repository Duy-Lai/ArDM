// Add different Wario files
// Which files to be process has to be declared under USER INPUT
// 
// root -L CombineWarios.C
// 
// Created by alex on 22.10.20.
// 
// Adapted by Duy Lai, 09.2022, requires param.h
// Added a combined HebingTree structure into the SuperWario.root

#include "param.h"

/////////////////////////////////////////////   USER INPUT   ///////////////////////////////////////////////////////////

string InputPath;
string OutputPath;

void input_output() {
  if (Source == "neutron") {
    if (bump == "both") {
      InputPath = "/mnt/raid/users/duylai/matched_neutron/";
      OutputPath = "/home/duylai/ardm/figures_neutron/";
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
      InputPath = "/mnt/raid/users/duylai/warios_DP2/";
      OutputPath = "/home/duylai/ardm/figures_dp/";
    } else {
      if (bump == "upper") {
        InputPath = "/mnt/raid/users/duylai/warios_dp_" + recoil_type + cut + "/1bump/";
      } else if (bump == "lower") {
        InputPath = "/mnt/raid/users/duylai/warios_dp_" + recoil_type + cut + "/2bump/";
      }
      OutputPath = "/home/duylai/ardm/figures_dp_" + recoil_type + cut + "/2bumps/";
    }
  }
  cout << InputPath << endl;
}

// Input details
string DataType = "*";  // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
// string DataType = "DP[12]";  // "DP1" or "DP2" or "unk" or "DP[12]" or "*"
int MinRunNr = 0;  // 63 // Minimum run number of runs to be processed; Neutrons - 13 and partially 62 is weak source
int MaxRunNr = 124;  // Maximum run number of runs to be processed

// Options
bool MergeTrees = true;  // Create new tree in newfile or not, may take a lot of space
bool ShowPlots = false;  // whether or not to show plots (in canvases)
bool SavePlots = false;  // to specified file type in OutputPath; only works if "ShowPlots = true"
bool SavePlotsThesis = false;  // to specified file type in OutputPath; only works if "ShowPlots = true"
string PlotFileType = "png";  // for non-"pdf"s, so "jpeg" or "png" this will mess up the GUI, but work...
bool FitBetaSpectrum = false;  // whether or not to do a (quick) beta fit on the S1L of all events
bool GaugeNDiscriminator = false;  // does a 2D gaussian fit on the logS2maxLoS1L vs F90 distribution

// Factor stretches (>1) or squeezes (<1) the plots, such that they are fitting the screen
float Factor = 0.80;

// X and Y offset determines starting position of first plot (default x=0, y=0)
int XOffset = 10;  // Alex's Home Zürich: 10, Alex's Home Basel: 0
int YOffset = 35;  // Alex's Home Zürich: 35, Alex's Home Basel: 880

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<string> GetStdoutFromCommand(string cmd);
double gauss_norm(double x, double* par);
double ar39_theo(double x, double* par);
double ar39_theo_smeared(double* t, double* par);
double ar39_theo_smeared_exp(double* t, double* par);
double gauss_2D(double* x, double* par);


// Main function
void CombineWarios() {
  input_output();
  string OutputFileName("SuperWario.root");

  // Extract list of filenames
  string command = "ls " + InputPath + "WarioRun*-" + DataType + "-" + Source + ".root";
  vector<string> pre_file_list = GetStdoutFromCommand(command);
  int start = InputPath.length() + 8; // getting the starting index of run number 00xx

  // Looping through the files and keeping the ones with appropriate run number
  // mainly to get the number of files to be processed for progress information...
  int files_to_be_processed = 0;
  string first_file;
  vector<string> consecutive_files;
  for (vector<string>::iterator t = pre_file_list.begin(); t != pre_file_list.end(); ++t) {
    const int run_nr_temp = stoi(t->substr(start, 4));
    if (run_nr_temp >= MinRunNr && run_nr_temp <= MaxRunNr) {
      if (files_to_be_processed == 0) {
        first_file = t->substr(0, t->length() - 1);
      } else {
        consecutive_files.push_back(t->c_str());
      }
      files_to_be_processed++;
    }
  }

  // Reading in first file and initialising histograms
  cout << "=========================================================================================" << endl;
  cout << "Initialising:" << first_file  << " (1/" << files_to_be_processed << ")" << endl;
  TFile* StandardFile = new TFile(first_file.c_str());

  auto aS1L_histo = (TH1I*)StandardFile->Get("aS1L");
  auto mS1L_histo = (TH1I*)StandardFile->Get("mS1L");
  auto matched_logS2maxLoS1L_vs_F90_histo = (TH2I*)StandardFile->Get("matched_logS2maxLoS1L_vs_F90");
  auto matched_ndisc_histo = (TH1I*)StandardFile->Get("matchedNdisc");
  auto pre_S1L_histo = (TH1I*)StandardFile->Get("pre_S1L");
  auto pre_F90_vs_S1L_histo = (TH2I*)StandardFile->Get("pre_F90_vs_S1L");
  TH2I* NDisc_cut_F90_vs_S1L_histo;
  TH2I* NDisc_cut_F90_vs_S1L_zoom_histo;
  if (SavePlotsThesis) {
    NDisc_cut_F90_vs_S1L_histo = (TH2I*)StandardFile->Get("NDisc_cut_F90_vs_S1L");
    NDisc_cut_F90_vs_S1L_zoom_histo = (TH2I*)StandardFile->Get("NDisc_cut_F90_vs_S1L_zoom");
  }
  auto pre_logS2maxLoS1L_vs_F90_histo = (TH2I*)StandardFile->Get("pre_logS2maxLoS1L_vs_F90");
  auto pre_ndisc_histo = (TH1I*)StandardFile->Get("preNdisc");
  auto S1maxFracPMT_histo = (TH1I*)StandardFile->Get("S1maxFracPMT");
  auto S2maxmaxFracPMT_histo = (TH1I*)StandardFile->Get("S2maxmaxFracPMT");
  auto S2maxFracS2L_histo = (TH1I*)StandardFile->Get("S2maxFracS2L");
  //auto QMatch_histo = (TH1I*)StandardFile->Get("QMatch");
  auto QMatch_vs_DT_histo = (TH2I*)StandardFile->Get("QMatch_vs_DT");
  auto S2maxLSCHI2_histo = (TH1I*)StandardFile->Get("S2maxLSCHI2");
  auto RoI_F90_vs_S1L_histo = (TH2I*)StandardFile->Get("RoI_F90_vs_S1L");
  auto RoI_logS2maxLoS1L_vs_F90_histo = (TH2I*)StandardFile->Get("RoI_logS2maxLoS1L_vs_F90");
  auto RoI_ndisc_histo = (TH1I*)StandardFile->Get("RoI_Ndisc");
  auto QMatch_histo = (TH1I*)StandardFile->Get("QMatch");
  auto S2maxLS_histo = (TH2I*)StandardFile->Get("S2max_LS_position");
  auto Radius_vs_DT_histo = (TH2I*)StandardFile->Get("Radius_vs_DT");
  auto cQMatch_histo = (TH1I*)StandardFile->Get("cQMatch");
  auto DT_vs_S1TTR_histo = (TH2I*)StandardFile->Get("DT_vs_TTR");
  auto DT_histo = (TH1I*)StandardFile->Get("S2maxDT");
  auto S1Nx_histo = (TH1I*)StandardFile->Get("S1Nx");
  auto nS2_histo = (TH1I*)StandardFile->Get("nS2");
  auto S2maxLY_histo = (TH1I*)StandardFile->Get("S2maxLY");
  auto S2LY_histo = (TH1I*)StandardFile->Get("S2LY");
  auto S1LY_histo = (TH1I*)StandardFile->Get("S1LY");
  // Average
  auto Average_S2maxFracS2L_vs_Position_histo = (TH2I*)StandardFile->Get("Average_S2maxFracS2L_vs_position");
  auto Sum_S2maxFracS2L = (TH2I*)StandardFile->Get("S2maxFracS2L_sums");
  auto Entries_S2maxFracS2L = (TH2I*)StandardFile->Get("S2maxFracS2L_entries");

  TList* treeList = new TList;
  treeList->Add(StandardFile->Get("HebingTree"));

  // StandardFile->ls();
  // KEY: TTree	HebingTree;1	HebingTree
  // KEY: TH1F	aS1L;1	all S1 L
  // KEY: TH1F	mS1L;1	matched S1 L
  // KEY: TH1F	pre_S1L;1	Preselected S1 L
  // KEY: TH2F	pre_F90_vs_S1L;1	Preselected F90 vs S1L
  // KEY: TH1F	S1maxFracPMT;1	S1maxFracPMT
  // KEY: TH1F	S2maxmaxFracPMT;1	S2maxmaxFracPMT
  // KEY: TH1F	S2maxFracS2L;1	S2maxFracS2L
  // KEY: TH1F	QMatch;1	QMatch
  // KEY: TH2F	QMatch_vs_DT;1	QMatch vs DT
  // KEY: TH1F	S2maxLSCHI2;1	S2maxLS CHI2
  // KEY: TH2F	RoI_F90_vs_S1L;1	F90 vs S1L
  // KEY: TH2F	RoI_logS2maxLoS1L_vs_F90;1	logS2maxLoS1L vs F90

  // Adding all consecutive files
  int processed_files = 1;
  for (vector<string>::iterator t = consecutive_files.begin(); t != consecutive_files.end(); ++t) {
    const int run_nr_temp = stoi(t->substr(start, 4));
    const string ending_temp = t->substr(start);  // extracting e. g. "0089-DP2-neutron.root"
    const string name_temp = t->substr(0, start + ending_temp.length() - 1);  // removing additional characters
    processed_files++;
    cout << "Adding:" << name_temp << " (" << processed_files << "/" << files_to_be_processed << ")" << endl;
    TFile* ToBeAdded = new TFile(name_temp.c_str());
    treeList->Add(ToBeAdded->Get("HebingTree"));

    // Adding histograms
    matched_logS2maxLoS1L_vs_F90_histo->Add((TH2I*)ToBeAdded->Get("matched_logS2maxLoS1L_vs_F90"));
    // matched_ndisc_histo->Add((TH1I*)ToBeAdded->Get("matchedNdisc"));
    pre_S1L_histo->Add((TH1I*)ToBeAdded->Get("pre_S1L"));
    pre_F90_vs_S1L_histo->Add((TH2I*)ToBeAdded->Get("pre_F90_vs_S1L"));
    if (SavePlotsThesis) {
      // NDoisc_cut_F90_vs_S1L_hist->Add((TH2I*)ToBeAdded->Get("NDisc_cut_F90_vs_S1L"));
      // NDisc_cut_F90_vs_S1L_zoom_histo->Add((TH2I*)ToBeAdded->Get("NDisc_cut_F90_vs_S1L_zoom"));
    }
    pre_logS2maxLoS1L_vs_F90_histo->Add((TH2I*)ToBeAdded->Get("pre_logS2maxLoS1L_vs_F90"));
    RoI_F90_vs_S1L_histo->Add((TH2I*)ToBeAdded->Get("RoI_F90_vs_S1L"));
    RoI_logS2maxLoS1L_vs_F90_histo->Add((TH2I*)ToBeAdded->Get("RoI_logS2maxLoS1L_vs_F90"));
    S1maxFracPMT_histo->Add((TH1I*)ToBeAdded->Get("S1maxFracPMT"));
    S2maxFracS2L_histo->Add((TH1I*)ToBeAdded->Get("S2maxFracS2L"));
    QMatch_histo->Add((TH1I*)ToBeAdded->Get("QMatch"));
    S2maxLS_histo->Add((TH2I*)ToBeAdded->Get("S2max_LS_position"));
    Radius_vs_DT_histo->Add((TH2I*)ToBeAdded->Get("Radius_vs_DT"));
    cQMatch_histo->Add((TH1I*)ToBeAdded->Get("cQMatch"));
    QMatch_vs_DT_histo->Add((TH2I*)ToBeAdded->Get("QMatch_vs_DT"));
    DT_vs_S1TTR_histo->Add((TH2I*)ToBeAdded->Get("DT_vs_TTR"));
    DT_histo->Add((TH1F*)ToBeAdded->Get("S2maxDT"));
    S1Nx_histo->Add((TH1I*)ToBeAdded->Get("S1Nx"));
    nS2_histo->Add((TH1I*)ToBeAdded->Get("nS2"));
    S2maxLY_histo->Add((TH1F*)ToBeAdded->Get("S2maxLY"));
    S2LY_histo->Add((TH1F*)ToBeAdded->Get("S2LY"));
    S1LY_histo->Add((TH1F*)ToBeAdded->Get("S1LY"));
    // Average
    Average_S2maxFracS2L_vs_Position_histo->Add((TH2I*)ToBeAdded->Get("Average_S2maxFracS2L_vs_position"));
    Sum_S2maxFracS2L->Add((TH1I*)ToBeAdded->Get("S2maxFracS2L_sums"));
    Entries_S2maxFracS2L->Add((TH1I*)ToBeAdded->Get("S2maxFracS2L_entries"));
  }
  
  // For the plots of average values of variables (vs positions), calculate the averages bin by bin
  int nr_bins = 100;
  int filled_bins = 0;
  int nr_global_bins = nr_bins * nr_bins;
  for (int bin_index = 1; bin_index <= nr_global_bins; bin_index++) {
    if (Entries_S2maxFracS2L->GetBinContent(bin_index) >= 1) {
      // std::cout << "Bin: " << bin_index << " out of " << nr_global_bins << std::endl;
      // std::cout << "Sum histo content: " << Sum_S2maxFracS2L->GetBinContent(bin_index) << std::endl;
      // std::cout << "Entries histo content: " << Entries_S2maxFracS2L->GetBinContent(bin_index) << std::endl;
      // std::cout << "Resulting average: " << Sum_S2maxFracS2L->GetBinContent(bin_index) / Entries_S2maxFracS2L->GetBinContent(bin_index) << std::endl;
      Average_S2maxFracS2L_vs_Position_histo->SetBinContent(bin_index, TMath::Max(Sum_S2maxFracS2L->
                GetBinContent(bin_index) / Entries_S2maxFracS2L->GetBinContent(bin_index), 0.901));
      filled_bins += 1;
    }
  }
  // cout << "=========================================================================================" << endl;
  // cout << "Non-empty bins: " << filled_bins << endl;
  // cout << "=========================================================================================" << endl;

  // Terminal information on scattering multiplicity
  int total_entries_fracS2L = S2maxFracS2L_histo->GetEntries();
  int entries_fracS2L_equal_one = S2maxFracS2L_histo->GetBinContent(S2maxFracS2L_histo->GetNbinsX() + 1);
  // cout << "Number of single scatters: " << entries_fracS2L_equal_one << endl;
  string ratio_fracS2L = ", \% (# single scatters) = " + to_string((float)entries_fracS2L_equal_one / (float)total_entries_fracS2L);
  S2maxFracS2L_histo->SetTitle((S2maxFracS2L_histo->GetTitle() + ratio_fracS2L).c_str());

  // Overflow bin correction
  S2maxFracS2L_histo->AddBinContent(S2maxFracS2L_histo->GetNbinsX(),
          S2maxFracS2L_histo->GetBinContent(S2maxFracS2L_histo->GetNbinsX() + 1));
  S2maxFracS2L_histo->SetBinContent(S2maxFracS2L_histo->GetNbinsX() + 1, 0);
  
  Average_S2maxFracS2L_vs_Position_histo->AddBinContent(Average_S2maxFracS2L_vs_Position_histo->GetNbinsZ(),
          Average_S2maxFracS2L_vs_Position_histo->GetBinContent(Average_S2maxFracS2L_vs_Position_histo->GetNbinsZ() + 1));
  Average_S2maxFracS2L_vs_Position_histo->SetBinContent(Average_S2maxFracS2L_vs_Position_histo->GetNbinsZ() + 1, 0);
  
  // cout << "entries overflow bin = " << entries_fracS2L_equal_one << endl;
  // cout << "entries total = " << total_entries_fracS2L << endl;
  // cout << "ratio = " << (float)entries_fracS2L_equal_one/(float)total_entries_fracS2L << endl;

  /*
  // Cut statistics (overall)
  int n_preselected = pre_S1L_histo->GetEntries();
  int nWarioCuts = 0;
  // Printing cut statistics
  cout << "Total number of S1 events " << n_aS1events << endl;
  cout << "Total number of matched S1 events " << n_mS1events << " ("  << 100.0 * (n_mS1events)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of unmatched S1 events " << n_umS1events << " ("  << 100.0 * (n_umS1events)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of preselected events " << n_preselected << " ("  << 100.0 * (n_preselected)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of S2FracL events " << n_S2FracL << " ("  << 100.0 * (n_S2FracL)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of QMatch events " << n_QMatch << " ("  << 100.0 * (n_QMatch)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of entries in RoI - I " << n_RoI_S1L << " ("  << 100.0 * (n_RoI_S1L)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of entries in RoI - II " << n_RoI_F90 << " ("  << 100.0 * (n_RoI_F90)/(n_aS1events) << " \%)" << endl;
  cout << "Total number of entries in Wario file " << nWarioCuts << " ("  << 100.0 * (nWarioCuts)/(n_aS1events) << " \%)" << endl;
  cout << "=========================================================================================" << endl;
  */
  //---------------------------------------------------
  //                      Plots
  //---------------------------------------------------
  // Plot positions
  int x_start = 80 + XOffset;
  int y_start = 0 + YOffset;
  float x_step = 400 * Factor;
  float y_step = 450 * Factor;
  float x_size = 400 * Factor;
  float y_size = 400 * Factor;

  // Adjust canvas size for "png"s and "jpeg"s because their resolution is the one of the initial canvas
  if ((SavePlots & (PlotFileType == "jpeg" | PlotFileType == "png")) | FitBetaSpectrum | GaugeNDiscriminator) {
    x_size = 1440;
    y_size = 1200;
    x_step = 0;
    y_step = 0;
    if (FitBetaSpectrum | GaugeNDiscriminator) {ShowPlots = false;}
  }

  // all other plots
  if (ShowPlots) {
    string bump_name;
    if (bump == "both") {
      bump_name = "";
    } else if (bump == "upper") {
      bump_name = "_upper";
    } else if (bump == "lower") {
      bump_name = "_lower";
    }
    /*
    // matched S1L
    auto c01 = new TCanvas("c01", "c01", x_start + 1 * x_step, y_start + 0 * y_step, x_size, y_size);
    c01->SetGridx();
    c01->SetGridy();
    mS1L_histo->DrawClone("");
    if (SavePlots) {c01->SaveAs((OutputPath + "mS1L" + "." + PlotFileType).c_str());}
    
    // matched logS2maxLoS1L vs S1F90
    auto c12 = new TCanvas("c12", "c12", x_start + 0 * x_step, y_start + 0 * y_step, x_size, y_size);
    c12->SetLogz();
    c12->SetGridx();
    c12->SetGridy();
    matched_logS2maxLoS1L_vs_F90_histo->DrawClone("colsz");
    if (SavePlots) {c12->SaveAs((OutputPath + "matched_logS2maxLoS1L_vs_F90" + "." + PlotFileType).c_str());}

    // matched NDiscriminator
    auto c17 = new TCanvas("c17", "c17", x_start + 1 * x_step, y_start + 0 * y_step, x_size, y_size);
    c17->SetLogy();
    c17->SetGridx();
    c17->SetGridy();
    matched_ndisc_histo->DrawClone("");
    if (SavePlots) {c17->SaveAs((OutputPath + "matched_NDisc" + "." + PlotFileType).c_str());}

    // preselected S1L
    auto c02 = new TCanvas("c02", "c02", x_start + 2 * x_step, y_start + 0 * y_step, x_size, y_size);
    c02->SetGridx();
    c02->SetGridy();
    pre_S1L_histo->DrawClone("");
    if (SavePlots) {c02->SaveAs((OutputPath + "pre_S1L" + "." + PlotFileType).c_str());}
    */
    // preselected S1F90 vs S1L
    auto c03 = new TCanvas("c03", "c03", x_start + 3 * x_step, y_start + 0 * y_step, x_size, y_size);
    c03->SetLogz();
    c03->SetGridx();
    c03->SetGridy();
    // pre_F90_vs_S1L_histo->SetStats(0);
    pre_F90_vs_S1L_histo->GetXaxis()->SetTitle("S1 [p. e.]");
    pre_F90_vs_S1L_histo->GetYaxis()->SetTitle("F90");
    pre_F90_vs_S1L_histo->DrawClone("colsz");
    if (SavePlots) {c03->SaveAs((OutputPath + "pre_F90_vs_S1L" + bump_name + "." + PlotFileType).c_str());}

    // preselected logS2maxLoS1L vs S1F90
    auto c13 = new TCanvas("c13", "c13", x_start + 0 * x_step, y_start + 0 * y_step, x_size, y_size);
    c13->SetLogz();
    c13->SetGridx();
    c13->SetGridy();
    pre_logS2maxLoS1L_vs_F90_histo->DrawClone("colsz");
    if (SavePlots) {c13->SaveAs((OutputPath + "pre_logS2maxLoS1L_vs_F90" + bump_name + "." + PlotFileType).c_str());}
    /*
    // preselected NDiscriminator
    auto c18 = new TCanvas("c18", "c18", x_start + 1 * x_step, y_start + 0 * y_step, x_size, y_size);
    c18->SetLogy();
    c18->SetGridx();
    c18->SetGridy();
    pre_ndisc_histo->DrawClone("");
    if (SavePlots) {c18->SaveAs((OutputPath + "pre_NDisc" + "." + PlotFileType).c_str());}
    */
    // S1maxFracPMT
    auto c04 = new TCanvas("c04", "c04", x_start + 0 * x_step, y_start + 1 * y_step, x_size, y_size);
    c04->SetLogy();
    c04->SetGridx();
    c04->SetGridy();
    S1maxFracPMT_histo->DrawClone("");
    if (SavePlots) {c04->SaveAs((OutputPath + "S1maxFracPMT" + bump_name + "." + PlotFileType).c_str());}
    /*
    // S2maxmaxFracPMT
    auto c05 = new TCanvas("c05", "c05", x_start + 1 * x_step, y_start + 1 * y_step, x_size, y_size);
    c05->SetLogy();
    c05->SetGridx();
    c05->SetGridy();
    S2maxmaxFracPMT_histo->DrawClone("");
    if (SavePlots) {c05->SaveAs((OutputPath + "S2maxmaxFracPMT" + "." + PlotFileType).c_str());}
    */
    // S2maxFracS2L
    auto c06 = new TCanvas("c06", "c06", x_start + 2 * x_step, y_start + 1 * y_step, x_size, y_size);
    c06->SetLogy();
    c06->SetGridx();
    c06->SetGridy();
    S2maxFracS2L_histo->DrawClone("");
    if (SavePlots) {c06->SaveAs((OutputPath + "S2maxFracS2L" + bump_name + "." + PlotFileType).c_str());}
    
    // QMatch
    auto c07 = new TCanvas("c07", "c07", x_start + 3 * x_step, y_start + 1 * y_step, x_size, y_size);
    c07->SetLogy();
    c07->SetGridx();
    c07->SetGridy();
    QMatch_histo->DrawClone("");
    if (SavePlots) {c07->SaveAs((OutputPath + "QMatch" + bump_name + "." + PlotFileType).c_str());}

    // auto c77 = new TCanvas("c77", "c77", x_start + 3 * x_step, y_start + 1 * y_step, x_size, y_size);
    // c77->SetLogy();
    // c77->SetGridx();
    // c77->SetGridy();
    // cQMatch_histo->DrawClone("");
    // if (SavePlots) {c77->SaveAs((OutputPath + "QMatch_current" + bump_name + "." + PlotFileType).c_str());}
    
    // QMatch vs S2maxDT
    auto c08 = new TCanvas("c08", "c08", x_start + 0 * x_step, y_start + 2 * y_step, x_size, y_size);
    c08->SetLogz();
    c08->SetGridx();
    c08->SetGridy();
    QMatch_vs_DT_histo->DrawClone("colsz");
    if (SavePlots) {c08->SaveAs((OutputPath + "QMatch_vs_DT" + bump_name + "." + PlotFileType).c_str());}
    /*
    // S2maxLSCHI2
    auto c09 = new TCanvas("c09", "c09", x_start + 1 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c09->SetLogy();
    c09->SetGridx();
    c09->SetGridy();
    S2maxLSCHI2_histo->DrawClone("");
    if (SavePlots) {c09->SaveAs((OutputPath + "S2maxLSCHI2" + "." + PlotFileType).c_str());}
    */
    // RoI: S1F90 vs S1L
    auto c10 = new TCanvas("c10", "c10", x_start + 2 * x_step, y_start + 2 * y_step, x_size, y_size);
    c10->SetLogz();
    c10->SetGridx();
    c10->SetGridy();
    RoI_F90_vs_S1L_histo->DrawClone("colsz");
    if (SavePlots) {c10->SaveAs((OutputPath + "RoI_F90_vs_S1L" + bump_name + "." + PlotFileType).c_str());}

    // RoI: logS2maxLoS1L vs S1F90
    // cout << "S2/S1 vs F90:" << endl;
    // cout << "Entries: " << RoI_logS2maxLoS1L_vs_F90_histo->GetEntries() << endl;
    // cout << "Integral: " << RoI_logS2maxLoS1L_vs_F90_histo->Integral(1, 200, 1, 200) << endl;
    // cout << "Underflow: " << RoI_logS2maxLoS1L_vs_F90_histo->GetBinContent(0) << endl;
    // cout << "Overflow: " << RoI_logS2maxLoS1L_vs_F90_histo->GetBinContent(40662) << endl;
    // cout << RoI_logS2maxLoS1L_vs_F90_histo->GetBin(201,201) << endl;
    // RoI_logS2maxLoS1L_vs_F90_histo->GetYaxis()->SetRangeUser(-11., 6.);
    auto c11 = new TCanvas("c11", "c11", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    c11->SetLogz();
    c11->SetGridx();
    c11->SetGridy();
    RoI_logS2maxLoS1L_vs_F90_histo->DrawClone("colsz");
    if (SavePlots) {c11->SaveAs((OutputPath + "RoI_logS2maxLoS1L_vs_F90" + bump_name + "." + PlotFileType).c_str());}
    /*
    // RoI: NDiscriminator
    auto c19 = new TCanvas("c19", "c19", x_start + 1 * x_step, y_start + 0 * y_step, x_size, y_size);
    c19->SetLogy();
    c19->SetGridx();
    c19->SetGridy();
    RoI_ndisc_histo->DrawClone("");
    if (SavePlots) {c19->SaveAs((OutputPath + "RoI_NDisc" + "." + PlotFileType).c_str());}
    */
    auto c21 = new TCanvas("c21", "c21", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    c21->SetLogz();
    c21->SetGridx();
    c21->SetGridy();
    S2maxLS_histo->DrawClone("colsz");
    if (SavePlots) {c21->SaveAs((OutputPath + "S2maxLSY_vs_S2maxLSX" + bump_name + "." + PlotFileType).c_str());}

    auto c22 = new TCanvas("c22", "c22", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    c22->SetLogz();
    c22->SetGridx();
    c22->SetGridy();
    // Radius_vs_DT_histo->SetStats(0);
    Radius_vs_DT_histo->DrawClone("colsz");
    if (SavePlots) {c22->SaveAs((OutputPath + "Radius_vs_DT" + bump_name + "." + PlotFileType).c_str());}

    auto c23 = new TCanvas("c23", "c23", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    c23->SetLogz();
    c23->SetGridx();
    c23->SetGridy();
    DT_vs_S1TTR_histo->DrawClone("colsz");
    if (SavePlots) {c23->SaveAs((OutputPath + "DT_vs_S1TTR" + bump_name + "." + PlotFileType).c_str());}

    auto c96 = new TCanvas("c96", "c96", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c96->SetLogy();
    c96->SetGridx();
    c96->SetGridy();
    DT_histo->DrawClone("");
    if (SavePlots) {c96->SaveAs((OutputPath + "DT" + bump_name + "." + PlotFileType).c_str());}

    auto c95 = new TCanvas("c95", "c95", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c95->SetLogy();
    c95->SetGridx();
    c95->SetGridy();
    S1Nx_histo->DrawClone("");
    if (SavePlots) {c95->SaveAs((OutputPath + "S1Nx" + bump_name + "." + PlotFileType).c_str());}

    auto c94 = new TCanvas("c94", "c94", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c94->SetLogy();
    c94->SetGridx();
    c94->SetGridy();
    nS2_histo->DrawClone("");
    if (SavePlots) {c94->SaveAs((OutputPath + "nS2" + bump_name + "." + PlotFileType).c_str());}

    auto c93 = new TCanvas("c93", "c93", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c93->SetLogy();
    c93->SetGridx();
    c93->SetGridy();
    S2maxLY_histo->DrawClone("");
    if (SavePlots) {c93->SaveAs((OutputPath + "S2maxLY" + bump_name + "." + PlotFileType).c_str());}

    auto c92 = new TCanvas("c92", "c92", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c92->SetLogy();
    c92->SetGridx();
    c92->SetGridy();
    S2LY_histo->DrawClone("");
    if (SavePlots) {c92->SaveAs((OutputPath + "S2LY" + bump_name + "." + PlotFileType).c_str());}

    auto c91 = new TCanvas("c91", "c91", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c91->SetLogy();
    c91->SetGridx();
    c91->SetGridy();
    S1LY_histo->DrawClone("");
    if (SavePlots) {c91->SaveAs((OutputPath + "S1LY" + bump_name + "." + PlotFileType).c_str());}

    auto c24 = new TCanvas("c24", "c24", x_start + 3 * x_step, y_start + 2 * y_step, x_size, y_size);
    // c24->SetLogz();
    c24->SetGridx();
    c24->SetGridy();
    Average_S2maxFracS2L_vs_Position_histo->GetZaxis()->SetRangeUser(0.9, 1.0);
    // Average_S2maxFracS2L_vs_Position_histo->SetStats(0);
    Average_S2maxFracS2L_vs_Position_histo->DrawClone("colsz");
    c24->SetRightMargin(0.16);
    if (SavePlots) {c24->SaveAs((OutputPath + "Average_S2maxFracS2L_vs_position" + bump_name + "." + PlotFileType).c_str());}
  }

  // Create a new file
  TFile newfile((InputPath + OutputFileName).c_str(), "recreate");

  if (MergeTrees) {
    TTree* newtree = TTree::MergeTrees(treeList);
    newtree->Write();
  }

  // Saving histograms to newfile
  matched_logS2maxLoS1L_vs_F90_histo->Write();
  pre_S1L_histo->Write();
  pre_F90_vs_S1L_histo->Write();
  pre_logS2maxLoS1L_vs_F90_histo->Write();
  RoI_F90_vs_S1L_histo->Write();
  RoI_logS2maxLoS1L_vs_F90_histo->Write();
  S1maxFracPMT_histo->Write();
  S2maxFracS2L_histo->Write();
  QMatch_histo->Write();
  S2maxLS_histo->Write();
  Radius_vs_DT_histo->Write();
  cQMatch_histo->Write();
  QMatch_vs_DT_histo->Write();
  DT_vs_S1TTR_histo->Write();
  Average_S2maxFracS2L_vs_Position_histo->Write();
  Sum_S2maxFracS2L->Write();
  Entries_S2maxFracS2L->Write();
  S1Nx_histo->Write();
  nS2_histo->Write();
  S2maxLY_histo->Write();

  // Closing newfile
  newfile.Close();

  // Gauge the n_discriminator variable by fitting the logS2LoS1L vs F90 distribution to a 2D gaussian
  if (GaugeNDiscriminator) {
    string bump_name;
    if (bump == "both") {
      bump_name = "";
    } else if (bump == "upper") {
      bump_name = "_upper";
    } else if (bump == "lower") {
      bump_name = "_lower";
    }
    // Rebinning the histogram
    // RoI_logS2maxLoS1L_vs_F90_histo->Rebin2D(2, 2);  // Rebining x and y axis (mergin n bins to one)
    // int n_bins_x = RoI_logS2maxLoS1L_vs_F90_histo->GetNbinsX();
    // int n_bins_y = RoI_logS2maxLoS1L_vs_F90_histo->GetNbinsY();
    // cout << "X_bins: " << n_bins_x << endl;
    // cout << "Y_bins: " << n_bins_y << endl;

    // Fit of 2D gaussian
    cout << "\n===================================================================================="
                 "=====" << endl;
    cout << "Fitting logS2maxLoS1L vs S1F90 distribution..." << endl;
    cout << "======================================================================================"
                 "===\n" << endl;
    // Initial parameter guesses and bounds
    // TF2 * f_n_disc = new TF2("f_n_disc", gauss_2D, 0.0, 1.0, -3.0, 5.0, 5);
    TF2* f_n_disc;
    if (Source == "neutron" && bump == "lower") {
      f_n_disc = new TF2("f_n_disc", gauss_2D, 0.50, 0.90, -1.4, 0.2, 6);
      f_n_disc->SetParameter(0, 5.0); // par[5] = amplitude
      f_n_disc->SetParLimits(0, 0.0, 10.0);
      f_n_disc->SetParameter(1, 0.75); // par[0] = mean F90
      f_n_disc->SetParLimits(1, 0.6, 0.9);
      f_n_disc->SetParameter(2, 0.10); // par[2] = sigma F90
      f_n_disc->SetParLimits(2, 0, 0.25);
      f_n_disc->SetParameter(3, -0.3); // par[1] = mean logS2maxLoS1L
      f_n_disc->SetParLimits(3, -1.2, 0.2);
      f_n_disc->SetParameter(4, 0.3); // par[3] = sigma logS2maxLoS1L
      f_n_disc->SetParLimits(4, 0.1, 0.5);
      f_n_disc->SetParameter(5, 0.0); // par[4] = correlation coeff
      f_n_disc->SetParLimits(5, -0.2, 0.2);
    } else {
      f_n_disc = new TF2("f_n_disc", gauss_2D, 0.50, 0.90, 0.2, 3.4, 6);
      f_n_disc->SetParameter(0, 5.0); // par[0] = amplitude
      f_n_disc->SetParLimits(0, 0.0, 100.0);
      f_n_disc->SetParameter(1, 0.73); // par[1] = mean F90
      f_n_disc->SetParLimits(1, 0.6, 0.9);
      f_n_disc->SetParameter(2, 0.10); // par[2] = sigma F90
      f_n_disc->SetParLimits(2, 0, 0.25);
      f_n_disc->SetParameter(3, 1.5); // par[3] = mean logS2maxLoS1L
      f_n_disc->SetParLimits(3, 0.2, 3.0);
      f_n_disc->SetParameter(4, 0.70); // par[4] = sigma logS2maxLoS1L
      f_n_disc->SetParLimits(4, 0.4, 1.0);
      f_n_disc->SetParameter(5, 0.0); // par[5] = correlation coeff
      f_n_disc->SetParLimits(5, -0.2, 0.2);
    }
    // Setting fit method and evaluating the fit
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    // TFitResultPtr r = pre_logS2maxLoS1L_vs_F90_histo->Fit("f_n_disc", "LR0S");  // Loglikelihood (chi^2)
    TFitResultPtr r = RoI_logS2maxLoS1L_vs_F90_histo->Fit("f_n_disc", "R0S");  // Least squares
    // string fit_method = "LogL";  // to be part of the final filenames...
    string fit_method = "LS";  // to be part of the final filenames...
    // Print fit results
    // r->Print();
    TMatrixDSym C = r->GetCorrelationMatrix();
    C.Print();
    // Extracting resulting center and half axes and drawing the resulting ellipse
    double mean_f90 = f_n_disc->GetParameter(1);
    double sigma_f90 = f_n_disc->GetParameter(2);
    double mean_logS2maxLoS1L = f_n_disc->GetParameter(3);
    double sigma_logS2maxLoS1L = f_n_disc->GetParameter(4);
    double correl_coeff = f_n_disc->GetParameter(5);
    double theta = 0.5 * TMath::ATan(2.0 * correl_coeff * sigma_f90 * sigma_logS2maxLoS1L /
          (sigma_f90 * sigma_f90 - sigma_logS2maxLoS1L * sigma_logS2maxLoS1L));
    double amp = f_n_disc->GetParameter(0);
    cout << "\n===================================================================================="
                 "=====" << endl;
    cout << "Resulting amplitude: " << amp << endl;
    cout << "Resulting S1F90 center: " << mean_f90 << endl;
    cout << "Resulting S1F90 half axis: " << sigma_f90 << endl;
    cout << "Resulting logS2maxLoS1L center: " << mean_logS2maxLoS1L << endl;
    cout << "Resulting logS2maxLoS1L half axis: " << sigma_logS2maxLoS1L << endl;
    cout << "======================================================================================"
                 "===\n" << endl;
    // Setting drawing options for fit results
    f_n_disc->SetRange(0.0, -3.0, 1.0, 5.0);
    f_n_disc->SetNpx(50);  // number of points when drawing the function
    f_n_disc->SetNpy(50);  // number of points when drawing the function
    // f_n_disc->SetLineWidth(2);
    // f_n_disc->SetLineColor(2);

    // Drawing 3D histogram
    auto c14 = new TCanvas("c14", "c14", x_start + 0 * x_step, y_start + 0 * y_step, x_size, y_size);
    // c14->SetLogz();
    c14->SetGridx();
    c14->SetGridy();
    gStyle->SetPalette(1);
    RoI_logS2maxLoS1L_vs_F90_histo->SetTitle("");
    // RoI_logS2maxLoS1L_vs_F90_histo->SetStats(0);
    RoI_logS2maxLoS1L_vs_F90_histo->GetXaxis()->SetTitle("F90");
    RoI_logS2maxLoS1L_vs_F90_histo->GetYaxis()->SetTitle("log_{10}#left(#frac{S2}{S1}#right)");
    RoI_logS2maxLoS1L_vs_F90_histo->SetTitleSize(0.06, "xyz");
    RoI_logS2maxLoS1L_vs_F90_histo->SetLabelSize(0.05, "xyz");
    RoI_logS2maxLoS1L_vs_F90_histo->DrawClone("lego2");
    f_n_disc->Draw("surf same");
    c14->SetTopMargin(0.01);
    c14->SetRightMargin(0.02);
    c14->SaveAs((OutputPath + "NDiscriminator/RoI_logS2maxLoS1L_vs_F90_with_" + fit_method + "fit_3D"
            + bump_name + "." + PlotFileType).c_str());
    // c14->GetListOfPrimitives()->Remove(f_n_disc);

    // Drawing 2D histogram
    auto c16 = new TCanvas("c16", "c16", x_start + 0 * x_step, y_start + 0 * y_step, x_size, y_size);
    // c16->SetLogz();
    c16->SetGridx();
    c16->SetGridy();
    // gStyle->SetPalette(1);
    RoI_logS2maxLoS1L_vs_F90_histo->SetLabelSize(0.06, "xyz");
    RoI_logS2maxLoS1L_vs_F90_histo->SetTitleOffset(0.9, "xy");
    RoI_logS2maxLoS1L_vs_F90_histo->DrawClone("");
    // f_n_disc->Draw("cont1 same");
    TEllipse* el16_1 = new TEllipse(mean_f90, mean_logS2maxLoS1L, sigma_f90, sigma_logS2maxLoS1L);
    el16_1->SetFillStyle(0);
    el16_1->SetLineColor(kRed);
    el16_1->SetLineWidth(6);
    // el16_1->SetFillColorAlpha(kRed, 0.5);
    el16_1->Draw("same");
    TEllipse* el16_2 = new TEllipse(mean_f90, mean_logS2maxLoS1L, 2 * sigma_f90, 2 * sigma_logS2maxLoS1L);
    el16_2->SetFillStyle(0);
    el16_2->SetLineColor(kGreen);
    el16_2->SetLineWidth(6);
    // el16_2->SetFillColorAlpha(kRed, 0.5);
    el16_2->Draw("same");
    TEllipse* el16_3 = new TEllipse(mean_f90, mean_logS2maxLoS1L, 3 * sigma_f90, 3 * sigma_logS2maxLoS1L);
    el16_3->SetFillStyle(0);
    el16_3->SetLineColor(kBlue);
    el16_3->SetLineWidth(6);
    // el16_3->SetFillColorAlpha(kRed, 0.5);
    el16_3->Draw("same");
    c16->SetTopMargin(0.03);
    c16->SetRightMargin(0.02);
    c16->SetLeftMargin(0.15);
    c16->SetBottomMargin(0.12);
    c16->SaveAs((OutputPath + "NDiscriminator/RoI_logS2maxLoS1L_vs_F90_with_" + fit_method + "fit_2D"
            + bump_name + "." + PlotFileType).c_str());
    // c16->GetListOfPrimitives()->Remove(f_n_disc);

    // Drawing 2D histogram with resulting 1-sigma ellipse
    auto c15 = new TCanvas("c15", "c15", x_start + 0 * x_step, y_start + 0 * y_step, x_size, y_size);
    c15->SetLogz();
    c15->SetGridx();
    c15->SetGridy();
    gStyle->SetPalette(kBird);
    RoI_logS2maxLoS1L_vs_F90_histo->DrawClone("colsz");
    TEllipse* el1 = new TEllipse(mean_f90, mean_logS2maxLoS1L, sigma_f90, sigma_logS2maxLoS1L);
    el1->SetFillStyle(3001);
    el1->SetFillColorAlpha(kRed, 0.5);
    el1->Draw("same");
    c15->SetTopMargin(0.03);
    c15->SetRightMargin(0.14);
    c15->SetLeftMargin(0.15);
    c15->SetBottomMargin(0.12);
    c15->SaveAs((OutputPath + "NDiscriminator/RoI_logS2maxLoS1L_vs_F90_with_" + fit_method + "ellipse"
            + bump_name + "." + PlotFileType).c_str());

    ofstream data_file;
    data_file.open(OutputPath + "NDiscriminator/fit_data.csv");
    data_file << "amplitude,mean_F90,sigma_F90,mean_logS2maxLoS1L,sigma_logS2maxLoS1L,correlation_coeff,theta" << endl;
    data_file << amp << "," << mean_f90 << "," << sigma_f90  << "," << mean_logS2maxLoS1L << ","
              << sigma_logS2maxLoS1L << "," << correl_coeff << "," << theta << endl;
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
  // par[0]  -  amplitude
  // par[1]  -  mean x
  // par[2]  -  sigma x
  // par[3]  -  mean y
  // par[4]  -  sigma y
  if (par[2] > 0 && par[4] > 0 && par[5] > 0) {
    double rxx = TMath::Power((x[0] - par[1]) / par[2], 2.0);
    double ryy = TMath::Power((x[1] - par[3]) / par[4], 2.0);
    double rxy = 2.0 * par[5] * (x[0] - par[1]) * (x[1] - par[3]) / par[2] / par[4];
    return TMath::InvPi() * par[0] / (2.0 * par[2] * par[4] * TMath::Sqrt(1 - par[5] * par[5]))
              * TMath::Exp(-(rxx + ryy - rxy) / 2.0 / (1 - par[5] * par[5]));
  }
  else {return 0.0;}
}


// Normed gauss distribution
double gauss_norm(double x, double* par) {
  if(par[0] < 0) {return 0;} //sigma must be non-negative
  else if(par[0] == 0) {
    if(x == 0) {return 1;}
    else {return 0;}
  }

  return 1. / par[0] / sqrt(2 * TMath::Pi()) * exp( - x*x / (2 * par[0]*par[0]));
}


// Theoretical Ar39 beta spectrum
double ar39_theo(double x, double* par) {
  // beta spectrum
  double me   = 511.0 ;   //number of photons detected if Edep = 511
  double Qval = 565.5 ;   //maximum kinetic energy in unit of me*c2

  double t0 = Qval / me;
  double t  = x / me;     //e- kinetic energy in unit of me


  if(t > t0 || t <= 0) {return 0;}

  double w  = 1 + t;
  double w0 = 1 + t0;
  double p  = sqrt(w*w -1 );

  double stat   = p * w * (w0 - w)*(w0-w);    //statistical part
  double alpha  = 1. / 137;                     //fine structure constant
  double z      = 19;                         //= 18(argon) + 1
  double pi     = TMath::Pi();

  double signum    = 1;
  double eta       = signum * alpha * z * w / p; // in "natural units" : alpha = e^2 / 4pi;
  // in SI units : alpha = 1/4pi/epsilon0 * e^2 / hbar / c
  double fermi     = 2 * pi * eta / (1 - exp(-2 * pi * eta) );
  double Enu       = w0 - w;
  double forbidden = p*p +  Enu*Enu; //=Pe^2 + Enu^2
  double n         = stat * fermi * forbidden;

  return n;
}


// Convolution of the theoretical Ar39 beta spectrum with a gaussian (with energy dependent sigma)
double ar39_theo_smeared(double* t, double* par) {
  // this function returns the value of the convolution function
  // conv between ar39_theo and gauss_norm evaluated at the energy t[0]
  // par[0] = normalization constant
  // par[1] = LY (npe / keV)
  // par[2] = resolution of single photon, in function  ar39_theo_smeared(double* , double*)
  // the resolution spSigma of the gaussian resolution function is :
  // sigma = sqrt( par[3]*par[3] + par[2]*E + (par[4]*E)^2  )
  // par[2] : from  gaussian reasoning, i.e. 1photon --> sigma, n photons --> sigma*sqrt(n)
  // par[3] : electronic noise
  // par[4] : non-uniform light collection
  // t[0] = number of detected photons

  double x = t[0] / par[1] ; //par[1] = npe / keV

  double sigma = 0;
  double min = 0, max = 700 * par[1], step=1; // min, max, integration step, convolution range
  double val = 0;

  //calculate the convolution between ar39_theo and resolution gauss-norm
  double y = 0;
  for(int dummyy=min; dummyy < max; dummyy += step) {
    y = dummyy / par[1];
    sigma = sqrt(par[3] * par[3] + par[2] * y + (par[4] * y) * (par[4] * y));

    //if (x-y) is too far (more than 5sigma) away from the gaussian mean --> ignore it
    if(sigma && fabs(x - y) > 5 * sigma) continue;

    double gausspart = gauss_norm(x - y, &sigma);
    double ar39part = ar39_theo(y, NULL); //argument (doube* par) in ar39_theo is not used at all
    val += gausspart * ar39part;
  }

  //val is now the conv(x),
  //applying overall normalization to conv-function
  val *= par[0];

  return val;
}

// Smeared theoretical beta spectrum of Ar39 plus exponential background
double ar39_theo_smeared_exp(double* t, double* par) {
  //par[0] = norm. const of ar39_theo_smeared
  //par[1] = LY, npe / keV
  //par[2] = single photon resolution
  //par[3] = electronic noise
  //par[4] = non-uniformity
  //par[5] = norm. const of exp. part
  //par[6] = tau
  return ar39_theo_smeared(t, par) + ((par[6]) ? par[5] * exp(-t[0] / par[6]) : 0 );
}