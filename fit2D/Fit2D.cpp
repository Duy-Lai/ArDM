#include <iostream>
#include <vector>
#include "SidebandUtils.h"
#include "SidebandUtils.cpp"
#include "PlotsManager.h"
#include "PlotsManager.cpp"

#define FMT_HEADER_ONLY
#include </home/duylai/cget/include/fmt/core.h>
// 1: "pip install cget"; 2: "cget install fmtlib/fmt"

const vector<double> vec_S1_lower_bound{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
const vector<double> vec_S1_upper_bound{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};
const int max_bin = vec_S1_upper_bound.size();

int min_bin_to_process = 0;
int max_bin_to_process = 10;
// vector<string> plots_to_draw{"CM_Neutron_ExtPdf", "CM_DP_ExtPdf", "CM_Neutron_error", "CM_DP_error"
//           , "ER_DP_ExtPdf", "ER_Neutron_ExtPdf", "NR_Neutron_ExtPdf", "NR_DP_ExtPdf"
//           , "histo_DP", "histo_DP_smooth", "histo_Neutron"
//           , "Neutron_discriminator", "histo_Neutron_NR", "histo_Neutron_NR_smooth"
//           , "Electron_discriminator", "histo_DP_ER", "histo_DP_ER_smooth"
//           };

vector<string> plots_to_draw{"CM_Neutron_ExtPdf", "CM_Neutron_error", "CM_DP_error", "CM_DP_ExtPdf"};

bool save_plots = true;
bool print_fit_result = false;
bool write_results_csv = false;
bool make_root_file = false;

string InputPath("/mnt/raid/users/duylai/warios_bkgs_template/");
string OutputPath("/home/duylai/ardm/fit2D/");

PlotsManager pm(plots_to_draw, OutputPath, save_plots);

const int MinRunNr = 0;
const int MaxRunNr = 106;

const int nbins = 50;
const int total_nbins_2D = Power(nbins+2, 2);
const double TH2_rescale = 8 * Power(nbins, -2);  // rescaling the TH2's amplitude when converting from a RooExtendPdf

// factor of the initial values for the DP template fit
const vector<double> N_ER_Neutron_init{0.5, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
// upper limit for the DP template fit to the Neutron data
const vector<double> F90_high_lim_fit_DP{0.90, 0.62, 0.54, 0.54, 0.54, 0.52, 0.52, 0.50, 0.50, 0.52};
// lower threshold on NDisc
const vector<double> NDisc_low_thresh{0.80, 0.75, 0.50, 0.75, 0.40, 0.50, 0.55, 0.55, 0.55, 0.50};
// lower threshold on F90 for the NR template 
const vector<double> F90_low_thresh_Neutron{0.0, 0.70, 0.56, 0.58, 0.52, 0.50, 0.50, 0.48, 0.46, 0.46};
// lower limit for the NR template fit to the DP data
const vector<double> F90_low_lim_fit_NR_DP{0.80, 0.78, 0.76, 0.76, 0.72, 0.74, 0.70, 0.72, 0.76, 0.76};
// lower thereshold on EDisc
const vector<double> EDisc_low_thresh{0.0, 0.0, 0.0, 0.90, 0.95, 0.97, 0.99, 0.99, 0.99, 0.99};
// upper threshold on F90 for the ER template
const vector<double> F90_high_thresh_DP{1.0, 0.78, 0.72, 0.66, 0.64, 0.62, 0.6, 0.6, 0.6, 0.6};

void Fit2D() {
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  if (max_bin_to_process > max_bin) max_bin_to_process = max_bin;

  // first line of the CM results file
  ofstream results_out;
  if (pm.in_plots_to_draw("CM_DP_ExtPdf") == false || (max_bin_to_process - min_bin_to_process) != max_bin)
            write_results_csv = false;
  if (write_results_csv) {
    results_out.open("neutron/CM_results.csv");
    results_out << "S1L_low,S1L_high,N_ER_DP,N_ER_DP_err,N_NR_DP,N_NR_DP_err,"
                << "N_ER_Neutron,N_ER_Neutron_err,N_NR_Neutron,N_NR_Neutron_err" << endl;
  }

  // variables of interest
  RooRealVar S1LY("S1LY", "S1LY", 0.0);
  RooRealVar F90("S1F90", "S1F90", 0.0, 1.0);
  RooRealVar S2maxLoS1Lcorr("S2maxLoS1Lcorr", "S2maxLoS1Lcorr", -3.0, 5.0);

  for (int bin_nr = min_bin_to_process; bin_nr < max_bin_to_process; bin_nr++) {
    if (!(bin_nr >= min_bin_to_process && bin_nr <= max_bin_to_process)) continue;
    cout << endl << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      S1L slice: " <<
              vec_S1_lower_bound[bin_nr] << " - " << vec_S1_upper_bound[bin_nr] <<
              " p.e.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;

    pm.set_S1L(vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]);

    // ############################ DP TEMPLATE -- DP DATA ############################

    // getting the DP data
    unique_ptr<RooDataSet> Data_DP = GetDataFromWorkspaces(InputPath, "none", MinRunNr, MaxRunNr, bin_nr);

    TH2F* histo_DP = (TH2F*)Data_DP->createHistogram("histo_DP", F90
              , Binning(nbins, 0.0, 1.0), YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    histo_DP->GetZaxis()->SetRangeUser(TH2_rescale, histo_DP->GetMaximum());
    histo_DP->SetTitle(fmt::format("histo_DP, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_DP->GetEntries()).c_str());
    if (pm.control("histo_DP", histo_DP)) continue;

    TH2F* histo_DP_smooth = (TH2F*)histo_DP->Clone("histo_DP_smooth");

    double DP_data_entries = histo_DP_smooth->GetEntries();
    RooRealVar N_ER_DP_RooVar("N_ER_DP", "N_ER_DP", 0.5*DP_data_entries, 1.5*DP_data_entries);
    N_ER_DP_RooVar.setVal(DP_data_entries);

    for (int i = 0; i < 1; ++i) histo_DP_smooth->Smooth(1);

    histo_DP_smooth->GetZaxis()->SetRangeUser(TH2_rescale, histo_DP_smooth->GetMaximum());
    histo_DP_smooth->SetTitle(fmt::format("histo_DP_smooth, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_DP_smooth->GetEntries()).c_str());
    if (pm.control("histo_DP_smooth", histo_DP_smooth)) continue;

    // DP template actually
    RooDataHist ER_template("ER_template", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_DP_smooth));
    RooHistPdf ER_Pdf("ER_Pdf", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), ER_template, 0);
    RooExtendPdf ER_DP_ExtPdf("ER_DP_ExtPdf", "logS2S1_vs_F90", ER_Pdf, N_ER_DP_RooVar);

    TH2F* ER_Pdf_ = (TH2F*)ER_Pdf.createHistogram("ER_Pdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    ER_Pdf_->GetZaxis()->SetRangeUser(ER_Pdf_->GetMinimum(), ER_Pdf_->GetMaximum());
    ER_Pdf_->SetTitle(fmt::format("ER_Pdf, S1L [{:.2f}, {:.2f}]", vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]).c_str());
    if (pm.control("ER_Pdf", ER_Pdf_)) continue;

    // control fit of the DP template
    RooFitResult* ER_DP_fit_result = ER_DP_ExtPdf.fitTo(*Data_DP, Range("FullRange"), Save(), PrintLevel(-5));
    // RooFitResult* ER_DP_fit_result = ER_DP_ExtPdf.fitTo(ER_template, Save(), PrintLevel(-5));
    if (print_fit_result) ER_DP_fit_result->Print();
    // cout << " chi²: " << ER_DP_ExtPdf.createChi2(ER_template, Range(""), Extended(true)
              // , DataError(RooAbsData::Poisson))->getVal() << endl << endl;

    TH2F* fit_ER_DP = (TH2F*)ER_DP_ExtPdf.createHistogram("ER_DP_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    fit_ER_DP->Scale(TH2_rescale);
    fit_ER_DP->GetZaxis()->SetRangeUser(TH2_rescale, fit_ER_DP->GetMaximum());
    fit_ER_DP->SetTitle(fmt::format("ER_DP_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_ER_DP_RooVar.getVal()).c_str());
    if (pm.control("ER_DP_ExtPdf", fit_ER_DP)) continue;

    // ############################ DP TEMPLATE -- NEUTRON DATA ############################

    // getting the Neutron data
    unique_ptr<RooDataSet> Data_Neutron = GetDataFromWorkspaces(InputPath, "neutron", MinRunNr, MaxRunNr, bin_nr);

    TH2F* histo_Neutron = (TH2F*)Data_Neutron->createHistogram("histo_Neutron", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    histo_Neutron->GetZaxis()->SetRangeUser(TH2_rescale, histo_Neutron->GetMaximum());
    histo_Neutron->SetTitle(fmt::format("histo_Neutron, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_Neutron->GetEntries()).c_str());
    if (pm.control("histo_Neutron", histo_Neutron)) continue;

    TH2F* histo_Neutron_tofit = (TH2F*) histo_Neutron->Clone("histo_Neutron_tofit");

    // constraing the fitting region of the DP template to the Neutron data
    int entries_removed = 0;
    int bins_removed = 0;
    for (int ix = 0; ix < nbins+2; ++ix) {
      for (int iy = 0; iy < nbins+2; ++iy) {
        if (fit_ER_DP->GetBinContent(ix, iy) == 0
            || histo_Neutron_tofit->GetXaxis()->GetBinCenter(ix) > F90_high_lim_fit_DP[bin_nr]) {
          entries_removed += histo_Neutron_tofit->GetBinContent(ix, iy);
          bins_removed += 1;
          histo_Neutron_tofit->SetBinContent(ix, iy, 0);
        }
      }
    }

    double Neutron_data_entries = histo_Neutron_tofit->GetEntries() - entries_removed - bins_removed;
    histo_Neutron_tofit->SetEntries(Neutron_data_entries);
    RooRealVar N_ER_Neutron_RooVar("N_ER_Neutron", "N_ER_Neutron", 0.0, 2.0*Neutron_data_entries);
    N_ER_Neutron_RooVar.setVal(Neutron_data_entries*N_ER_Neutron_init[bin_nr]);

    RooDataHist Roohisto_Neutron("Roohisto_Neutron", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_Neutron_tofit));
    RooExtendPdf ER_Neutron_ExtPdf("ER_Neutron_ExtPdf", "logS2S1_vs_F90", ER_Pdf, N_ER_Neutron_RooVar);

    // DP template fit to the Neutron data
    RooFitResult* ER_Neutron_fit_result = ER_Neutron_ExtPdf.fitTo(Roohisto_Neutron, Save(), PrintLevel(-1));
    // RooFitResult* ER_Neutron_fit_result = ER_Neutron_ExtPdf.fitTo(*Data_Neutron, Range("FullRange"), Save(), PrintLevel(-1));
    if (print_fit_result) ER_Neutron_fit_result->Print();
    // cout << " chi²: " << ER_Neutron_ExtPdf.createChi2(Roohisto_Neutron, Range(""), Extended(true)
              // , DataError(RooAbsData::Poisson))->getVal() << endl << endl;

    TH2F* fit_ER_Neutron = (TH2F*)ER_Neutron_ExtPdf.createHistogram("ER_Neutron_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    // fit_ER_Neutron->Scale(Power(N_ER_Neutron_RooVar.getVal(), -1));
    fit_ER_Neutron->Scale(TH2_rescale);
    fit_ER_Neutron->GetZaxis()->SetRangeUser(TH2_rescale, fit_ER_Neutron->GetMaximum());
    fit_ER_Neutron->SetTitle(fmt::format("ER_Neutron_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_ER_Neutron_RooVar.getVal()).c_str());
    if (pm.control("ER_Neutron_ExtPdf", fit_ER_Neutron)) continue;

    // ############################ DP TEMPLATE SUBTRACTION ############################

    // extracted NR events to be put in this TH2F
    TH2F* histo_Neutron_NR = new TH2F("histo_Neutron_NR", Data_Neutron->GetTitle()
              , nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    histo_Neutron_NR->GetXaxis()->SetTitle(F90.GetTitle());
    histo_Neutron_NR->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    TH2F* Neutron_discriminator = new TH2F("Neutron_discriminator"
              , fmt::format("Neutron_discriminator, S1L [{:.2f}, {:.2f}]"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]).c_str(), nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    Neutron_discriminator->GetXaxis()->SetTitle(F90.GetTitle());
    Neutron_discriminator->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    double N_NR_Neutron_estimated = histo_Neutron->GetEntries()-N_ER_Neutron_RooVar.getVal();
    double N_NR_Neutron_normalized = N_NR_Neutron_estimated / histo_Neutron->GetEntries();

    double fit_ER_Neutron_integral = 0;
    double Neutron_NR_subtracted_entries = 0;
    double fit_DP_template_Neutron_Chi2 = 0.0;
    double fit_DP_template_Neutron_dof = -1.0;  // 1 fit parameter
    double fit_DP_template_Neutron_Chi2_o_dof = 0.0;
    double fit_DP_template_Neutron_p_value = 0.0;
    
    // filling histo_Neutron_NR and Neutron_descriminator
    for (int ix = 0; ix < nbins+2; ++ ix) {
      for (int iy = 0; iy < nbins+2; ++ iy) {
          // cout << histo_Neutron->GetXaxis()->GetBinLowEdge(ix) << endl;
        double data_bin_content = histo_Neutron->GetBinContent(ix, iy);
        double fit_bin_content = fit_ER_Neutron->GetBinContent(ix, iy);
        fit_ER_Neutron_integral += fit_bin_content;
        if (data_bin_content != 0) {
          double fit_error = data_bin_content - fit_bin_content;
          double fit_relative_error = fit_error / data_bin_content;   // small for ERs, large for NRs
          if (fit_ER_DP->GetBinContent(ix, iy) != 0
            && histo_Neutron_tofit->GetXaxis()->GetBinCenter(ix) <= F90_high_lim_fit_DP[bin_nr]) {
            fit_DP_template_Neutron_Chi2 += fit_error * fit_error / fit_bin_content;
            fit_DP_template_Neutron_dof += 1;
          }
          Neutron_discriminator->SetBinContent(ix, iy, fit_relative_error);
          // if (fit_relative_error > 1.0 - N_NR_Neutron_normalized && fit_relative_error <= 1.0) {
          // if (fit_relative_error >= 0.55 && fit_relative_error <= 1.0
          if (fit_relative_error >= NDisc_low_thresh[bin_nr] && fit_relative_error <= 1.0 
          && histo_Neutron->GetXaxis()->GetBinCenter(ix) >= F90_low_thresh_Neutron[bin_nr]) {
            double bin_val = fit_error*Power(fit_relative_error, 0);
            histo_Neutron_NR->SetBinContent(ix, iy, bin_val);
                      // , Erf(histo_Neutron_NR->GetXaxis()->GetBinCenter(ix)));
                      // , TanH(4.0*histo_Neutron_NR->GetXaxis()->GetBinCenter(ix))*fit_error);
            Neutron_NR_subtracted_entries += bin_val;
          }
        }
      }
    }
    // cout << "Fit integral: " <<  fit_ER_Neutron_integral << endl;
    fit_DP_template_Neutron_Chi2_o_dof = fit_DP_template_Neutron_Chi2 / fit_DP_template_Neutron_dof;
    fit_DP_template_Neutron_p_value = ROOT::Math::chisquared_cdf_c(fit_DP_template_Neutron_Chi2, fit_DP_template_Neutron_dof);
    cout << "DP template fit to Neutron data, chi2: " << fit_DP_template_Neutron_Chi2 << endl;
    cout << "DP template fit to Neutron data, dof: " << fit_DP_template_Neutron_dof << endl;
    cout << "DP template fit to Neutron data, chi2_o_dof: " << fit_DP_template_Neutron_Chi2_o_dof << endl;
    cout << "DP template fit to Neutron data, p_value: " << fit_DP_template_Neutron_p_value << endl;

    double Neutron_ER_subtracted_entries = histo_Neutron->GetEntries() - Neutron_NR_subtracted_entries;

    Neutron_discriminator->GetZaxis()->SetRangeUser(-1.0, 1.0);

    histo_Neutron_NR->SetEntries(Neutron_NR_subtracted_entries);
    histo_Neutron_NR->GetZaxis()->SetRangeUser(TH2_rescale, histo_Neutron_NR->GetMaximum());
    histo_Neutron_NR->SetTitle(fmt::format("histo_Neutron_NR, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_Neutron_NR->GetEntries()).c_str());

    // overriding C++ logical short-circuit, requires to "purge" manually
    bool Neutron_discriminator_fin = pm.control("Neutron_discriminator", Neutron_discriminator, false, false,
            fmt::format("chi2: {:.2f}, dof: {:.2f}, chi2_o_dof: {:.2f}, p_value: {:.2f}", fit_DP_template_Neutron_Chi2,
            fit_DP_template_Neutron_dof, fit_DP_template_Neutron_Chi2_o_dof, fit_DP_template_Neutron_p_value));
    bool histo_Neutron_NR_subtracted_fin = pm.control("histo_Neutron_NR", histo_Neutron_NR, true, false);

    if (Neutron_discriminator_fin || histo_Neutron_NR_subtracted_fin) {
      pm.purge();
      continue;
    }

    // smoothen extracted NR events
    TH2F* histo_Neutron_NR_smooth = (TH2F*)histo_Neutron_NR->Clone("histo_Neutron_NR_smooth");
    histo_Neutron_NR_smooth->Smooth(1);
    histo_Neutron_NR_smooth->SetTitle(fmt::format("histo_Neutron_NR_smooth, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_Neutron_NR_smooth->GetEntries()).c_str());
    histo_Neutron_NR_smooth->GetZaxis()->SetRangeUser(TH2_rescale, histo_Neutron_NR_smooth->GetMaximum());
    if (pm.control("histo_Neutron_NR_smooth", histo_Neutron_NR_smooth)) continue;

    // ############################ NR MODEL - SUBTRACTED NEUTRON DATA ############################

    RooRealVar N_NR_Neutron_RooVar("N_NR_Neutron", "N_NR_Neutron", 0.5*Neutron_NR_subtracted_entries, 1.5*Neutron_NR_subtracted_entries);
    N_NR_Neutron_RooVar.setVal(Neutron_NR_subtracted_entries);

    RooDataHist NR_template("NR_template", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr)
              , Import(*histo_Neutron_NR_smooth));
    RooHistPdf NR_Pdf("ER_Pdf", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), NR_template, 0);
    RooExtendPdf NR_Neutron_ExtPdf("NR_Neutron_ExtPdf", "logS2S1_vs_F90", NR_Pdf, N_NR_Neutron_RooVar);

    // control fit of the NR template
    RooFitResult* NR_Neutron_fit_result = NR_Neutron_ExtPdf.fitTo(NR_template, SumW2Error(true), Save(), PrintLevel(-5));
    if (print_fit_result) NR_Neutron_fit_result->Print();

    TH2F* fit_NR_Neutron = (TH2F*)NR_Neutron_ExtPdf.createHistogram("NR_Neutron_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    fit_NR_Neutron->Scale(TH2_rescale);
    fit_NR_Neutron->GetZaxis()->SetRangeUser(TH2_rescale, fit_NR_Neutron->GetMaximum());
    fit_NR_Neutron->SetTitle(fmt::format("NR_Neutron_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_NR_Neutron_RooVar.getVal()).c_str());
    if (pm.control("NR_Neutron_ExtPdf", fit_NR_Neutron)) continue;

    double fit_NR_Neutron_integral = 0;
    for (int i = 0; i < total_nbins_2D; ++i) {
      fit_NR_Neutron_integral += fit_NR_Neutron->GetBinContent(i);
    }
    // cout << fit_NR_Neutron_integral << endl;

    // ############################ NR MODEL - DP DATA ############################

    TH2F* histo_DP_tofit = (TH2F*)histo_DP_smooth->Clone("histo_DP_tofit");
    // TH2F* histo_DP_tofit = (TH2F*)histo_DP->Clone("histo_DP_tofit");

    //  constraing fitting region of the NR template to the DP data
    entries_removed = 0;
    bins_removed = 0;
    for (int ix = 0; ix < total_nbins_2D; ++ix) {
      for (int iy = 0; iy < total_nbins_2D; ++iy) {
        if (fit_NR_Neutron->GetBinContent(ix, iy) == 0
                  || histo_DP_tofit->GetXaxis()->GetBinCenter(ix) < F90_low_lim_fit_NR_DP[bin_nr]) {
          entries_removed += histo_DP_tofit->GetBinContent(ix, iy);
          bins_removed += 1;
          histo_DP_tofit->SetBinContent(ix, iy, 0);
        }
      }
    }

    double NR_instr_entries = histo_DP_tofit->GetEntries() - entries_removed - bins_removed;
    histo_DP_tofit->SetEntries(NR_instr_entries);
    RooRealVar N_NR_DP_RooVar("N_NR_DP_RooVar", "N_NR_DP_RooVar", 0.0, 2.0*NR_instr_entries);
    N_NR_DP_RooVar.setVal(NR_instr_entries);

    RooDataHist Roohisto_DP("Roohisto_DP", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_DP_tofit));
    RooExtendPdf NR_DP_ExtPdf("NR_DP_ExtPdf", "logS2S1_vs_F90", NR_Pdf, N_NR_DP_RooVar);

    // NR template fit to the DP data
    RooFitResult* NR_DP_fit_result = NR_DP_ExtPdf.fitTo(Roohisto_DP, SumW2Error(true), Save(), PrintLevel(-1));
    if (print_fit_result) NR_DP_fit_result->Print();

    TH2F* fit_NR_DP = (TH2F*)NR_DP_ExtPdf.createHistogram("NR_DP_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    fit_NR_DP->Scale(TH2_rescale);
    fit_NR_DP->GetZaxis()->SetRangeUser(TH2_rescale, fit_NR_DP->GetMaximum());
    fit_NR_DP->SetTitle(fmt::format("NR_DP_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_NR_DP_RooVar.getVal()).c_str());
    if (pm.control("NR_DP_ExtPdf", fit_NR_DP)) continue;

    // ############################ ER MODEL - DP DATA ############################

    // extracted ER events to be put in the TH2F
    TH2F* histo_DP_ER = new TH2F("histo_DP_ER", Data_DP->GetTitle()
              , nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    histo_DP_ER->GetXaxis()->SetTitle(F90.GetTitle());
    histo_DP_ER->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    TH2F* Electron_discriminator = new TH2F("Electron_discriminator"
              , fmt::format("Electron_discriminator, S1L [{:.2f}, {:.2f}]"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]).c_str(), nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    Electron_discriminator->GetXaxis()->SetTitle(F90.GetTitle());
    Electron_discriminator->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    double DP_ER_subtracted_entries = histo_DP->GetEntries();
    double fit_NR_template_DP_Chi2 = 0.0;
    double fit_NR_template_DP_dof = -1.0;  // 1 fit parameter
    double fit_NR_template_DP_Chi2_o_dof = 0.0;
    double fit_NR_template_DP_p_value = 0.0;
    double fit_NR_template_DP_EDisc_sum = 0.0;

    // filling histo_DP_ER and Electron_discriminator
    for (int ix = 0; ix < nbins+2; ++ ix) {
      for (int iy = 0; iy < nbins+2; ++ iy) {
        double data_bin_content = histo_DP_smooth->GetBinContent(ix, iy);
        // double data_bin_content = histo_DP->GetBinContent(ix, iy);
        double fit_bin_content = fit_NR_DP->GetBinContent(ix, iy);
        if (data_bin_content != 0) {
          double fit_error = data_bin_content - fit_bin_content;
          // double fit_error = Max(data_bin_content - fit_bin_content, 0.001);
          double fit_relative_error = fit_error / data_bin_content;
          if (fit_NR_DP->GetBinContent(ix, iy) != 0
            && histo_DP_tofit->GetXaxis()->GetBinCenter(ix) >= F90_low_lim_fit_NR_DP[bin_nr]) {
            fit_NR_template_DP_Chi2 += fit_error * fit_error / fit_bin_content;
            fit_NR_template_DP_dof += 1;
            fit_NR_template_DP_EDisc_sum += Abs(fit_relative_error);
          }
          Electron_discriminator->SetBinContent(ix, iy, fit_relative_error);
          if (fit_relative_error > EDisc_low_thresh[bin_nr] && fit_relative_error <= 1.0
                    && histo_DP_smooth->GetXaxis()->GetBinCenter(ix) < F90_high_thresh_DP[bin_nr]) {
                    // && histo_DP->GetXaxis()->GetBinCenter(ix) < F90_high_thresh_DP[bin_nr]) {
            histo_DP_ER->SetBinContent(ix, iy, fit_error);
            DP_ER_subtracted_entries -= fit_error;
          }
        }
      }
    }
    fit_NR_template_DP_Chi2_o_dof = fit_NR_template_DP_Chi2 / fit_NR_template_DP_dof;
    fit_NR_template_DP_p_value = ROOT::Math::chisquared_cdf_c(fit_NR_template_DP_Chi2, fit_NR_template_DP_dof);
    cout << "NR template fit to DP data, chi2: " << fit_NR_template_DP_Chi2 << endl;
    cout << "NR template fit to DP data, dof: " << fit_NR_template_DP_dof << endl;
    cout << "NR template fit to DP data, chi2_o_dof: " << fit_NR_template_DP_Chi2_o_dof << endl;
    cout << "NR template fit to DP data, p_value: " << fit_NR_template_DP_p_value << endl;
    cout << "NR template fit to DP data, Avg. EDisc: " << fit_NR_template_DP_EDisc_sum / (fit_NR_template_DP_dof + 1) << endl;

    histo_DP_ER->SetEntries(DP_ER_subtracted_entries);
    histo_DP_ER->GetZaxis()->SetRangeUser(TH2_rescale, histo_DP_ER->GetMaximum());
    histo_DP_ER->SetTitle(fmt::format("histo_DP_ER, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], DP_ER_subtracted_entries).c_str());
    // if (pm.control("histo_DP_ER", histo_DP_ER)) continue;

    // cout <<"Edisc min: " << Electron_discriminator->GetMinimum() << endl;
    // Electron_discriminator->GetZaxis()->SetRangeUser(0.9, 1.0);
    Electron_discriminator->GetZaxis()->SetRangeUser(-1.0, 1.0);

    // overriding C++ logical short-circuit, requires to "purge" manually
    bool Electron_discriminator_fin = pm.control("Electron_discriminator", Electron_discriminator, false, false,
            fmt::format("chi2: {:.2f}, dof: {:.2f}, chi2_o_dof: {:.2f}, p_value: {:.2f}", fit_NR_template_DP_Chi2,
            fit_NR_template_DP_dof, fit_NR_template_DP_Chi2_o_dof, fit_NR_template_DP_p_value));
    bool histo_DP_ER_subtracted_fin = pm.control("histo_DP_ER", histo_DP_ER, true, false);

    if (Electron_discriminator_fin || histo_DP_ER_subtracted_fin) {
      pm.purge();
      continue;
    }

    // "smoothened" ER template, actually we took the smoothened DP data so don't smooth it again
    TH2F* histo_DP_ER_smooth = (TH2F*)histo_DP_ER->Clone("histo_DP_ER_smooth");
    // histo_DP_ER_smooth->Smooth(1);
    histo_DP_ER_smooth->SetTitle(fmt::format("histo_DP_ER_smooth, S1L [{:.2f}, {:.2f}], entries {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], histo_DP_ER_smooth->GetEntries()).c_str());
    histo_DP_ER_smooth->GetZaxis()->SetRangeUser(TH2_rescale, histo_DP_ER_smooth->GetMaximum());
    if (pm.control("histo_DP_ER_smooth", histo_DP_ER_smooth)) continue;

    // ############################ COMPLETE MODEL - NEUTRON DATA ############################

    // the CM should fit to everywhere, the zero bins result in fitting problems of the CM
    // so we need to set them to arbitrarily small values
    for (int i = 0; i < total_nbins_2D; ++i) {
      if (histo_DP_ER_smooth->GetBinContent(i) == 0)
                histo_DP_ER_smooth->SetBinContent(i,1e-29);
      if (histo_Neutron_NR_smooth->GetBinContent(i) == 0)
                histo_Neutron_NR_smooth->SetBinContent(i,1e-29);
    }

    RooDataHist ER_template_new("ER_template", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr)
              , Import(*histo_DP_ER_smooth));
    RooHistPdf ER_Pdf_new("ER_Pdf", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), ER_template_new, 0);
    RooExtendPdf ER_DP_ExtPdf_new("ER_DP_ExtPdf", "logS2S1_vs_F90", ER_Pdf_new, N_ER_DP_RooVar);

    RooDataHist NR_template_new = RooDataHist("NR_template", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr)
              , Import(*histo_Neutron_NR_smooth));
    RooHistPdf NR_Pdf_new("NR_Pdf", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), NR_template_new, 0);
    RooExtendPdf NR_Neutron_ExtPdf_new("NR_Neutron_ExtPdf", "logS2S1_vs_F90", NR_Pdf_new, N_NR_Neutron_RooVar);

    // combining the ER and NR templates
    RooAddPdf Complete_ExtPdf("Complete_ExtPdf", "Complete_ExtPdf", RooArgList(ER_DP_ExtPdf_new, NR_Neutron_ExtPdf_new));

    TH2F* CM_template = (TH2F*)Complete_ExtPdf.createHistogram("CM_template", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    CM_template->Scale(TH2_rescale);
    CM_template->GetZaxis()->SetRangeUser(TH2_rescale, CM_template->GetMaximum());

    if (pm.control("CM_template", CM_template)) continue;

    // histo_Neutron_tofit = (TH2F*) histo_Neutron->Clone("histo_Neutron_tofit");

    // for (int ix = 0; ix < nbins+2; ++ix) {
    //   for (int iy = 0; iy < nbins+2; ++iy) {
    //     if (CM_template->GetBinContent(ix, iy) == 0) {
    //       histo_Neutron_tofit->SetBinContent(ix, iy, 0);
    //     }
    //   }
    // }

    N_ER_DP_RooVar.setRange(0.5*Neutron_ER_subtracted_entries, 1.5*Neutron_ER_subtracted_entries);
    N_ER_DP_RooVar.setVal(Neutron_ER_subtracted_entries);
    N_NR_Neutron_RooVar.setRange(0.5*Neutron_NR_subtracted_entries, 1.5*Neutron_NR_subtracted_entries);
    N_NR_Neutron_RooVar.setVal(Neutron_NR_subtracted_entries);

    // RooDataHist Neutron_Data("Neutron_Data", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_Neutron_tofit));
    RooDataHist Neutron_Data("Neutron_Data", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_Neutron));

    // CM fit to the Neutron data
    RooFitResult* CM_Neutron_fit_result = Complete_ExtPdf.fitTo(Neutron_Data, Range("FullRange"), SumW2Error(true), Save(), PrintLevel(-5));
    if (print_fit_result) CM_Neutron_fit_result->Print();

    double N_ER_Neutron = N_ER_DP_RooVar.getVal();
    double N_NR_Neutron = N_NR_Neutron_RooVar.getVal();

    double N_ER_Neutron_err = N_ER_DP_RooVar.getError();
    double N_NR_Neutron_err = N_NR_Neutron_RooVar.getError();

    TH2F* fit_CM_Neutron = (TH2F*)Complete_ExtPdf.createHistogram("CM_Neutron_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    fit_CM_Neutron->Scale(TH2_rescale);
    fit_CM_Neutron->GetZaxis()->SetRangeUser(TH2_rescale, fit_CM_Neutron->GetMaximum());
    fit_CM_Neutron->SetTitle(fmt::format("CM_Neutron_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N_ER = {:.2e}, N_NR = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_ER_DP_RooVar.getVal(), N_NR_Neutron_RooVar.getVal()).c_str());
    if (pm.control("CM_Neutron_ExtPdf", fit_CM_Neutron)) continue;

    // relative errors of the CM fit to the Neutron data
    TH2F* CM_Neutron_error = new TH2F("CM_Neutron_error", Data_DP->GetTitle()
              , nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    CM_Neutron_error->GetXaxis()->SetTitle(F90.GetTitle());
    CM_Neutron_error->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    double fit_CM_template_Neutron_Chi2 = 0.0;
    double fit_CM_template_Neutron_dof = -1.0;  // 1 fit parameter
    double fit_CM_template_Neutron_Chi2_o_dof = 0.0;
    double fit_CM_template_Neutron_p_value = 0.0;
    double fit_CM_template_Neutron_EDisc_sum = 0.0;

    for (int ix = 0; ix < nbins + 2; ++ix) {
      for (int iy = 0; iy < nbins + 2; ++iy) {
        double data_bin_content = histo_Neutron->GetBinContent(ix, iy);
        double fit_bin_content = fit_CM_Neutron->GetBinContent(ix, iy);
        if (data_bin_content != 0) {
          double fit_error = data_bin_content - fit_bin_content;
          double fit_relative_error = fit_error / data_bin_content;
          CM_Neutron_error->SetBinContent(ix, iy, fit_relative_error);
          fit_CM_template_Neutron_Chi2 += fit_error * fit_error / fit_bin_content;
          fit_CM_template_Neutron_dof += 1;
          fit_CM_template_Neutron_EDisc_sum += Abs(fit_relative_error);
        }
      }
    }
    fit_CM_template_Neutron_Chi2_o_dof = fit_CM_template_Neutron_Chi2 / fit_CM_template_Neutron_dof;
    fit_CM_template_Neutron_p_value = ROOT::Math::chisquared_cdf_c(fit_CM_template_Neutron_Chi2, fit_CM_template_Neutron_dof);
    cout << "CM template fit to Neutron data, chi2: " << fit_CM_template_Neutron_Chi2 << endl;
    cout << "CM template fit to Neutron data, dof: " << fit_CM_template_Neutron_dof << endl;
    cout << "CM template fit to Neutron data, chi2_o_dof: " << fit_CM_template_Neutron_Chi2_o_dof << endl;
    cout << "CM template fit to Neutron data, p_value: " << fit_CM_template_Neutron_p_value << endl;
    cout << "CM template fit to Neutron data, Avg. Error: " << fit_CM_template_Neutron_EDisc_sum / (fit_CM_template_Neutron_dof + 1) << endl;

    CM_Neutron_error->GetZaxis()->SetRangeUser(-1.0, 1.0);
    CM_Neutron_error->SetTitle(fmt::format("CM_Neutron_error, S1L [{:.2f}, {:.2f}]"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]).c_str());
    if (pm.control("CM_Neutron_error", CM_Neutron_error, false, false,
            fmt::format("chi2: {:.2e}, dof: {:.2f}, chi2_o_dof: {:.2e}, p_value: {:.2f}", fit_CM_template_Neutron_Chi2,
            fit_CM_template_Neutron_dof, fit_CM_template_Neutron_Chi2_o_dof, fit_CM_template_Neutron_p_value))) continue;

    // ############################ COMPLETE MODEL - DP DATA ############################

    // histo_DP_tofit = (TH2F*) histo_DP->Clone("histo_DP_tofit");

    // for (int ix = 0; ix < nbins+2; ++ix) {
    //   for (int iy = 0; iy < nbins+2; ++iy) {
    //     if (CM_template->GetBinContent(ix, iy) == 0) {
    //       histo_DP_tofit->SetBinContent(ix, iy, 0);
    //     }
    //   }
    // }

    N_ER_DP_RooVar.setRange(0.9*histo_DP->GetEntries(), 1.1*histo_DP->GetEntries());
    N_ER_DP_RooVar.setVal(0.8*histo_DP->GetEntries());
    N_NR_Neutron_RooVar.setRange(0, 1e4);
    N_NR_Neutron_RooVar.setVal(10);

    // RooDataHist DP_Data("DP_Data", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_DP_tofit));
    // RooDataHist DP_Data("DP_Data", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_DP_smooth));
    RooDataHist DP_Data("DP_Data", "logS2S1_vs_F90", RooArgSet(F90, S2maxLoS1Lcorr), Import(*histo_DP));

    // CM fit to the DP data
    RooFitResult* CM_DP_fit_result = Complete_ExtPdf.fitTo(DP_Data, Range("FullRange"), SumW2Error(true), Save(), PrintLevel(-5));
    if (print_fit_result) CM_DP_fit_result->Print();

    double N_ER_DP = N_ER_DP_RooVar.getVal();
    double N_NR_DP = N_NR_Neutron_RooVar.getVal();  // instrumental neutrons in DP data

    double N_ER_DP_err = N_ER_DP_RooVar.getError();
    double N_NR_DP_err = N_NR_Neutron_RooVar.getError();

    // fitting results of the CM in the current S1 bin
    string one_line = fmt::format("{:.1f},{:.1f},{},{},{},{},{},{},{},{}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]
              , N_ER_DP, N_ER_DP_err, N_NR_DP, N_NR_DP_err
              , N_ER_Neutron, N_ER_Neutron_err, N_NR_Neutron, N_NR_Neutron_err);
    
    if (write_results_csv) results_out << one_line << endl;
    else cout << endl << "------ CM results: ------" << endl << one_line << endl << endl;

    // relative errors of the CM fit to the DP data
    TH2F* fit_CM_DP = (TH2F*)Complete_ExtPdf.createHistogram("CM_DP_ExtPdf", F90, Binning(nbins, 0.0, 1.0)
              , YVar(S2maxLoS1Lcorr, Binning(nbins, -3.0, 5.0)));
    fit_CM_DP->Scale(TH2_rescale);
    fit_CM_DP->GetZaxis()->SetRangeUser(TH2_rescale, fit_CM_DP->GetMaximum());
    fit_CM_DP->SetTitle(fmt::format("CM_DP_ExtPdf, S1L [{:.2f}, {:.2f}], fit param. N_ER = {:.2e}, N_NR = {:.2e}"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr], N_ER_DP_RooVar.getVal(), N_NR_Neutron_RooVar.getVal()).c_str());
    if (pm.control("CM_DP_ExtPdf", fit_CM_DP)) continue;

    TH2F* CM_DP_error = new TH2F("CM_DP_error", Data_DP->GetTitle()
              , nbins, 0.0, 1.0, nbins, -3.0, 5.0);
    CM_DP_error->GetXaxis()->SetTitle(F90.GetTitle());
    CM_DP_error->GetYaxis()->SetTitle(S2maxLoS1Lcorr.GetTitle());

    double fit_CM_template_DP_Chi2 = 0.0;
    double fit_CM_template_DP_dof = -1.0;  // 1 fit parameter
    double fit_CM_template_DP_Chi2_o_dof = 0.0;
    double fit_CM_template_DP_p_value = 0.0;
    double fit_CM_template_DP_EDisc_sum = 0.0;

    for (int ix = 0; ix < nbins + 2; ++ix) {
      for (int iy = 0; iy < nbins + 2; ++iy) {
        double data_bin_content = histo_DP->GetBinContent(ix, iy);
        // double data_bin_content = histo_DP_smooth->GetBinContent(ix, iy);
        double fit_bin_content = fit_CM_DP->GetBinContent(ix, iy);
        if (data_bin_content != 0) {
          double fit_error = data_bin_content - fit_bin_content;
          double fit_relative_error = fit_error / data_bin_content;
          CM_DP_error->SetBinContent(ix, iy, fit_relative_error);
          fit_CM_template_DP_Chi2 += fit_error * fit_error / fit_bin_content;
          fit_CM_template_DP_dof += 1;
          fit_CM_template_DP_EDisc_sum += Abs(fit_relative_error);
        }
      }
    }
    fit_CM_template_DP_Chi2_o_dof = fit_CM_template_DP_Chi2 / fit_CM_template_DP_dof;
    fit_CM_template_DP_p_value = ROOT::Math::chisquared_cdf_c(fit_CM_template_DP_Chi2, fit_CM_template_DP_dof);
    cout << "CM template fit to Neutron data, chi2: " << fit_CM_template_DP_Chi2 << endl;
    cout << "CM template fit to Neutron data, dof: " << fit_CM_template_DP_dof << endl;
    cout << "CM template fit to Neutron data, chi2_o_dof: " << fit_CM_template_DP_Chi2_o_dof << endl;
    cout << "CM template fit to Neutron data, p_value: " << fit_CM_template_DP_p_value << endl;
    cout << "CM template fit to Neutron data, Avg. Error: " << fit_CM_template_DP_EDisc_sum / (fit_CM_template_DP_dof + 1) << endl;

    CM_DP_error->GetZaxis()->SetRangeUser(-1.0, 1.0);
    CM_DP_error->SetTitle(fmt::format("CM_DP_error, S1L [{:.2f}, {:.2f}]"
              , vec_S1_lower_bound[bin_nr], vec_S1_upper_bound[bin_nr]).c_str());
    if (pm.control("CM_DP_error", CM_DP_error, false, true,
            fmt::format("chi2: {:.2e}, dof: {:.2f}, chi2_o_dof: {:.2e}, p_value: {:.2f}", fit_CM_template_DP_Chi2,
            fit_CM_template_DP_dof, fit_CM_template_DP_Chi2_o_dof, fit_CM_template_DP_p_value))) continue;
    
    // the following lines are not executed unless there are "fake" plots in plots_to_draw,
    // these plots will be removed in check_plots_to_draw(), then the current loop repeats
    pm.check_plots_to_draw();
    --bin_nr;
  }
  if (make_root_file) pm.make_root_file();
}
