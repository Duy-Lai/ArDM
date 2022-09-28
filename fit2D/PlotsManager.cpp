#include "PlotsManager.h"
using namespace std;

void PlotsManager::check_plots_to_draw() {
  std::cout << std::endl << "Nothing has been drawn. Checking plots_to_draw:" << std::endl;
  int k = 0;
  bool kill;
  while (k < nbr_plots_to_draw) {
    kill = true;
    for (int i = 0; i < all_plots.size(); ++i) {
      if (plots_to_draw[k] == all_plots[i]) {
        kill = false;
        break;
      }
    }
    if (kill) {
      std::cout << "- " << plots_to_draw[k] << " removed." << std::endl;
      plots_to_draw.erase(plots_to_draw.begin() + k);
    }
    else ++k;
    nbr_plots_to_draw = plots_to_draw.size();
  }
  purge();
}

bool PlotsManager::control(std::string plot, TH1* histo, bool log_scale = true, bool purge_after_draw = true,
          std::string annotation = "") {
  all_plots.push_back(plot);
  all_histos.push_back(histo);
  vec_log_scale.push_back(log_scale);
  all_annots.push_back(annotation);
  if (in_plots_to_draw(plot)) pd_match_count++;
  if (pd_match_count == nbr_plots_to_draw) {
    draw_all();
    if (purge_after_draw) purge();
    return true;
  }
  return false;
}

void PlotsManager::draw_all() {
  for (int i = 0; i < canvases.size(); ++i) canvases[i]->Close();
  canvases.clear();
  for (int i = 0; i < all_plots.size(); ++i) {
    if (in_plots_to_draw(all_plots[i])) {
      canvases.push_back(new TCanvas(all_plots[i].c_str(), all_plots[i].c_str(), x_start + 3 * x_step, y_start + 2 * y_step
                , x_size, y_size));
      int j = canvases.size() - 1;
      if (vec_log_scale[i] && all_log_scale) canvases[j]->SetLogz();
      else gStyle->SetPalette(kLightTemperature);
      canvases[j]->SetGridx();
      canvases[j]->SetGridy();
      canvases[j]->SetTopMargin(0.08);
      canvases[j]->SetLeftMargin(0.09);
      canvases[j]->SetRightMargin(0.13);
      // gStyle->SetPalette(kRedBlue);
      // TColor::InvertPalette();
      all_histos[i]->SetStats(0);
      // all_histos[i]->GetXaxis()->SetTitleSize(0.05);
      // all_histos[i]->GetYaxis()->SetTitleSize(0.05);
      // all_histos[i]->GetXaxis()->SetLabelSize(0.05);
      // all_histos[i]->GetYaxis()->SetLabelSize(0.05);
      // all_histos[i]->GetZaxis()->SetLabelSize(0.05);
      histo_list->Add(all_histos[i]->Clone(all_histos[i]->GetTitle()));
      if (save_plots) {
        all_histos[i]->Draw("colsz");
        if (all_annots[i] != "") {
          TText *ann = new TText(0.01, 4.5, all_annots[i].c_str());
          // TText *ann = new TText(0.45, 0.76, "S1");
          ann->SetTextSize(0.05);
          ann->SetTextColor(kBlack);
          ann->Draw("Same");
        }
        canvases[j]->SaveAs(fmt::format("{}{}/{:.0f}_{:.0f}.png", OutputPath, all_plots[i], S1_low, S1_high).c_str());
      }
      else all_histos[i]->DrawClone("colsz");
      gStyle->SetPalette(kBird);
    }
  }
}

bool PlotsManager::in_plots_to_draw(std::string plot) {
  for (int i = 0; i < nbr_plots_to_draw; ++i) {
    if (plot == plots_to_draw[i]) return true;
  }
  return false;
}

void PlotsManager::make_root_file() {
  TFile* newfile = new TFile((OutputPath + "histos.root").c_str(), "recreate");
  histo_list->Write();
  newfile->Close();
}

void PlotsManager::purge() {
  for (int i = 0; i < all_histos.size(); ++i) delete all_histos[i];
  all_histos.clear();
  all_plots.clear();
  all_annots.clear();
  pd_match_count = 0;
}

void PlotsManager::set_S1L(double low, double high) {
  S1_low = low;
  S1_high = high;
}
