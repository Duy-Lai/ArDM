// managing necessary plots to produce, with lazy-evaluation approach
// inspired by SidebandFit1D.cpp of Alex Stauffer
// written by Duy Lai, 09.2022

#ifndef PLOTSMANAGER_H
#define PLOTSMANAGER_H

#include<iostream>
#include<vector>
#include<string>
#include<sys/stat.h>

#define FMT_HEADER_ONLY
#include </home/duylai/cget/include/fmt/core.h>
// 1: "pip install cget"; 2: "cget install fmtlib/fmt"

class PlotsManager {
  std::vector<std::string> plots_to_draw;
  std::string OutputPath;
  bool all_log_scale;
  bool save_plots;
  double S1_low;
  double S1_high;

  // names of the plots
  std::vector<std::string> all_plots;
  // the plots themselves
  std::vector<TH1*> all_histos;
  std::vector<bool> vec_log_scale;
  std::vector<std::string> all_annots;
  // canvases on which the plots are drawn
  std::vector<TCanvas*> canvases;
  // list of the plots to put in histos.root, which are to be reproduced by plot_fit2D.py
  TList* histo_list;

  // current number of plots of interest (i.e. to draw)
  int pd_match_count;
  // total number of plots to draw
  int nbr_plots_to_draw;
  void draw_all();

  int XOffset = 10;
  int YOffset = 35;
  int x_start = 80 + XOffset;
  int y_start = 0 + YOffset;
  float x_size = 2100;
  float y_size = 1500;
  float x_step = 0;
  float y_step = 0;

  public:
  PlotsManager(std::vector<std::string> pd, std::string op, bool sp, bool ls = true)
  : OutputPath(op), save_plots(sp), all_log_scale(ls) {
    histo_list = new TList;
    pd_match_count = 0;
    // only adds the plot once, in case plots_to_draw contains two or more identical plot names
    std::cout << "Plots to draw:" << std::endl;
    for (int i = 0; i < pd.size(); ++i) {
      if (in_plots_to_draw(pd[i]) == false) {
        std::cout << "- " << pd[i] << std::endl;
        plots_to_draw.push_back(pd[i]);
        ++nbr_plots_to_draw;
        // creating plot folders within OutputPath if they don't exist yet
        struct stat buffer;
        if (stat(pd[i].c_str(), &buffer) != 0) mkdir((OutputPath + pd[i]).c_str(), 0777);
      }
    }
  }
  // check if plots_to_draw contains bs names, e.g. "Human_desciminator", and delete them in plots_to_draw
  // to be put at the end of Fit2D, if check_plots_to_draw() is provoked, the current S1 bin is treated again
  // after removing the bs plot names
  void check_plots_to_draw();
  // control(..) adds a plot to collection, then checks if all the plots to draw are added or not
  // will return a true if it is the case so that Fit2D.cpp continues to a new S1 bins without having going to
  // the end, e.g. stops at Neutron_descriminator if we don't want/need the Complete Model
  bool control(std::string, TH1*, bool, bool, std::string);
  // check if a plot (its name) is in plots_to_draw or not
  bool in_plots_to_draw(std::string);
  // make a .root file with histo_list, can be used by plot_fit2D.py
  void make_root_file();
  // remove all plots, their names, canvases, etc. when an S1 bin is done
  void purge();
  // set S1 lower and upper bounds
  void set_S1L(double, double);
};

#endif