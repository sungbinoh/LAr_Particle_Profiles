#ifndef canvas_margins_h
#define canvas_margins_h

#include "TStyle.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "TH2D.h"

using namespace std;

class canvas_margins{

 public:
  void fixOverlay();
  void setTDRStyle();
  void canvas_margin(TCanvas *c1);
  void canvas_margin_limit(TCanvas *c1);
  void canvas_margin(TCanvas *c1, TPad *c1_up, TPad *c1_down);
  void canvas_margin_twoplots(TCanvas *c1, TPad *c1_up, TPad *c1_down);
  void hist_axis(TH1D *hist);
  void hist_axis_limit(TH1D *hist);
  void hist_axis(THStack *hist);
  void hist_axis(TGraph *hist);
  void hist_axis(TGraphAsymmErrors *hist);
  void hist_axis(TH2D *hist);
  void hist_axis(TH2F *hist);
  void hist_axis(THStack *hist, TH1D *hist_compare);
  void hist_axis(TH1D *hist, TH1D *hist_compare);
  void hist_axis_twoplots(TH1D *hist, TH1D *hist_compare);
  void hist_axis(TGraph *hist, TGraph *hist_compare);

  canvas_margins();
  virtual ~canvas_margins();
};

#endif
