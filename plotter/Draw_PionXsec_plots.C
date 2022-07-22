#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"

Double_t langaufun(Double_t *x, Double_t *par) {
  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298;
  Double_t np = 500.0;
  Double_t sc = 5.0;// convolution extends to +-sc Gaussian sigmas 
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;

  mpc=par[1];
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow)/np;

  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status, TString FunName)
{
  Int_t i;
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");

  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  TFitResultPtr fitres = his->Fit(FunName,"RBOSQ"); // fit within specified range, use ParLimits, do not plot
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }

  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  Status[0] = fitres->CovMatrixStatus();

  return (ffit);              // return fit function                        
}

void Draw_MC_vs_FakeData(TString filename, TString histname, TString TitleX, double xmin, double xmax, double rebin){

  

}

void Draw_MC_vs_Data(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  double mc_scale = 1.;
  //mc_scale = 1.5440743 * 2.;

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    //maphist[this_hist_name] -> Scale(mc_scale);
    maphist[this_hist_name] -> Rebin(rebin);
  }
  double N_mc_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  TH1D *hist_data = (TH1D*)gDirectory -> Get("htrack_" + histname + "_0") -> Clone();
  hist_data -> Rebin(rebin);
  double data_max = hist_data -> GetMaximum();
  double N_data_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);

  mc_scale = 2.0 * N_data_isSelectedPart / N_mc_isSelectedPart;
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    maphist[this_hist_name] -> Scale(mc_scale);
  }


  TString title_y = "Events";
  TString nameofhistogram = histname + "Draw_MC_vs_Data" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);

  ////////////////////////////////////
  // == Pad 1
  ////////////////////////////////////
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.05 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();
  
  TH1D *pad1_template = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  pad1_template -> SetStats(0);
  pad1_template -> GetXaxis() -> SetTitle(nameofhistogram);
  pad1_template -> GetXaxis() -> SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitleSize(0);
  pad1_template -> GetYaxis() -> SetLabelSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleSize(0.07);
  pad1_template -> GetYaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetTitle("Events");
  pad1_template -> GetYaxis() -> SetRangeUser(0.1, data_max * 10000.); // == logy
  //pad1_template -> GetYaxis() -> SetRangeUser(0., data_max * 1.5);
  pad1_template -> Draw("hist");

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);
 
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  TH1D * mc_sum = (TH1D*)maphist["htrack_" + histname + "_1"] -> Clone();
  for(int i = 1; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    TString this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
    TString this_legend_str = pi_type_str[i] + " " + this_N_event;
    if(i != 1) mc_sum -> Add(maphist[this_hist_name]);
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetFillColor(colour_array[i]);
    maphstack[hstack] -> Add(maphist[this_hist_name]);
    maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "f");
  }
  TString mc_sum_N_event = Form("%.1f", mc_sum -> Integral());
  mc_sum -> SetLineColor(kWhite);
  maplegend[legend]->AddEntry(mc_sum, "MC Sum " + mc_sum_N_event, "l");
  double mc_max = mc_sum -> GetMaximum(); 
  if(mc_max > data_max){
    pad1_template -> GetYaxis() -> SetRangeUser(0.1, mc_max * 10000.); // == logy
    //pad1_template -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
    pad1_template -> Draw("hist");
  }

  hist_data -> SetLineColor(kBlack);
  hist_data -> SetMarkerColor(kBlack);
  hist_data -> SetMarkerStyle(20);
  hist_data -> SetMarkerSize(1);
  hist_data -> SetLineWidth(2);
  TString data_N_event = Form("%.1f", hist_data -> Integral());
  maplegend[legend]->AddEntry(hist_data, "Observed " + data_N_event, "lp");
  maphstack[hstack] -> Draw("histsame");
  hist_data -> Draw("epsame");
  maplegend[legend] -> Draw("same");
  gPad->RedrawAxis();

  ////////////////////////////////////
  // == Pad 2
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.05 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();

  TH1D * pad2_template = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  pad2_template -> Draw("hist");
  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(TitleX);
  pad2_template -> GetXaxis() -> SetTitleSize(0.15);
  pad2_template -> GetXaxis() -> SetLabelSize(0.125);
  pad2_template -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.4);
  pad2_template -> GetYaxis() -> SetLabelSize(0.09);
  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.0, 3.0);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");

  TH1D * data_mc_ratio = (TH1D*)hist_data -> Clone();
  data_mc_ratio -> Divide(mc_sum);
  data_mc_ratio -> Draw("psame");

  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.6, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(data_mc_ratio, "Obs./Pred.", "lp");

  TLine *pad2_line = new TLine(xmin, 1, xmax, 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

  maplegend["bottom" + legend] -> Draw("same");

  gPad->RedrawAxis();

  ////////////////////////////////////
  // == Latex
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.63, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_vs_Data_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
  f_data -> Close();
}

void Draw_2D_MC_and_Data(TString filename, TString histname, TString TitleX, TString TitleY, TString beam_P, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TH2D* mc_sum = (TH2D*)gDirectory -> Get("htrack_" + histname + "_0") -> Clone();
  for(int i = 1; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    TH2D* this_hist = (TH2D*)gDirectory -> Get(this_hist_name);
    mc_sum -> Add(this_hist);
  }
  double N_mc_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  TH2D *hist_data = (TH2D*)gDirectory -> Get("htrack_" + histname + "_0") -> Clone();
  double N_data_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);

  double mc_scale = N_data_isSelectedPart / N_mc_isSelectedPart;
  mc_sum -> Scale(mc_scale);

  mc_sum -> RebinX(rebin_x);
  mc_sum -> RebinY(rebin_y);
  hist_data -> RebinX(rebin_x);
  hist_data -> RebinY(rebin_y);

  double data_max = hist_data -> GetMaximum();
  double mc_max = mc_sum -> GetMaximum();
  double z_max = max(data_max, mc_max);

  TCanvas *c = new TCanvas("", "", 800, 800);
  c->SetRightMargin(0.15);
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleY);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");
  
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.56, 0.96, "Run 1, " + beam_P + " GeV/c Beam");

  mc_sum -> Draw("colzsame");
  //gPad->RedrawAxis();

  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Draw_2D_MC_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  template_h -> Draw("colz");
  latex_ArgoNeuT.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.56, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  hist_data -> Draw("colzsame");
  //gPad->RedrawAxis();
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Draw_2D_Data_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  f_mc -> Close();
  f_data -> Close();
  
}

void Draw_MC_shape_comparison(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double mc_max = -1.;
  for(int i = 0; i < N_pi_type; i++){
    if(i == 1 || i == 2 || i == 3 || i == 5 || i == 7){
      TString this_hist_name = Form("htrack_" + histname + "_%d", i);
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      maphist[this_hist_name] -> Scale(1. / maphist[this_hist_name] -> Integral());
      double this_max = maphist[this_hist_name] -> GetMaximum();
      if(mc_max < this_max) mc_max = this_max;
    }
  }

  TString nameofhistogram = histname + "Draw_MC_shape_comparison" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> GetYaxis() -> SetRangeUser(0., 0.4);
  template_h -> Draw();
  
  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  for(int i = 1; i < N_pi_type; i++){
    if(i == 1 || i == 2 || i == 3 || i == 5 || i == 7){
      TString this_hist_name = Form("htrack_" + histname + "_%d", i);
      TString this_legend_str = pi_type_str[i];
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetLineStyle(i);
      maphist[this_hist_name] -> SetLineWidth(4);
      maphist[this_hist_name] -> Draw("histsame");
      maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "l");
    }
  }
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();
  
  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_MC_daughter_shape_comparison(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double mc_max = -1.;
  const int N_PDGs = 3;
  TString PDGs_str[N_PDGs] = {"proton", "pion", "other"};
  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    maphist[this_hist_name] -> Rebin(rebin);
    maphist[this_hist_name] -> Scale(1. / maphist[this_hist_name] -> Integral());
    double this_max = maphist[this_hist_name] -> GetMaximum();
    if(mc_max < this_max) mc_max = this_max;
  }

  TString nameofhistogram = histname + "Draw_MC_shape_comparison" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 600, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    TString this_legend_str = PDGs_str[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(i + 1);
    maphist[this_hist_name] -> SetLineWidth(4);
    maphist[this_hist_name] -> Draw("histsame");
    maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();

  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_daughter_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_PionXsec_plots(){

  setTDRStyle();
  TString file_suffix = "_noBeamXY.root";
  // ==== Data vs MC
  // == After beam cuts
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "BeamP", "P_{Beam Inst.} (MeV/c)", "2.0", 1000., 3000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "BeamKE", "KE_{Beam Inst.} (MeV)", "2.0", 1000., 3000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "PandoraSlice", "Pass PandoraSlice", "2.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "CaloSize", "Pass CaloSize", "2.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dx", "dx_{Beam}^{Calo Start} / #sigma_{x}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dy", "dy_{Beam}^{Calo Start} / #sigma_{y}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dz", "dz_{Beam}^{Calo Start} / #sigma_{z}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dxy", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "2.0", 0., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_costh", "cos#theta_{Beam Track}", "2.0", -1., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "daughter_michel_score", "Daughter Michel Score", "2.0", 0., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "2.0", 0., 100., 100.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "Beam_alt_length", "Track Length (cm)", "2.0", 0., 600., 10.);

 /*
  Draw_MC_vs_Data("_PionXsec_2.0GeV.root", "BeamP", "P_{Beam Inst.} (MeV/c)", "2.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV.root", "BeamKE", "KE_{Beam Inst.} (MeV)", "2.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV.root", "Beam_alt_length", "Track Length (cm)", "2.0", 0., 600., 10.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV.root", "daughter_michel_score", "Michel Score", "2.0", 0., 1., 1.);
  */
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "BeamP", "P_{Beam Inst.} (MeV/c)", "1.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "BeamKE", "KE_{Beam Inst.} (MeV)", "1.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "PandoraSlice", "Pass PandoraSlice", "1.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "CaloSize", "Pass CaloSize", "1.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dx", "dx_{Beam}^{Calo Start} / #sigma_{x}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dy", "dy_{Beam}^{Calo Start} / #sigma_{y}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dz", "dz_{Beam}^{Calo Start} / #sigma_{z}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dxy", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "1.0", 0., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_costh", "cos#theta_{Beam Track}", "1.0", -1., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "daughter_michel_score", "Daughter Michel Score", "1.0", 0., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "1.0", 0., 100., 100.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "Beam_alt_length", "Track Length (cm)", "1.0", 0., 600., 10.);

  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "BeamP", "P_{Beam Inst.} (MeV/c)", "0.5", 0., 1000., 20.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "BeamKE", "KE_{Beam Inst.} (MeV)", "0.5", 0., 1000., 20.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "PandoraSlice", "Pass PandoraSlice", "0.5", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "CaloSize", "Pass CaloSize", "0.5", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dx", "dx_{Beam}^{Calo Start} / #sigma_{x}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dy", "dy_{Beam}^{Calo Start} / #sigma_{y}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dz", "dz_{Beam}^{Calo Start} / #sigma_{z}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dxy", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "0.5", 0., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_costh", "cos#theta_{Beam Track}", "0.5", -1., 1., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "daughter_michel_score", "Daughter Michel Score", "0.5", 0., 1., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "0.5", 0., 100., 100.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "Beam_alt_length", "Track Length (cm)", "0.5", 0., 600., 10.);

  //Draw_MC_vs_Data("_PionXsec_0.3GeV.root", "BeamP", "P_{Beam Inst.} (MeV/c)", "0.3", 0., 1000., 20.);

  // == Before beam cuts
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "BeamP_precut", "P_{Beam Inst.} (MeV/c)", "2.0", 1000., 3000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "BeamKE_precut", "KE_{Beam Inst.} (MeV)", "2.0", 1000., 3000., 20.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "PandoraSlice_precut", "Pass PandoraSlice", "2.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "CaloSize_precut", "Pass CaloSize", "2.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dx_precut", "dx_{Beam}^{Calo Start} / #sigma_{x}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dy_precut", "dy_{Beam}^{Calo Start} / #sigma_{y}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dz_precut", "dz_{Beam}^{Calo Start} / #sigma_{z}", "2.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_dxy_precut", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "2.0", 0., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "beam_costh_precut", "cos#theta_{Beam Track}", "2.0", -1., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "daughter_michel_score_precut", "Daughter Michel Score", "2.0", 0., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_2.0GeV" + file_suffix, "chi2_proton_precut", "#chi^{2}_{proton}", "2.0", 0., 100., 100.);

  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "BeamP_precut", "P_{Beam Inst.} (MeV/c)", "1.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "BeamKE_precut", "KE_{Beam Inst.} (MeV)", "1.0", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "PandoraSlice_precut", "Pass PandoraSlice", "1.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "CaloSize_precut", "Pass CaloSize", "1.0", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dx_precut", "dx_{Beam}^{Calo Start} / #sigma_{x}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dy_precut", "dy_{Beam}^{Calo Start} / #sigma_{y}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dz_precut", "dz_{Beam}^{Calo Start} / #sigma_{z}", "1.0", -10., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_dxy_precut", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "1.0", 0., 10., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "beam_costh_precut", "cos#theta_{Beam Track}", "1.0", -1., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "daughter_michel_score_precut", "Daughter Michel Score", "1.0", 0., 1., 1.);
  Draw_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "chi2_proton_precut", "#chi^{2}_{proton}", "1.0", 0., 100., 100.);

  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "BeamP_precut", "P_{Beam Inst.} (MeV/c)", "0.5", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "BeamKE_precut", "KE_{Beam Inst.} (MeV)", "0.5", 0., 2000., 20.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "PandoraSlice_precut", "Pass PandoraSlice", "0.5", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "CaloSize_precut", "Pass CaloSize", "0.5", 0., 2., 1.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dx_precut", "dx_{Beam}^{Calo Start} / #sigma_{x}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dy_precut", "dy_{Beam}^{Calo Start} / #sigma_{y}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dz_precut", "dz_{Beam}^{Calo Start} / #sigma_{z}", "0.5", -10., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_dxy_precut", "dxy_{Beam}^{Calo Start} / #sigma_{xy}", "0.5", 0., 10., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "beam_costh_precut", "cos#theta_{Beam Track}", "0.5", -1., 1., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "daughter_michel_score_precut", "Daughter Michel Score", "0.5", 0., 1., 10.);
  Draw_MC_vs_Data("_PionXsec_0.5GeV" + file_suffix, "chi2_proton_precut", "#chi^{2}_{proton}", "0.5", 0., 100., 100.);

  // ==== Draw 2D plots
  Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "beam_inst_XY_precut", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "1.0", -50, 10, 5., 400., 440., 5.);
  Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "beam_inst_XY", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "1.0", -50, 10, 5., 400., 440., 5.);

  Draw_2D_MC_and_Data("_PionXsec_0.5GeV" + file_suffix, "beam_inst_XY_precut", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "0.5", -50, 10, 20., 400., 440., 20.);
  Draw_2D_MC_and_Data("_PionXsec_0.5GeV" + file_suffix, "beam_inst_XY", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "0.5", -50, 10, 20., 400., 440., 20.);

  // ==== MC shape comparison
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "2.0", -100., 100., 5.);
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "2.0", -100., 100., 5.);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "1.0", -100., 100., 5.);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "1.0", -100., 100., 5.);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "0.5", -100., 100., 5.);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "0.5", -100., 100., 5.);

  // ==== MC daughter shape comparison
  Draw_MC_daughter_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "trackScore", "Track Score", "2.0", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "emScore", "e/#gamma Score", "2.0", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "2.0", 0., 100., 100.);
  Draw_MC_daughter_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "trackScore", "Track Score", "1.0", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "emScore", "e/#gamma Score", "1.0", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "1.0", 0., 100., 100.);
  Draw_MC_daughter_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "trackScore", "Track Score", "0.5", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "emScore", "e/#gamma Score", "0.5", 0., 1., 5.);
  Draw_MC_daughter_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "chi2_proton", "#chi^{2}_{proton}", "0.5", 0., 100., 100.);
}
