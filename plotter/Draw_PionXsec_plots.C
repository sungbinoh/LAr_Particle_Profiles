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

void Draw_MC_vs_FakeData(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool logy){

  TString suffix = "hdaughter_reco_";
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double max_fake = -1.;
  gDirectory -> Cd(dir);
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form(suffix + histname + "_%d", i);
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      if(i == 0){
	max_fake = maphist[this_hist_name] -> GetMaximum();
      }
    }
    else maphist[this_hist_name] = nullptr;
  }

  TString title_y = "Events";
  TString nameofhistogram = histname + "Draw_MC_vs_FakeData" + beam_P;
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
  if(logy) mappad[pad1] -> SetLogy();

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
  pad1_template -> GetYaxis() -> SetRangeUser(0., max_fake * 1.8);
  if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, max_fake * 10000.);
  pad1_template -> Draw("hist");

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  TH1D * mc_sum = (TH1D*)maphist[suffix + histname + "_1"] -> Clone();
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form(suffix + histname + "_%d", i);
    //cout << "[Draw_MC_vs_FakeData] this_hist_name : " << this_hist_name << endl;
    if(maphist[this_hist_name] == nullptr) continue;
    TString this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
    TString this_legend_str = pi_type_str[i] + " " + this_N_event;
    Int_t this_color = colour_array[i];
    if(logy){
      this_hist_name = Form(suffix + histname + "_%d", N_pi_type - i - 1);
      this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
      this_legend_str = pi_type_str[N_pi_type - i - 1] + " " + this_N_event;
      this_color = colour_array[N_pi_type - i - 1];
    }
    if(this_legend_str.Contains("Fake")){
      maphist[this_hist_name] -> SetLineColor(kBlack);
      maphist[this_hist_name] -> SetMarkerColor(kBlack);
      maphist[this_hist_name] -> SetMarkerStyle(20);
      maphist[this_hist_name] -> SetMarkerSize(1);
      maphist[this_hist_name] -> SetLineWidth(2);
      maplegend[legend]->AddEntry(maphist[this_hist_name], "Fake Data", "lp");
    }
    else{
      if(!this_legend_str.Contains("Inel")) mc_sum -> Add(maphist[this_hist_name]);
      maphist[this_hist_name] -> SetLineColor(this_color);
      maphist[this_hist_name] -> SetFillColor(this_color);
      maphstack[hstack] -> Add(maphist[this_hist_name]);
      maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "f");
    }
  }
  TString this_fakedata_name = Form(suffix + histname + "_%d", 0);
  maphstack[hstack] -> Draw("histsame");
  maphist[this_fakedata_name] -> Draw("epsame");
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
  pad2_template -> GetYaxis() -> SetTitle("#frac{Fake Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.4);
  pad2_template -> GetYaxis() -> SetLabelSize(0.09);
  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.0, 3.0);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");

  TH1D * data_mc_ratio = (TH1D*)maphist[this_fakedata_name] -> Clone();
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_vs_FakeData/MC_vs_FakeData_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
  mapcanvas[canvas] -> Close();
}

void Draw_MC_vs_Data(TString filename, TString dir, TString suffix, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool logy){
  
  double mc_scale = 1.;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double N_mc_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);
  gDirectory -> Cd(dir);
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form(suffix + histname + "_%d", i);
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  double N_data_isSelectedPart = ((TH1D*)gDirectory -> Get("Cutflow")) -> GetBinContent(1);
  gDirectory -> cd(dir);
  TH1D *hist_data = (TH1D*)gDirectory -> Get(suffix + histname + "_0") -> Clone();
  hist_data -> Rebin(rebin);
  double data_max = hist_data -> GetMaximum();

  mc_scale = 2.0 * N_data_isSelectedPart / N_mc_isSelectedPart;
  //mc_scale = 2.0 * 119545.00 / 156701.00;
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form(suffix + histname + "_%d", i);
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> Scale(mc_scale);
    }
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
  if(logy) mappad[pad1] -> SetLogy();
  
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
  //pad1_template -> GetYaxis() -> SetRangeUser(0., 1500.);
  pad1_template -> GetYaxis() -> SetRangeUser(0., data_max * 1.5);
  if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, data_max * 10000.); // == logy
  pad1_template -> Draw("hist");

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);
 
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  TH1D * mc_sum = (TH1D*)maphist[suffix + histname + "_1"] -> Clone();
  for(int i = 1; i < N_pi_type; i++){
    TString this_hist_name = Form(suffix + histname + "_%d", i);
    if(maphist[this_hist_name] != nullptr){    
      TString this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
      TString this_legend_str = pi_type_str[i] + " " + this_N_event;
      if(i != 1) mc_sum -> Add(maphist[this_hist_name]);
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetFillColor(colour_array[i]);
      maphstack[hstack] -> Add(maphist[this_hist_name]);
      maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "f");
    }
  }
  TString mc_sum_N_event = Form("%.1f", mc_sum -> Integral());
  mc_sum -> SetLineColor(kWhite);
  maplegend[legend]->AddEntry(mc_sum, "MC Sum " + mc_sum_N_event, "l");
  double mc_max = mc_sum -> GetMaximum(); 
  if(mc_max > data_max){
    pad1_template -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
    if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, mc_max * 10000.); // == logy
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_vs_Data/MC_vs_Data_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
  f_data -> Close();
  mapcanvas[canvas] -> Close();
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/2D/Draw_2D_MC_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  template_h -> Draw("colz");
  latex_ArgoNeuT.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.56, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  hist_data -> Draw("colzsame");
  //gPad->RedrawAxis();
  pdfname = WORKING_DIR + "/output/plots/PionXsec/2D/Draw_2D_Data_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  f_mc -> Close();
  f_data -> Close();
  
}

void Draw_2D_MC(TString filename, TString dir, TString histname, TString TitleX, TString TitleY, TString beam_P, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y, bool logz){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  TString this_hist_name = histname;
  TH2D *this_hist = (TH2D*)gDirectory -> Get(this_hist_name) -> Clone();
  
  this_hist -> RebinX(rebin_x);
  this_hist -> RebinY(rebin_y);
  //this_hist -> Scale(1. / this_hist -> Integral() );
  double z_max = this_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 800);
  c->SetRightMargin(0.15);
  if(logz) c->SetLogz();
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
  template_h -> GetZaxis() -> SetRangeUser(0.0, z_max * 1.1 );
  if(logz) template_h -> GetZaxis() -> SetRangeUser(1.0, 1000. );
  template_h -> Draw("colz");

  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.56, 0.96, "Run 1, " + beam_P + " GeV/c Beam");

  this_hist -> Draw("colzsame");
  gPad->RedrawAxis();
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/2D/Draw_2D_MC_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  f_mc -> Close();
}

void Draw_MC_shape_comparison(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, double ymax){

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
      int this_N_bin = maphist[this_hist_name] -> GetNbinsX();
      double last_content = maphist[this_hist_name] -> GetBinContent(this_N_bin);
      double overflow = maphist[this_hist_name] -> GetBinContent(this_N_bin + 1);
      maphist[this_hist_name] -> SetBinContent(this_N_bin, last_content + overflow);

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
  template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> Draw();
  
  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_comparison/MC_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_MC_daughter_shape_comparison(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool logy){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;
  const int N_PDGs = 4;
  //TString PDGs_str[N_PDGs] = {"proton", "pion", "muon", "other"};
  TString PDGs_str[N_PDGs] = {"muon", "other", "proton", "pion"};
  //Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  Int_t colour_array[] = {867, 416, 632, 800};

  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    //cout << "[Draw_MC_daughter_shape_comparison] " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
	maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
	maphist[this_hist_name] -> Rebin(rebin);
	int this_N_bin = maphist[this_hist_name] -> GetNbinsX();
	double overflow = maphist[this_hist_name] -> GetBinContent(this_N_bin + 1);
	double integral= maphist[this_hist_name] -> Integral();
	maphist[this_hist_name] -> Scale(1. / (overflow + integral) );
	double last_content = maphist[this_hist_name] -> GetBinContent(this_N_bin);
	overflow = maphist[this_hist_name] -> GetBinContent(this_N_bin + 1);
	maphist[this_hist_name] -> SetBinContent(this_N_bin, last_content + overflow);
	double this_max = maphist[this_hist_name] -> GetMaximum();
	if(mc_max < this_max) mc_max = this_max;
    }
    else maphist[this_hist_name] = nullptr;
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
  if(logy) mapcanvas[canvas] -> SetLogy();

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  //template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> GetYaxis() -> SetRangeUser(0., 0.08);
  if(logy) template_h -> GetYaxis() -> SetRangeUser(0.0001, 1000.);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    TString this_legend_str = PDGs_str[i];
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      //maphist[this_hist_name] -> SetLineStyle(i + 1);
      maphist[this_hist_name] -> SetLineStyle(1);
      if(histname.Contains("chi2") || histname.Contains("nHits")) maphist[this_hist_name] -> SetLineStyle(1);
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_MC_daughter_cutflow(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool draw_eff){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;

  const int N_cuts = 8;
  TString cutflow_str[N_cuts] = {"nocut", "cos", "dist", "Nhits", "emScore", "TrackScore", "chi2_proton", "mean_dEdx"};
  TString cutflow_lgd_str[N_cuts] = {"No cut", "cos#theta_{Beam, Daughter}", "Distance(beam end, daughter start)", "N_{Hits}", "e#gamma Score", "Track Score", "#chi^{2}_{proton}", "<dE/dx>_{Truncated}"};

  if(histname.Contains("pion")){
    cutflow_str[6] = cutflow_str[6] + "_pionID";
    cutflow_str[7] = cutflow_str[7] + "_pionID";
  }
  else if(histname.Contains("proton")){
    //N_cuts = 7;
    cutflow_str[6] = cutflow_str[6] + "_protonID";
  }
  
  for(int i = 0; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    //cout << "[Draw_MC_daughter_cutflow] this_hist_name : " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }

  // == Dividing by nocut hit
  for(int i = 1; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    TString nocut_hist_name = "hdaughter_" + histname + "_" + cutflow_str[0] + "_" + dir;
    if(draw_eff && maphist[this_hist_name] != nullptr){
      //cout << "[Draw_MC_daughter_cutflow] Dividing " << this_hist_name << endl;
      maphist[this_hist_name] -> Divide(maphist[nocut_hist_name]);
    }
  }
  TString nocut_hist_name = "hdaughter_" + histname + "_" + cutflow_str[0] + "_" + dir;
  if(draw_eff) maphist[nocut_hist_name] -> Divide(maphist[nocut_hist_name]);

  mc_max = maphist["hdaughter_" + histname + "_" + cutflow_str[1] + "_" + dir] -> GetMaximum();
  TString nameofhistogram = histname + "Draw_MC_daughter_cutflow" + beam_P + dir;
  if(draw_eff) nameofhistogram = nameofhistogram + "true";
  else nameofhistogram = nameofhistogram + "false";
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
  if(draw_eff) template_h -> GetYaxis() -> SetTitle("Eff.");
  else template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    TString this_legend_str = cutflow_lgd_str[i];
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      //maphist[this_hist_name] -> SetLineStyle(i + 1);
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
  if(draw_eff) pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_Cutflow_" + histname + "_" + beam_P + "GeV_eff.pdf";
  else pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_Cutflow_" + histname + "_" + beam_P + "GeV_hist.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_MC_daughter_purity(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;

  const int N_particles = 5;
  TString particle_str[N_particles] = {"all", "pion", "proton", "muon", "other"};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + particle_str[i] + "_" + dir;
    //cout << "[Draw_MC_daughter_purity] this_hist_name : " << this_hist_name << endl;
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    maphist[this_hist_name] -> Rebin(rebin);
  }

  for(int i = 1; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + particle_str[i] + "_" + dir;
    TString nocut_hist_name = "hdaughter_" + histname + "_" + particle_str[0] + "_" + dir;
    //cout << "[Draw_MC_daughter_purity] Dividing " << this_hist_name << endl;
    maphist[this_hist_name] -> Divide(maphist[nocut_hist_name]);
  }
  TString nocut_hist_name = "hdaughter_" + histname + "_" + particle_str[0] + "_" + dir;
  maphist[nocut_hist_name] -> Divide(maphist[nocut_hist_name]);

  TString nameofhistogram = histname + "Draw_MC_daughter_cutflow" + beam_P;
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
  template_h -> GetYaxis() -> SetTitle("Fraction");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 1.5);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + particle_str[i] + "_" + dir;
    TString this_legend_str = particle_str[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    //maphist[this_hist_name] -> SetLineStyle(i + 1);                                                                                                                                                    
    maphist[this_hist_name] -> SetLineWidth(4);
    maphist[this_hist_name] -> Draw("histsame");
    maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }
  maphist[nocut_hist_name] -> Draw("histsame");

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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_Purity_" + histname + "_" + beam_P + "GeV_overlap.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Draw_MC_daughter_efficiency(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;

  const int N_particles = 4;
  TString particle_str[N_particles] = {"pion", "proton", "muon", "other"};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + particle_str[i] + "_" + histname + "_" + dir;
    TString nocut_hist_name = "hdaughter_" + particle_str[i] + "_TrueP_nocut_" + dir;

    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    TH1D * this_nocut_hist = (TH1D*)gDirectory -> Get(nocut_hist_name) -> Clone();
    
    maphist[this_hist_name] -> Rebin(rebin);
    this_nocut_hist -> Rebin(rebin);
    maphist[this_hist_name] -> Divide(this_nocut_hist);
  }

  TString nameofhistogram = histname + "Draw_MC_daughter_efficiency" + beam_P;
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
  template_h -> GetYaxis() -> SetTitle("Efficiency");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 1.5);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + particle_str[i] + "_" + histname + "_" + dir;
    TString this_legend_str = particle_str[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_Efficiency_" + histname + "_" + beam_P + "GeV_overlap.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Stack_MC_daughter(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;
  const int N_PDGs = 4;
  TString PDGs_str[N_PDGs] = {"proton", "pion", "muon", "other"};
  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    //cout << "[Stack_MC_daughter] " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }

  TString nameofhistogram = histname + "Stack_MC_daughter" + beam_P + dir;
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
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(4000);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  TH1D *mc_sum = (TH1D*)maphist["hdaughter_" + PDGs_str[1] + "_" + histname] -> Clone();
  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    TString this_legend_str = PDGs_str[i];
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetFillColor(colour_array[i]);
      maphstack[hstack] -> Add(maphist[this_hist_name]);
      maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "l");
      if(i != 1) mc_sum -> Add(maphist[this_hist_name]);
    }
  }
  mc_max = mc_sum -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> Draw();

  maphstack[hstack] -> Draw("histsame");
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
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/MC_daughter_stack_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();
}

void Run_Draw_MC_vs_FakeData(TString file_prefix, TString file_suffix, TString P_str, TString dir){

  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "QE_N_pion_michel_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "QE_N_pion_gaus_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "QE_N_pion_likelihood_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);

  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 10., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);
 
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "gaus_fit_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 10., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "gaus_fit_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "gaus_fit_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "likelihood_fit_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 10., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "likelihood_fit_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "likelihood_fit_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_NC_40_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);

  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_NC_40_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "stop_pion_EQEmE_NC_4_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "likelihood_fit_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_FakeData("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "gaus_fit_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);

}

void Run_Draw_MC_vs_Data(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  if(P_str == "0.5") Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "BeamP_" + dir, "P_{Beam Inst.} (MeV/c)", P_str, 0., 1500., 20., true);
  if(P_str == "1.0") Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "BeamP_" + dir, "P_{Beam Inst.} (MeV/c)", P_str, 0., 2000., 20., true);
  if(P_str == "2.0") Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "BeamP_" + dir, "P_{Beam Inst.} (MeV/c)", P_str, 1000., 3000., 20., true);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "PandoraSlice_" + dir, "Pass PandoraSlice", P_str, 0., 2., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "CaloSize_" + dir, "Pass CaloSize", P_str, 0., 2., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_dx_" + dir, "dx_{Beam}^{Calo Start} / #sigma_{x}", P_str, -10., 10., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_dy_" + dir, "dy_{Beam}^{Calo Start} / #sigma_{y}", P_str, -10., 10., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_dz_" + dir, "dz_{Beam}^{Calo Start} / #sigma_{z}", P_str, -10., 10., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_dxy_" + dir, "dxy_{Beam}^{Calo Start} / #sigma_{xy}", P_str, 0., 10., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_costh_" + dir, "cos#theta_{Beam Track}", P_str, -1., 1., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "daughter_michel_score_" + dir, "Daughter Michel Score", P_str, 0., 1., 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "chi2_proton_" + dir, "#chi^{2}_{proton}", P_str, 0., 100., 100., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "htrack_", "beam_alt_len_" + dir, "Beam Track Length (cm)", P_str, 0., 600., 10., true);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_michel_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_gaus_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_likelihood_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "stop_pion_EQEmE_NC_4_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "stop_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "likelihood_fit_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "gaus_fit_pion_EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "stop_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 20., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "stop_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "stop_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "gaus_fit_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 20., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "gaus_fit_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "gaus_fit_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "likelihood_fit_pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 20., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "likelihood_fit_pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "likelihood_fit_pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_michel_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_gaus_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_likelihood_fitted_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);


  /*
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "QE_N_pion_michel_" + dir, "N_{#pi^{+}}(Michel > 0.5)", P_str, -0.5, 9.5, 1., true);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_KE_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 10., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_cos_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_length_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_KE_KE400_" + dir, "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, 0., 1000., 10., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_cos_KE400_" + dir, "cos_{Beam, #pi^{#pm}}", P_str, -1., 1., 5., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "pion_length_KE400_" + dir, "L_{Daughter, #pi^{#pm}} (cm)", P_str, 0., 300., 10., false);

  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 20., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_NC_40_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_NC_40_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  Draw_MC_vs_Data("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_", "EQEmE_NC_4_KE400_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", P_str, -3000., 3000., 40., false);
  */
}

void Run_Draw_2D_MC(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_pion_TrueP_vs_proton_chi2_" + dir, "P_{True}^{Start} (MeV)", "#chi^{2}_{proton}", P_str, 0., 1400., 20., 0., 500., 10., true);
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_proton_TrueP_vs_proton_chi2_" + dir, "P_{True}^{Start} (MeV)", "#chi^{2}_{proton}", P_str, 0., 1400., 20., 0., 500., 10., true);
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_pion_TrueP_vs_proton_chi2_offset_" + dir, "P_{True}^{Start} (MeV)", "#chi^{2}_{proton}", P_str, 0., 1400., 20., 0., 500., 10., true);
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_proton_TrueP_vs_proton_chi2_offset_" + dir, "P_{True}^{Start} (MeV)", "#chi^{2}_{proton}", P_str, 0., 1400., 20., 0., 500., 10., true);

  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_EQEmE_NC_4_vs_KE_daughter_" + dir + "_1", "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, -1000., 500., 40., 0., 1000., 20., false);
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_EQEmE_NC_4_vs_KE_daughter_" + dir + "_2", "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, -1000., 500., 40., 0., 1000., 20., false);
  Draw_2D_MC("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "hdaughter_reco_EQEmE_NC_4_vs_KE_daughter_" + dir + "_3", "E_{QE, #pi^{#pm}} - E_{Beam,Reco} (MeV)", "KE_{Daughter, #pi^{#pm}}(MeV)", P_str, -1000., 500., 40., 0., 1000., 20., false);

}

void Run_Draw_MC_daughter_shape_comparison(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_not_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_dist_" + dir, "Distance(beam end, daughter start) (cm)", P_str, 0., 100., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "trackScore_" + dir, "Track Score", P_str, 0., 1., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "emScore_" + dir, "e/#gamma Score", P_str, 0., 1., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "chi2_proton_" + dir, "#chi^{2}_{proton}", P_str, 0., 100., 100., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "PFP_nHits_" + dir, "N_{Hits}", P_str, 0., 1000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "mean_dEdx_" + dir, "<dE/dx>_{Truncated}^{Last 10 cm} (MeV/cm)", P_str, 0., 20., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_proton_chi2_" + dir, "#chi^{2}_{proton}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_pion_chi2_" + dir, "#chi^{2}_{pion}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_muon_chi2_" + dir, "#chi^{2}_{muon}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_proton_chi2_offset_" + dir, "#chi^{2}_{proton}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_pion_chi2_offset_" + dir, "#chi^{2}_{pion}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_muon_chi2_offset_" + dir, "#chi^{2}_{muon}", P_str, 0., 200., 50., false);

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "fitted_chi2_proton_" + dir, "#chi^{2}_{proton, fit}", P_str, 0., 10., 100., false);

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "cos_theta_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 1., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "cos_theta_not_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 1., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "dist_beam_end_" + dir, "Distance(beam end, daughter start) (cm)", P_str, 0., 200., 5., false);

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, 0., 2000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_NC_40_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, 0., 2000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_NC_4_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, 0., 2000., 10., false);

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -3000., 3000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_NC_40_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -3000., 3000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -3000., 3000., 10., false);

}

void Run_Draw_MC_daughter_cutflow(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "pion_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50., false);
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "proton_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1400., 50., false);

  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "pion_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50., true);
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "proton_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1400., 50., true);

  Draw_MC_daughter_efficiency("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "TrueP_chi2_proton_protonID", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50.);
  Draw_MC_daughter_efficiency("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "TrueP_mean_dEdx_pionID", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50.);

}

void Run_Draw_MC_daughter_purity(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Draw_MC_daughter_purity("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "pion_purity_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50.);
  Draw_MC_daughter_purity("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "proton_purity_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1400., 50.);
}

void Run_Stack_MC_daughter(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, -2000., 2000., 10.);
  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_NC_40_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, -2000., 2000., 10.);
  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_NC_4_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, -2000., 2000., 10.);

  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -2000., 2000., 10.);
  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_NC_40_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -2000., 2000., 10.);
  Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_NC_4_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -2000., 2000., 10.);

  //Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_anglecut_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, -2000., 2000., 10.);
  //Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_anglecut_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -2000., 2000., 10.);
  //Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQE_distcut_" + dir, "E_{QE, #pi^{#pm}} (MeV)", P_str, -2000., 2000., 10.);
  //Stack_MC_daughter("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "EQEmE_distcut_" + dir, "E_{QE, #pi^{#pm}} - E_{Beam,true} (MeV)", P_str, -2000., 2000., 10.);
}

void Draw_PionXsec_plots(){

  setTDRStyle();
  TString file_suffix = "_noBeamXY.root";
  file_suffix = "_EQE_RecoOnly.root";

  // ==== Fake Data vs MC
  Run_Draw_MC_vs_FakeData("PionXsec", file_suffix, "1.0", "beam_window_Preweight_piOnly");
  //Run_Draw_MC_vs_FakeData("PionXsec", file_suffix, "2.0", "Preweight_piOnly");

  // ==== Data vs MC
  /*
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "Preweight_piOnly");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "precut_noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "precut_Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "precut_Preweight_piOnly");


  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "Preweight_piOnly");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "precut_noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "precut_Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "precut_Preweight_piOnly");
  */
  //Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "Preweight_piOnly");
  //Run_Draw_MC_vs_Data("PionXsec", file_suffix, "2.0", "Preweight_piOnly");
  //Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "1.0", "beam_window_Preweight_piOnly");
  /*
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "Preweight_piOnly");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "precut_noweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "precut_Preweight");
  Run_Draw_MC_vs_Data("PionXsec", file_suffix, "0.5", "precut_Preweight_piOnly");
  */

  //file_suffix = "_norm_chi2_sigma.root";
  //file_suffix = "_this_chi2_sigma.root";
  // ==== Stacked MC

  //Run_Stack_MC_daughter("PionXsec", file_suffix, "1.0", "noweight");
  //Run_Stack_MC_daughter("PionXsec", file_suffix, "1.0", "Preweight_piOnly");

  // ==== MC daughter shape comparison
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "2.0", "noweight");
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "2.0", "Preweight_piOnly");
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "1.0", "noweight");
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "1.0", "Preweight_piOnly");
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "0.5", "noweight");
  //Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "0.5", "Preweight_piOnly");

  // ==== Daughter Selection
  //Run_Draw_MC_daughter_purity("PionXsec", file_suffix, "1.0", "Preweight_piOnly");
  //Run_Draw_MC_daughter_cutflow("PionXsec", file_suffix, "1.0", "Preweight_piOnly");

  // ==== Draw 2D plots
  //Run_Draw_2D_MC("PionXsec", file_suffix, "1.0", "Preweight_piOnly");
  /*
  Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "beam_inst_XY_precut", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "1.0", -50, 10, 5., 400., 440., 5.);
  Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "beam_inst_XY", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "1.0", -50, 10, 5., 400., 440., 5.);

  Draw_2D_MC_and_Data("_PionXsec_0.5GeV" + file_suffix, "beam_inst_XY_precut", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "0.5", -50, 10, 20., 400., 440., 20.);
  Draw_2D_MC_and_Data("_PionXsec_0.5GeV" + file_suffix, "beam_inst_XY", "X_{Beam Inst} (cm)", "Y_{Beam Inst} (cm)", "0.5", -50, 10, 20., 400., 440., 20.);
  */

  // ==== MC shape comparison
  /*
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "2.0", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "2.0", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "KEffTruth_precut", "KE_{True}^{FF}  (MeV)", "2.0", 0., 2500., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "PffTruth_precut", "P_{True}^{FF}  (MeV/c)", "2.0", 0., 2500., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_2.0GeV" + file_suffix, "BeamP_true", "P_{True}^{Start} (MeV/c)", "2.0", 0., 2500., 10., 0.2);

  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "1.0", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "1.0", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "KEffTruth_precut", "KE_{True}^{FF}  (MeV)", "1.0", 0., 1500., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "PffTruth_precut", "KE (P) _{True}^{FF}  (MeV/c)", "1.0", 0., 1500., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_1.0GeV" + file_suffix, "BeamP_true", "P_{True}^{Start} (MeV/c)", "1.0", 0., 1500., 10., 0.2);

  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "BeamKE_loss", "KE_{loss}^{True}  (MeV)", "0.5", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "BeamKE_loss_true_mass", "KE_{loss}^{True}  (MeV)", "0.5", -100., 100., 5., 0.4);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "KEffTruth_precut", "KE_{True}^{FF}  (MeV)", "0.5", 0., 1000., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "PffTruth_precut", "P_{True}^{FF}  (MeV/c)", "0.5", 0., 1000., 10., 0.2);
  Draw_MC_shape_comparison("_PionXsec_0.5GeV" + file_suffix, "BeamP_true", "P_{True}^{Start} (MeV/c)", "0.5", 0., 1000., 10., 0.2);
  */
}
