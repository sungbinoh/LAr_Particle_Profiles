#include "canvas_margin.h"
#include "mylib.h"

void Draw_MC_vs_Data(TString filename, TString dir, TString P_str, TString histname, TString title_X, double x_min, double x_max, double rebin, bool logy){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";

  // == Call MC Histograms
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);
  int N_type = 0;
  if(dir.Contains("proton")) N_type = N_p_type;
  else if(dir.Contains("pion")) N_type = N_pi_type;
  TH1D * mc_sum;
  for(int i = 0; i < N_type; i++){
    TString this_hist_name = Form(histname + "_" + dir + "_%d", i);

    if(i == 1) mc_sum = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();

    if((TH1D*)gDirectory -> Get(this_hist_name)){
      if(i == 1) mc_sum =(TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      else if (i != 0){
	mc_sum -> Add((TH1D*)gDirectory -> Get(this_hist_name));
      }
      
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }
  mc_sum -> Rebin(rebin);

  // == Call Data Histogram
  TFile *f_data = new TFile(root_file_path + "data" + filename);
  gDirectory -> cd(dir);
  TH1D *hist_data = (TH1D*)gDirectory -> Get(histname + "_" + dir + "_0") -> Clone();
  hist_data -> Rebin(rebin);

  // == Get Scale and apply
  double data_integ = hist_data -> Integral();
  double mc_integ = mc_sum -> Integral();
  double mc_scale = data_integ / mc_integ;
  for(int i = 0; i < N_type; i++){
    TString this_hist_name = Form(histname + "_" + dir + "_%d", i);
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> Scale(mc_scale);
    }
  }
  mc_sum -> Scale(mc_scale);

  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  // == Top pad
  TPad *pad_1 = new TPad("", "", 0, 0.25, 1, 1);
  pad_1 -> SetTopMargin( 0.07 );
  pad_1 -> SetBottomMargin( 0.05 );
  pad_1 -> SetLeftMargin( 0.15 );
  pad_1 -> SetRightMargin( 0.03 );
  pad_1 -> Draw();
  pad_1 -> cd();
  if(logy) pad_1 -> SetLogy();

  double y_max = hist_data -> GetMaximum();
  if(y_max < mc_sum -> GetMaximum()) y_max = mc_sum ->GetMaximum();
  TH1D *pad1_template = new TH1D("", "", 1, x_min, x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  pad1_template -> SetStats(0);
  pad1_template -> GetXaxis() -> SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitleSize(0);
  pad1_template -> GetYaxis() -> SetLabelSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleSize(0.07);
  pad1_template -> GetYaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetTitle("Events");
  pad1_template -> GetYaxis() -> SetRangeUser(0., y_max * 1.8);
  if(logy) pad1_template -> GetYaxis() -> SetRangeUser(0.1, y_max * 10000.);
  pad1_template -> Draw("hist");

  TLegend *l_pad1 = new TLegend(0.20, 0.70, 0.90, 0.90);
  l_pad1 -> SetFillColor(kWhite);
  l_pad1 -> SetLineColor(kWhite);
  l_pad1 -> SetBorderSize(1);
  l_pad1 -> SetFillStyle(4000);
  l_pad1 -> SetShadowColor(0);
  l_pad1 -> SetEntrySeparation(0.3);
  l_pad1 -> SetNColumns(3);

  THStack * hist_mc_stack = new THStack("", "");
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  for(int i = 1; i < N_type; i++){
    TString this_hist_name = Form(histname + "_" + dir + "_%d", i);
    if(maphist[this_hist_name] != nullptr){
      TString this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
      TString this_legend_str = "";
      if(dir.Contains("proton")) this_legend_str = p_type_str[i] + " " + this_N_event;
      else if(dir.Contains("pion")) this_legend_str = pi_type_str[i] + " " + this_N_event;
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetFillColor(colour_array[i]);
      hist_mc_stack -> Add(maphist[this_hist_name]);
      l_pad1 -> AddEntry(maphist[this_hist_name], this_legend_str, "f");
    }
  }
  TString mc_sum_N_event = Form("%.1f", mc_sum -> Integral());
  mc_sum -> SetLineColor(kWhite);
  l_pad1 -> AddEntry(mc_sum, "MC Sum " + mc_sum_N_event, "l");

  hist_data -> SetLineColor(kBlack);
  hist_data -> SetMarkerColor(kBlack);
  hist_data -> SetMarkerStyle(20);
  hist_data -> SetMarkerSize(1);
  hist_data -> SetLineWidth(2);
  TString data_N_event = Form("%.1f", hist_data -> Integral());
  l_pad1 -> AddEntry(hist_data, "Observed " + data_N_event, "lp");
  hist_mc_stack -> Draw("histsame");
  hist_data -> Draw("epsame");
  l_pad1 -> Draw("same");
  gPad->RedrawAxis();

  // == Bottom pad
  c -> cd();
  TPad *pad_2 = new TPad("", "", 0, 0, 1, 0.25);
  pad_2 -> SetTopMargin( 0.05 );
  pad_2 -> SetBottomMargin( 0.4 );
  pad_2 -> SetLeftMargin( 0.15 );
  pad_2 -> SetRightMargin( 0.03 );
  pad_2 -> Draw();
  pad_2 -> cd();
  TH1D * pad2_template = new TH1D("", "", 1, x_min, x_max);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  pad2_template -> Draw("hist");
  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(title_X);
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
  
  TLegend *l_pad2 = new TLegend(0.2, 0.85, 0.6, 0.95);
  l_pad2 -> SetBorderSize(0);
  l_pad2 -> SetNColumns(3);
  l_pad2 -> AddEntry(data_mc_ratio, "Obs./Pred.", "lp");

  TLine *pad2_line = new TLine(x_min, 1, x_max, 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

  l_pad2 -> Draw("same");

  gPad->RedrawAxis();

  // == Latex
  c -> cd();
  TLatex latex_ProtoDUNE, latex_data;
  latex_ProtoDUNE.SetNDC();
  latex_data.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data.DrawLatex(0.63, 0.96, "Run 1, " + P_str + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/MC_vs_Data/MC_vs_Data_" + histname + "_" + P_str + "GeV_" + dir + ".pdf";
  c -> SaveAs(pdfname);

  f_mc -> Close();
  f_data -> Close();
  c -> Close();
}

void Run_MC_vs_Data(TString filename, TString P_str, TString dir){

  Draw_MC_vs_Data(filename, "proton_BeamScraper_" + dir, P_str, "htrack_KE_ff_reco", "KE_{ff}^{reco} [MeV]", 0., 700., 5., true);
}

void Draw_MC_vs_Data(){
  setTDRStyle();
  TString file_suffix = ".root";
  Run_MC_vs_Data("_Beam_Study_1.0GeV" + file_suffix, "1.0", "noweight");

  //file_suffix = "_beam.root";
  //Run_MC_vs_Data("_PionXsec_1.0GeV" + file_suffix, "1.0", "noweight");
  //mc_PionXsec_1.0GeV_beam.root
}
