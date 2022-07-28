#include "canvas_margin.h"
#include "mylib.h"

void Fit_Beam_P_Data_and_MC(TString filename, TString beam_P, double xmin, double xmax, double rebin){

  TString histname_BeamP = "BeamP_precut";
  TString histname_BeamP_true = "BeamP_true";

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname_BeamP + "_%d", i);
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    maphist[this_hist_name] -> Rebin(rebin);
  }
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname_BeamP_true + "_%d", i);
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    maphist[this_hist_name] -> Rebin(rebin);
  }
  TH1D *hist_mc_sum_beam_inst = (TH1D*)maphist["htrack_" + histname_BeamP + "_1"] -> Clone();
  TH1D *hist_mc_sum_beam_true = (TH1D*)maphist["htrack_" + histname_BeamP_true + "_1"] -> Clone();
  for(int i = 2; i < N_pi_type; i++){
    TString this_hist_name_beam_inst = Form("htrack_" + histname_BeamP + "_%d", i);
    TString this_hist_name_beam_true = Form("htrack_" + histname_BeamP_true + "_%d", i);
    hist_mc_sum_beam_inst -> Add(maphist[this_hist_name_beam_inst]);
    hist_mc_sum_beam_true -> Add(maphist[this_hist_name_beam_true]);
  }
  

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  TH1D *hist_data = (TH1D*)gDirectory -> Get("htrack_" + histname_BeamP + "_0") -> Clone();
  hist_data -> Rebin(rebin);
  hist_data -> Scale(1. / hist_data -> Integral());
  hist_mc_sum_beam_inst -> Scale(0.7 / hist_mc_sum_beam_inst -> Integral());
  hist_mc_sum_beam_true -> Scale(0.5 / hist_mc_sum_beam_true -> Integral());
  double data_max = hist_data -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  // == Fit total mc and data P Beam Inst. and true
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("P_{Beam} (MeV/c)");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., data_max * 1.5);
  template_h -> Draw();

  hist_data -> SetLineColor(kBlack);
  hist_mc_sum_beam_inst -> SetLineColor(kGreen);
  hist_mc_sum_beam_true -> SetLineColor(kBlue);
  hist_data -> Draw("histsame");
  hist_mc_sum_beam_inst -> Draw("histsame");
  hist_mc_sum_beam_true -> Draw("histsame");
  double data_fit_x_min = hist_data -> GetMean() - 2.0 * hist_data -> GetRMS();
  double data_fit_x_max = hist_data -> GetMean() + 2.0 * hist_data -> GetRMS();
  double mc_fit_x_min = hist_mc_sum_beam_inst -> GetMean() - 2.0 * hist_mc_sum_beam_inst -> GetRMS();
  double mc_fit_x_max = hist_mc_sum_beam_inst -> GetMean() + 2.0 * hist_mc_sum_beam_inst -> GetRMS();
  double mc_true_fit_x_min = hist_mc_sum_beam_true -> GetMean() - 1.5 * hist_mc_sum_beam_true -> GetRMS();
  double mc_true_fit_x_max = hist_mc_sum_beam_true -> GetMean() + 1.5 * hist_mc_sum_beam_true -> GetRMS();

  TF1 *data_gaus = new TF1("data_gaus", "gaus", data_fit_x_min, data_fit_x_max);
  TF1 *mc_gaus = new TF1("mc_gaus", "gaus", mc_fit_x_min, mc_fit_x_max);
  TF1 *mc_true_gaus = new TF1("mc_true_gaus", "gaus", mc_true_fit_x_min, mc_true_fit_x_max);
  data_gaus -> SetLineColor(kGray);
  mc_gaus -> SetLineColor(kSpring+3);
  mc_true_gaus -> SetLineColor(kCyan);
  data_gaus -> SetLineStyle(7);
  mc_gaus -> SetLineStyle(7);
  mc_true_gaus -> SetLineStyle(7);
  data_gaus -> SetLineWidth(3);
  mc_gaus -> SetLineWidth(3);
  mc_true_gaus -> SetLineWidth(3);
  hist_data -> Fit(data_gaus, "R", "", data_fit_x_min, data_fit_x_max);
  hist_mc_sum_beam_inst -> Fit(mc_gaus, "R", "", mc_fit_x_min, mc_fit_x_max);
  hist_mc_sum_beam_true -> Fit(mc_true_gaus, "R", "", mc_true_fit_x_min, mc_true_fit_x_max);
  data_gaus -> Draw("lsame");
  mc_gaus -> Draw("lsame");
  mc_true_gaus -> Draw("lsame");
  
  TLegend *l = new TLegend(0.20, 0.72, 0.90, 0.92);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> SetNColumns(3);
  l -> AddEntry(hist_data, "Data_{Beam Inst.}", "l");
  l -> AddEntry(hist_mc_sum_beam_inst, "MC_{Beam Inst.}", "l");
  l -> AddEntry(hist_mc_sum_beam_true, "MC_{true}", "l");
  l -> AddEntry(data_gaus, Form("#mu : %.3f, #sigma : %.3f", data_gaus -> GetParameter(1), data_gaus -> GetParameter(2)), "l");
  l -> AddEntry(mc_gaus, Form("#mu : %.3f, #sigma : %.3f", mc_gaus -> GetParameter(1), mc_gaus -> GetParameter(2)), "l");
  l -> AddEntry(mc_true_gaus, Form("#mu : %.3f, #sigma : %.3f", mc_true_gaus -> GetParameter(1), mc_true_gaus -> GetParameter(2)), "l");
  l -> Draw("same");

  double data_conv_mu = data_gaus -> GetParameter(1) - mc_true_gaus -> GetParameter(1);
  double data_conv_sigma = sqrt(pow(data_gaus -> GetParameter(2), 2) - pow(mc_true_gaus -> GetParameter(2), 2));
  double mc_conv_mu = mc_gaus -> GetParameter(1) - mc_true_gaus -> GetParameter(1);
  double mc_conv_sigma = sqrt(pow(mc_gaus -> GetParameter(2), 2) - pow(mc_true_gaus -> GetParameter(2), 2));
  cout << "data conv (mu, sigma) = () : (" << data_conv_mu << ", " << data_conv_sigma << "), sigma / <P> = " << data_conv_sigma / mc_true_gaus -> GetParameter(1)<< endl;
  cout << "mc conv (mu, sigma) = () : (" << mc_conv_mu << ", " << mc_conv_sigma << "), sigma / <P> = " << mc_conv_sigma / mc_true_gaus -> GetParameter(1)<< endl;

  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/Beam_" + beam_P + "GeV_P_fit_Data_and_MC.pdf");

  // ==== Draw pion and muon beam
  // == Call histograms
  TH1D *hist_muon = (TH1D*)maphist["htrack_" + histname_BeamP + "_3"] -> Clone();
  TH1D *hist_muon_true = (TH1D*)maphist["htrack_" + histname_BeamP_true + "_3"] -> Clone();
  TH1D *hist_pion = (TH1D*)maphist["htrack_" + histname_BeamP + "_1"] -> Clone();
  TH1D *hist_pion_true = (TH1D*)maphist["htrack_" + histname_BeamP_true + "_1"] -> Clone();
  hist_pion -> Add(maphist["htrack_" + histname_BeamP + "_2"]);
  hist_pion_true -> Add(maphist["htrack_" + histname_BeamP_true + "_2"]);

  hist_muon -> Scale(1. / hist_muon -> Integral());
  hist_muon_true -> Scale(0.7 / hist_muon_true -> Integral());
  hist_pion -> Scale(1.0 / hist_pion -> Integral());
  hist_pion_true -> Scale(0.7 / hist_pion_true -> Integral());

  // == muon
  double muon_max = hist_muon -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., muon_max * 1.5);
  template_h -> Draw();
  hist_muon -> SetLineColor(kGreen);
  hist_muon_true -> SetLineColor(kBlue);
  hist_muon -> Draw("histsame");
  hist_muon_true -> Draw("histsame");
  double muon_fit_x_min = hist_muon -> GetMean() - 1.5 * hist_muon -> GetRMS();
  double muon_fit_x_max = hist_muon -> GetMean() + 1.5 * hist_muon -> GetRMS();
  double muon_true_fit_x_min = hist_muon_true -> GetMean() - 1.5 * hist_muon_true -> GetRMS();
  double muon_true_fit_x_max = hist_muon_true -> GetMean() + 1.5 * hist_muon_true -> GetRMS();

  TF1 *muon_gaus = new TF1("muon_gaus", "gaus", muon_fit_x_min, muon_fit_x_max);
  TF1 *muon_true_gaus = new TF1("muon_true_gaus", "gaus", muon_true_fit_x_min, muon_true_fit_x_max);
  muon_gaus -> SetLineColor(kSpring+3);
  muon_true_gaus -> SetLineColor(kCyan);
  muon_gaus -> SetLineStyle(7);
  muon_true_gaus -> SetLineStyle(7);
  muon_gaus -> SetLineWidth(3);
  muon_true_gaus -> SetLineWidth(3);
  hist_muon -> Fit(muon_gaus, "R", "", muon_fit_x_min, muon_fit_x_max);
  hist_muon_true -> Fit(muon_true_gaus, "R", "", muon_true_fit_x_min, muon_true_fit_x_max);
  muon_gaus -> Draw("lsame");
  muon_true_gaus -> Draw("lsame");

  TLegend *l_muon = new TLegend(0.20, 0.72, 0.90, 0.92);
  l_muon -> SetFillColor(kWhite);
  l_muon -> SetLineColor(kWhite);
  l_muon -> SetBorderSize(1);
  l_muon -> SetFillStyle(4000);
  l_muon -> SetShadowColor(0);
  l_muon -> SetNColumns(2);
  l_muon -> AddEntry(hist_muon, "MC_{Beam Inst.}", "l");
  l_muon -> AddEntry(hist_muon_true, "MC_{true}", "l");
  l_muon -> AddEntry(muon_gaus, Form("#mu : %.3f, #sigma : %.3f", muon_gaus -> GetParameter(1), muon_gaus -> GetParameter(2)), "l");
  l_muon -> AddEntry(muon_true_gaus, Form("#mu : %.3f, #sigma : %.3f", muon_true_gaus -> GetParameter(1), muon_true_gaus -> GetParameter(2)), "l");
  l_muon -> Draw("same");

  double muon_conv_mu = muon_gaus -> GetParameter(1) - muon_true_gaus -> GetParameter(1);
  double muon_conv_sigma = sqrt(pow(muon_gaus -> GetParameter(2), 2) - pow(muon_true_gaus -> GetParameter(2), 2));
  cout << "muon conv (mu, sigma) = () : (" << muon_conv_mu << ", " << muon_conv_sigma << "), sigma / <P> = " << muon_conv_sigma / muon_true_gaus -> GetParameter(1)<< endl;

  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/Beam_" + beam_P + "GeV_P_fit_MC_muon.pdf");

  // == pion
  double pion_max = hist_pion -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., pion_max * 1.5);
  template_h -> Draw();
  hist_pion -> SetLineColor(kGreen);
  hist_pion_true -> SetLineColor(kBlue);
  hist_pion -> Draw("histsame");
  hist_pion_true -> Draw("histsame");
  double pion_fit_x_min = hist_pion -> GetMean() - 2.0 * hist_pion -> GetRMS();
  double pion_fit_x_max = hist_pion -> GetMean() + 2.0 * hist_pion -> GetRMS();
  double pion_true_fit_x_min = hist_pion_true -> GetMean() - 2.0 * hist_pion_true -> GetRMS();
  double pion_true_fit_x_max = hist_pion_true -> GetMean() + 2.0 * hist_pion_true -> GetRMS();

  TF1 *pion_gaus = new TF1("pion_gaus", "gaus", pion_fit_x_min, pion_fit_x_max);
  TF1 *pion_true_gaus = new TF1("pion_true_gaus", "gaus", pion_true_fit_x_min, pion_true_fit_x_max);
  pion_gaus -> SetLineColor(kSpring+3);
  pion_true_gaus -> SetLineColor(kCyan);
  pion_gaus -> SetLineStyle(7);
  pion_true_gaus -> SetLineStyle(7);
  pion_gaus -> SetLineWidth(3);
  pion_true_gaus -> SetLineWidth(3);
  hist_pion -> Fit(pion_gaus, "R", "", pion_fit_x_min, pion_fit_x_max);
  hist_pion_true -> Fit(pion_true_gaus, "R", "", pion_true_fit_x_min, pion_true_fit_x_max);
  pion_gaus -> Draw("lsame");
  pion_true_gaus -> Draw("lsame");

  TLegend *l_pion = new TLegend(0.20, 0.72, 0.90, 0.92);
  l_pion -> SetFillColor(kWhite);
  l_pion -> SetLineColor(kWhite);
  l_pion -> SetBorderSize(1);
  l_pion -> SetFillStyle(4000);
  l_pion -> SetShadowColor(0);
  l_pion -> SetNColumns(2);
  l_pion -> AddEntry(hist_pion, "MC_{Beam Inst.}", "l");
  l_pion -> AddEntry(hist_pion_true, "MC_{true}", "l");
  l_pion -> AddEntry(pion_gaus, Form("#mu : %.3f, #sigma : %.3f", pion_gaus -> GetParameter(1), pion_gaus -> GetParameter(2)), "l");
  l_pion -> AddEntry(pion_true_gaus, Form("#mu : %.3f, #sigma : %.3f", pion_true_gaus -> GetParameter(1), pion_true_gaus -> GetParameter(2)), "l");
  l_pion -> Draw("same");

  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/Beam_" + beam_P + "GeV_P_fit_MC_pion.pdf");

  double pion_conv_mu = pion_gaus -> GetParameter(1) - pion_true_gaus -> GetParameter(1);
  double pion_conv_sigma = sqrt(pow(pion_gaus -> GetParameter(2), 2) - pow(pion_true_gaus -> GetParameter(2), 2));
  cout << "pion conv (mu, sigma) = () : (" << pion_conv_mu << ", " << pion_conv_sigma << "), sigma / <P> = " << pion_conv_sigma / pion_true_gaus -> GetParameter(1)<< endl;

  c -> Close();

}

void Fit_BeamP_Res_MC(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double fit_x_min, double fit_x_max, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  for(int i = 0; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    maphist[this_hist_name] -> Rebin(rebin);
  }
  TH1D *hist_mc_sum_beam_res = (TH1D*)maphist["htrack_" + histname + "_1"] -> Clone();
  for(int i = 2; i < N_pi_type; i++){
    TString this_hist_name = Form("htrack_" + histname + "_%d", i);
    hist_mc_sum_beam_res -> Add(maphist[this_hist_name]);
  }
  hist_mc_sum_beam_res -> Scale(1. / hist_mc_sum_beam_res -> Integral());
  double this_max = hist_mc_sum_beam_res -> GetMaximum();
  cout << "[Fit_BeamP_Res_MC] this_max : " << this_max << endl;
  TCanvas *c = new TCanvas("", "", 600, 800);
  canvas_margin(c);
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
  template_h -> GetYaxis() -> SetRangeUser(0., this_max * 1.5);
  template_h -> Draw();

  hist_mc_sum_beam_res -> SetLineColor(kBlue);
  hist_mc_sum_beam_res -> Draw("histsame");

  TF1 *this_gaus = new TF1("this_gaus", "gaus", fit_x_min, fit_x_max);
  this_gaus -> SetLineColor(kCyan);
  this_gaus -> SetLineStyle(7);
  this_gaus -> SetLineWidth(3);
  hist_mc_sum_beam_res -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
  this_gaus -> Draw("lsame");

  TLegend *l = new TLegend(0.20, 0.72, 0.90, 0.92);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> AddEntry(hist_mc_sum_beam_res, "#Delta p / p_{True}", "l");
  l -> AddEntry(this_gaus, Form("#mu : %.3f, #sigma : %.4f", this_gaus -> GetParameter(1), this_gaus -> GetParameter(2)), "l");
  l -> Draw("same");
  
  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/BeamP_Res_" + beam_P + "GeV.pdf");

  // ==== Draw pion and muon beam
  // == Call histograms
  TH1D *hist_muon = (TH1D*)maphist["htrack_" + histname + "_3"] -> Clone();
  TH1D *hist_pion = (TH1D*)maphist["htrack_" + histname + "_1"] -> Clone();
  hist_pion -> Add(maphist["htrack_" + histname + "_2"]);

  hist_muon -> Scale(1.0 / hist_muon -> Integral());
  hist_pion -> Scale(1.0 / hist_pion -> Integral());

  // == muon
  double muon_max = hist_muon -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., muon_max * 1.5);
  template_h -> Draw();
  hist_muon -> SetLineColor(kBlue);
  hist_muon -> Draw("histsame");
  double muon_fit_x_min = hist_muon -> GetMean() - 1.5 * hist_muon -> GetRMS();
  double muon_fit_x_max = hist_muon -> GetMean() + 1.5 * hist_muon -> GetRMS();

  TF1 *muon_gaus = new TF1("muon_gaus", "gaus", muon_fit_x_min, muon_fit_x_max);
  muon_gaus -> SetLineColor(kCyan);
  muon_gaus -> SetLineStyle(7);
  muon_gaus -> SetLineWidth(3);
  hist_muon -> Fit(muon_gaus, "R", "", muon_fit_x_min, muon_fit_x_max);
  muon_gaus -> Draw("lsame");

  TLegend *l_muon = new TLegend(0.20, 0.72, 0.90, 0.92);
  l_muon -> SetFillColor(kWhite);
  l_muon -> SetLineColor(kWhite);
  l_muon -> SetBorderSize(1);
  l_muon -> SetFillStyle(4000);
  l_muon -> SetShadowColor(0);
  l_muon -> AddEntry(hist_muon, "#Delta p / p_{True}", "l");
  l_muon -> AddEntry(muon_gaus, Form("#mu : %.3f, #sigma : %.4f", muon_gaus -> GetParameter(1), muon_gaus -> GetParameter(2)), "l");
  l_muon -> Draw("same");

  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/BeamP_Res_" + beam_P + "GeV_muon.pdf");

  // == pion
  double pion_max = hist_pion -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., pion_max * 1.5);
  template_h -> Draw();
  hist_pion -> SetLineColor(kBlue);
  hist_pion -> Draw("histsame");
  double pion_fit_x_min = hist_pion -> GetMean() - 1.5 * hist_pion -> GetRMS();
  double pion_fit_x_max = hist_pion -> GetMean() + 1.5 * hist_pion -> GetRMS();

  TF1 *pion_gaus = new TF1("pion_gaus", "gaus", pion_fit_x_min, pion_fit_x_max);
  pion_gaus -> SetLineColor(kCyan);
  pion_gaus -> SetLineStyle(7);
  pion_gaus -> SetLineWidth(3);
  hist_pion -> Fit(pion_gaus, "R", "", pion_fit_x_min, pion_fit_x_max);
  pion_gaus -> Draw("lsame");

  TLegend *l_pion = new TLegend(0.20, 0.72, 0.90, 0.92);
  l_pion -> SetFillColor(kWhite);
  l_pion -> SetLineColor(kWhite);
  l_pion -> SetBorderSize(1);
  l_pion -> SetFillStyle(4000);
  l_pion -> SetShadowColor(0);
  l_pion -> AddEntry(hist_pion, "#Delta p / p_{True}", "l");
  l_pion -> AddEntry(pion_gaus, Form("#mu : %.3f, #sigma : %.4f", pion_gaus -> GetParameter(1), pion_gaus -> GetParameter(2)), "l");
  l_pion -> Draw("same");

  c -> SaveAs("./output/plots/PionXsec/Beam_P_fitting/BeamP_Res_" + beam_P + "GeV_pion.pdf");

  c -> Close();  

}


void Fit_Beam_P(){

  cout << "============================" << endl;
  cout << "============START===========" << endl;
  cout << "============================" << endl;


  setTDRStyle();
  TString file_suffix = "_noBeamXY.root";

  Fit_Beam_P_Data_and_MC("_PionXsec_0.5GeV" + file_suffix, "0.5", 0., 1000., 10.);
  Fit_Beam_P_Data_and_MC("_PionXsec_1.0GeV" + file_suffix, "1.0", 0., 1500., 10.);
  Fit_Beam_P_Data_and_MC("_PionXsec_2.0GeV" + file_suffix, "2.0", 0., 2500., 10.);

  Fit_BeamP_Res_MC("_PionXsec_0.5GeV" + file_suffix, "BeamP_Res_precut", "(P_{True} - P_{Beam Inst.}) / P_{True}", "0.5", -0.2, 0.2, -0.05, 0.05, 10.);
  Fit_BeamP_Res_MC("_PionXsec_1.0GeV" + file_suffix, "BeamP_Res_precut", "(P_{True} - P_{Beam Inst.}) / P_{True}", "1.0", -0.2, 0.2, -0.05, 0.05, 1.); 
  Fit_BeamP_Res_MC("_PionXsec_2.0GeV" + file_suffix, "BeamP_Res_precut", "(P_{True} - P_{Beam Inst.}) / P_{True}", "2.0", -0.2, 0.2, -0.05, 0.05, 10.);

}

