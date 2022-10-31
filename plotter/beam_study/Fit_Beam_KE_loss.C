#include "canvas_margin.h"
#include "mylib.h"

void Fit_Beam_DeltaKE(TString data_or_mc, TString filename, TString KE_diff_str, TString dir, TString beam_P, TString this_KE_str, TString this_KE_legend, TString title_X, double xmin, double xmax, double rebin){

  TString histname = "htrack_KE_diff_" + KE_diff_str + "_KE_beam_inst" + this_KE_str + "_" + dir;
  TString output_name = "KE_diff_" + KE_diff_str + "_KE" + this_KE_str;

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> Cd(dir);
  for(int i = 0; i < 9; i++){
    TString this_hist_name = Form(histname + "_%d", i);
    cout << "[Fit_Beam_DeltaKE] this_hist_name : " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }
  
  TH1D *h_KE_diff;
  if(data_or_mc == "data"){
    h_KE_diff = (TH1D*)maphist[histname + "_0"] -> Clone();
  }
  else{
    if(dir.Contains("pion")){
      h_KE_diff = (TH1D*)maphist[histname + "_1"] -> Clone();
      h_KE_diff -> Add(maphist[histname+ "_2"]);
    }
    if(dir.Contains("proton")){
      h_KE_diff = (TH1D*)maphist[histname + "_2"] -> Clone();
      //if(!KE_diff_str.Contains("FittedFF")) h_KE_diff -> Add(maphist[histname+ "_1"]);

      if(maphist[histname+ "_1"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_1"]);
      if(maphist[histname+ "_2"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_2"]);
      if(maphist[histname+ "_3"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_3"]);
      if(maphist[histname+ "_4"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_4"]);
      if(maphist[histname+ "_5"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_5"]);
      if(maphist[histname+ "_6"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_6"]);
      if(maphist[histname+ "_7"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_7"]);
      if(maphist[histname+ "_8"] != nullptr) h_KE_diff -> Add(maphist[histname+ "_8"]);
    }
  }
  double y_max = h_KE_diff -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  // == Fit total mc and data P Beam Inst. and true
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(title_X + " [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  h_KE_diff -> SetLineColor(kBlack);
  h_KE_diff -> SetLineWidth(2);
  h_KE_diff -> Draw("histsame");

  double max_x = h_KE_diff -> GetBinCenter(h_KE_diff -> GetMaximumBin());
  double fit_x_min = max_x - 35.;
  double fit_x_max = max_x + 35.;

  TF1 *f_gaus = new TF1("f_gaus", "gaus", fit_x_min, fit_x_max);
  f_gaus -> SetLineColor(kRed);
  f_gaus -> SetLineWidth(3);
  h_KE_diff -> Fit(f_gaus, "R", "", fit_x_min, fit_x_max);
  f_gaus -> Draw("lsame");
  
  TF1 *f_gaus_extended = new TF1("f_gaus_extended", "gaus", xmin, xmax);
  f_gaus_extended -> SetParameter(0, f_gaus -> GetParameter(0));
  f_gaus_extended -> SetParameter(1, f_gaus -> GetParameter(1));
  f_gaus_extended -> SetParameter(2, f_gaus -> GetParameter(2));
  f_gaus_extended -> SetLineColor(kRed);
  f_gaus_extended -> SetLineStyle(7);
  f_gaus_extended -> SetLineWidth(2);
  f_gaus_extended -> Draw("lsame");

  double x_3sigma = f_gaus -> GetParameter(1) + 3.0 * f_gaus -> GetParameter(2);
  TLine * line_3sigma = new TLine(x_3sigma, 0., x_3sigma, y_max);
  line_3sigma -> SetLineStyle(7);
  line_3sigma -> SetLineColor(kBlue);
  line_3sigma -> Draw("lsame");
  TString x_3sigma_str = Form("3#sigma : %.3f", x_3sigma);

  TLegend *l = new TLegend(0.20, 0.72, 0.90, 0.92);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> AddEntry(h_KE_diff, title_X, "l");
  l -> AddEntry(f_gaus, Form("#mu : %.3f, #sigma : %.3f", f_gaus -> GetParameter(1), f_gaus -> GetParameter(2)), "l");
  l -> AddEntry(f_gaus_extended, "Extended", "l");
  l -> AddEntry(line_3sigma, x_3sigma_str, "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_KE;
  latex_ProtoDUNE.SetNDC();
  latex_KE.SetNDC();
  latex_KE.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_KE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, "KE_{Beam Inst.} = " + this_KE_legend);

  TString particle_str = "";
  if(dir.Contains("pion")) particle_str = "pion";
  if(dir.Contains("proton")) particle_str = "proton";
  TString pdf_name = "./output/plots/BeamStudy/KE_fit/" + particle_str + "/Fit_" + data_or_mc + "_" + particle_str + "_" + beam_P + "GeV_" + output_name + ".pdf";
  c -> SaveAs(pdf_name);

  c -> Close();
}

void Run_Fit_Beam_DeltaKE(TString filename, TString dir, TString beam_P, double xmin, double xmax, double rebin){

  int KE_step = 100;

  int pion_KE_min = 700;
  int pion_KE_max = 1100;
  int pion_N_step = (pion_KE_max - pion_KE_min) / KE_step;
  for(int i = 0; i < pion_N_step; i++){
    TString this_KE_str = Form("%dto%dMeV", pion_KE_min + i * KE_step, pion_KE_min + (i + 1) * KE_step);
    TString this_KE_legend = Form("%d - %dMeV", pion_KE_min + i * KE_step, pion_KE_min + (i + 1) * KE_step);
    Fit_Beam_DeltaKE("mc", filename, "RecoBeam_TrueFF", "pion_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{true} [MeV]", -300., 300., rebin);
  }

  int proton_KE_min = 300;
  int proton_KE_max = 600;
  int proton_N_step = (proton_KE_max - proton_KE_min) / KE_step;
  for(int i = 0; i < proton_N_step; i++){
    TString this_KE_str = Form("%dto%dMeV", proton_KE_min + i * KE_step, proton_KE_min + (i + 1) * KE_step);
    TString this_KE_legend = Form("%d - %dMeV", proton_KE_min + i * KE_step, proton_KE_min + (i + 1) * KE_step);
    Fit_Beam_DeltaKE("mc", filename, "RecoBeam_TrueFF", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{true}", -450., 450., rebin);
    Fit_Beam_DeltaKE("mc", filename, "RecoBeam_FittedFF", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{fitted}", -450., 450., rebin);
    Fit_Beam_DeltaKE("mc", filename, "TrueFF_FittedFF", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{ff}^{true} - KE_{ff}^{fitted}", -50., 50., 1.);
    Fit_Beam_DeltaKE("data", filename, "RecoBeam_FittedFF", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{fitted}", -450., 450., rebin);

    Fit_Beam_DeltaKE("mc", filename, "RecoBeam_TrueFF_nonscraper", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{true}", -450., 450., rebin);
    Fit_Beam_DeltaKE("mc", filename, "RecoBeam_FittedFF_nonscraper", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{fitted}", -450., 450., rebin);
    Fit_Beam_DeltaKE("data", filename, "RecoBeam_FittedFF_nonscraper", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, "KE_{Beam Inst.} - KE_{ff}^{fitted}", -450., 450., rebin);
  }
}

void Fit_Beam_KE_loss(){

  cout << "============================" << endl;
  cout << "=========== START ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = ".root";
  //data_Beam_Study_1.0GeV.root
  Run_Fit_Beam_DeltaKE("_Beam_Study_1.0GeV" + file_suffix, "noweight", "1.0", -300., 300., 4.);

}
