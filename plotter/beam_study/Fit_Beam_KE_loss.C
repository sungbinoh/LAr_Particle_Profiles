#include "canvas_margin.h"
#include "mylib.h"

void Fit_Beam_DeltaKE(TString filename, TString dir, TString beam_P, TString this_KE_str, double xmin, double xmax, double rebin){

  TString histname = "htrack_KE_diff_nonscraper_Beam" + this_KE_str + "_" + dir;

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);
  for(int i = 1; i < 3; i++){
    TString this_hist_name = Form(histname + "_%d", i);
    cout << "[Fit_Beam_DeltaKE] this_hist_name : " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }

  TH1D *h_pi = (TH1D*)maphist[histname + "_1"] -> Clone();
  if(maphist[histname+ "_1"] != nullptr) h_pi -> Add(maphist[histname+ "_2"]);
  double y_max = h_pi -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  // == Fit total mc and data P Beam Inst. and true
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE_{Beam Inst.} - KE_{ff}^{true} [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  h_pi -> SetLineColor(kBlack);
  h_pi -> SetLineWidth(2);
  h_pi -> Draw("histsame");

  double max_x = h_pi -> GetBinCenter(h_pi -> GetMaximumBin());
  //double fit_x_min = h_pi -> GetMean() - 1.5 * h_pi -> GetRMS();
  //double fit_x_max = h_pi -> GetMean() + 1.5 * h_pi -> GetRMS();
  double fit_x_min = max_x - 40.;
  double fit_x_max = max_x + 40.;


  TF1 *f_gaus = new TF1("f_gaus", "gaus", fit_x_min, fit_x_max);
  f_gaus -> SetLineColor(kRed);
  //f_gaus -> SetLineStyle(7);
  f_gaus -> SetLineWidth(3);
  h_pi -> Fit(f_gaus, "R", "", fit_x_min, fit_x_max);
  f_gaus -> Draw("lsame");
  
  TF1 *f_gaus_extended = new TF1("f_gaus_extended", "gaus", xmin, xmax);
  f_gaus_extended -> SetParameter(0, f_gaus -> GetParameter(0));
  f_gaus_extended -> SetParameter(1, f_gaus -> GetParameter(1));
  f_gaus_extended -> SetParameter(2, f_gaus -> GetParameter(2));
  f_gaus_extended -> SetLineColor(kRed);
  f_gaus_extended -> SetLineStyle(7);
  f_gaus_extended -> SetLineWidth(2);
  f_gaus_extended -> Draw("lsame");

  TLegend *l = new TLegend(0.20, 0.72, 0.90, 0.92);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> AddEntry(h_pi, "KE_{Beam Inst.} - KE_{ff}^{ture}", "l");
  l -> AddEntry(f_gaus, Form("#mu : %.3f, #sigma : %.3f", f_gaus -> GetParameter(1), f_gaus -> GetParameter(2)), "l");
  l -> AddEntry(f_gaus_extended, "Extended", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_KE;
  latex_ProtoDUNE.SetNDC();
  latex_KE.SetNDC();
  latex_KE.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_KE.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, this_KE_str);

  c -> SaveAs("./output/plots/BeamStudy/KE_fit/Fit_" + beam_P + "GeV_" + histname + "_nonscraper.pdf");

  c -> Close();
}

void Run_Fit_Beam_DeltaKE(TString filename, TString dir, TString beam_P, double xmin, double xmax, double rebin){

  int KE_min = 600;
  int KE_max = 1200;
  int KE_step = 100;
  int N_step = (KE_max - KE_min) / KE_step;
  vector<TString> KE_str_vec;
  for(int i = 0; i < N_step; i++){
    //TString this_KE_str = Form("htrack_KE_diff_BeamKE%dto%dMeV_", KE_min + i * KE_step, KE_min + (i + 1) * KE_step);
    TString this_KE_str = Form("KE%dto%dMeV", KE_min + i * KE_step, KE_min + (i + 1) * KE_step);

    cout << "this_KE_str : " << this_KE_str << endl;
    
    Fit_Beam_DeltaKE(filename, dir, beam_P, this_KE_str, xmin, xmax, rebin);
  }
}

void Fit_Beam_KE_loss(){

  cout << "============================" << endl;
  cout << "=========== START ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = "_beam_study.root";

  Run_Fit_Beam_DeltaKE("_PionXsec_1.0GeV" + file_suffix, "precut_noweight", "1.0", -300., 300., 4.);

}

