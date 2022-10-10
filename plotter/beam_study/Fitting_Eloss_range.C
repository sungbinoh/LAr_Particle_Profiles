#include "canvas_margin.h"
#include "mylib.h"

void Fit_and_Draw(TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  //gDirectory -> ls();
  vector<double> KE_vec;
  vector<double> KE_err_vec;
  vector<double> mean_vec;
  vector<double> mean_err_vec;
  vector<double> std_vec;
  
  int KE_low = 700;
  int KE_high = 1100;
  int KE_step = 50;
  int N_KE_steps = (KE_high - KE_low) / KE_step;
  for(int i = 0; i < N_KE_steps; i++){
    //KE1000to1100MeV
    TString KE_range = Form("KE%dto%dMeV", KE_low + KE_step * i, KE_low + KE_step * (i + 1) );
    TString this_histname = "htrack_" + histname + "_Beam" + KE_range + "_" + dir;
    cout << "[Fit_and_Draw] Fitting " << this_histname  + "_1" << endl;
    if(!((TH1D*)gDirectory -> Get(this_histname + "_1"))) continue;
    cout << "[Fit_and_Draw] Found the histogram" << endl;
    TH1D *this_hist = (TH1D*)gDirectory -> Get(this_histname + "_1") -> Clone();
    if((TH1D*)gDirectory -> Get(this_histname + "_2")) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_2")); 
    double max_x = this_hist -> GetBinCenter(this_hist -> GetMaximumBin());
    double fit_x_min = max_x - width;
    double fit_x_max = max_x + width;
    TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
    this_hist -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
    double this_mean = this_gaus -> GetParameter(1);
    double this_std = this_gaus -> GetParameter(2);
    double this_mean_err = this_gaus -> GetParError(1);
    mean_vec.push_back(this_mean);
    std_vec.push_back(this_std);
    mean_err_vec.push_back(this_mean_err);

    double this_KE = (KE_low + 0.) + (KE_step + 0.) * (i + 0.) + (KE_step + 0.) / 2.0;
    KE_vec.push_back(this_KE);
    KE_err_vec.push_back((KE_step + 0.) / 2.0);
    cout << "[Fit_and_Draw] this_KE : " << this_KE << " +- " << (KE_step + 0.) / 2.0 << endl;
  }

  TGraphErrors *Eloss_ratio_gr = new TGraphErrors(50, &KE_vec[0], &mean_vec[0], &KE_err_vec[0], &mean_err_vec[0]);

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH1D* template_h = new TH1D("", "", 1., 600., 1200.);
  template_h -> SetStats(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle(title_X);
  template_h -> GetYaxis() -> SetTitle(title_Y);
  template_h -> GetYaxis() -> SetRangeUser(-50., 100.);
  template_h -> Draw();

  Eloss_ratio_gr -> SetLineColor(kGreen);
  Eloss_ratio_gr -> SetLineWidth(3);
  Eloss_ratio_gr -> Draw("epsame");

  TF1 * f_pol_1 = new TF1("f_pol_1", "pol1", KE_low + 0., KE_high + 0.);
  TF1 * f_pol_2 = new TF1("f_pol_2", "pol2" , KE_low + 0., KE_high + 0.);
  f_pol_1 -> SetStats(0);
  f_pol_2 -> SetStats(0);
  Eloss_ratio_gr -> Fit(f_pol_1, "R0", "", KE_low + 0., KE_high + 0.);
  Eloss_ratio_gr -> Fit(f_pol_2, "R0", "", KE_low + 0., KE_high + 0.);

  f_pol_1 -> SetLineColor(kRed);
  f_pol_2 -> SetLineColor(kRed);
  f_pol_1 -> SetLineStyle(7);
  f_pol_2 -> SetLineStyle(7);
  if(histname.Contains("diff_true")) f_pol_1 -> Draw("lsame");
  else f_pol_2 -> Draw("lsame");

  TF1 * f_comparison_pol_2 = new TF1("f_comparison_pol_2", "pol2" , KE_low + 0., KE_high + 0.);
  f_comparison_pol_2 -> SetParameter(0, 1.52626e+02);
  f_comparison_pol_2 -> SetParameter(1, -4.68246e-01);
  f_comparison_pol_2 -> SetParameter(2, 3.48254e-04);
  f_comparison_pol_2 -> SetLineColor(kBlue);
  f_comparison_pol_2 -> SetLineStyle(7);
  f_comparison_pol_2 -> Draw("lsame");

  TString fit_fuction_str = "Fitted p_{1}x + p_{0}";
  if(!histname.Contains("diff_true")) fit_fuction_str = "Fitted p_{2} x^{2} + p_{1}x + p_{0}";

  TString fit_result_str = "";
  if(histname.Contains("diff_true")) fit_result_str = Form("%.2e + %.2e x", f_pol_1 -> GetParameter(0), f_pol_1 -> GetParameter(1));
  else fit_result_str = Form("%.2e + %.2e x + %.2e x^{2}", f_pol_2 -> GetParameter(0), f_pol_2 -> GetParameter(1), f_pol_2 -> GetParameter(2));

  TLegend *l = new TLegend(0.18, 0.7, 0.55, 0.85);
  l -> AddEntry(f_pol_1, fit_fuction_str, "l");
  l -> AddEntry(f_comparison_pol_2, "#mu + 9.48", "l");
  l -> SetLineColor(kWhite); 
  l -> Draw("same");

  double p0, p1, p2, p0_err, p1_err, p2_err;
  if(histname.Contains("diff_true")){
    p0 = f_pol_1 -> GetParameter(0);
    p1 = f_pol_1 -> GetParameter(1);
    p2 = f_pol_1 -> GetParameter(2);
    p0_err = f_pol_1 -> GetParError(0);
    p1_err = f_pol_1 -> GetParError(1);
    p2_err = f_pol_1 -> GetParError(2);
  }
  else{
    p0 = f_pol_2 -> GetParameter(0);
    p1 = f_pol_2 -> GetParameter(1);
    p2 = f_pol_2 -> GetParameter(2);
    p0_err = f_pol_2 -> GetParError(0);
    p1_err = f_pol_2 -> GetParError(1);
    p2_err = f_pol_2 -> GetParError(2);
  }
  TString p0_str, p1_str, p2_str;
  p0_str = Form("p_{0} = %.2e #pm %.2e", p0, p0_err);
  p1_str = Form("p_{1} = %.2e #pm %.2e", p1, p1_err);
  p2_str = Form("p_{2} = %.2e #pm %.2e", p2, p2_err);
  TLatex latex_p0, latex_p1, latex_p2;
  latex_p0.SetNDC();
  latex_p1.SetNDC();
  latex_p2.SetNDC();
  latex_p0.SetTextSize(0.035);
  latex_p1.SetTextSize(0.035);
  latex_p2.SetTextSize(0.035);
  latex_p0.DrawLatex(0.20, 0.66, p0_str);
  latex_p1.DrawLatex(0.20, 0.62, p1_str);
  if(!histname.Contains("diff_true")) latex_p2.DrawLatex(0.20, 0.58, p2_str);

  TLatex latex_ProtoDUNE, latex_sample;
  latex_ProtoDUNE.SetNDC();
  latex_sample.SetNDC();
  latex_sample.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_sample.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_sample.DrawLatex(0.95, 0.96, "MC ( #pi^{#pm} Elas. & Inel.)");

  TLine *this_line = new TLine(0., 1, 50., 1.);
  this_line -> SetLineStyle(5);
  this_line -> SetLineColor(kRed);
  this_line -> Draw("same");

  TLine *best_line = new TLine(23., 0.8, 23., 1.2); 
  best_line -> SetLineStyle(5);
  best_line -> SetLineColor(kRed);
  best_line -> Draw("same");

  c -> SaveAs("./output/plots/BeamStudy/Upstream_Eloss/Eloss_" + histname + "_" + dir + ".pdf");

  c -> Close();
}


void Fitting_Eloss_range(){
  setTDRStyle();
  
  Fit_and_Draw("_PionXsec_1.0GeV_beam_study.root", "noweight", "KE_diff_nonscraper", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.} - KE_{ff}^{true}) [MeV]", 40.);
  Fit_and_Draw("_PionXsec_1.0GeV_beam_study.root", "noweight", "KE_diff_true_nonscraper", "KE_{Beam Inst.} [MeV]", "#mu (KE_{true} - KE_{ff}^{true}) [MeV]", 40.);
  Fit_and_Draw("_PionXsec_1.0GeV_beam_study.root", "noweight", "KE_diff_beam_true_nonscraper", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.} - KE^{true}) [MeV]", 40.);

}
