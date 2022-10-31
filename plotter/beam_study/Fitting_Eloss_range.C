#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"
Double_t langaufun(Double_t *x, Double_t *par) {
  Double_t invsq2pi = 0.398942280401;// Control constants
  //Double_t mpshift = -0.22278298l
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


void Fit_and_Draw(TString data_or_mc, TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_out = new TFile(root_file_path + "Fitting_Eloss_" + data_or_mc + "_" + histname + "_" + dir + filename, "RECREATE");
  TFile *f_mc = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> cd(dir);
  vector<double> KE_vec;
  vector<double> KE_err_vec;
  vector<double> mean_vec;
  vector<double> mean_err_vec;
  vector<double> std_vec;
  
  int KE_low = 700;
  int KE_high = 1100;
  if(dir.Contains("proton")){
    KE_low = 300.;
    KE_high = 600.;
  }
  int KE_step = 50;
  int N_KE_steps = (KE_high - KE_low) / KE_step;
  for(int i = 0; i < N_KE_steps; i++){
    TString KE_range = Form("%dto%dMeV", KE_low + KE_step * i, KE_low + KE_step * (i + 1) );
    TString this_histname = "htrack_KE_diff_" + histname + "_nonscraper_KE_beam_inst" + KE_range + "_" + dir;
    if(histname.Contains("KE_ff_reco")) this_histname = "htrack_diff_" + histname + "_KE_beam_inst" + KE_range + "_" + dir;
    cout << "[Fit_and_Draw] Fitting " << this_histname  + "_1" << endl;

    cout << "[Fit_and_Draw] Found the histogram" << endl;
    TH1D *this_hist = nullptr;
    if(data_or_mc == "mc"){
      if(dir.Contains("pion")){
	if(!((TH1D*)gDirectory -> Get(this_histname + "_1"))) continue;
	this_hist = (TH1D*)gDirectory -> Get(this_histname + "_1") -> Clone();

	if(((TH1D*)gDirectory -> Get(this_histname + "_2"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_2"));

	if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
 	if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));

      }
      if(dir.Contains("proton")){
	if(!((TH1D*)gDirectory -> Get(this_histname + "_2"))) continue;
	this_hist = (TH1D*)gDirectory -> Get(this_histname + "_2") -> Clone();

	if(((TH1D*)gDirectory -> Get(this_histname + "_1"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_1"));

	if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));

      }
    }
    else if(data_or_mc == "data"){
      if(!((TH1D*)gDirectory -> Get(this_histname + "_0"))) continue;
      this_hist = (TH1D*)gDirectory -> Get(this_histname + "_0") -> Clone();
    }
    else continue;

    if(this_hist != nullptr) cout << "[Fit_and_Draw] Found the histogram" << endl;
    this_hist -> Rebin(rebin);
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

    f_out -> cd();
    this_hist -> SetName(this_histname);
    this_hist -> Write();
    f_mc -> cd();
    gDirectory -> cd(dir);
  }

  TGraphErrors *Eloss_ratio_gr = new TGraphErrors(N_KE_steps, &KE_vec[0], &mean_vec[0], &KE_err_vec[0], &mean_err_vec[0]);

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);
  TH1D* template_h = new TH1D("", "", 1., KE_low - 2.0 * KE_step, KE_high + 2.0 * KE_step);
  template_h -> SetStats(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle(title_X);
  template_h -> GetYaxis() -> SetTitle(title_Y);
  template_h -> GetYaxis() -> SetRangeUser(-50., 100.);
  template_h -> Draw();

  Eloss_ratio_gr -> SetLineColor(kGreen);
  Eloss_ratio_gr -> SetLineWidth(3);
  Eloss_ratio_gr -> Draw("epsame ");

  TF1 * f_pol_1 = new TF1("f_pol_1", "pol1", KE_low + 0., KE_high + 0.);
  f_pol_1 -> SetParameter(0,  30.);
  f_pol_1 -> SetParameter(1,  0.);
  TF1 * f_pol_2 = new TF1("f_pol_2", "pol2" , KE_low + 0., KE_high + 0.);
  //f_pol_1 -> SetStats(0);
  //f_pol_2 -> SetStats(0);
  Eloss_ratio_gr -> Fit(f_pol_1, "R0", "", KE_low + 0., KE_high + 0.);
  Eloss_ratio_gr -> Fit(f_pol_2, "R0", "", KE_low + 0., KE_high + 0.);

  f_pol_1 -> SetLineColor(kRed);
  f_pol_2 -> SetLineColor(kRed);
  f_pol_1 -> SetLineStyle(7);
  f_pol_2 -> SetLineStyle(7);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) f_pol_1 -> Draw("lsame");
  else f_pol_2 -> Draw("lsame");

  TF1 * f_comparison_pol_2 = new TF1("f_comparison_pol_2", "pol2" , KE_low + 0., KE_high + 0.);
  f_comparison_pol_2 -> SetParameter(0, 3.97800e+01);
  f_comparison_pol_2 -> SetParameter(1, -2.39599e-01);
  f_comparison_pol_2 -> SetParameter(2, 4.49754e-04);
  f_comparison_pol_2 -> SetLineColor(kBlue);
  f_comparison_pol_2 -> SetLineStyle(7);
  f_comparison_pol_2 -> Draw("lsame");

  TString fit_fuction_str = "";
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) fit_fuction_str = "Fitted p_{1}x + p_{0}";
  else fit_fuction_str = "Fitted p_{2} x^{2} + p_{1}x + p_{0}";

  TString fit_result_str = "";
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) fit_result_str = Form("%.2e + %.2e x", f_pol_1 -> GetParameter(0), f_pol_1 -> GetParameter(1));
  else fit_result_str = Form("%.2e + %.2e x + %.2e x^{2}", f_pol_2 -> GetParameter(0), f_pol_2 -> GetParameter(1), f_pol_2 -> GetParameter(2));

  TLegend *l = new TLegend(0.18, 0.7, 0.55, 0.85);
  l -> AddEntry(f_pol_1, fit_fuction_str, "l");
  //l -> AddEntry(f_comparison_pol_2, "#mu + 9.48", "l");
  l -> SetLineColor(kWhite); 
  l -> Draw("same");

  double p0, p1, p2, p0_err, p1_err, p2_err;
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")){
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
  p0_str = Form("p_{0} = %.3e #pm %.2e", p0, p0_err);
  p1_str = Form("p_{1} = %.3e #pm %.2e", p1, p1_err);
  p2_str = Form("p_{2} = %.3e #pm %.2e", p2, p2_err);
  TLatex latex_p0, latex_p1, latex_p2;
  latex_p0.SetNDC();
  latex_p1.SetNDC();
  latex_p2.SetNDC();
  latex_p0.SetTextSize(0.035);
  latex_p1.SetTextSize(0.035);
  latex_p2.SetTextSize(0.035);
  latex_p0.DrawLatex(0.20, 0.66, p0_str);
  latex_p1.DrawLatex(0.20, 0.62, p1_str);
  //if(!histname.Contains("diff_true")) latex_p2.DrawLatex(0.20, 0.58, p2_str);
  if(!(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF"))) latex_p2.DrawLatex(0.20, 0.58, p2_str);


  // ========== For fit up/down shapes ========== //
  TH1D *hint = new TH1D("hint", "Fitted Gaussian with .95 conf.band", 1000., KE_low, KE_high);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
  hint->SetStats(false);
  hint->SetFillColorAlpha(kRed, 0.3);
  hint -> SetMarkerSize(0);
  hint -> SetLineColor(kRed);
  hint->Draw("e3 same");

  TH1D *h_up = new TH1D("h_up", "", 1000., KE_low, KE_high);
  TH1D *h_down = new TH1D("h_down","", 1000., KE_low, KE_high);
  for(int i = 1; i < 1001; i++){
    double this_content = hint -> GetBinContent(i);
    double this_err = hint -> GetBinError(i);
    h_up -> SetBinContent(i, this_content + this_err);
    h_down -> SetBinContent(i, this_content - this_err);
    h_up -> SetBinError(i, (this_content + this_err) / 20.);
    h_down ->  SetBinError(i, (this_content - this_err) / 20.);
  }

  TF1 * h_up_pol_2 = new TF1("h_up_pol_2", "pol2", KE_low + 0., KE_high + 0.);
  h_up -> Fit(h_up_pol_2, "R0", "", KE_low + 0., KE_high + 0.);
  h_up_pol_2 -> SetLineColor(kBlue);
  h_up_pol_2 -> Draw("lsame");

  TF1 * h_down_pol_2 = new TF1("h_down_pol_2", "pol2", KE_low + 0., KE_high + 0.);
  h_down -> Fit(h_down_pol_2, "R0", "", KE_low + 0., KE_high + 0.);
  h_down_pol_2 -> SetLineColor(kCyan);
  h_down_pol_2 ->Draw("lsame");

  TLegend * l_err = new TLegend(0.2, 0.30, 0.92, 0.35);
  l_err -> SetNColumns(2);
  l_err -> AddEntry(h_up_pol_2, "Fitted Up", "l");
  l_err -> AddEntry(h_down_pol_2, "Fitted Down", "l");
  l_err -> Draw("same");

  TString p0_up, p1_up, p2_up;
  p0_up = Form("p_{0} = %.3e", h_up_pol_2 -> GetParameter(0));
  p1_up = Form("p_{1} = %.3e", h_up_pol_2 -> GetParameter(1));
  p2_up = Form("p_{2} = %.3e", h_up_pol_2 -> GetParameter(2));
  TLatex latex_p0_up, latex_p1_up, latex_p2_up;
  latex_p0_up.SetNDC();
  latex_p1_up.SetNDC();
  latex_p2_up.SetNDC();
  latex_p0_up.SetTextSize(0.035);
  latex_p1_up.SetTextSize(0.035);
  latex_p2_up.SetTextSize(0.035);
  latex_p0_up.DrawLatex(0.20, 0.26, p0_up);
  latex_p1_up.DrawLatex(0.20, 0.22, p1_up);
  latex_p2_up.DrawLatex(0.20, 0.18, p2_up);

  TString p0_down, p1_down, p2_down;
  p0_down = Form("p_{0} = %.3e", h_down_pol_2 -> GetParameter(0));
  p1_down = Form("p_{1} = %.3e", h_down_pol_2 -> GetParameter(1));
  p2_down = Form("p_{2} = %.3e", h_down_pol_2 -> GetParameter(2));
  TLatex latex_p0_down, latex_p1_down, latex_p2_down;
  latex_p0_down.SetNDC();
  latex_p1_down.SetNDC();
  latex_p2_down.SetNDC();
  latex_p0_down.SetTextSize(0.035);
  latex_p1_down.SetTextSize(0.035);
  latex_p2_down.SetTextSize(0.035);
  latex_p0_down.DrawLatex(0.60, 0.26, p0_down);
  latex_p1_down.DrawLatex(0.60, 0.22, p1_down);
  latex_p2_down.DrawLatex(0.60, 0.18, p2_down);

  // ================================================== //

  TString latex_sample_str = "";
  if(dir.Contains("pion")) latex_sample_str = "(#pi^{#pm} Elas. & Inel.)";
  if(dir.Contains("proton")) latex_sample_str = "(proton #chi^{2}_{proton} < 10)";
  if(data_or_mc == "mc") latex_sample_str = "MC " + latex_sample_str;
  if(data_or_mc == "data") latex_sample_str = "Data " + latex_sample_str;
  TLatex latex_ProtoDUNE, latex_sample;
  latex_ProtoDUNE.SetNDC();
  latex_sample.SetNDC();
  latex_sample.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_sample.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_sample.DrawLatex(0.95, 0.97, latex_sample_str);

  TLine *this_line = new TLine(0., 1, 50., 1.);
  this_line -> SetLineStyle(5);
  this_line -> SetLineColor(kRed);
  this_line -> Draw("same");

  // == For pol_1 difference
  double pol_1_low  = f_pol_1 -> Eval(KE_low);
  double pol_1_high = f_pol_1 -> Eval(KE_high);
  TLine *pol_1_line_low = new TLine(KE_low - 2.0 * KE_step, pol_1_low, KE_high + 2.0 * KE_step, pol_1_low);
  pol_1_line_low -> SetLineWidth(1);
  pol_1_line_low -> SetLineStyle(5);
  pol_1_line_low -> SetLineColor(kBlack);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) pol_1_line_low -> Draw("same");
  TLine *pol_1_line_high = new TLine(KE_low - 2.0 * KE_step, pol_1_high, KE_high + 2.0 * KE_step, pol_1_high);
  pol_1_line_high -> SetLineWidth(1);
  pol_1_line_high -> SetLineStyle(5);
  pol_1_line_high -> SetLineColor(kBlack);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) pol_1_line_high -> Draw("same");
  double pol_1_diff = fabs(pol_1_high - pol_1_low);

  TString pol_1_diff_str = Form("#Delta = %.2f MeV", pol_1_diff);
  TLatex latex_pol_1_diff;
  latex_pol_1_diff.SetNDC();
  latex_pol_1_diff.SetTextSize(0.03);
  
  // == For pol_2 difference
  double pol_2_diff_low = fabs(f_pol_2 -> Eval(KE_low) - f_comparison_pol_2 -> Eval(KE_low));
  double pol_2_diff_high = fabs(f_pol_2 -> Eval(KE_high) - f_comparison_pol_2 -> Eval(KE_high));
  TString pol_2_diff_low_str = Form("#Delta = %.2f MeV", pol_2_diff_low);
  TString pol_2_diff_high_str = Form("#Delta = %.2f MeV", pol_2_diff_high);
  TLatex latex_pol_2_diff_low, latex_pol_2_diff_high;
  latex_pol_2_diff_low.SetNDC();
  latex_pol_2_diff_low.SetTextSize(0.03);
  latex_pol_2_diff_high.SetNDC();
  latex_pol_2_diff_high.SetTextSize(0.03);
  

  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) latex_pol_1_diff.DrawLatex(0.2, 0.45, pol_1_diff_str);
  else{
    latex_pol_2_diff_low.DrawLatex(0.3, 0.41, pol_2_diff_low_str);
    latex_pol_2_diff_high.DrawLatex(0.8,0.67, pol_2_diff_high_str);
  }

  TString pdfname = "";
  if(dir.Contains("pion")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/pion/Eloss_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
   if(dir.Contains("proton")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/proton/Eloss_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
  c -> SaveAs(pdfname);

  c -> Close();

  f_mc -> Close();
  f_out -> Close();
  KE_vec.clear();
  KE_err_vec.clear();
  mean_vec.clear();
  mean_err_vec.clear();
  std_vec.clear();
}

void Fit_and_Draw_Landau(TString data_or_mc, TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_out = new TFile(root_file_path + "Fitting_Eloss_Landau_" + data_or_mc + "_" + histname + "_" + dir + filename, "RECREATE");
  TFile *f_mc = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> cd(dir);
  vector<double> KE_vec;
  vector<double> KE_err_vec;
  vector<double> mean_vec;
  vector<double> mean_err_vec;
  vector<double> std_vec;

  int KE_low = 700;
  int KE_high = 1100;
  if(dir.Contains("proton")){
    KE_low = 300.;
    KE_high = 600.;
  }
  int KE_step = 50;
  int N_KE_steps = (KE_high - KE_low) / KE_step;

  for(int i = 0; i < N_KE_steps; i++){
    TString KE_range = Form("%dto%dMeV", KE_low + KE_step * i, KE_low + KE_step * (i + 1) );
    TString this_histname = "htrack_KE_diff_" + histname + "_nonscraper_KE_beam_inst" + KE_range + "_" + dir;
    
    cout << "[Fit_and_Draw] Fitting " << this_histname  + "_1" << endl;

    cout << "[Fit_and_Draw] Found the histogram" << endl;
    TH1D *this_hist = nullptr;
    if(data_or_mc == "mc"){
      if(dir.Contains("pion")){
        if(!((TH1D*)gDirectory -> Get(this_histname + "_1"))) continue;
        this_hist = (TH1D*)gDirectory -> Get(this_histname + "_1") -> Clone();
        if(((TH1D*)gDirectory -> Get(this_histname + "_2"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_2"));

        if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));

      }
      if(dir.Contains("proton")){
        if(!((TH1D*)gDirectory -> Get(this_histname + "_2"))) continue;
        this_hist = (TH1D*)gDirectory -> Get(this_histname + "_2") -> Clone();
        //if(((TH1D*)gDirectory -> Get(this_histname + "_1"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_1"));
	/*
        if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
	if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));
	*/
      }
    }
    else if(data_or_mc == "data"){
      if(!((TH1D*)gDirectory -> Get(this_histname + "_0"))) continue;
      this_hist = (TH1D*)gDirectory -> Get(this_histname + "_0") -> Clone();
    }
    else continue;

    if(this_hist != nullptr) cout << "[Fit_and_Draw] Found the histogram" << endl;
    this_hist -> Rebin(rebin);
    double max_x = this_hist -> GetBinCenter(this_hist -> GetMaximumBin());
    double fit_x_min = max_x - width;
    double fit_x_max = max_x + width;
    //  ffit->SetParNames("Width","MPV","Area","GSigma");
    Double_t fitting_range[2];
    fitting_range[0] = 0.;
    fitting_range[1] = 50.;
    Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
    sv[0] = 3.;
    sv[1] = max_x;
    sv[2] = this_hist -> Integral() * 0.05;
    sv[3] = 1.;
    for(int j=0; j<4; ++j){
      pllo[j] = 0.01*sv[j];
      plhi[j] = 100*sv[j];
    }
    Double_t chisqr;
    Int_t    ndf;
    Int_t    status;
    TF1 *this_Langau = langaufit(this_hist, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_data");
    this_Langau -> SetLineColor(kCyan);
    this_Langau -> SetLineStyle(7);
    this_Langau -> SetLineWidth(3);
    
    double this_MPV = this_Langau -> GetParameter(1);
    double this_MPV_err = this_Langau -> GetParError(1);
    mean_vec.push_back(this_MPV);
    std_vec.push_back(this_MPV_err);
    mean_err_vec.push_back(this_MPV_err);

    double this_KE = (KE_low + 0.) + (KE_step + 0.) * (i + 0.) + (KE_step + 0.) / 2.0;
    KE_vec.push_back(this_KE);
    KE_err_vec.push_back((KE_step + 0.) / 2.0);
    cout << "[Fit_and_Draw] this_KE : " << this_KE << " +- " << (KE_step + 0.) / 2.0 << endl;

    f_out -> cd();
    this_hist -> SetName(this_histname);
    this_hist -> Write();
    f_mc -> cd();
    gDirectory -> cd(dir);
  }

  TGraphErrors *Eloss_ratio_gr = new TGraphErrors(N_KE_steps, &KE_vec[0], &mean_vec[0], &KE_err_vec[0], &mean_err_vec[0]);

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);
  TH1D* template_h = new TH1D("", "", 1., KE_low - 2.0 * KE_step, KE_high + 2.0 * KE_step);
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
  f_pol_1 -> SetParameter(0,  30.);
  f_pol_1 -> SetParameter(1,  0.);

  Eloss_ratio_gr -> Fit(f_pol_1, "R0", "", KE_low + 0., KE_high + 0.);
  f_pol_1 -> SetLineColor(kRed);
  f_pol_1 -> SetLineStyle(7);
  f_pol_1 -> Draw("lsame");


  // ========== For fit up/down shapes ========== //
  TH1D *hint = new TH1D("hint", "Fitted Gaussian with .95 conf.band", 1000., KE_low, KE_high);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
  hint->SetStats(false);
  hint->SetFillColorAlpha(kRed, 0.3);
  hint -> SetMarkerSize(0);
  hint -> SetLineColor(kRed);
  hint -> Draw("e3 same"); 

  TString fit_fuction_str = "";
  fit_fuction_str = "Fitted p_{1}x + p_{0}";
  
  TString fit_result_str = "";
  fit_result_str = Form("%.2e + %.2e x", f_pol_1 -> GetParameter(0), f_pol_1 -> GetParameter(1));

  TLegend *l = new TLegend(0.18, 0.7, 0.55, 0.85);
  l -> AddEntry(f_pol_1, fit_fuction_str, "l");
  l -> SetLineColor(kWhite);
  l -> Draw("same");

  double p0, p1, p2, p0_err, p1_err, p2_err;
  p0 = f_pol_1 -> GetParameter(0);
  p1 = f_pol_1 -> GetParameter(1);
  //p2 = f_pol_1 -> GetParameter(2);
  p0_err = f_pol_1 -> GetParError(0);
  p1_err = f_pol_1 -> GetParError(1);
  //p2_err = f_pol_1 -> GetParError(2);
  
  TString p0_str, p1_str, p2_str;
  p0_str = Form("p_{0} = %.3e #pm %.2e", p0, p0_err);
  p1_str = Form("p_{1} = %.3e #pm %.2e", p1, p1_err);
  p2_str = Form("p_{2} = %.3e #pm %.2e", p2, p2_err);
  TLatex latex_p0, latex_p1, latex_p2;
  latex_p0.SetNDC();
  latex_p1.SetNDC();
  latex_p2.SetNDC();
  latex_p0.SetTextSize(0.035);
  latex_p1.SetTextSize(0.035);
  latex_p2.SetTextSize(0.035);
  latex_p0.DrawLatex(0.20, 0.66, p0_str);
  latex_p1.DrawLatex(0.20, 0.62, p1_str);

  TString latex_sample_str = "";
  if(dir.Contains("pion")) latex_sample_str = "(#pi^{#pm} Elas. & Inel.)";
  if(dir.Contains("proton")) latex_sample_str = "(proton #chi^{2}_{proton} < 10)";
  if(data_or_mc == "mc") latex_sample_str = "MC " + latex_sample_str;
  if(data_or_mc == "data") latex_sample_str = "Data " + latex_sample_str;
  TLatex latex_ProtoDUNE, latex_sample;
  latex_ProtoDUNE.SetNDC();
  latex_sample.SetNDC();
  latex_sample.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_sample.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_sample.DrawLatex(0.95, 0.97, latex_sample_str);

  TLine *this_line = new TLine(0., 1, 50., 1.);
  this_line -> SetLineStyle(5);
  this_line -> SetLineColor(kRed);
  this_line -> Draw("same");

  // == For pol_1 difference
  double pol_1_low  = f_pol_1 -> Eval(KE_low);
  double pol_1_high = f_pol_1 -> Eval(KE_high);
  TLine *pol_1_line_low = new TLine(KE_low - 2.0 * KE_step, pol_1_low, KE_high + 2.0 * KE_step, pol_1_low);
  pol_1_line_low -> SetLineWidth(1);
  pol_1_line_low -> SetLineStyle(5);
  pol_1_line_low -> SetLineColor(kBlack);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) pol_1_line_low -> Draw("same");
  TLine *pol_1_line_high = new TLine(KE_low - 2.0 * KE_step, pol_1_high, KE_high + 2.0 * KE_step, pol_1_high);
  pol_1_line_high -> SetLineWidth(1);
  pol_1_line_high -> SetLineStyle(5);
  pol_1_line_high -> SetLineColor(kBlack);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) pol_1_line_high -> Draw("same");
  double pol_1_diff = fabs(pol_1_high - pol_1_low);

  TString pol_1_diff_str = Form("#Delta = %.2f MeV", pol_1_diff);
  TLatex latex_pol_1_diff;
  latex_pol_1_diff.SetNDC();
  latex_pol_1_diff.SetTextSize(0.03);
  if(histname.Contains("TrueBeam_TrueFF") || histname.Contains("TrueBeam_FittedFF")) latex_pol_1_diff.DrawLatex(0.2, 0.45, pol_1_diff_str); 

  TString pdfname = "";
  if(dir.Contains("pion")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/pion/Eloss_Landau_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
  if(dir.Contains("proton")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/proton/Eloss_Landau_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
  c -> SaveAs(pdfname);

  c -> Close();

  f_mc -> Close();
  f_out -> Close();
  KE_vec.clear();
  KE_err_vec.clear();
  mean_vec.clear();
  mean_err_vec.clear();
  std_vec.clear();

}

void Verify_and_Draw(TString data_or_mc, TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_out = new TFile(root_file_path + "Fitting_Eloss_" + data_or_mc + "_" + histname + "_" + dir + filename, "RECREATE");
  TFile *f_mc = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> cd(dir);
  vector<double> KE_vec;
  vector<double> KE_err_vec;
  vector<double> mean_vec;
  vector<double> mean_err_vec;
  vector<double> std_vec;

  int KE_low = 700;
  int KE_high = 1100;
  if(dir.Contains("proton")){
    KE_low = 300.;
    KE_high = 600.;
  }
  int KE_step = 50;
  int N_KE_steps = (KE_high - KE_low) / KE_step;
  for(int i = 0; i < N_KE_steps; i++){
    TString KE_range = Form("%dto%dMeV", KE_low + KE_step * i, KE_low + KE_step * (i + 1) );
    TString this_histname = "htrack_KE_diff_" + histname + "_nonscraper_KE_beam_inst" + KE_range + "_" + dir;
    if(histname.Contains("KE_ff_reco")) this_histname = "htrack_diff_" + histname + "_KE_beam_inst" + KE_range + "_" + dir;
    cout << "[Fit_and_Draw] Fitting " << this_histname  + "_1" << endl;

    cout << "[Fit_and_Draw] Found the histogram" << endl;
    TH1D *this_hist = nullptr;
    if(data_or_mc == "mc"){
      if(dir.Contains("pion")){
        if(!((TH1D*)gDirectory -> Get(this_histname + "_1"))) continue;
        this_hist = (TH1D*)gDirectory -> Get(this_histname + "_1") -> Clone();

        if(((TH1D*)gDirectory -> Get(this_histname + "_2"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_2"));

        if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));

      }
      if(dir.Contains("proton")){
        if(!((TH1D*)gDirectory -> Get(this_histname + "_2"))) continue;
        this_hist = (TH1D*)gDirectory -> Get(this_histname + "_2") -> Clone();

        //if(((TH1D*)gDirectory -> Get(this_histname + "_1"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_1"));
	/*
        if(((TH1D*)gDirectory -> Get(this_histname + "_3"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_3"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_4"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_4"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_5"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_5"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_6"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_6"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_7"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_7"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_8"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_8"));
        if(((TH1D*)gDirectory -> Get(this_histname + "_9"))) this_hist -> Add((TH1D*)gDirectory -> Get(this_histname + "_9"));
	*/
      }
    }
    else if(data_or_mc == "data"){
      if(!((TH1D*)gDirectory -> Get(this_histname + "_0"))) continue;
      this_hist = (TH1D*)gDirectory -> Get(this_histname + "_0") -> Clone();
    }
    else continue;

    if(this_hist != nullptr) cout << "[Fit_and_Draw] Found the histogram" << endl;
    this_hist -> Rebin(rebin);
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

    f_out -> cd();
    this_hist -> SetName(this_histname);
    this_hist -> Write();
    f_mc -> cd();
    gDirectory -> cd(dir);
  }

  TGraphErrors *Eloss_ratio_gr = new TGraphErrors(N_KE_steps, &KE_vec[0], &mean_vec[0], &KE_err_vec[0], &mean_err_vec[0]);

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);
  TH1D* template_h = new TH1D("", "", 1., KE_low - 2.0 * KE_step, KE_high + 2.0 * KE_step);
  template_h -> SetStats(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle(title_X);
  template_h -> GetYaxis() -> SetTitle(title_Y);
  template_h -> GetYaxis() -> SetRangeUser(-50., 100.);
  template_h -> Draw();

  Eloss_ratio_gr -> SetLineColor(kGreen);
  Eloss_ratio_gr -> SetLineWidth(3);
  Eloss_ratio_gr -> Draw("epsame ");


  TString latex_sample_str = "";
  if(dir.Contains("pion")) latex_sample_str = "(#pi^{#pm} Elas. & Inel.)";
  if(dir.Contains("proton")) latex_sample_str = "(proton #chi^{2}_{proton} < 10)";
  if(data_or_mc == "mc") latex_sample_str = "MC " + latex_sample_str;
  if(data_or_mc == "data") latex_sample_str = "Data " + latex_sample_str;
  TLatex latex_ProtoDUNE, latex_sample;
  latex_ProtoDUNE.SetNDC();
  latex_sample.SetNDC();
  latex_sample.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_sample.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_sample.DrawLatex(0.95, 0.97, latex_sample_str);

  TLine *this_line = new TLine(0., 1, 50., 1.);
  this_line -> SetLineStyle(5);
  this_line -> SetLineColor(kRed);
  this_line -> Draw("same");

  
  TString pdfname = "";
  if(dir.Contains("pion")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/pion/Eloss_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
  if(dir.Contains("proton")) pdfname = "./output/plots/BeamStudy/Upstream_Eloss/proton/Eloss_" + data_or_mc + "_" + histname + "_" + dir + ".pdf";
  c -> SaveAs(pdfname);

  c -> Close();

  f_mc -> Close();
  f_out -> Close();
  KE_vec.clear();
  KE_err_vec.clear();
  mean_vec.clear();
  mean_err_vec.clear();
  std_vec.clear();



}





void Run_Fit_and_Draw(TString data_or_mc, TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString particle_str_arr[2] = {"proton", "pion"};
  int N_particles = 2;
  if(histname.Contains("Fitted")) N_particles = 1;
  TString cutflow_arr[1] = {"BeamWindow"};
  for(int i = 0; i < N_particles; i++){
    for(int j = 0; j < 1; j++){
      TString this_dir = particle_str_arr[i] + "_" + cutflow_arr[j] + "_" + dir;
      Fit_and_Draw(data_or_mc, filename, this_dir, histname, title_X, title_Y, width, rebin);
      if(histname == "TrueBeam_TrueFF") Fit_and_Draw_Landau(data_or_mc, filename, this_dir, histname, title_X, title_Y, width, rebin);
    }
  }
}

void Run_Verify_and_Draw(TString data_or_mc, TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString particle_str_arr[2] = {"proton", "pion"};
  int N_particles = 2;
  if(histname.Contains("Fitted")) N_particles = 1;
  TString cutflow_arr[2] = {"BeamWindow", "BeamScraper"};
  for(int i = 0; i < N_particles; i++){
    for(int j = 0; j < 2; j++){
      TString this_dir = particle_str_arr[i] + "_" + cutflow_arr[j] + "_" + dir;
      Verify_and_Draw(data_or_mc, filename, this_dir, histname, title_X, title_Y, width, rebin);
      //if(histname == "TrueBeam_TrueFF") Fit_and_Draw_Landau(data_or_mc, filename, this_dir, histname, title_X, title_Y, width, rebin);
    }
  }
}

void Fitting_Eloss_range(){
  setTDRStyle();

  // === For study
  /*
  Run_Fit_and_Draw("mc", "_Beam_Study_1.0GeV.root", "noweight", "RecoBeam_TrueFF", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.} - KE_{ff}^{true}) [MeV]", 35., 4.);
  Run_Fit_and_Draw("mc", "_Beam_Study_1.0GeV.root", "noweight", "RecoBeam_FittedFF", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.} - KE_{ff}^{fitted}) [MeV]", 35., 4.);
  Run_Fit_and_Draw("mc", "_Beam_Study_1.0GeV.root", "noweight", "TrueBeam_TrueFF", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.}^{true} - KE_{ff}^{true}) [MeV]", 35., 5.);
  Run_Fit_and_Draw("data", "_Beam_Study_1.0GeV.root", "noweight", "RecoBeam_FittedFF", "KE_{Beam Inst.} [MeV]", "#mu (KE_{Beam Inst.} - KE_{ff}^{fitted}) [MeV]", 35., 4.);
  */
  // == After study, for performance validation
  Run_Verify_and_Draw("mc", "_Beam_Study_1.0GeV.root", "noweight", "KE_ff_reco_ElasTrue_KE_ff_true", "KE_{Beam Inst.} [MeV]", "#mu (KE_{ff}^{reco} - KE_{ff}^{true}) [MeV]", 35., 5.);
  Run_Verify_and_Draw("mc", "_Beam_Study_1.0GeV.root", "noweight", "KE_ff_reco_AllTrue_KE_ff_true", "KE_{Beam Inst.} [MeV]", "#mu (KE_{ff}^{reco} - KE_{ff}^{true}) [MeV]", 35., 4.);
  
}
