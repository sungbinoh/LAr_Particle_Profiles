#include "canvas_margin.h"
#include "mylib.h"

Double_t SB_Dist(Double_t *x, Double_t *par){
  Double_t mean = par[0];
  Double_t height = par[1];
  Double_t sigma = par[2];
  Double_t beta = par[3];
  Double_t a = par[4];
  Double_t b = par[5];

  Double_t arg = x[0]-mean;
  Double_t numer = height;
  Double_t denom = 1. + pow( fabs(arg) / sigma, beta);

  Double_t func = numer / denom + a * x[0] + b;
  return func;
}

Double_t SB_Dist_test(Double_t *x, Double_t *par){
  Double_t mean = par[0];
  Double_t height = par[1];
  Double_t sigma = par[2];
  Double_t beta = par[3];
  Double_t a0 = par[4];
  Double_t a1 = par[5];
  Double_t a2 = par[6];

  Double_t arg = x[0]-mean;
  Double_t numer = height;
  Double_t denom = 1. + pow( fabs(arg) / sigma, beta);

  Double_t func = numer / denom + a0 * exp( (-1.) * a1 * (x[0] - a2));
  return func;
}

Double_t SB_Dist_signal(Double_t *x, Double_t *par){
  Double_t mean = par[0];
  Double_t height = par[1];
  Double_t sigma = par[2];
  Double_t beta = par[3];
  Double_t a0 = par[4];
  Double_t a1 = par[5];
  Double_t a2 = par[6];

  Double_t arg = x[0]-mean;
  Double_t numer = height;
  Double_t denom = 1. + pow( fabs(arg) / sigma, beta);


  Double_t func = numer / denom;
   return func;
}

Double_t SB_Dist_bkg(Double_t *x, Double_t *par){
  Double_t mean = par[0];
  Double_t height = par[1];
  Double_t sigma = par[2];
  Double_t beta = par[3];
  Double_t a = par[4];
  Double_t b = par[5];
  Double_t a2 = par[6];

  Double_t func = a * x[0] + b;
  return func;
}

void Fit_true_deltaEQE(TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  //TFile *f_out = new TFile(root_file_path + "/Fitting_results/Fitting_" + KE_or_P + "loss_" + data_or_mc + "_" + histname + "_" + dir + filename, "RECREATE");
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  ofstream file_integ(input_file_dir + "/output/txt/True_study/EQE_fit/Signal_yield_" + dir + ".txt");
  cout << "[Fit_true_deltaEQE] Open mc" + filename << endl;
  gDirectory -> cd(dir);

  int KE_low = 100;
  int KE_high = 1100;
  int KE_step = 50.;
  int N_KE_steps = (KE_high - KE_low) / KE_step;
  for(int i = 0; i < N_KE_steps; i++){
    TString KE_range = Form("%dto%dMeV", KE_low + KE_step * i, KE_low + KE_step * (i + 1) );
    TString this_histname = "hdaughter_" + histname + "_KE_in" + KE_range + "_" + dir;
    double this_KE_low = KE_low + KE_step * i;
    double this_KE_high = KE_low + KE_step * (i + 1);

    TH1D *this_hist = nullptr;
    if(!((TH1D*)gDirectory -> Get(this_histname + "_0"))) continue;
    this_hist = (TH1D*)gDirectory -> Get(this_histname + "_0") -> Clone();
  
    if(this_hist == nullptr) continue;
    else{
      cout << "[Fit_true_deltaEQE] Found " << this_histname << endl;
    }
    this_hist -> Rebin(rebin);

    int max_bin = this_hist -> GetMaximumBin();
    double max_x = this_hist -> GetBinCenter(max_bin);
    double fit_x_min = 0.3 * (0. - this_KE_low);
    double fit_x_max = 0.5 * this_KE_low;
    //TF1 *this_f = new TF1("this_f", SB_Dist, fit_x_min, fit_x_max, 9);
    TF1 *this_f = new TF1("this_f", SB_Dist, fit_x_min, fit_x_max, 6);

    this_f->SetParameters(max_x, this_hist -> GetMaximum(), 2., 1.);
    this_f->SetParNames("Mean","Height", "Sigma", "Beta", "a0", "a1", "a2", "a3", "a4");
    this_f -> SetParLimits(0, -5., 5.);
    this_f -> SetParLimits(3, 0., 10.);
    //this_f -> SetParLimits(4, 0., 100.);
    //this_f -> SetParLimits(6, 0.1, 10.);

    this_hist -> Fit(this_f, "R0", "", fit_x_min, fit_x_max);

    TCanvas *c_temp = new TCanvas("", "", 600, 600);
    canvas_margin(c_temp);
    gStyle -> SetOptStat(1111);
    gStyle->SetOptFit(0);
    TH1D *h_temp = new TH1D("", "", 1., -1000., 1000.);
    h_temp -> SetStats(0);
    gStyle -> SetOptTitle(0);
    gStyle -> SetLineWidth(2);
    h_temp -> GetYaxis() -> SetRangeUser(this_hist -> GetMaximum() * (-0.2), this_hist -> GetMaximum() * 1.5);
    h_temp -> GetYaxis() -> SetLabelSize(0.035);
    h_temp -> GetXaxis() -> SetLabelSize(0.035);
    h_temp -> GetXaxis() -> SetTitle("#Delta E_{QE} [MeV]");
    h_temp -> Draw();
    this_hist -> SetLineColor(kBlack);
    this_hist -> SetMarkerColor(kBlack);
    this_hist -> Draw("e0same");
    
    this_f -> SetLineColor(kOrange-3);
    //this_f -> SetLineStyle(2);
    this_f -> SetLineWidth(3);
    this_f -> Draw("lsame");
    
    double par_out[6];
    this_f -> GetParameters(par_out);
    TF1 *f_signal = new TF1("f_signal", SB_Dist_signal, -500., 500., 6);
    TF1 *f_bkg = new TF1("f_bkg", SB_Dist_bkg, fit_x_min, fit_x_max, 6);
    f_signal -> SetParameters(par_out);
    f_bkg -> SetParameters(par_out);

    f_signal -> SetLineColor(kBlue);
    f_bkg -> SetLineColor(kRed);
    f_signal -> SetLineStyle(2);
    f_bkg -> SetLineStyle(2);
    f_signal -> Draw("lsame");
    f_bkg -> Draw("lsame");
    this_hist -> Draw("e0same");

    TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9);
    l -> AddEntry(this_f, "Global Fit", "l");
    l -> AddEntry(f_signal, "Signal Fit", "l");
    l -> AddEntry(f_bkg, "Bkg. Fit", "l");
    l -> Draw("same");
    
    vector<TString> latex_str;
    latex_str.push_back("f(x) = #frac{a}{1 + ( |x - #mu| / #sigma)^{#beta} } + a_{1}x + a_{0}");
    latex_str.push_back(Form("#mu = %.2f", par_out[0]));
    latex_str.push_back(Form("a = %.2f", par_out[1]));
    latex_str.push_back(Form("#sigma = %.2f", par_out[2]));
    latex_str.push_back(Form("#beta = %.2f", par_out[3]));
    latex_str.push_back(Form("a_{1} = %.2f", par_out[4]));
    latex_str.push_back(Form("a_{0} = %.2f", par_out[5]));
    //latex_str.push_back(Form("a2 = %.2f", par_out[6]));

    vector<TLatex> Latex_vec;
    TLatex tex_1, tex_2, tex_3, tex_4, tex_5, tex_6, tex_7, tex_8;
    Latex_vec.push_back(tex_1);
    Latex_vec.push_back(tex_2);
    Latex_vec.push_back(tex_3);
    Latex_vec.push_back(tex_4);
    Latex_vec.push_back(tex_5);
    Latex_vec.push_back(tex_6);
    Latex_vec.push_back(tex_7);
    //Latex_vec.push_back(tex_8);

    for(unsigned int i = 0; i < Latex_vec.size(); i++){
      Latex_vec.at(i).SetNDC();
      Latex_vec.at(i).SetTextSize(0.035);
      double this_y = 0.80 - (i + 0.) * 0.04;
      if(i == 0) this_y = 0.85; 
      Latex_vec.at(i).DrawLatex(0.20, this_y, latex_str.at(i));
    }

    TLatex latex_ProtoDUNE, latex_sample;
    latex_ProtoDUNE.SetNDC();
    latex_sample.SetNDC();
    latex_sample.SetTextAlign(31);
    latex_ProtoDUNE.SetTextSize(0.03);
    latex_sample.SetTextSize(0.03);
    latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_sample.DrawLatex(0.95, 0.96, "KE_{in} " + KE_range);

    TString pdfname_temp = "./output/plots/PionXsec/True_study/EQE_fit/EQE_fit_" + dir + "_KE_in" + KE_range + ".pdf";
    c_temp -> SaveAs(pdfname_temp);
    c_temp -> Close();

    // == Integrals
    TAxis *axis = this_hist->GetXaxis();
    double bin_width = axis->GetBinWidth(10);
    //double signal_integ = f_signal -> Integral(-500., 500.);
    double global_integ = 0.;//this_f -> Integral(fit_x_min, fit_x_max) / bin_width;
    double signal_integ = 0.;//f_signal -> Integral(fit_x_min, fit_x_max) / bin_width;

    double bkg_integ = 0.;
    int integ_steps = 100000;
    double step = (fit_x_max - fit_x_min) / (integ_steps + 0.);
    for(int j = 0; j < integ_steps; j++){
      double this_x = fit_x_min + step * (i + 0.);
      double global_y = this_f -> Eval(this_x);
      double signal_y = f_signal -> Eval(this_x);
      double bkg_y = f_bkg -> Eval(this_x);
      double global_part = global_y * step;
      double signal_part = signal_y * step;
      double bkg_part = bkg_y * step;
      global_integ += global_part;
      signal_integ += signal_part;
      if(bkg_y > 0.) bkg_integ += bkg_part;
    }
    global_integ = global_integ / bin_width;
    signal_integ = signal_integ / bin_width;
    bkg_integ = bkg_integ / bin_width;
    double hist_integ = 0.;
    int bin_min = axis->FindBin(fit_x_min);
    int bin_max = axis->FindBin(fit_x_max);
    hist_integ = this_hist -> Integral(bin_min, bin_max);
    hist_integ -= this_hist -> GetBinContent(bin_min) * (fit_x_min - axis->GetBinLowEdge(bin_min) ) / axis->GetBinWidth(bin_min);
    hist_integ -= this_hist -> GetBinContent(bin_max) * (axis->GetBinUpEdge(bin_min) - fit_x_max) / axis->GetBinWidth(bin_min);

    file_integ << KE_range << "\tbin_width\t" << bin_width << "\tglobal_integ\t" << global_integ << "\tsignal_integ\t" << signal_integ << "\tbkg_integ\t" << bkg_integ << "\thist_integ\t" << hist_integ << "\thist_integ-bkg_integ\t" << hist_integ - bkg_integ << endl;  

  }

  file_integ.close();
}

void Run_Fit_true_deltaEQE(TString filename, TString dir, TString histname, TString title_X, TString title_Y, double width, double rebin){

  TString more_dir[14] = {"", "nocut",
			"pion_APA3", "pion_BeamScraper", "pion_BeamWindow", "pion_MichelScore", "pion_ProtonVeto",
			"true_gaus_fit_pion_mEQEcut", "true_gaus_fit_pion", "true_likelihood_fit_pion_mEQEcut", "true_likelihood_fit_pion",
			"true_stop_pion_KE400", "true_stop_pion_mEQEcut", "true_stop_pion"};
  for(int i = 1; i < 14; i++){
    TString this_dir = more_dir[i] + "_" + dir;
    Fit_true_deltaEQE(filename, this_dir, histname, title_X, title_Y, width, rebin);
  }
}

void Fitting_deltaEQE(){
  setTDRStyle();

  Run_Fit_true_deltaEQE("_PionXsec_1.0GeV_true_test.root", "noweight", "best_pion_EQE_NC_4", "P_{spec.} [MeV]", "#mu (P_{ff}^{reco} - P_{ff}^{true}) [MeV]", 100., 20.);
  
}
