void Fit_and_Draw(TString filename, double width){
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/";

  TFile *f_output = new TFile(root_file_path + "Fit" + filename, "RECREATE");

  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  vector<double> Eloss;
  vector<double> Eloss_err;
  vector<double> mean_vec_mc;
  vector<double> std_vec_mc;
  vector<double> mean_vec_data;
  vector<double> std_vec_data;
  for(int i = 0; i < 50; i++){
    TString this_hist_name = Form("htrack_length_ratio_eloss%dMeV", i);
    TH1D *this_hist = (TH1D*)gDirectory -> Get(this_hist_name);

    int max_bin = this_hist -> GetMaximumBin();
    double max_x = this_hist -> GetBinCenter(max_bin);
    double fit_x_min = max_x - width;
    double fit_x_max = max_x + width;
    TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
    this_hist -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
    double this_mean = this_gaus -> GetParameter(1);
    double this_std = this_gaus -> GetParameter(2);
    mean_vec_mc.push_back(this_mean);
    std_vec_mc.push_back(this_std);

    f_output -> cd();
    this_hist -> SetName(this_hist_name + "_MC");
    this_hist -> Write();
    f_mc -> cd();

    double this_Eloss = i + 0.;
    Eloss.push_back(this_Eloss);
    Eloss_err.push_back(0.5);
  }

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  //f_data -> cd();
 for(int i = 0; i < 50; i++){
    TString this_hist_name = Form("htrack_length_ratio_eloss%dMeV", i);
    TH1D *this_hist = (TH1D*)gDirectory -> Get(this_hist_name);

    int max_bin = this_hist -> GetMaximumBin();
    double max_x = this_hist -> GetBinCenter(max_bin);
    double fit_x_min = max_x - width;
    double fit_x_max = max_x + width;
    TF1 *this_gaus = new TF1("fit_gaus", "gaus", fit_x_min, fit_x_max);
    this_hist -> Fit(this_gaus, "R", "", fit_x_min, fit_x_max);
    double this_mean = this_gaus -> GetParameter(1);
    double this_std = this_gaus -> GetParameter(2);
    mean_vec_data.push_back(this_mean);
    std_vec_data.push_back(this_std);
 
    f_output ->cd();
    this_hist -> SetName(this_hist_name + "_Data");
    this_hist -> Write();
    f_data -> cd();
 }

  TGraphErrors *Eloss_ratio_mc_gr = new TGraphErrors(50, &Eloss[0], &mean_vec_mc[0], &Eloss_err[0], &std_vec_mc[0]);
  TGraphErrors *Eloss_ratio_data_gr = new TGraphErrors(50, &Eloss[0], &mean_vec_data[0], &Eloss_err[0], &std_vec_data[0]);

  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  TH1D* template_h = new TH1D("", "", 1., 0., 50.);
  template_h ->GetXaxis() -> SetTitle("E_{loss} [MeV]");
  template_h ->GetYaxis() -> SetTitle("Fitted <L_{Reco} / KEtoRange(E_{beam} - E_{loss})>");
  template_h -> GetYaxis() -> SetRangeUser(0.8, 1.2);
  template_h -> Draw();

  Eloss_ratio_mc_gr -> SetLineColor(kGreen);
  Eloss_ratio_mc_gr -> SetLineWidth(3);
  Eloss_ratio_data_gr -> SetLineColor(kBlue);
  Eloss_ratio_data_gr -> SetLineWidth(1);
  //Eloss_ratio_mc_gr -> SetLineStyle(5);
  //Eloss_ratio_data_gr -> SetLineStyle(3);
  Eloss_ratio_mc_gr -> Draw("epsame");
  Eloss_ratio_data_gr -> Draw("epsame");

  TLegend *l = new TLegend(0.1, 0.7, 0.4, 0.9);
  l -> AddEntry(Eloss_ratio_mc_gr, "1 GeV MC #pm #sigma", "l");
  //l -> SetLineColor(kWhite);
  l -> AddEntry(Eloss_ratio_data_gr, "1 GeV Data #pm #sigma", "l");
  l -> Draw("same");

  TLine *this_line = new TLine(0., 1, 50., 1.);
  this_line -> SetLineStyle(5);
  this_line -> SetLineColor(kRed);
  this_line -> Draw("same");

  TLine *best_line = new TLine(23., 0.8, 23., 1.2); 
  best_line -> SetLineStyle(5);
  best_line -> SetLineColor(kRed);
  best_line -> Draw("same");

  c -> SaveAs("./output/plots/Eloss_scan_proton_1GeV.pdf");
}


void Fitting_Eloss_range(){
  
  Fit_and_Draw("penergy_Eloss_scan.root", 0.08);

}
