void Draw_fits_and_BetheBloch(TString filename, TString event, double r_min, double r_max, double max_dEdx){
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/";

  TFile *f_profile = new TFile(root_file_path + "hists.root");
  TGraph *this_profile = (TGraph*)gDirectory->Get("pion_range_vs_dEdx");
  
  TFile *f_fit = new TFile(root_file_path + filename);
  
  TString event_str = event;
  TGraph *dEdx_gr = (TGraph*) gDirectory -> Get("dEdx_" + event_str);
  TGraph *dEdx_bestfit_gr = (TGraph*) gDirectory -> Get("dEdx_bestfit_"+ event_str);

  TCanvas *c = new TCanvas("", "", 1000, 600);
  gStyle->SetOptStat(0);
  TH1D* template_h = new TH1D("", "", 1., r_min, r_max);
  template_h ->GetXaxis() -> SetTitle("Range [cm]");
  template_h ->GetYaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetYaxis() -> SetRangeUser(0., max_dEdx);
  template_h -> Draw();

  this_profile -> Draw("same");

  dEdx_gr -> SetLineColor(kBlue - 8);
  dEdx_gr -> SetLineStyle(7);
  dEdx_gr -> Draw("same");
  
  dEdx_bestfit_gr -> SetLineColor(kBlue);
  dEdx_bestfit_gr -> Draw("same");

  TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.9);
  l -> AddEntry(dEdx_gr, "Prefit dE/dx", "l");
  l -> AddEntry(dEdx_bestfit_gr, "Best-fit dE/dx", "l");
  l -> AddEntry(this_profile, "Bethe-Bloch dE/dx", "l");
  l -> Draw("same");

  c -> SaveAs("./output/plots/" + event_str + ".pdf");
}


void Draw_bestfit(){
  
  Draw_fits_and_BetheBloch("effval_0p5_nofit.root", "Run46728044_Evt14_PID211", 0., 200., 10.);
  Draw_fits_and_BetheBloch("effval_0p5_nofit.root", "Run46700441_Evt215_PID211", 0., 200., 10.);
  Draw_fits_and_BetheBloch("effval_0p5_nofit.root", "Run22612262_Evt443_PID211", 0., 200., 10.);
  Draw_fits_and_BetheBloch("effval_0p3.root", "Run46732160_Evt191_PID211", 0., 30., 20.);
  Draw_fits_and_BetheBloch("effval_0p5_michell0p5_daughter_fit.root", "Run46726942_Evt123_Nhit20", 0., 50., 10.);

}
