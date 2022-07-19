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

void Fit_Gaus_MC_only(TString filename, TString histname, TString TitleX, double xmin, double xmax, double width){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TH1D *hist_mc = (TH1D*)gDirectory -> Get("htrack_" + histname + "_0");

  TCanvas *c = new TCanvas("", "", 600, 800);
  c -> SetRightMargin( 0.03 );
  gStyle->SetOptStat(0);

  TH1D* template_h = new TH1D("", "", 1., xmin, xmax);
  double max_x = hist_mc -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h ->GetXaxis() -> SetTitle(TitleX);
  template_h ->GetXaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetTitle("Events");
  template_h -> Draw();

  hist_mc -> SetLineColor(kGreen);
  hist_mc -> Draw("histsame");
  int max_bin_x_mc = hist_mc -> GetMaximumBin();
  double max_x_mc = hist_mc -> GetBinCenter(max_bin_x_mc);
  TF1 *mc_gaus = new TF1("mc_gaus", "gaus", max_x_mc - width, max_x_mc + width);
  mc_gaus -> SetLineColor(kSpring+3);
  mc_gaus -> SetLineStyle(7);
  mc_gaus -> SetLineWidth(3);
  hist_mc -> Fit(mc_gaus, "R", "", max_x_mc - width, max_x_mc + width);
  mc_gaus -> Draw("lsame");
  
  TLegend *l = new TLegend(0.2, 0.75, 0.5, 0.9);
  l -> AddEntry(hist_mc, "MC", "l");
  l -> Draw("same");

  TPaveText *pt_mc = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  pt_mc -> AddText(Form("#mu : %.3f, #sigma : %.3f", mc_gaus -> GetParameter(1), mc_gaus -> GetParameter(2)));
  ((TText*)pt_mc->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_mc -> SetFillStyle(4000);
  pt_mc -> SetFillColor(kWhite);
  
  pt_mc -> Draw("same");
  c -> SaveAs("./output/" + histname + "_mc_only.pdf");

}

void Fit_Gaus_MC_and_Data(TString filename, TString histname, TString TitleX, double xmin, double xmax, double width){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TH1D *hist_mc = (TH1D*)gDirectory -> Get("htrack_" + histname + "_0");

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  TH1D *hist_data = (TH1D*)gDirectory -> Get("htrack_" + histname + "_0");

  TCanvas *c = new TCanvas("", "", 800, 600);
  c -> SetRightMargin( 0.03 );
  gStyle->SetOptStat(0);

  TH1D* template_h = new TH1D("", "", 1., xmin, xmax);
  double max_x = hist_mc -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h ->GetXaxis() -> SetTitle(TitleX);
  template_h ->GetXaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetTitle("Events");
  template_h -> Draw();

  hist_mc -> SetLineColor(kGreen);
  hist_data -> SetLineColor(kBlue);
  hist_mc -> Draw("histsame");
  hist_data -> Draw("histsame");
  int max_bin_x_mc = hist_mc -> GetMaximumBin();
  double max_x_mc = hist_mc -> GetBinCenter(max_bin_x_mc);
  int max_bin_x_data = hist_data -> GetMaximumBin();
  double max_x_data = hist_data -> GetBinCenter(max_bin_x_data);
  TF1 *mc_gaus = new TF1("mc_gaus", "gaus", max_x_mc - width, max_x_mc + width);
  TF1 *data_gaus = new TF1("data_gaus", "gaus", max_x_data - width, max_x_data + width);
  mc_gaus -> SetLineColor(kSpring+3);
  data_gaus ->SetLineColor(kCyan);
  mc_gaus -> SetLineStyle(7);
  data_gaus-> SetLineStyle(7);
  mc_gaus -> SetLineWidth(3);
  data_gaus-> SetLineWidth(3);
  hist_mc -> Fit(mc_gaus, "R", "", max_x_mc - width, max_x_mc + width);
  hist_data -> Fit(data_gaus, "R", "", max_x_data - width, max_x_data + width);
  mc_gaus -> Draw("lsame");
  data_gaus -> Draw("lsame");
  
  TLegend *l = new TLegend(0.2, 0.7, 0.4, 0.9);
  l -> AddEntry(hist_data, "Data", "l");
  l -> AddEntry(hist_mc, "MC", "l");
  l -> Draw("same");

  TPaveText *pt_mc = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  TPaveText *pt_data = new TPaveText(0.17, 0.37, 0.50, 0.53,"NDC");
  pt_mc -> AddText(Form("#mu : %.3f, #sigma : %.3f", mc_gaus -> GetParameter(1), mc_gaus -> GetParameter(2)));
  ((TText*)pt_mc->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_data -> AddText(Form("#mu : %.3f, #sigma : %.3f", data_gaus -> GetParameter(1), data_gaus -> GetParameter(2)));
  ((TText*)pt_data->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_mc -> SetFillStyle(4000);
  pt_data -> SetFillStyle(4000);
  pt_mc -> SetFillColor(kWhite);
  pt_data -> SetFillColor(kWhite);

  pt_mc -> Draw("same");
  pt_data -> Draw("same");
  c -> SaveAs("./output/" + histname + "_mc_vs_data.pdf");

}

void Fit_and_Draw(TString filename, double width){
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";

  //TFile *f_output = new TFile(root_file_path + "Fit" + filename, "RECREATE");

  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TH1D *BeamKEtoRange_mc = (TH1D*)gDirectory -> Get("htrack_length_BeamKEtoRange_0");
  TH1D *fitted_range_mc = (TH1D*)gDirectory ->Get("htrack_length_fitted_0");
  TH1D *BeamP_mc = (TH1D*)gDirectory ->Get("htrack_BeamKE_0");
  TH1D *fittedP_mc = (TH1D*)gDirectory ->Get("htrack_fittedKE_0");
  TH1D *fitted_dP_mc = (TH1D*)gDirectory ->Get("htrack_fitted_dKE_0");

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  TH1D *BeamKEtoRange_data = (TH1D*)gDirectory ->Get("htrack_length_BeamKEtoRange_0");
  TH1D *fitted_range_data = (TH1D*)gDirectory ->Get("htrack_length_fitted_0");
  TH1D *BeamP_data = (TH1D*)gDirectory ->Get("htrack_BeamKE_0");
  TH1D *fittedP_data = (TH1D*)gDirectory ->Get("htrack_fittedKE_0");
  TH1D *fitted_dP_data = (TH1D*)gDirectory ->Get("htrack_fitted_dKE_0");

  TCanvas *c = new TCanvas("", "", 600, 800);
  c -> SetRightMargin( 0.03 );
  gStyle->SetOptStat(0);
  
  TH1D* template_h = new TH1D("", "", 1., 0., 200.);
  double max_x = BeamKEtoRange_mc -> GetMaximum();
  //template_h -> GetXaxis() -> SetRangeUser(0., 200.);
  template_h -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h ->GetXaxis() -> SetTitle("Range [cm]");
  template_h ->GetXaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetLabelSize(0.03);
  template_h ->GetYaxis() -> SetTitle("Events");
  template_h -> Draw();

  BeamKEtoRange_mc -> SetLineColor(kBlue);
  fitted_range_mc -> SetLineColor(kGreen);
  BeamKEtoRange_mc -> Draw("histsame");
  fitted_range_mc -> Draw("histsame");
 
  TLegend *l = new TLegend(0.2, 0.75, 0.6, 0.95);
  l -> AddEntry(BeamKEtoRange_mc, "Beam KE_{inst} to Range", "l");
  l -> AddEntry(fitted_range_mc, "Best fit range", "l");
  l -> SetFillStyle(4000);
  l -> Draw("same");
  
  c -> SaveAs("./output/plots/Range_comparison_MC.pdf");

  max_x = BeamKEtoRange_data -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h -> Draw();
  BeamKEtoRange_data -> SetLineColor(kBlue);
  fitted_range_data -> SetLineColor(kGreen);
  BeamKEtoRange_data -> Draw("histsame");
  fitted_range_data -> Draw("histsame");

  l -> Draw("same");
  c -> SaveAs("./output/Range_comparison_Data.pdf");

  max_x = BeamP_mc -> GetMaximum();
  TH1D* template_h2 = new TH1D("", "", 1., 0., 800.);
  template_h2 -> GetXaxis() -> SetTitle("KE [MeV]");
  template_h2 -> GetXaxis() -> SetLabelSize(0.03);
  template_h2 -> GetXaxis() -> SetRangeUser(0., 800.);  
  template_h2 -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h2 -> GetYaxis() -> SetLabelSize(0.03);
  template_h2 -> GetYaxis() -> SetTitle("Events");
  template_h2 -> Draw();
  BeamP_mc -> SetLineColor(kBlue);
  fittedP_mc -> SetLineColor(kGreen);
  BeamP_mc -> Draw("histsame");
  fittedP_mc -> Draw("histsame");
  int max_bin_beamP_mc = BeamP_mc -> GetMaximumBin();
  double max_x_beamP_mc = BeamP_mc -> GetBinCenter(max_bin_beamP_mc);
  int max_bin_fittedP_mc = fittedP_mc -> GetMaximumBin();
  double max_x_fittedP_mc = fittedP_mc -> GetBinCenter(max_bin_fittedP_mc);
  TF1 *beamP_mc_gaus = new TF1("beamP_mc_gaus", "gaus", max_x_beamP_mc - width, max_x_beamP_mc + width);
  TF1 *fittedP_mc_gaus = new TF1("fittedP_mc_gaus", "gaus", max_x_fittedP_mc - width, max_x_fittedP_mc + width);
  cout << "max_x_fittedP_mc : " << max_x_fittedP_mc << endl;
  beamP_mc_gaus -> SetLineColor(kCyan);
  fittedP_mc_gaus ->SetLineColor(kSpring+3);
  beamP_mc_gaus -> SetLineStyle(7);
  fittedP_mc_gaus-> SetLineStyle(7);
  beamP_mc_gaus -> SetLineWidth(3);
  fittedP_mc_gaus-> SetLineWidth(3);
  BeamP_mc -> Fit(beamP_mc_gaus, "R", "", max_x_beamP_mc - width, max_x_beamP_mc + width);
  fittedP_mc -> Fit(fittedP_mc_gaus, "R", "", max_x_fittedP_mc - width, max_x_fittedP_mc + width);
  beamP_mc_gaus -> Draw("lsame");
  fittedP_mc_gaus -> Draw("lsame");

  TLegend *l2 = new TLegend(0.2, 0.7, 0.4, 0.9);
  l2 -> AddEntry(BeamP_mc, "Beam KE_{inst}", "l");
  l2 -> AddEntry(fittedP_mc, "Fitted beam KE", "l");
  l2 -> Draw("same");

  TPaveText *pt_BeamP = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  TPaveText *pt_fittedP = new TPaveText(0.17, 0.37, 0.50, 0.53,"NDC");
  pt_BeamP -> AddText(Form("#mu : %.3f, #sigma : %.3f", beamP_mc_gaus -> GetParameter(1), beamP_mc_gaus -> GetParameter(2)));
  ((TText*)pt_BeamP->GetListOfLines()->Last())->SetTextColor(kCyan);
  pt_fittedP -> AddText(Form("#mu : %.3f, #sigma : %.3f", fittedP_mc_gaus -> GetParameter(1), fittedP_mc_gaus -> GetParameter(2)));
  ((TText*)pt_fittedP->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_BeamP -> SetFillStyle(4000);
  pt_fittedP -> SetFillStyle(4000);
  pt_BeamP -> SetFillColor(kWhite);
  pt_fittedP -> SetFillColor(kWhite);

  pt_BeamP -> Draw("same");
  pt_fittedP -> Draw("same");
  c -> SaveAs("./output/P_comparison_mc.pdf");

  max_x = BeamP_data -> GetMaximum();
  template_h2 -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h2 -> Draw();
  BeamP_data -> SetLineColor(kBlue);
  fittedP_data -> SetLineColor(kGreen);
  BeamP_data -> Draw("histsame");
  fittedP_data -> Draw("histsame");
  int max_bin_beamP_data = BeamP_data -> GetMaximumBin();
  double max_x_beamP_data = BeamP_data -> GetBinCenter(max_bin_beamP_data);
  int max_bin_fittedP_data = fittedP_data -> GetMaximumBin();
  double max_x_fittedP_data = fittedP_data -> GetBinCenter(max_bin_fittedP_data);
  TF1 *beamP_data_gaus = new TF1("beamP_mdata_gaus", "gaus", max_x_beamP_data - width, max_x_beamP_data + width);
  TF1 *fittedP_data_gaus = new TF1("fittedP_data_gaus", "gaus", max_x_fittedP_data - width, max_x_fittedP_data + width);
  cout << "max_x_fittedP_data : " << max_x_fittedP_data << endl;
  beamP_data_gaus -> SetLineColor(kCyan);
  fittedP_data_gaus ->SetLineColor(kSpring+3);
  beamP_data_gaus -> SetLineStyle(7);
  fittedP_data_gaus-> SetLineStyle(7);
  beamP_data_gaus -> SetLineWidth(3);
  fittedP_data_gaus-> SetLineWidth(3);
  BeamP_data -> Fit(beamP_data_gaus, "R", "", max_x_beamP_data - width, max_x_beamP_data + width);
  fittedP_data -> Fit(fittedP_data_gaus, "R", "", max_x_fittedP_data - width, max_x_fittedP_data + width);
  beamP_data_gaus -> Draw("lsame");
  fittedP_data_gaus -> Draw("lsame");

  l2 -> Draw("same");

  TPaveText *pt_BeamP2 = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  TPaveText *pt_fittedP2 = new TPaveText(0.17, 0.37, 0.50, 0.53,"NDC");
  pt_BeamP2 -> AddText(Form("#mu : %.3f, #sigma : %.3f", beamP_data_gaus -> GetParameter(1), beamP_data_gaus -> GetParameter(2)));
  ((TText*)pt_BeamP2->GetListOfLines()->Last())->SetTextColor(kCyan);
  pt_fittedP2 -> AddText(Form("#mu : %.3f, #sigma : %.3f", fittedP_data_gaus -> GetParameter(1), fittedP_data_gaus -> GetParameter(2)));
  ((TText*)pt_fittedP2->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_BeamP2 -> SetFillStyle(4000);
  pt_fittedP2 -> SetFillStyle(4000);
  pt_BeamP2 -> SetFillColor(kWhite);
  pt_fittedP2 -> SetFillColor(kWhite);

  pt_BeamP2 -> Draw("same");
  pt_fittedP2 -> Draw("same");
  c -> SaveAs("./output/P_comparison_data.pdf");

  max_x = fitted_dP_data -> GetMaximum();
  TH1D* template_h3 = new TH1D("", "", 1., -400., 400.);
  template_h3 -> GetXaxis() -> SetTitle("KE_{inst} - KE_{fit} [MeV]");
  template_h3 -> GetXaxis() -> SetLabelSize(0.03);
  template_h3 -> GetYaxis() -> SetRangeUser(0., max_x * 1.2);
  template_h3 -> GetYaxis() -> SetLabelSize(0.03);
  template_h3 -> GetYaxis() -> SetTitle("Events");
  template_h3 -> Draw();
  fitted_dP_data -> SetLineColor(kBlue);
  fitted_dP_mc -> SetLineColor(kGreen);
  fitted_dP_data -> Draw("histsame");
  fitted_dP_mc -> Draw("histsame");
  int max_bin_fitted_dP_data = fitted_dP_data -> GetMaximumBin();
  double max_x_fitted_dP_data = fitted_dP_data -> GetBinCenter(max_bin_fitted_dP_data);
  int max_bin_fitted_dP_mc = fitted_dP_mc -> GetMaximumBin();
  double max_x_fitted_dP_mc = fitted_dP_mc -> GetBinCenter(max_bin_fitted_dP_mc);
  TF1 *fitted_dP_data_gaus = new TF1("fitted_dP_data_gaus", "gaus", max_x_fitted_dP_data - width, max_x_fitted_dP_data + width);
  fitted_dP_data_gaus -> SetLineColor(kCyan);
  fitted_dP_data_gaus -> SetLineStyle(7);
  fitted_dP_data_gaus -> SetLineWidth(3);
  fitted_dP_data -> Fit(fitted_dP_data_gaus, "R", "", max_x_fitted_dP_data - width, max_x_fitted_dP_data + width);
  fitted_dP_data -> Draw("lsame");

  cout << "max_bin_fitted_dP_data : " << max_bin_fitted_dP_data << ", max_x_fitted_dP_data : " << max_x_fitted_dP_data << endl;

  TF1 *fitted_dP_mc_gaus = new TF1("fitted_dP_mc_gaus", "gaus", max_x_fitted_dP_mc - width, max_x_fitted_dP_mc + width);
  fitted_dP_mc_gaus -> SetLineColor(kSpring+3);
  fitted_dP_mc_gaus -> SetLineStyle(7);
  fitted_dP_mc_gaus -> SetLineWidth(3);
  fitted_dP_mc -> Fit(fitted_dP_mc_gaus, "R", "", max_x_fitted_dP_mc - width, max_x_fitted_dP_mc + width);
  fitted_dP_mc -> Draw("lsame");

  cout << "max_bin_fitted_dP_mc : " << max_bin_fitted_dP_mc << ", max_x_fitted_dP_mc : " << max_x_fitted_dP_mc << endl;

  TLegend *l3 = new TLegend(0.2, 0.7, 0.4, 0.9);
  l3 -> AddEntry(fitted_dP_data, "Data", "l");
  l3 -> AddEntry(fitted_dP_mc, "MC", "l");
  l3 -> Draw("same");

  TPaveText *pt_fitted_dP_data = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  pt_fitted_dP_data -> AddText(Form("#mu : %.3f, #sigma : %.3f", fitted_dP_data_gaus -> GetParameter(1), fitted_dP_data_gaus -> GetParameter(2)));
  ((TText*)pt_fitted_dP_data->GetListOfLines()->Last())->SetTextColor(kCyan);
  pt_fitted_dP_data -> SetFillStyle(4000);
  pt_fitted_dP_data -> SetFillColor(kWhite);
  pt_fitted_dP_data -> Draw("same");

  TPaveText *pt_fitted_dP_mc = new TPaveText(0.17, 0.37, 0.50, 0.53,"NDC");
  pt_fitted_dP_mc -> AddText(Form("#mu : %.3f, #sigma : %.3f", fitted_dP_mc_gaus -> GetParameter(1), fitted_dP_mc_gaus -> GetParameter(2)));
  ((TText*)pt_fitted_dP_mc->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_fitted_dP_mc -> SetFillStyle(4000);
  pt_fitted_dP_mc -> SetFillColor(kWhite);
  pt_fitted_dP_mc -> Draw("same");

  c -> SaveAs("./output/dP_comparison_data.pdf");

  template_h3 -> Draw();
  fitted_dP_data -> Draw("histsame");
  fitted_dP_mc -> Draw("histsame");
  l3 -> Draw("same");
  
  Double_t fitting_range[2];
  fitting_range[0] = -10.;
  fitting_range[1] = 200.;
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  sv[0] = 3.;
  sv[1] = 20.;
  sv[2] = fitted_dP_data -> Integral() * 0.05;
  sv[3] = 6.;
  for(int j=0; j<4; ++j){
    pllo[j] = 0.01*sv[j];
    plhi[j] = 100*sv[j];
  }
  Double_t chisqr;
  Int_t    ndf;
  Int_t    status;  
  TF1 *Langau_data = langaufit(fitted_dP_data, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_data");
  Langau_data -> SetLineColor(kCyan);
  Langau_data -> SetLineStyle(7);
  Langau_data -> SetLineWidth(3);
  Langau_data -> Draw("lsame");
  //Langau_data -> Draw();
  //Langau_data_clone  -> Draw();

  std::cout << "************** MPV : " << Langau_data->GetParameter(1) << " +/- " << Langau_data->GetParError(1) << std::endl;
  std::cout << "************** Chi^2/NDF : " << Langau_data->GetChisquare()/Langau_data->GetNDF() << std::endl;

  TF1 *Langau_mc = langaufit(fitted_dP_mc, fitting_range, sv,pllo,plhi,fp,fpe,&chisqr,&ndf,&status, "Langau_mc");
  Langau_mc -> SetLineColor(kSpring+3);
  Langau_mc -> SetLineStyle(7);
  Langau_mc -> SetLineWidth(3);
  Langau_mc -> Draw("lsame");
  std::cout << "************** MPV : " << Langau_mc->GetParameter(1) << " +/- " << Langau_mc->GetParError(1) << std::endl;
  std::cout << "************** Chi^2/NDF : " << Langau_mc->GetChisquare()/Langau_mc->GetNDF() << std::endl;

  TPaveText *pt_fitted_dP_data_Langau = new TPaveText(0.61, 0.67, 0.94, 0.83, "NDC");
  pt_fitted_dP_data_Langau -> AddText(Form("MPV : %.3f, #sigma_{Gaus} : %.3f", Langau_data -> GetParameter(1), Langau_data -> GetParameter(3)));
  ((TText*)pt_fitted_dP_data_Langau->GetListOfLines()->Last())->SetTextColor(kCyan);
  pt_fitted_dP_data_Langau -> SetFillStyle(4000);
  pt_fitted_dP_data_Langau -> SetFillColor(kWhite);
  pt_fitted_dP_data_Langau -> Draw("same");

  TPaveText *pt_fitted_dP_mc_Langau = new TPaveText(0.17, 0.37, 0.50, 0.53,"NDC");
  pt_fitted_dP_mc_Langau -> AddText(Form("MPV : %.3f, #sigma_{Gaus} : %.3f", Langau_mc -> GetParameter(1), Langau_mc -> GetParameter(3)));
  ((TText*)pt_fitted_dP_mc_Langau->GetListOfLines()->Last())->SetTextColor(kSpring+3);
  pt_fitted_dP_mc_Langau -> SetFillStyle(4000);
  pt_fitted_dP_mc_Langau -> SetFillColor(kWhite);
  pt_fitted_dP_mc_Langau -> Draw("same");


  c -> SaveAs("./output/dP_comparison_Langau.pdf");

}


void Draw_range_fit_result(){

  setTDRStyle();
  Fit_and_Draw("penergy_range_fit_scraper_veto.root", 100.);
  //Fit_Gaus_MC_and_Data("penergy_range_fit_scraper_veto.root", "dKE_fitted_vs_Truth", "KE_{true} - KE_{fitted}", -400., 400., 40.);
  Fit_Gaus_MC_only("penergy_range_fit_scraper_veto.root", "dKE_fitted_vs_Truth", "KE_{true,ff} - KE_{fitted}", -100., 100., 8.);

}
