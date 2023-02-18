#include "canvas_margin.h"
#include "mylib.h"

void Draw_2D_KE_each_Nhit(TString filename, TString method, TString this_Nhit_str, TString this_Nhit_legend, TString particle, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  TString data_or_MC = "";
  bool isMC = false;
  if(filename.Contains("MC")){
    data_or_MC = "MC";
    isMC = true;
  }  
  else data_or_MC = "Data";

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/PionKEScale/";
  TFile *f_input = new TFile(root_file_path + filename);
  gDirectory -> Cd(method);

  TString this_histname = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_" + particle;

  TH2D *this_hist = nullptr;
  if((TH2D*)gDirectory-> Get(this_histname)) this_hist = (TH2D*)gDirectory -> Get(this_histname) -> Clone();

  if(this_hist == nullptr) return;
  this_hist -> RebinX(rebin_x);
  this_hist -> RebinY(rebin_y);

  double z_max = this_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleY);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  this_hist -> SetLineWidth(1);
  this_hist -> SetLineColor(kGreen);
  this_hist -> Draw("colsame");

  TF1 *yequlx = new TF1("yequlx", "x", xmin, xmax);
  yequlx -> SetLineColor(kRed);
  yequlx -> Draw("lsame");

  TLegend *l = new TLegend(0.3, 0.7, 0.5, 0.9);
  l -> AddEntry(yequlx, "y = x", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_Nhit;
  latex_ProtoDUNE.SetNDC();
  latex_Nhit.SetNDC();
  latex_Nhit.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_Nhit.SetTextSize(0.025);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  if(isMC) latex_Nhit.DrawLatex(0.95, 0.97, data_or_MC + ", " + particle + ", " + this_Nhit_legend);
  else latex_Nhit.DrawLatex(0.95, 0.97, data_or_MC + ", " + this_Nhit_legend);
  c -> SaveAs("./output/plots/PionKE/2D/" + data_or_MC + "_" + method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_" + particle + ".pdf");
  c -> Close();
}

void Draw_2D_KE_each_Nhit_all_MC(TString filename, TString method, TString this_Nhit_str, TString this_Nhit_legend, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  TString data_or_MC = "";
  bool isMC = false;
  if(filename.Contains("MC")){
    data_or_MC = "MC";
    isMC = true;
  }
  else data_or_MC = "Data";

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/PionKEScale/";
  TFile *f_input = new TFile(root_file_path + filename);
  gDirectory -> Cd(method);

  TString this_histname = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_";

  TH2D *pion_hist = nullptr;
  TH2D *muon_hist = nullptr;
  TH2D *proton_hist = nullptr;
  TH2D *other_hist = nullptr;

  if((TH2D*)gDirectory-> Get(this_histname + "pion")) pion_hist = (TH2D*)gDirectory -> Get(this_histname + "pion") -> Clone();
  if((TH2D*)gDirectory-> Get(this_histname + "muon")) muon_hist = (TH2D*)gDirectory -> Get(this_histname + "muon") -> Clone();
  if((TH2D*)gDirectory-> Get(this_histname + "proton")) proton_hist = (TH2D*)gDirectory -> Get(this_histname + "proton") -> Clone();
  if((TH2D*)gDirectory-> Get(this_histname + "other")) other_hist = (TH2D*)gDirectory -> Get(this_histname + "other") -> Clone();

  if(pion_hist == nullptr) return;
  if(muon_hist != nullptr) pion_hist -> Add(muon_hist);
  if(proton_hist != nullptr) pion_hist -> Add(proton_hist);
  if(other_hist != nullptr) pion_hist -> Add(other_hist);

  pion_hist -> RebinX(rebin_x);
  pion_hist -> RebinY(rebin_y);

  double z_max = pion_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH2D *template_h = new TH2D("", "", 1, xmin, xmax, 1, ymin, ymax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleY);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetZaxis() -> SetTitle("Events");
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");

  pion_hist -> SetLineWidth(1);
  pion_hist -> SetLineColor(kGreen);
  pion_hist -> Draw("colsame");

  TF1 *yequlx = new TF1("yequlx", "x", xmin, xmax);
  yequlx -> SetLineColor(kRed);
  yequlx -> Draw("lsame");

  TLegend *l = new TLegend(0.3, 0.7, 0.5, 0.9);
  l -> AddEntry(yequlx, "y = x", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_Nhit;
  latex_ProtoDUNE.SetNDC();
  latex_Nhit.SetNDC();
  latex_Nhit.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_Nhit.SetTextSize(0.025);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  if(isMC) latex_Nhit.DrawLatex(0.95, 0.97, data_or_MC + ", all, " + this_Nhit_legend);
  c -> SaveAs("./output/plots/PionKE/2D/" + data_or_MC + "_" + method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_all.pdf");
  c -> Close();
}

void Draw_efficiency_KE_range(TString filename, TString MC_or_Data, TString method, TString this_Nhit_str, TString this_Nhit_legend, TString particle, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  TString data_or_MC = "";
  bool isMC = false;
  if(filename.Contains("MC")){
    data_or_MC = "MC";
    isMC = true;
  }
  else data_or_MC = "Data";

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/PionKEScale/";
  TFile *f_input = new TFile(root_file_path + "PionKEScale_1.0_" + MC_or_Data + filename);

  gDirectory -> Cd(method);
  TString this_2D_histname = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_" + particle;
  TH2D *this_2D_hist = nullptr;
  if((TH2D*)gDirectory -> Get(this_2D_histname)) this_2D_hist = (TH2D*)gDirectory -> Get(this_2D_histname) -> Clone();
  gDirectory -> Cd("../Denom");
  TString this_1D_histname = "KE_beam_" + this_Nhit_str + "_" + particle;
  TH1D *this_1D_hist = nullptr;
  if((TH1D*)gDirectory -> Get(this_1D_histname)) this_1D_hist = (TH1D*)gDirectory -> Get(this_1D_histname) -> Clone();

  if(this_2D_hist == nullptr || this_1D_hist == nullptr) return;

  // == Bin size : 1D (1 MeV), 2D (5 MeV)
  this_1D_hist -> Rebin(rebin_y);
  this_1D_hist -> Rebin(5);
  this_2D_hist -> RebinX(rebin_x);
  this_2D_hist -> RebinY(rebin_y);

  int Nbins_X = this_2D_hist -> GetNbinsX();
  cout << "[Draw_efficiency_KE_range] Nbins_X : " << Nbins_X << endl;
  TH1D * projected_2D_hist = this_2D_hist -> ProjectionY("projected_2D_hist", 1, Nbins_X, "");

  TH1D * eff_hist = (TH1D*)projected_2D_hist -> Clone();
  eff_hist -> Divide(this_1D_hist);

  double max_y = this_1D_hist -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleY);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 1.2);
  //template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.5);
  template_h -> Draw("hist");

  eff_hist -> SetLineWidth(1);
  eff_hist -> SetLineColor(kGreen);
  eff_hist -> Draw("histsame");

  /*
  this_1D_hist -> SetLineColor(kBlue);
  projected_2D_hist -> SetLineColor(kCyan);
  this_1D_hist -> Draw("histsame");
  projected_2D_hist -> Draw("histsame");
  */
  TF1 *yequlx = new TF1("yequlx", "1", xmin, xmax);
  yequlx -> SetLineColor(kRed);
  yequlx -> Draw("lsame");

  TLegend *l = new TLegend(0.7, 0.4, 0.9, 0.5);
  l -> AddEntry(yequlx, "y = 1", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_Nhit;
  latex_ProtoDUNE.SetNDC();
  latex_Nhit.SetNDC();
  latex_Nhit.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_Nhit.SetTextSize(0.025);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  if(isMC) latex_Nhit.DrawLatex(0.95, 0.97, data_or_MC + ", " + particle + ", " + this_Nhit_legend);
  else latex_Nhit.DrawLatex(0.95, 0.97, data_or_MC + ", " + this_Nhit_legend);
  c -> SaveAs("./output/plots/PionKE/Eff/" + data_or_MC + "_" + method + "_Eff_" + this_Nhit_str + "_" + particle + ".pdf");
  c -> Close(); 

}

void Draw_efficiency_KE_range_MC_vs_Data(TString filename, TString method, TString this_Nhit_str, TString this_Nhit_legend, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){
  

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/PionKEScale/";

  // == Call MC
  TFile *f_MC = new TFile(root_file_path + "PionKEScale_1.0_MC" + filename);
  gDirectory -> Cd(method);
  TString MC_2D_histname_pion = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_pion";
  TString MC_2D_histname_muon = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_muon";
  TString MC_2D_histname_proton = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_proton";
  TString MC_2D_histname_other = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_other";

  TH2D *MC_2D_hist_pion = nullptr;
  TH2D *MC_2D_hist_muon = nullptr;
  TH2D *MC_2D_hist_proton = nullptr;
  TH2D *MC_2D_hist_other = nullptr;

  if((TH2D*)gDirectory -> Get(MC_2D_histname_pion)) MC_2D_hist_pion = (TH2D*)gDirectory -> Get(MC_2D_histname_pion) -> Clone();
  if((TH2D*)gDirectory -> Get(MC_2D_histname_muon)) MC_2D_hist_muon = (TH2D*)gDirectory -> Get(MC_2D_histname_muon) -> Clone();
  if((TH2D*)gDirectory -> Get(MC_2D_histname_proton)) MC_2D_hist_proton = (TH2D*)gDirectory -> Get(MC_2D_histname_proton) -> Clone();
  if((TH2D*)gDirectory -> Get(MC_2D_histname_other)) MC_2D_hist_other = (TH2D*)gDirectory -> Get(MC_2D_histname_other) -> Clone();

  if(MC_2D_hist_pion == nullptr) return;

  TH2D *MC_2D_hist = (TH2D*)MC_2D_hist_pion -> Clone();
  if(MC_2D_hist_muon != nullptr) MC_2D_hist -> Add(MC_2D_hist_muon);
  if(MC_2D_hist_proton != nullptr) MC_2D_hist -> Add(MC_2D_hist_proton);
  if(MC_2D_hist_other != nullptr) MC_2D_hist -> Add(MC_2D_hist_other);

  gDirectory -> Cd("../Denom");
  TString MC_1D_histname_pion = "KE_beam_" + this_Nhit_str + "_pion";
  TString MC_1D_histname_muon = "KE_beam_" + this_Nhit_str + "_muon";
  TString MC_1D_histname_proton = "KE_beam_" + this_Nhit_str + "_proton";
  TString MC_1D_histname_other = "KE_beam_" + this_Nhit_str + "_other";

  TH1D *MC_1D_hist_pion = nullptr;
  TH1D *MC_1D_hist_muon = nullptr;
  TH1D *MC_1D_hist_proton = nullptr;
  TH1D *MC_1D_hist_other = nullptr;

  if((TH1D*)gDirectory -> Get(MC_1D_histname_pion)) MC_1D_hist_pion = (TH1D*)gDirectory -> Get(MC_1D_histname_pion) -> Clone();
  if((TH1D*)gDirectory -> Get(MC_1D_histname_muon)) MC_1D_hist_muon = (TH1D*)gDirectory -> Get(MC_1D_histname_muon) -> Clone();
  if((TH1D*)gDirectory -> Get(MC_1D_histname_proton)) MC_1D_hist_proton = (TH1D*)gDirectory -> Get(MC_1D_histname_proton) -> Clone();
  if((TH1D*)gDirectory -> Get(MC_1D_histname_other)) MC_1D_hist_other = (TH1D*)gDirectory -> Get(MC_1D_histname_other) -> Clone();

  if(MC_1D_hist_pion == nullptr) return;

  TH1D *MC_1D_hist = (TH1D*)MC_1D_hist_pion -> Clone();
  if(MC_1D_hist_muon != nullptr) MC_1D_hist -> Add(MC_1D_hist_muon);
  if(MC_1D_hist_proton != nullptr) MC_1D_hist -> Add(MC_1D_hist_proton);
  if(MC_1D_hist_other != nullptr) MC_1D_hist -> Add(MC_1D_hist_other);

  // == Call Data
  TFile *f_Data = new TFile(root_file_path + "PionKEScale_1.0_Data" + filename);
  gDirectory -> Cd(method);
  TString Data_2D_histname = method + "_KE_fit_vs_KE_BB_" + this_Nhit_str + "_Data";
  TH2D *Data_2D_hist = nullptr;
  if((TH2D*)gDirectory -> Get(Data_2D_histname)) Data_2D_hist = (TH2D*)gDirectory -> Get(Data_2D_histname) -> Clone();
  if(Data_2D_hist == nullptr) return;
  gDirectory -> Cd("../Denom");
  TString Data_1D_histname = "KE_beam_" + this_Nhit_str + "_Data";
  TH1D *Data_1D_hist = nullptr;
  if((TH1D*)gDirectory -> Get(Data_1D_histname)) Data_1D_hist = (TH1D*)gDirectory -> Get(Data_1D_histname) -> Clone();
  if(Data_1D_hist == nullptr) return;

  // == Rebinning, Bin size : 1D (1 MeV), 2D (5 MeV)
  MC_2D_hist -> RebinX(rebin_x);
  MC_2D_hist -> RebinY(rebin_y);
  MC_1D_hist -> Rebin(rebin_y);
  MC_1D_hist -> Rebin(5);
  Data_2D_hist -> RebinX(rebin_x);
  Data_2D_hist -> RebinY(rebin_y);
  Data_1D_hist -> Rebin(rebin_y);
  Data_1D_hist -> Rebin(5);

  // == Projection
  int Nbins_X = MC_2D_hist -> GetNbinsX();
  TH1D * MC_projected_2D_hist = MC_2D_hist -> ProjectionY("MC_projected_2D_hist", 1, Nbins_X, "");
  TH1D * Data_projected_2D_hist = Data_2D_hist -> ProjectionY("Data_projected_2D_hist", 1, Nbins_X, "");

  TH1D * MC_eff_hist = (TH1D*)MC_projected_2D_hist -> Clone();
  MC_eff_hist -> Divide(MC_1D_hist);

  TH1D * Data_eff_hist = (TH1D*)Data_projected_2D_hist -> Clone();
  Data_eff_hist -> Divide(Data_1D_hist);

  // == Efficiency graphes
  vector<double> eff_MC, eff_MC_eyl, eff_MC_eyh, eff_Data, eff_Data_eyl, eff_Data_eyh, KE, KE_exl, KE_exh;
  for(int i = 1; i < Nbins_X + 1; i++){
    double denom_MC = MC_1D_hist -> GetBinContent(i);
    double numer_MC = MC_projected_2D_hist -> GetBinContent(i);
    double this_eff_MC = numer_MC / denom_MC;
    double this_err_MC_low = this_eff_MC - TEfficiency::ClopperPearson(denom_MC, numer_MC, 0.683, false);
    double this_err_MC_high = TEfficiency::ClopperPearson(denom_MC, numer_MC, 0.683, true) - this_eff_MC;

    double denom_Data = Data_1D_hist -> GetBinContent(i);
    double numer_Data = Data_projected_2D_hist -> GetBinContent(i);
    double this_eff_Data = numer_Data / denom_Data;
    double this_err_Data_low = this_eff_Data - TEfficiency::ClopperPearson(denom_Data, numer_Data, 0.683, false);
    double this_err_Data_high = TEfficiency::ClopperPearson(denom_Data, numer_Data, 0.683, true) - this_eff_Data;

    double this_KE = MC_1D_hist -> GetBinCenter(i);
    double this_KE_err = MC_1D_hist -> GetBinWidth(i) / 2.;

    eff_MC.push_back(this_eff_MC);
    eff_MC_eyl.push_back(this_err_MC_low);
    eff_MC_eyh.push_back(this_err_MC_high);
    eff_Data.push_back(this_eff_Data);
    eff_Data_eyl.push_back(this_err_Data_low);
    eff_Data_eyh.push_back(this_err_Data_high);
    KE.push_back(this_KE);
    KE_exl.push_back(this_KE_err);
    KE_exh.push_back(this_KE_err);

    cout << i << ", this_eff_MC : " << this_eff_MC << " + " << this_err_MC_high << " - " << this_err_MC_low << ", for " << numer_MC << " / " << denom_MC << endl;
    cout << i << ", this_eff_Data : " << this_eff_Data << " + " << this_err_Data_high << " - " << this_err_Data_low << ", for " << numer_Data << " / " << denom_Data << endl;
  }

  TGraphAsymmErrors *gr_MC = new TGraphAsymmErrors(Nbins_X, &KE[0], &eff_MC[0], &KE_exl[0], &KE_exh[0], &eff_MC_eyl[0], &eff_MC_eyh[0]);
  gr_MC -> SetFillColor(kRed);
  gr_MC -> SetFillStyle(3004);

  TGraphAsymmErrors *gr_Data = new TGraphAsymmErrors(Nbins_X, &KE[0], &eff_Data[0], &KE_exl[0], &KE_exh[0], &eff_Data_eyl[0], &eff_Data_eyh[0]);
  gr_Data-> SetFillColor(kBlue);
  gr_Data-> SetFillStyle(3005);


  // == Draw
  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle(TitleY);
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 1.2);
  template_h -> Draw("hist");

  MC_eff_hist -> SetLineWidth(2);
  MC_eff_hist -> SetLineColor(kRed);
  MC_eff_hist -> Draw("histsame");

  Data_eff_hist -> SetLineWidth(2);
  Data_eff_hist -> SetLineColor(kBlue);
  Data_eff_hist -> Draw("histsame");

  gr_MC -> Draw("e2same");
  gr_Data -> Draw("e2same");

  TF1 *yequlx = new TF1("yequlx", "1", xmin, xmax);
  yequlx -> SetLineColor(kBlack);
  yequlx -> SetLineWidth(2);
  yequlx -> SetLineStyle(5);
  yequlx -> Draw("lsame");  

  TLegend *l = new TLegend(0.7, 0.83, 0.9, 0.93);
  l -> AddEntry(yequlx, "y = 1", "l");
  l -> AddEntry(MC_eff_hist, "MC Eff.", "l");
  l -> AddEntry(Data_eff_hist, "Data Eff.","l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_Nhit;
  latex_ProtoDUNE.SetNDC();
  latex_Nhit.SetNDC();
  latex_Nhit.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_Nhit.SetTextSize(0.025);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Nhit.DrawLatex(0.95, 0.97, "All, " + this_Nhit_legend);
  c -> SaveAs("./output/plots/PionKE/Eff/Comparison_" + method + "_Eff_" + this_Nhit_str + "_all.pdf");
  c -> Close();
  
  f_MC -> Close();
  f_Data -> Close();
}

void Run_Draw_2D_MC_and_Data(TString filename, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  int Nhit_low = 0;
  int Nhit_high = 300;
  int Nhit_step = 30;
  int Nhit_steps = (Nhit_high - Nhit_low) / Nhit_step;

  TString particles[4] = {"pion", "proton", "muon", "other"};
  int N_particles = 4;

  for(int i = 0; i < Nhit_steps; i++){

    TString this_Nhit_str = Form("Nhits%dto%d", Nhit_step * i, Nhit_step * (i + 1));
    TString this_Nhit_legend = Form("N_{hits} : %d - %d", Nhit_step * i, Nhit_step * (i + 1));
    if(filename.Contains("MC")){

      for(int j = 0; j < N_particles; j++){
	TString this_particle = particles[j];
	Draw_2D_KE_each_Nhit(filename, "Gaussian", this_Nhit_str, this_Nhit_legend, this_particle, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
	Draw_2D_KE_each_Nhit(filename, "Likelihood", this_Nhit_str, this_Nhit_legend, this_particle, TitleX, TitleY, xmin,xmax, rebin_x, ymin, ymax, rebin_y);
      }
      Draw_2D_KE_each_Nhit_all_MC(filename, "Likelihood", this_Nhit_str, this_Nhit_legend, TitleX, TitleY, xmin,xmax, rebin_x, ymin, ymax, rebin_y);
    }
    else{
      Draw_2D_KE_each_Nhit(filename, "Gaussian", this_Nhit_str, this_Nhit_legend, "Data", TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
      Draw_2D_KE_each_Nhit(filename, "Likelihood", this_Nhit_str, this_Nhit_legend, "Data", TitleX, TitleY, xmin,xmax, rebin_x, ymin, ymax, rebin_y);
    }
  }
}

void Run_Draw_Eff(TString filename, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){
  int Nhit_low = 0;
  int Nhit_high = 300;
  int Nhit_step = 30;
  int Nhit_steps = (Nhit_high - Nhit_low) / Nhit_step;

  TString MC_or_Data[2] = {"MC", "Data"};
  TString particles[4] = {"pion", "proton", "muon", "other"};
  int N_particles = 4;

  for(int i = 0; i < Nhit_steps; i++){

    TString this_Nhit_str = Form("Nhits%dto%d", Nhit_step * i, Nhit_step * (i + 1));
    TString this_Nhit_legend = Form("N_{hits} : %d - %d", Nhit_step * i, Nhit_step * (i + 1));


    for(int j = 0; j < 2; j++){
      TString this_MC_or_Data= MC_or_Data[j];

      if(this_MC_or_Data.Contains("MC")){
	for(int k = 0; k < N_particles; k++){
	  TString this_particle = particles[k];
	  Draw_efficiency_KE_range(filename, this_MC_or_Data, "Likelihood", this_Nhit_str, this_Nhit_legend, this_particle, "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
	  Draw_efficiency_KE_range(filename, this_MC_or_Data, "Gaussian", this_Nhit_str, this_Nhit_legend, this_particle, "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
	}
      }
      else{
	Draw_efficiency_KE_range(filename, this_MC_or_Data, "Likelihood", this_Nhit_str, this_Nhit_legend, "Data", "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
	Draw_efficiency_KE_range(filename, this_MC_or_Data, "Gaussian", this_Nhit_str, this_Nhit_legend, "Data", "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
      }
    }

    Draw_efficiency_KE_range_MC_vs_Data(filename, "Likelihood", this_Nhit_str, this_Nhit_legend, "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    Draw_efficiency_KE_range_MC_vs_Data(filename, "Gaussian", this_Nhit_str, this_Nhit_legend, "KE_{range} [MeV]", "Eff.", xmin, xmax, rebin_x, ymin, ymax, rebin_y);
  }

}

void KEScale(){

  cout << "============================" << endl;
  cout << "=========== START ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = ".root";

  Run_Draw_2D_MC_and_Data("PionKEScale_1.0_Data_1GeV_test.root", "KE_{fitted} [MeV]", "KE_{range} [MeV]", 0., 1500., 5., 0., 1500., 5.);
  Run_Draw_2D_MC_and_Data("PionKEScale_1.0_MC_1GeV_test.root", "KE_{fitted} [MeV]", "KE_{range} [MeV]", 0., 1500., 5., 0., 1500., 5.);
  Run_Draw_Eff("_1GeV_test_Fit_MinNhit_15.root", "KE_{fitted} [MeV]", "KE_{range} [MeV]", 0., 500., 20., 0., 1500., 20.);
}
