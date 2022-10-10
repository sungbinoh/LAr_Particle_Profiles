#include "canvas_margin.h"
#include "mylib.h"

void Draw_XY_each_KE(TString filename, TString dir, TString this_KE_str, TString TitleX, TString TitleY, TString beam_P, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  //htrack_beam_inst_XY_nonscraper_BeamKE600to700MeV_precut_noweight_1;1
  TString this_histname_all = "htrack_beam_inst_XY_Beam" + this_KE_str + "_" + dir + "_0";
  TString histname_scraper = "htrack_beam_inst_XY_scraper_Beam" + this_KE_str + "_" + dir;
  TString histname_nonscraper = "htrack_beam_inst_XY_nonscraper_Beam" + this_KE_str + "_" + dir;

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);
  if((TH2D*)gDirectory -> Get(this_histname_all)){
    maphist2D[this_histname_all + "mc"] = (TH2D*)gDirectory -> Get(this_histname_all) -> Clone();
    maphist2D[this_histname_all + "mc"] -> RebinX(rebin_x);
    maphist2D[this_histname_all + "mc"] -> RebinY(rebin_y);
  }
  for(int i = 1; i < 3; i++){
    TString this_histname_scraper = Form(histname_scraper + "_%d", i);
    TString this_histname_nonscraper = Form(histname_nonscraper + "_%d", i);
    if((TH2D*)gDirectory -> Get(this_histname_scraper)){
      maphist2D[this_histname_scraper] = (TH2D*)gDirectory -> Get(this_histname_scraper) -> Clone();
      maphist2D[this_histname_scraper] -> RebinX(rebin_x);
      maphist2D[this_histname_scraper] -> RebinY(rebin_y);
    }
    else  maphist2D[this_histname_scraper] = nullptr;
    if((TH2D*)gDirectory-> Get(this_histname_nonscraper)){
      maphist2D[this_histname_nonscraper] = (TH2D*)gDirectory -> Get(this_histname_nonscraper) -> Clone();
      maphist2D[this_histname_nonscraper] -> RebinX(rebin_x);
      maphist2D[this_histname_nonscraper] -> RebinY(rebin_y);
    }
    else  maphist2D[this_histname_nonscraper] = nullptr;
  }

  TFile *f_data = new TFile(root_file_path + "data" + filename);
  gDirectory -> Cd(dir);
  if((TH2D*)gDirectory -> Get(this_histname_all)){
    maphist2D[this_histname_all + "data"] = (TH2D*)gDirectory -> Get(this_histname_all) -> Clone();
    maphist2D[this_histname_all + "data"] -> RebinX(rebin_x);
    maphist2D[this_histname_all + "data"] -> RebinY(rebin_y);
  }
  else maphist2D[this_histname_all + "data"] = nullptr;

  if(maphist2D[histname_scraper + "_1"] == nullptr) return;

  cout << histname_scraper + "_1" << endl;
  TH2D * h_scraper = (TH2D*)maphist2D[histname_scraper + "_1"] -> Clone();
  if(maphist2D[histname_scraper + "_2"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_2"]);
  TH2D * h_nonscraper = (TH2D*)maphist2D[histname_nonscraper + "_1"] -> Clone();
  if(maphist2D[histname_nonscraper + "_2"] != nullptr) h_nonscraper -> Add(maphist2D[histname_nonscraper + "_2"]);
  double z_max = h_scraper -> GetMaximum();
  if(z_max < h_nonscraper -> GetMaximum()) z_max = h_nonscraper -> GetMaximum();

  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  //c->SetRightMargin(0.15);
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

  h_nonscraper -> SetMarkerColor(kGreen);
  h_scraper -> SetMarkerColor(kRed);
  h_nonscraper -> SetLineColor(kGreen);
  h_scraper -> SetLineColor(kRed);
  h_nonscraper -> SetFillColor(kWhite);
  h_scraper -> SetFillColor(kWhite);

  h_nonscraper -> Draw("boxsame");
  h_scraper -> Draw("boxsame");
 
  TEllipse *circle = new TEllipse(-29.6, 422, 1.5 * 4.8);
  circle -> SetLineColor(kBlue);
  circle -> SetLineStyle(7);
  circle -> SetFillColorAlpha(kWhite, 0.);
  circle -> Draw("lsame");

  TEllipse *circle_HY = new TEllipse(-29.1637, 421.76, 4.50311 * 1.5, 3.83908 * 1.5);
  circle_HY -> SetLineColor(kOrange);
  circle_HY -> SetLineStyle(7);
  circle_HY -> SetFillColorAlpha(kWhite, 0.);
  circle_HY -> Draw("lsame");
 
  TLegend *l = new TLegend(0.20, 0.6, 0.94, 0.94);
  l -> SetNColumns(2);
  l -> AddEntry(h_nonscraper, "Non-scraper", "lf");
  l -> AddEntry(h_scraper, "Scraper", "lf");
  l -> AddEntry(circle, "1.5 #sigma", "l");
  l -> AddEntry(circle_HY, "Heng-Ye's ellipse", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_KE;
  latex_ProtoDUNE.SetNDC();
  latex_KE.SetNDC();
  latex_KE.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_KE.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, this_KE_str);

  c -> SaveAs("./output/plots/BeamStudy/Beam_scraper/BeamXY_scraper" + beam_P + "GeV_" + this_KE_str + "_" + dir + ".pdf");

  // == Draw mc all
  if(maphist2D[this_histname_all + "mc"] == nullptr) return;
  z_max = maphist2D[this_histname_all + "mc"] -> GetMaximum();
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");
  TH2D * h_mc = (TH2D*)maphist2D[this_histname_all + "mc"] -> Clone();
  h_mc -> SetMarkerColor(kGreen);
  h_mc -> SetLineColor(kGreen);
  h_mc -> SetFillColor(kWhite);
  h_mc -> Draw("boxsame");
  circle -> Draw("lsame");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, this_KE_str);

  c -> SaveAs("./output/plots/BeamStudy/Beam_scraper/BeamXY_mc_scraper" + beam_P + "GeV_" + this_KE_str + "_" + dir + ".pdf");

  // == Draw data all
  if(maphist2D[this_histname_all + "data"] == nullptr) return;
  z_max = maphist2D[this_histname_all + "data"] -> GetMaximum();
  template_h -> GetZaxis() -> SetRangeUser(0., z_max * 1.5);
  template_h -> Draw("colz");
  TH2D * h_data = (TH2D*)maphist2D[this_histname_all + "data"] -> Clone();
  h_data -> SetMarkerColor(kGreen);
  h_data -> SetLineColor(kGreen);
  h_data -> SetFillColor(kWhite);
  h_data -> Draw("boxsame");
  circle -> Draw("lsame");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, this_KE_str);

  c -> SaveAs("./output/plots/BeamStudy/Beam_scraper/BeamXY_data_scraper" + beam_P + "GeV_" + this_KE_str + "_" + dir + ".pdf");

  c -> Close();

}

void Run_Draw_2D_MC_and_Data(TString filename, TString dir, TString TitleX, TString TitleY, TString beam_P, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  int KE_min = 600;
  int KE_max = 1200;
  int KE_step = 100;
  int N_step = (KE_max - KE_min) / KE_step;
  vector<TString> KE_str_vec;
  for(int i = 0; i < N_step; i++){
    TString this_KE_str = Form("KE%dto%dMeV", KE_min + i * KE_step, KE_min + (i + 1) * KE_step);
    cout << "this_KE_str : " << this_KE_str << endl;
    Draw_XY_each_KE(filename, dir, this_KE_str, TitleX, TitleY, beam_P, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
  }

}

void Draw_Beam_Scraper(){

  cout << "============================" << endl;
  cout << "=========== START ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = "_beam_study.root";

  Run_Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "precut_noweight", "X_{ff} [cm]", "Y_{ff} [cm]", "1.0", -50, 10, 5., 400., 460., 5.);
  Run_Draw_2D_MC_and_Data("_PionXsec_1.0GeV" + file_suffix, "noweight", "X_{ff} [cm]", "Y_{ff} [cm]", "1.0", -50, 10, 5., 400., 460., 5.);

}
