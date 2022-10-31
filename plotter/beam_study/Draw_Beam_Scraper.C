#include "canvas_margin.h"
#include "mylib.h"

void Draw_XY_each_KE(TString data_or_mc, TString filename, TString key, TString dir, TString beam_P, TString this_KE_str, TString this_KE_legend, TString TitleX, TString TitleY, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  TString histname_scraper = "htrack_beam_inst_XY_" + key + "_scraper_KE_beam_inst" + this_KE_str + "_" + dir;
  TString histname_nonscraper = "htrack_beam_inst_XY_" + key + "_nonscraper_KE_beam_inst" + this_KE_str + "_" + dir;

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_input = new TFile(root_file_path + data_or_mc + filename);
  if(data_or_mc == "data") f_input = new TFile(root_file_path + "data_Beam_Study_1.0GeV_wider_window.root");
  gDirectory -> Cd(dir);
  for(int i = 0; i < 9; i++){
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

  if(maphist2D[histname_scraper + "_0"] == nullptr) return;

  cout << histname_scraper + "_1" << endl;
  TH2D * h_scraper;
  TH2D * h_nonscraper;
  if(data_or_mc == "data"){
    h_scraper = (TH2D*)maphist2D[histname_scraper + "_0"] -> Clone();
    h_nonscraper = (TH2D*)maphist2D[histname_nonscraper + "_0"] -> Clone();
  }
  else{
    h_scraper = (TH2D*)maphist2D[histname_scraper + "_1"] -> Clone();
    if(maphist2D[histname_scraper + "_2"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_2"]);
    
    if(maphist2D[histname_scraper + "_3"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_3"]);
    if(maphist2D[histname_scraper + "_4"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_4"]);
    if(maphist2D[histname_scraper + "_5"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_5"]);
    if(maphist2D[histname_scraper + "_6"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_6"]);
    if(maphist2D[histname_scraper + "_7"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_7"]);
    if(maphist2D[histname_scraper + "_8"] != nullptr) h_scraper -> Add(maphist2D[histname_scraper + "_8"]);

    h_nonscraper = (TH2D*)maphist2D[histname_nonscraper + "_1"] -> Clone();
    if(maphist2D[histname_nonscraper + "_2"] != nullptr) h_nonscraper -> Add(maphist2D[histname_nonscraper + "_2"]);
    
  }

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
  //h_scraper -> Draw("boxsame");
 
  TEllipse *circle;
  if(data_or_mc == "mc")  circle = new TEllipse(-29.6, 422, 1.4 * 4.8);
  if(data_or_mc == "data") circle = new TEllipse(-32.16, 422.7, 1.2 * 4.8);
  circle -> SetLineColor(kBlue);
  circle -> SetLineStyle(7);
  circle -> SetFillColorAlpha(kWhite, 0.);
  circle -> Draw("lsame");

  TEllipse *circle_HY;
  if(data_or_mc == "mc") circle_HY = new TEllipse(-29.1637, 421.76, 4.50311 * 1.5, 3.83908 * 1.5);
  if(data_or_mc == "data") circle_HY = new TEllipse(-31.3139, 422.116, 3.79366 * 1.5, 3.48005 * 1.5);
  circle_HY -> SetLineColor(kOrange);
  circle_HY -> SetLineStyle(7);
  circle_HY -> SetFillColorAlpha(kWhite, 0.);
  circle_HY -> Draw("lsame");
 
  TLegend *l = new TLegend(0.20, 0.6, 0.94, 0.94);
  l -> SetNColumns(2);
  l -> AddEntry(h_nonscraper, "Non-scraper", "lf");
  l -> AddEntry(h_scraper, "Scraper", "lf");
  l -> AddEntry(circle, "1.2 #sigma", "l");
  l -> AddEntry(circle_HY, "Heng-Ye's ellipse", "l");
  l -> Draw("same");

  TLatex latex_ProtoDUNE, latex_KE;
  latex_ProtoDUNE.SetNDC();
  latex_KE.SetNDC();
  latex_KE.SetTextAlign(31);
  latex_ProtoDUNE.SetTextSize(0.03);
  latex_KE.SetTextSize(0.03);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_KE.DrawLatex(0.95, 0.96, "KE_{Beam Inst.} = " + this_KE_legend);

  c -> SaveAs("./output/plots/BeamStudy/Beam_scraper/BeamXY_" + data_or_mc + "_scraper" + beam_P + "GeV_" + key + "_" + this_KE_str + "_" + dir + ".pdf");

  c -> Close();

}

void Run_Draw_2D_MC_and_Data(TString filename, TString dir, TString TitleX, TString TitleY, TString beam_P, double xmin, double xmax, double rebin_x, double ymin, double ymax, double rebin_y){

  int KE_step = 100;

  int pion_KE_min = 700;
  int pion_KE_max = 1100;
  int pion_N_step = (pion_KE_max - pion_KE_min) / KE_step;
  for(int i = 0; i < pion_N_step; i++){
    TString this_KE_str = Form("%dto%dMeV", pion_KE_min + i * KE_step, pion_KE_min + (i + 1) * KE_step);
    TString this_KE_legend = Form("%d - %dMeV", pion_KE_min + i * KE_step, pion_KE_min + (i + 1) * KE_step);
    Draw_XY_each_KE("mc", filename, "True", "pion_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
  }

  int proton_KE_min = 300;
  int proton_KE_max = 600;
  int proton_N_step = (proton_KE_max - proton_KE_min) / KE_step;
  for(int i = 0; i < proton_N_step; i++){
    TString this_KE_str = Form("%dto%dMeV", proton_KE_min + i * KE_step, proton_KE_min + (i + 1) * KE_step);
    TString this_KE_legend = Form("%d - %dMeV", proton_KE_min + i * KE_step, proton_KE_min + (i + 1) * KE_step);
    Draw_XY_each_KE("mc", filename, "True", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    Draw_XY_each_KE("mc", filename, "FittedElas", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    Draw_XY_each_KE("mc", filename, "FittedFakeData", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    Draw_XY_each_KE("mc", filename, "FittedData", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    
    Draw_XY_each_KE("data", filename, "FittedFakeData", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
    Draw_XY_each_KE("data", filename, "FittedData", "proton_BeamWindow_" + dir, beam_P, this_KE_str, this_KE_legend, TitleX, TitleY, xmin, xmax, rebin_x, ymin, ymax, rebin_y);
  }
}

void Draw_Beam_Scraper(){

  cout << "============================" << endl;
  cout << "=========== START ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = ".root";

  Run_Draw_2D_MC_and_Data("_Beam_Study_1.0GeV" + file_suffix, "noweight", "X_{ff} [cm]", "Y_{ff} [cm]", "1.0", -50, 10, 1., 400., 460., 1.);
}
