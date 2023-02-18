#include "canvas_margin.h"
#include "mylib.h"

void Draw_MC_daughter_shape_comparison(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool logy = false){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;
  const int N_PDGs = 4;
  TString PDGs_str[N_PDGs] = {"muon", "other", "proton", "pion"};
  Int_t colour_array[] = {867, 416, 632, 800};

  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    cout << "[Draw_MC_daughter_shape_comparison] " << this_hist_name << endl;

    if((TH1D*)gDirectory -> Get(this_hist_name)){
      cout << "[Draw_MC_daughter_shape_comparison] Found " << this_hist_name << ", " << ((TH1D*)gDirectory -> Get(this_hist_name)) -> Integral() << endl;
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      int this_N_bin = maphist[this_hist_name] -> GetNbinsX();
      double overflow = maphist[this_hist_name] -> GetBinContent(this_N_bin + 1);
      double integral= maphist[this_hist_name] -> Integral();
      maphist[this_hist_name] -> Scale(1. / (overflow + integral) );
      double last_content = maphist[this_hist_name] -> GetBinContent(this_N_bin);
      overflow = maphist[this_hist_name] -> GetBinContent(this_N_bin + 1);
      maphist[this_hist_name] -> SetBinContent(this_N_bin, last_content + overflow);
      double this_max = maphist[this_hist_name] -> GetMaximum();
      if(mc_max < this_max) mc_max = this_max;
    }
    else maphist[this_hist_name] = nullptr;
  }

  TCanvas *c = new TCanvas("","",800,800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  if(logy) c -> SetLogy();

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);

  double max_y = -1.;
  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    TString this_legend_str = PDGs_str[i];
    if(maphist[this_hist_name] != nullptr){
      double this_y = maphist[this_hist_name] -> GetMaximum();
      if(this_y > max_y) max_y = this_y;
    }
  }
  template_h -> GetYaxis() -> SetRangeUser(0., max_y * 1.8);
  if(logy) template_h -> GetYaxis() -> SetRangeUser(0.001, max_y * 100.);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  for(int i = 0; i < N_PDGs; i++){
    TString this_hist_name = "hdaughter_" + PDGs_str[i] + "_" + histname;
    TString this_legend_str = PDGs_str[i];
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetLineStyle(1);
      if(histname.Contains("chi2") || histname.Contains("nHits")) maphist[this_hist_name] -> SetLineStyle(1);
      maphist[this_hist_name] -> SetLineWidth(4);
      maphist[this_hist_name] -> Draw("histsame");
      l->AddEntry(maphist[this_hist_name], this_legend_str, "l");
    }
  }
  l -> Draw("same");

  gPad->RedrawAxis();

  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/Beam_1GeV/MC_daughter_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  if(logy) pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/Beam_1GeV/Logy_MC_daughter_comparison_" + histname + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  f_mc -> Close();
  c -> Close();

}

void Draw_MC_daughter_efficiency(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;

  const int N_particles = 4;
  TString particle_str[N_particles] = {"pion", "proton", "muon", "other"};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + particle_str[i] + "_" + histname + "_" + dir;
    TString nocut_hist_name = "hdaughter_" + particle_str[i] + "_TrueP_nocut_" + dir;

    maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
    TH1D * this_nocut_hist = (TH1D*)gDirectory -> Get(nocut_hist_name) -> Clone();

    maphist[this_hist_name] -> Rebin(rebin);
    this_nocut_hist -> Rebin(rebin);
    maphist[this_hist_name] -> Divide(this_nocut_hist);
  }


  TCanvas *c = new TCanvas("","",800,800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Efficiency");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_particles; i++){
    TString this_hist_name = "hdaughter_" + particle_str[i] + "_" + histname + "_" + dir;
    TString this_legend_str = particle_str[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineWidth(4);
    maphist[this_hist_name] -> Draw("histsame");
    l->AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> Draw("same");

  gPad->RedrawAxis();

  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/Beam_1GeV/MC_daughter_Efficiency_" + histname + "_" + beam_P + "GeV_overlap.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  f_mc -> Close();

}

void Draw_MC_daughter_cutflow(TString filename, TString dir, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, bool draw_eff){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> cd(dir);
  double mc_max = -1.;

  const int N_cuts = 7;
  TString cutflow_str[99] = {"nocut", "cos", "dist", "Nhits", "emScore", "TrackScore", "chi2_proton", "mean_dEdx"};
  TString cutflow_lgd_str[99] = {"No cut", "cos#theta_{Beam, Daughter}", "Distance(beam end, daughter start)", "N_{Hits}", "e#gamma Score", "Track Score", "#chi^{2}_{proton}", "<dE/dx>_{Truncated}"};

  if(histname.Contains("pion")){
    cutflow_str[6] = cutflow_str[6] + "_pionID";
    cutflow_str[7] = cutflow_str[7] + "_pionID";
  }
  else if(histname.Contains("proton")){
    cutflow_str[6] = cutflow_str[6] + "_protonID";
  }

  for(int i = 0; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
    }
    else maphist[this_hist_name] = nullptr;
  }

  // == Dividing by nocut hit 
  for(int i = 1; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    TString nocut_hist_name = "hdaughter_" + histname + "_" + cutflow_str[0] + "_" + dir;
    if(draw_eff && maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> Divide(maphist[nocut_hist_name]);
    }
  }
  TString nocut_hist_name = "hdaughter_" + histname + "_" + cutflow_str[0] + "_" + dir;
  if(draw_eff) maphist[nocut_hist_name] -> Divide(maphist[nocut_hist_name]);

  mc_max = maphist["hdaughter_" + histname + "_" + cutflow_str[1] + "_" + dir] -> GetMaximum();
  TString nameofhistogram = histname + "Draw_MC_daughter_cutflow" + beam_P + dir;
  if(draw_eff) nameofhistogram = nameofhistogram + "true";
  else nameofhistogram = nameofhistogram + "false";
 
  TCanvas *c = new TCanvas("","",800,800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  if(draw_eff) template_h -> GetYaxis() -> SetTitle("Eff.");
  else template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(0);
  l -> SetFillStyle(4000);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  Int_t colour_array[] = {632, 800, 867, 416, 901, 432, 400, 920};
  for(int i = 0; i < N_cuts; i++){
    TString this_hist_name = "hdaughter_" + histname + "_" + cutflow_str[i] + "_" + dir;
    TString this_legend_str = cutflow_lgd_str[i];
    if(maphist[this_hist_name] != nullptr){
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetLineWidth(4);
      maphist[this_hist_name] -> Draw("histsame");
      l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
    }
  }
  l -> Draw("same");

  gPad->RedrawAxis();

  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  if(draw_eff) pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/Beam_1GeV/MC_daughter_Cutflow_" + histname + "_" + beam_P + "GeV_eff.pdf";
  else pdfname = WORKING_DIR + "/output/plots/PionXsec/Daughters/Beam_1GeV/MC_daughter_Cutflow_" + histname + "_" + beam_P + "GeV_hist.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  f_mc -> Close();

}

void Run_Draw_MC_daughter_shape_comparison(TString file_prefix, TString file_suffix, TString P_str, TString dir){

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_not_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_dist_" + dir, "Distance(beam end, daughter start) (cm)", P_str, 0., 100., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "trackScore_" + dir, "Track Score", P_str, 0., 1., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "emScore_" + dir, "e/#gamma Score", P_str, 0., 1., 5., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "chi2_proton_" + dir, "#chi^{2}_{proton}", P_str, 0., 100., 100., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "PFP_nHits_" + dir, "N_{Hits}", P_str, 0., 1000., 10., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_pion_chi2_" + dir, "#chi^{2}_{pion}", P_str, 0., 200., 50., false);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "this_muon_chi2_" + dir, "#chi^{2}_{muon}", P_str, 0., 200., 50., false);

  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., true);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_not_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., true);
  Draw_MC_daughter_shape_comparison("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "beam_cos_broken_track_" + dir, "cos#theta(beam, daughter)", P_str, -1., 1., 2., true);
}

void Run_Draw_MC_daughter_cutflow(TString file_prefix, TString file_suffix, TString P_str, TString dir){
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "pion_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50., false);
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "proton_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1400., 50., false);

  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "pion_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50., true);
  Draw_MC_daughter_cutflow("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "proton_TrueP", "P_{True}^{Start} (MeV)", P_str, 0., 1400., 50., true);

  Draw_MC_daughter_efficiency("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "TrueP_chi2_proton_protonID", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50.);
  Draw_MC_daughter_efficiency("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, "TrueP_chi2_proton_pionID", "P_{True}^{Start} (MeV)", P_str, 0., 1200., 50.);

  
}

void Draw_Daughter_Plots(){
  TString file_suffix = "_IsQE.root";
  Run_Draw_MC_daughter_shape_comparison("PionXsec", file_suffix, "1.0", "noweight");
  Run_Draw_MC_daughter_cutflow("PionXsec", file_suffix, "1.0", "noweight"); 
}
