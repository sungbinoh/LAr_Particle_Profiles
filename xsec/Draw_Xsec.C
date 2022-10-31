#include "canvas_margin.h"
#include "mylib.h"

void Draw_Eslice_distribution(TString filename, TString dir, TString interaction, TString TitleX, TString beam_P, double xmin, double xmax, bool logy){

  TString suffix = "hdaughter_KE_pion_";
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path = input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double max_fake = -1.;
  gDirectory -> Cd(dir);
  TString this_hist_name_init = suffix + "init_" + dir;
  TString this_hist_name_int = suffix + "end_" + dir;
  if(interaction.Contains("QE")) this_hist_name_int = suffix + "QE_int_" + dir;
  TString this_hist_name_end = suffix + "end_" + dir;
  cout << "[Draw_Eslice_distribution] Hist init : " << this_hist_name_init << endl;
  cout << "[Draw_Eslice_distribution] Hist int : " << this_hist_name_int << endl;
  cout << "[Draw_Eslice_distribution] Hist end : " << this_hist_name_end << endl;
  if((TH1D*)gDirectory -> Get(this_hist_name_init + "_1")){
    maphist[this_hist_name_init] = (TH1D*)gDirectory -> Get(this_hist_name_init + "_1") -> Clone();
    maphist[this_hist_name_init] -> Add((TH1D*)gDirectory -> Get(this_hist_name_init + "_2"));
    Rebin_with_overflow(this_hist_name_init, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_init] = nullptr;
  if((TH1D*)gDirectory -> Get(this_hist_name_int + "_1")){
    if(!interaction.Contains("QE")) this_hist_name_int = suffix + "int_" + dir;
    maphist[this_hist_name_int] = (TH1D*)gDirectory -> Get(this_hist_name_end + "_1") -> Clone();
    Rebin_with_overflow(this_hist_name_int, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_int] = nullptr;
  if((TH1D*)gDirectory -> Get(this_hist_name_end + "_1")){
    maphist[this_hist_name_end] = (TH1D*)gDirectory -> Get(this_hist_name_end + "_1") -> Clone();
    maphist[this_hist_name_end] -> Add((TH1D*)gDirectory -> Get(this_hist_name_end + "_2"));
    Rebin_with_overflow(this_hist_name_end, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_end] = nullptr;

  if(maphist[this_hist_name_init] == nullptr || maphist[this_hist_name_int] == nullptr || maphist[this_hist_name_end] == nullptr) return; 

  // == Make Incident KE distribution
  TH1D *h_incident = (TH1D*)maphist[this_hist_name_end + "rebin"] -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    double sum_N_end = 0.;
    double sum_N_init = 0.;
    for(int j = 1; j < i + 1; j++){
      cout << "[sum_N_end] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_end += maphist[this_hist_name_end + "rebin"] -> GetBinContent(j);
    }
    for(int j = 1; j < i + 0; j++){
      cout << "[sum_N_init] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_init += maphist[this_hist_name_init + "rebin"] -> GetBinContent(j);
    }
    double this_content = sum_N_end - sum_N_init;
    double this_err = sqrt(this_content);
    h_incident -> SetBinContent(i, this_content);
    h_incident -> SetBinError(i, this_err);
    cout << "[Draw_Eslice_distribution] " << interaction << ", " << i << " : sum_N_end : " << sum_N_end << ", sum_N_init : " << sum_N_init << ", this_content : " << this_content << endl;
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  if(logy) c -> SetLogy();

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  template_h -> SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events / MeV");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> Draw("hist");

  // == Draw Initial KE distribution
  TH1D *h_init = (TH1D*) maphist[this_hist_name_init + "rebin"] -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_init -> GetBinContent(i);
    double this_error = h_init -> GetBinError(i);
    double this_width = KE_binning[i] - KE_binning[i - 1];
    double normalized_content = this_content / this_width;
    double normalized_err = this_error / this_width;
    h_init -> SetBinContent(i, normalized_content);
    h_init -> SetBinError(i, normalized_err);
  }
  double y_max = h_init -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  h_init -> SetLineColor(kGreen);
  h_init -> SetMarkerStyle(9);
  h_init -> SetMarkerColor(kGreen);
  h_init -> Draw("epsame");
  gPad->RedrawAxis();

  TLatex latex_ProtoDUNE, latex_data_POT, latex_legend;
  latex_ProtoDUNE.SetNDC();
  latex_data_POT.SetNDC();
  latex_legend.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_legend.SetTextSize(0.05);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.60, "Initial #pi^{+}");

  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Initial_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  // == Draw End KE distribution
  TH1D *h_end = (TH1D*) maphist[this_hist_name_end + "rebin"] -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_end -> GetBinContent(i);
    double this_error = h_end -> GetBinError(i);
    double this_width =  KE_binning[i] - KE_binning[i - 1];
    double normalized_content = this_content / this_width;
    double normalized_err = this_error / this_width;
    h_end -> SetBinContent(i, normalized_content);
    h_end -> SetBinError(i, normalized_err);
  }
  y_max = h_end -> GetMaximum();
  template_h ->GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();
  h_end -> SetLineColor(kGreen);
  h_end -> SetMarkerStyle(9);
  h_end -> SetMarkerColor(kGreen);
  h_end -> Draw("epsame");
  gPad->RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.60, "End #pi^{+}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_End_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  // == Draw Int KE distribution
  TH1D *h_int = (TH1D*) maphist[this_hist_name_int + "rebin"] -> Clone();
  for(int i = 0; i < N_KE_bins ; i++){
    double this_content = maphist[this_hist_name_int + "rebin"] -> GetBinContent(i);
    cout << "[Draw_Eslice_distribution] int content, " << i << ", " << this_content << endl;
  }
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_int -> GetBinContent(i);
    double this_error = h_int -> GetBinError(i);
    double this_width =  KE_binning[i] - KE_binning[i - 1];
    double normalized_content = this_content / this_width;
    double normalized_err = this_error / this_width;
    h_int -> SetBinContent(i, normalized_content);
    h_int -> SetBinError(i, normalized_err);
    //cout << "[Draw_Eslice_distribution] int content, " << i << ", " << this_content << ", " << normalized_content << endl;
  }
  y_max = h_int -> GetMaximum();
  template_h ->GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();
  h_int -> SetLineColor(kGreen);
  h_int -> SetMarkerStyle(9);
  h_int -> SetMarkerColor(kGreen);
  h_int -> Draw("epsame");
  gPad->RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.60, "Int. #pi^{+}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Interaction_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  // == Draw Inc KE distribution
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_incident -> GetBinContent(i);
    double this_error = h_incident -> GetBinError(i);
    double this_width =  KE_binning[i] - KE_binning[i - 1];
    double normalized_content = this_content / this_width;
    double normalized_err = this_error / this_width;
    h_incident -> SetBinContent(i, normalized_content);
    h_incident -> SetBinError(i, normalized_err);
  }
  y_max = h_incident -> GetMaximum();
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();
  h_incident -> SetLineColor(kGreen);
  h_incident -> SetMarkerStyle(9);
  h_incident -> SetMarkerColor(kGreen);
  h_incident -> Draw("epsame");
  gPad->RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.60, "Inc. #pi^{+}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Incident_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  // == Make Log ratio distribution
  TH1D *h_ratio = (TH1D*) h_incident -> Clone();
  TH1D *h_diff = (TH1D*) h_incident -> Clone();
  h_diff -> Add(h_int, -1);
  h_ratio -> Divide(h_diff);
  int N_bin = h_ratio -> GetNbinsX();
  for(int i = 1; i < N_KE_bins; i++){
    if(h_ratio -> GetBinContent(i) > 0.){
      double this_content = h_ratio -> GetBinContent(i);
      double this_error = h_ratio -> GetBinError(i);
      double log_content = log(this_content);
      double log_error = log(1. + this_error / this_content );
      double this_width =  KE_binning[i] - KE_binning[i - 1];
      h_ratio -> SetBinContent(i, log_content / this_width);
      h_ratio -> SetBinError(i, log_error / this_width);
    }
  }
  y_max = h_ratio -> GetMaximum();
  template_h -> GetYaxis() -> SetTitle("#frac{1}{#DeltaE}log#frac{N_{inc.}}{N_{inc.} - N_{int.}}");
  template_h ->GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h ->GetYaxis() -> SetRangeUser(0., 0.02);
  template_h -> Draw();
  h_ratio -> SetLineColor(kGreen);
  h_ratio -> SetMarkerStyle(9);
  h_ratio -> SetMarkerColor(kGreen);
  h_ratio -> Draw("epsame");
  gPad->RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.60, "#frac{1}{#DeltaE}log#frac{N_{inc.}}{N_{inc.} - N_{int.}}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Log_ratio_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  // == Make Cross section distribution
  TH1D *h_xsec = (TH1D*) h_ratio -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    if(h_xsec -> GetBinContent(i) > 0.){
      double this_content = h_xsec -> GetBinContent(i);
      double this_err = h_xsec -> GetBinError(i);
      double this_KE = (KE_binning[i] + KE_binning[i - 1]) / 2.0;
      double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
      h_xsec -> SetBinContent(i, this_content * this_dEdx);
      h_xsec -> SetBinError(i, this_err * this_dEdx);
    }
  }
  h_xsec -> Scale(10000. * xsec_unit);
  y_max = h_xsec -> GetMaximum();
  template_h -> GetYaxis() -> SetTitle("#sigma_{Inelastic} [mb]");
  if(interaction.Contains("QE")) template_h -> GetYaxis() -> SetTitle("#sigma_{QE} [mb]");
  template_h -> GetYaxis() -> SetRangeUser(0., 1000.);
  template_h -> Draw();
  
  TFile *f_xsec_template = new TFile(input_file_dir + "/xsec/exclusive_xsec.root");
  TGraph *g_xsec_template = (TGraph*) gDirectory -> Get("total_inel_KE") -> Clone();
  if(interaction.Contains("QE")) g_xsec_template = (TGraph*) gDirectory -> Get("inel_KE") -> Clone();
  g_xsec_template -> SetLineColor(kRed);
  g_xsec_template -> SetLineWidth(2);
  g_xsec_template -> Draw("lsame");
  h_xsec -> SetLineColor(kGreen);
  h_xsec -> SetMarkerStyle(9);
  h_xsec -> SetMarkerColor(kGreen);
  h_xsec -> Draw("epsame");
  gPad->RedrawAxis();
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  //latex_legend.DrawLatex(0.75, 0.75, "#sigma");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Xsec_" + interaction + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
}

void Draw_Eslice_distribution_inc_true(TString filename, TString dir, TString interaction, TString true_gr_name, TString TitleX, TString beam_P, double xmin, double xmax, double ymax, bool logy){

  TString suffix = "hdaughter_KE_pion_";
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path = input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double max_fake = -1.;
  gDirectory -> Cd(dir);

  TString this_hist_name_init = suffix + "init_" + dir;
  TString this_hist_name_inc = suffix + "inc_" + dir;
  TString this_hist_name_int = suffix + "int_" + interaction + "_" + dir;
  TString this_hist_name_end =  suffix + "end_" + dir;

  cout << "[Draw_Eslice_distribution_inc_true] Hist int : " << this_hist_name_int << endl;

  if((TH1D*)gDirectory -> Get(this_hist_name_init + "_211")){
    maphist[this_hist_name_init] = (TH1D*)gDirectory -> Get(this_hist_name_init + "_211") -> Clone();
    //maphist[this_hist_name_init] -> Add((TH1D*)gDirectory -> Get(this_hist_name_init + "_2"));
    Rebin_with_overflow( this_hist_name_init, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_init] = nullptr;
  if((TH1D*)gDirectory -> Get(this_hist_name_inc + "_211")){
    maphist[this_hist_name_inc] = (TH1D*)gDirectory -> Get(this_hist_name_inc + "_211") -> Clone();
    //maphist[this_hist_name_inc] -> Add((TH1D*)gDirectory -> Get(this_hist_name_inc + "_2"));
    Rebin_with_overflow(this_hist_name_inc, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_inc] = nullptr;
  if((TH1D*)gDirectory -> Get(this_hist_name_int + "_211")){
    maphist[this_hist_name_int] = (TH1D*)gDirectory -> Get(this_hist_name_int + "_211") -> Clone();
    //maphist[this_hist_name_int] -> Add((TH1D*)gDirectory -> Get(this_hist_name_int + "_2"), -1);
    Rebin_with_overflow(this_hist_name_int, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_int] = nullptr;
  if((TH1D*)gDirectory -> Get(this_hist_name_end + "_211")){
    maphist[this_hist_name_end] = (TH1D*)gDirectory -> Get(this_hist_name_end + "_211") -> Clone();
    Rebin_with_overflow(this_hist_name_end, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_end] = nullptr;

  if(maphist[this_hist_name_init] == nullptr || maphist[this_hist_name_inc] == nullptr || maphist[this_hist_name_int] == nullptr || maphist[this_hist_name_end] == nullptr) return;
  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  if(logy) c -> SetLogy();

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  template_h -> SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events / MeV");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> Draw("hist");

  //if(!interaction.Contains("InElas")) maphist[this_hist_name_int + "rebin"] -> Scale(0.5); // FIXME : after runnign with if(pitype_str != "0") cuts
  TH1D* h_init = (TH1D*)maphist[this_hist_name_init + "rebin"] -> Clone();
  TH1D* h_inc = (TH1D*)maphist[this_hist_name_inc + "rebin"] -> Clone();
  TH1D* h_int = (TH1D*)maphist[this_hist_name_int + "rebin"] -> Clone();

  // == Draw initial
  TH1D *template_h_KE = new TH1D("", "", 1, xmin, xmax + 200.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  template_h_KE -> SetLineWidth(2);
  template_h_KE -> SetStats(0);
  template_h_KE -> GetXaxis() -> SetTitle(TitleX);
  template_h_KE -> GetXaxis() -> SetTitleSize(0.05);
  template_h_KE -> GetXaxis() -> SetLabelSize(0.035);
  template_h_KE -> GetYaxis() -> SetTitle("Events / MeV");
  template_h_KE -> GetYaxis() -> SetTitleSize(0.05);
  template_h_KE -> GetYaxis() -> SetLabelSize(0.035);
  template_h_KE -> GetYaxis() -> SetRangeUser(0., 1200.);
  template_h_KE -> Draw("hist");
  double y_max = h_init -> GetMaximum();
  template_h_KE -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h_KE -> Draw("hist");
  h_init -> SetLineColor(kGreen);
  h_init -> SetMarkerStyle(9);
  h_init -> SetMarkerColor(kGreen);
  h_init -> Draw("epsame");
  TLatex latex_ProtoDUNE, latex_data_POT, latex_textbox;
  latex_ProtoDUNE.SetNDC();
  latex_data_POT.SetNDC();
  latex_textbox.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_textbox.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_textbox.DrawLatex(0.30, 0.65, "Initial #pi^{+}");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Initial_" + beam_P + "GeV_true_all.pdf";
  c -> SaveAs(pdfname);

  // == Draw incident
  y_max = h_inc -> GetMaximum();
  template_h_KE -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h_KE -> Draw("hist");
  h_inc -> SetLineColor(kGreen);
  h_inc -> SetMarkerStyle(9);
  h_inc -> SetMarkerColor(kGreen);
  h_inc -> Draw("epsame");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_textbox.DrawLatex(0.30, 0.65, "Incident #pi^{+}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Incident_" + beam_P + "GeV_true_all.pdf";
  c -> SaveAs(pdfname);

  // == Draw interaction
  y_max = h_int -> GetMaximum();
  template_h_KE -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h_KE -> Draw("hist");
  h_int -> SetLineColor(kGreen);
  h_int -> SetMarkerStyle(9);
  h_int -> SetMarkerColor(kGreen);
  h_int -> Draw("epsame");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_textbox.DrawLatex(0.20, 0.65, "Interaction #pi^{+} (" + interaction + ")");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Interaction_" + interaction + "_" + beam_P + "GeV_true_all.pdf";
  c -> SaveAs(pdfname);


  // == Draw cross section
  template_h -> Draw();
  TH1D * h_xsec_all = Make_cross_section_histogram(this_hist_name_inc + "rebin", this_hist_name_int + "rebin");
  TH1D * h_maden_inc = Make_incident_histogram(this_hist_name_init + "rebin", this_hist_name_end + "rebin");

  TFile *f_xsec_template = new TFile(input_file_dir + "/xsec/exclusive_xsec.root");
  TString this_gr_name = true_gr_name;
  TGraph *g_xsec_template;
  /*
  if(this_gr_name == ""){
    this_gr_name = "total_inel_KE";
    vector<double> x_qe;
    vector<double> y_qe;
    int N_gr_points = (TGraph*) gDirectory -> Get(this_gr_name) -> GetN();
    for(int i = 0; i < N_gr_points; i++){
      double out_x, out_y;
      double this_y = 0.;
      (TGraph*) gDirectory -> Get(this_gr_name) -> GetPoint(i,out_x, out_y); 
      (TGraph*) gDirectory -> Get("abs_KE") -> GetPoint(i,out_x, this_y);
      out_y = out_y - this_y;
      (TGraph*) gDirectory -> Get("") -> GetPoint(i,out_x, this_y);
      out_y = out_y - this_y;
      (TGraph*) gDirectory -> Get("") -> GetPoint(i,out_x, this_y);
      out_y = out_y - this_y;
      (TGraph*) gDirectory -> Get("") -> GetPoint(i,out_x, this_y);
      out_y = out_y - this_y;
      (TGraph*) gDirectory -> Get("") -> GetPoint(i,out_x, this_y);
      out_y = out_y - this_y;

    }
  }
  */

  g_xsec_template = (TGraph*) gDirectory -> Get(this_gr_name) -> Clone();
  g_xsec_template -> SetLineColor(kRed);
  g_xsec_template -> SetLineWidth(2);
  g_xsec_template -> Draw("lsame");
  h_xsec_all -> SetLineColor(kGreen);
  h_xsec_all -> SetMarkerStyle(9);
  h_xsec_all -> SetMarkerColor(kGreen);
  h_xsec_all -> Draw("epsame");
  gPad->RedrawAxis();

  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_textbox.DrawLatex(0.60, 0.75, "#sigma_{#pi^{+}} (" + interaction + ")");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Xsec_" + interaction + "_" + beam_P + "GeV_true_all.pdf";
  c -> SaveAs(pdfname);


  // == Compare inc histograms
  y_max = h_inc -> GetMaximum(); 
  template_h_KE -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h_KE -> Draw("hist");
  h_inc -> SetLineColor(kGreen);
  h_inc -> SetMarkerStyle(9);
  h_inc -> SetMarkerColor(kGreen);
  //h_inc -> SetLineWidth(2);
  h_inc -> Draw("epsame");
  h_maden_inc -> SetLineColor(kRed);
  h_maden_inc -> SetMarkerStyle(9);
  h_maden_inc -> SetMarkerColor(kGreen);
  h_maden_inc -> Draw("histsame");
  h_maden_inc -> SetLineWidth(2);
  h_inc -> Draw("epsame");
  TLegend *l = new TLegend(0.2, 0.7, 0.6, 0.9);
  l -> AddEntry(h_inc, "Original Inc.", "lp");
  l -> AddEntry(h_maden_inc, "Produced Inc. Using Init. and End", "l");
  l -> Draw("same");
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_textbox.DrawLatex(0.30, 0.65, "Incident #pi^{+}");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Incident_Comparison_" + beam_P + "GeV_true_all.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
}

void test_Eslice_assumptions(TString graph_str, TString legend_str, TString interaction_def, TString beam_P){
  TString input_file_dir = getenv("LArProf_WD");
  TFile *f_xsec_template = new TFile(input_file_dir + "/xsec/exclusive_xsec.root");
  TGraph *g_xsec_template = (TGraph*) gDirectory -> Get(graph_str) -> Clone();
  double KE_low = 0.;
  double KE_high = 1200.;
  double KE_step = 5.;
  double N_steps = (KE_high - KE_low) / KE_step;
  vector<double> KE_vec;
  vector<double> y_vec;
  double y_max = 0.;
  for(int i = 0; i < N_steps; i++){
    double this_KE = KE_low + KE_step * (i + 0.);
    double this_xsec = g_xsec_template -> Eval(this_KE);
    double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
    double this_y = this_xsec / this_dEdx;
    KE_vec.push_back(this_KE);
    y_vec.push_back(this_y);
    if(y_max < this_y) y_max = this_y;
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, 0., 1000.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  template_h -> SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#sigma #frac{1}{#frac{dE}{dx}} (mb cm / MeV)");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw("hist");

  TGraph *g_xsec_over_dEdx = new TGraph(N_steps, &KE_vec[0], &y_vec[0]);
  g_xsec_over_dEdx -> SetLineColor(kRed);
  g_xsec_over_dEdx -> SetLineWidth(2);
  g_xsec_over_dEdx -> Draw("lsame");
  gPad->RedrawAxis();

  TLatex latex_ProtoDUNE, latex_data_POT, latex_legend, latex_interaction;
  latex_ProtoDUNE.SetNDC();
  latex_data_POT.SetNDC();
  latex_legend.SetNDC();
  latex_interaction.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_legend.SetTextSize(0.05);
  latex_interaction.SetTextSize(0.05);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.80, "#pi^{+} " + legend_str);
  latex_interaction.DrawLatex(0.25, 0.75, interaction_def);

  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");  
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Sigma_over_dEdx_pion_" + graph_str + ".pdf";
  c -> SaveAs(pdfname);

  c -> Close(); 
}

void test_impact_of_assumption(TString graph_str, TString legend_str, TString interaction_def, double bin_size, TString bin_size_str, TString beam_P){

  TString input_file_dir = getenv("LArProf_WD");
  TFile *f_xsec_template = new TFile(input_file_dir + "/xsec/exclusive_xsec.root");
  TGraph *g_xsec_template = (TGraph*) gDirectory -> Get(graph_str) -> Clone();
  double KE_low = 0.;
  double KE_high = 1200.;
  double KE_step = bin_size;
  int N_steps = (KE_high - KE_low) / KE_step;
  vector<double> KE_vec;
  vector<double> KE_err_vec;  
  vector<double> y_vec;
  vector<double> y_err_vec;
  vector<double> approx_xsec_vec;
  double y_max = 0.;
  for(int i = 0; i < N_steps; i++){
    double this_KE = KE_low + KE_step * (i + 0.);
    
    double this_integral = 0.;
    double integ_step = 0.1;
    int N_integ_steps = KE_step / integ_step;
    for(int j = 0; j < N_integ_steps; j++){
      double this_integ_KE = this_KE + integ_step * (j + 0.);
      double this_xsec = g_xsec_template -> Eval(this_integ_KE);
      double this_dEdx = dEdx_Bethe_Bloch(this_integ_KE, mass_pion);
      double this_area = integ_step * this_xsec / this_dEdx;
      this_integral += this_area;
    }
    
    double this_KE_center = this_KE + KE_step / (2.0);
    double this_dEdx_center = dEdx_Bethe_Bloch(this_KE_center, mass_pion);
    double this_approx_xsec = this_integral * this_dEdx_center / KE_step;
    KE_vec.push_back(this_KE_center);
    approx_xsec_vec.push_back(this_approx_xsec);
    KE_err_vec.push_back(KE_step / 2.0); 
    y_err_vec.push_back(0.);
    if(y_max < this_approx_xsec) y_max = this_approx_xsec;
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, 0., 1000.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  template_h -> SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle("KE [MeV]");
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("#sigma (mb)");
  template_h -> GetYaxis() -> SetTitleSize(0.05);
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw("hist");

  g_xsec_template -> SetLineColor(kRed);
  g_xsec_template -> SetLineWidth(2);
  g_xsec_template -> Draw("lsame");

  TGraphErrors *g_approx_xsec = new TGraphErrors(N_steps, &KE_vec[0], &approx_xsec_vec[0], &KE_err_vec[0], &y_err_vec[0]);
  g_approx_xsec -> SetMarkerStyle(9);
  g_approx_xsec -> SetMarkerColor(kBlue);
  g_approx_xsec -> SetLineColor(kBlue);
  //g_approx_xsec -> SetLineStyle(7);
  g_approx_xsec -> SetLineWidth(1);
  g_approx_xsec -> Draw("epsame");
    
  TLegend *l = new TLegend(0.60, 0.70, 0.93, 0.93);
  l -> AddEntry(g_xsec_template, "GEANT4 Expectation", "l");
  l -> AddEntry(g_approx_xsec, "Approx #sigma", "l");
  l -> Draw("same");
  gPad->RedrawAxis();

  TLatex latex_ProtoDUNE, latex_data_POT, latex_legend, latex_interaction;
  latex_ProtoDUNE.SetNDC();
  latex_data_POT.SetNDC();
  latex_legend.SetNDC();
  latex_interaction.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_legend.SetTextSize(0.05);
  latex_interaction.SetTextSize(0.05);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_legend.DrawLatex(0.25, 0.80, "#pi^{+} " + legend_str);
  latex_interaction.DrawLatex(0.25, 0.75, interaction_def);

  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Approx_sigma_pion_" + graph_str + "_" + bin_size_str + ".pdf";
  c -> SaveAs(pdfname);

  c -> Close();
  delete l;
}

void Run_Draw_Eslice(TString file_prefix, TString file_suffix, TString interaction, TString P_str, TString dir){
  //Draw_Eslice_distribution("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, interaction, "KE [MeV]", P_str, 0., 1200., false);
  //Draw_Eslice_distribution_inc_true("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, interaction, "KE [MeV]", P_str, 0., 1000., false);
  TString interactions_str[7] = {"InElas", "Absorption", "ChargeEx", "DoubleChargeEx", "PiProd", "QuasiElas", "EQE_pass"};
  TString true_gr_names[7] = {"total_inel_KE", "abs_KE", "cex_KE", "dcex_KE", "prod_KE", "inel_KE", "inel_KE"};
  double y_maxs[7] = {1200., 1000., 500., 400., 500., 600., 600.};
  for(int i = 0; i < 7; i++){
    Draw_Eslice_distribution_inc_true("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, interactions_str[i], true_gr_names[i], "KE [MeV]", P_str, 0., 1000., y_maxs[i], false);
  }
}

void Draw_Xsec(){

  setTDRStyle();

  TString filename = "PionXsec_1.0GeV_QE_xsec.root";
  TString file_suffix = "_QE_xsec.root";
  //Run_Draw_Slice_distribution("PionXsec", file_suffix, "QE_", "1.0", "beam_window_Preweight_piOnly");
  //Run_Draw_Slice_distribution("PionXsec", file_suffix, "", "1.0", "beam_window_Preweight_piOnly");


  // == E slice
  //mc_PionXsec_1.0GeV_Eslice_test.root
  //mc_PionXsec_1.0GeV_beam_study.root
  file_suffix = "_Eslice_test.root";
  //Run_Draw_Eslice("PionXsec", file_suffix, "QE_", "1.0", "beam_window_Preweight_piOnly");
  //Run_Draw_Eslice("PionXsec", file_suffix, "", "1.0", "beam_window_Preweight_piOnly");
  Run_Draw_Eslice("PionXsec", file_suffix, "", "1.0", "nocut_noweight");
    


  TString xsec_template_strings[6] = {"abs_KE", "inel_KE", "cex_KE", "dcex_KE", "prod_KE", "total_inel_KE"};
  TString xsec_legend_strings[6] = {"Absorption",
				    "Quasi Elastic",
				    "Charge Exchange",
				    "Double Charge Exchange",
				    "Pion Production",
				    "Total Inelastic"};
  TString interaction_definition_strings[6] = {"(#pi^{+} + N #rightarrow N' + X)",
					       "(#pi^{+} + N #rightarrow #pi^{+} + N' + X)",
					       "(#pi^{+} + N #rightarrow #pi^{0} + N' + X)",
					       "(#pi^{+} + N #rightarrow #pi^{-} + N' + X)",
					       "(#pi^{+} + N #rightarrow n#pi + N' + X)",
					       ""
  };
  for(int i = 0; i < 6; i++){
    //test_Eslice_assumptions(xsec_template_strings[i], xsec_legend_strings[i], interaction_definition_strings[i], "1.0");
    //test_impact_of_assumption(xsec_template_strings[i], xsec_legend_strings[i], interaction_definition_strings[i], 20., "20MeV", "1.0");
    //test_impact_of_assumption(xsec_template_strings[i], xsec_legend_strings[i], interaction_definition_strings[i], 50., "50MeV", "1.0");
    //test_impact_of_assumption(xsec_template_strings[i], xsec_legend_strings[i], interaction_definition_strings[i], 100., "100MeV", "1.0");
    //test_impact_of_assumption(xsec_template_strings[i], xsec_legend_strings[i], interaction_definition_strings[i], 200., "200MeV", "1.0");
  }
}
