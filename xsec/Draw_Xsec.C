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

void Draw_Eslice_distribution_inc_true(TString filename, TString dir, TString interaction, TString TitleX, TString beam_P, double xmin, double xmax, bool logy){

  TString suffix = "hdaughter_KE_pion_";
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path = input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double max_fake = -1.;
  gDirectory -> Cd(dir);
  TString this_hist_name_int = suffix + "end_" + dir;
  TString this_hist_name_init = suffix + "init_" + dir;
  if(interaction.Contains("QE")) this_hist_name_int = suffix + "QE_int_" + dir;
  TString this_hist_name_end = suffix + "end_" + dir;
  TString this_hist_name_inc = suffix + "inc_" + dir;
  cout << "[Draw_Eslice_distribution_inc_true] Hist int : " << this_hist_name_int << endl;
  cout << "[Draw_Eslice_distribution_inc_true] Hist inc : " << this_hist_name_inc << endl;
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
  if((TH1D*)gDirectory -> Get(this_hist_name_inc + "_1")){
    maphist[this_hist_name_inc] = (TH1D*)gDirectory -> Get(this_hist_name_inc + "_1") -> Clone();
    maphist[this_hist_name_inc] -> Add((TH1D*)gDirectory -> Get(this_hist_name_inc + "_2"));
    Rebin_with_overflow(this_hist_name_inc, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_inc] = nullptr;

  if(maphist[this_hist_name_int] == nullptr || maphist[this_hist_name_inc] == nullptr) return;

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

  // == Draw Initial KE log ratio
  TH1D* h_inc = (TH1D*)maphist[this_hist_name_inc + "rebin"] -> Clone();
  TH1D* h_init = (TH1D*)maphist[this_hist_name_init + "rebin"] -> Clone();
  TH1D* h_ratio = (TH1D*)maphist[this_hist_name_inc + "rebin"] -> Clone();
  TH1D* h_diff = (TH1D*)maphist[this_hist_name_inc + "rebin"] -> Clone();
  h_diff -> Add(maphist[this_hist_name_int + "rebin"], -1);
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_ratio -> GetBinContent(i);
    double this_error = h_ratio -> GetBinError(i);
    cout <<"h_ratio before divide, " << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
  }
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_diff -> GetBinContent(i);
    double this_error = h_diff -> GetBinError(i);
    cout << "h_diff, " << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
  }
  for(int i = 1; i < N_KE_bins; i++){
    double this_ratio_content = h_ratio -> GetBinContent(i);
    double this_ratio_error = h_ratio -> GetBinError(i);
    double this_diff_content = h_diff -> GetBinContent(i);
    double this_diff_error = h_diff -> GetBinError(i);
    double this_ratio = this_ratio_content / this_diff_content;
    double approx_err = this_ratio * sqrt(pow(this_ratio_error/this_ratio_content, 2) + pow(this_diff_error / this_diff_content, 2));
    cout << "Approx divided err, " << i << ", this_ratio : " << this_ratio << ", approx_err : " << approx_err << endl;
    h_ratio -> SetBinContent(i, this_ratio);
    h_ratio -> SetBinError(i, approx_err);
  }
  //h_ratio -> Divide(h_diff);
  for(int i = 1; i < N_KE_bins - 1; i++){
    if(h_ratio -> GetBinContent(i) > 0.){
      double this_content = h_ratio -> GetBinContent(i);
      double this_error = h_ratio -> GetBinError(i);
      cout << "h_ratio after divide, " << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
      double log_content = log(this_content);
      double log_error = log(1. + this_error / this_content );
      double this_width =  KE_binning[i] - KE_binning[i - 1];
      h_ratio -> SetBinContent(i, log_content / this_width);
      h_ratio -> SetBinError(i, log_error / this_width);
    }
  }

  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_ratio -> GetBinContent(i);
    double this_error = h_ratio -> GetBinError(i);
    cout << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
  }  
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
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_xsec -> GetBinContent(i);
    cout << i << ", this_content : " << this_content << endl;
  }

  double y_max = h_xsec -> GetMaximum();
  template_h -> GetYaxis() -> SetTitle("#sigma_{Inelastic} [mb]");
  if(interaction.Contains("QE")) template_h -> GetYaxis() -> SetTitle("#sigma_{QE} [mb]");
  template_h -> GetYaxis() -> SetRangeUser(0., 1000.);
  template_h -> Draw();

  h_init -> Scale(800. / h_init -> GetMaximum() );
  h_inc -> Scale(600. / h_inc -> GetMaximum() );
  h_init -> SetLineWidth(2);
  h_inc -> SetLineWidth(2);
  h_init -> SetLineColor(kBlue);
  h_inc -> SetLineColor(kOrange);
  h_init -> Draw("histsame");
  h_inc -> Draw("histsame");

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

  TLatex latex_ProtoDUNE, latex_data_POT;
  latex_ProtoDUNE.SetNDC();
  latex_data_POT.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Xsec_" + interaction + beam_P + "GeV_true_inc.pdf";
  c -> SaveAs(pdfname);

  // == Histogram for Ninit correction
  TH1D *h_init_corr = (TH1D*)maphist[this_hist_name_init + "rebin"] -> Clone();
  //h_init_corr -> Divide(maphist[this_hist_name_inc + "rebin"]);
  for(int i = 1; i < N_KE_bins; i++){
    double this_numer = maphist[this_hist_name_init + "rebin"] -> GetBinContent(i + 1);
    double this_denom = maphist[this_hist_name_inc + "rebin"] -> GetBinContent(i);
    double this_numer_err = maphist[this_hist_name_init + "rebin"] -> GetBinError(i + 1);
    double this_denom_err = maphist[this_hist_name_inc + "rebin"] -> GetBinError(i);

    double this_ratio = 0.;
    double approx_err = 0.;
    if(this_denom > 0. && this_numer > 0.){
      this_ratio = this_numer / this_denom;
      approx_err = this_ratio * sqrt(pow(this_numer_err/this_numer, 2) + pow(this_denom_err / this_denom, 2));
    }
    double this_KE = (KE_binning[i] + KE_binning[i - 1]) / 2.0;
    double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
    double this_width =  KE_binning[i] - KE_binning[i - 1];
    h_init_corr -> SetBinContent(i, this_ratio* this_dEdx / this_width);
    h_init_corr -> SetBinError(i, approx_err* this_dEdx / this_width);
    cout << "h_init_corr, " << i << ", this_numer : " << this_numer << ", this_denom : " << this_denom << ", this_ratio : " << this_ratio << ", approx_err : " << approx_err << endl;
  }
  h_init_corr -> Scale(10000. * xsec_unit);
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_init_corr -> GetBinContent(i);
    double this_error = h_init_corr -> GetBinError(i);
    cout << "h_init_corr\t" << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
  }

  TH1D *h_log = (TH1D*)maphist[this_hist_name_inc + "rebin"] -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    double this_inc = maphist[this_hist_name_inc + "rebin"]  -> GetBinContent(i);
    double this_inc_err = maphist[this_hist_name_inc + "rebin"]  ->GetBinError(i);
    double this_int = maphist[this_hist_name_int + "rebin"]  -> GetBinContent(i);
    double this_int_err= maphist[this_hist_name_int + "rebin"]  ->GetBinError(i);
    double this_init = maphist[this_hist_name_init + "rebin"]  -> GetBinContent(i + 1);
    double this_init_err= maphist[this_hist_name_init + "rebin"]  ->GetBinError(i);
    
    double this_numer = this_inc;
    double this_denom = this_inc - this_int + this_init;
    double this_ratio = this_numer / this_denom;
    double this_err = this_ratio * sqrt( 1. / this_numer + 1./this_denom );

    if(this_ratio > 0.){
      double log_content = log(this_ratio);
      double log_error = fabs(log(1. + this_err / this_ratio ));
      double this_width =  KE_binning[i] - KE_binning[i - 1];
      double this_KE = (KE_binning[i] + KE_binning[i - 1]) / 2.0;
      double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
      h_log -> SetBinContent(i, log_content * this_dEdx / this_width);
      h_log -> SetBinError(i, log_error * this_dEdx / this_width);
    }
    else{
      h_log -> SetBinContent(i, 0.);
      h_log -> SetBinError(i, 0.);
    }
  }
  h_log -> Scale(10000. * xsec_unit);
  for(int i = 1; i < N_KE_bins; i++){
    double this_content = h_log -> GetBinContent(i);
    double this_error = h_log -> GetBinError(i);
    cout << "h_log\t" << i << ", this_content : " << this_content << ", this_error : " << this_error << endl;
  }

  h_log -> Add(h_init_corr);
  template_h -> GetYaxis() -> SetTitle("#sigma_{Inelastic} [mb]");
  if(interaction.Contains("QE")) template_h -> GetYaxis() -> SetTitle("#sigma_{QE} [mb]");
  template_h -> GetYaxis() -> SetRangeUser(0., 1000.);
  template_h -> Draw();
  g_xsec_template -> Draw("lsame");
  h_log -> SetLineColor(kGreen);
  h_log -> SetMarkerStyle(9);
  h_log -> SetMarkerColor(kGreen);
  h_log -> Draw("epsame");
  gPad->RedrawAxis();

  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Xsec_" + interaction + beam_P + "GeV_true_inc_corr.pdf";
  c -> SaveAs(pdfname);


  // == Just ratio
  h_ratio = (TH1D*)maphist[this_hist_name_int + "rebin"] -> Clone();
  h_ratio -> Divide(maphist[this_hist_name_inc + "rebin"]);
  h_xsec = (TH1D*) h_ratio -> Clone();
  for(int i = 1; i < N_KE_bins; i++){
    if(h_xsec -> GetBinContent(i) > 0.){
      double this_content = h_xsec -> GetBinContent(i);
      double this_err = h_xsec -> GetBinError(i);
      double this_KE = (KE_binning[i] + KE_binning[i - 1]) / 2.0;
      double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
      h_xsec -> SetBinContent(i, this_content * this_dEdx / 50.);
      h_xsec -> SetBinError(i, this_err * this_dEdx / 50.);
      //cout << i << ", this_content * this_dEdx : " << this_content * this_dEdx << endl;
    }
  }
  h_xsec -> Scale(10000. * xsec_unit);
  y_max = h_xsec -> GetMaximum();
  template_h -> GetYaxis() -> SetTitle("#sigma_{Inelastic} [mb]");
  if(interaction.Contains("QE")) template_h -> GetYaxis() -> SetTitle("#sigma_{QE} [mb]");
  template_h -> GetYaxis() -> SetRangeUser(0., 1000.);
  template_h -> Draw();
  g_xsec_template -> Draw("lsame");
  h_xsec -> SetLineColor(kGreen);
  h_xsec -> SetMarkerStyle(9);
  h_xsec -> SetMarkerColor(kGreen);
  h_xsec -> Draw("epsame");
  gPad->RedrawAxis();

  latex_ProtoDUNE.DrawLatex(0.16, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.69, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/Xsec/Eslice_pion_Xsec_" + interaction + beam_P + "GeV_true_inc_ratio.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
}

void Draw_thisEslice_inc_true(TString filename, TString dir, TString interaction, TString TitleX, TString beam_P, double xmin, double xmax, bool logy){

  TString suffix = "hdaughter_KE_pion_";
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path = input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double max_fake = -1.;
  gDirectory -> Cd(dir);
  TString this_hist_name_int = suffix + "end_" + dir;
  TString this_hist_name_init = suffix + "init_" + dir;
  if(interaction.Contains("QE")) this_hist_name_int = suffix + "QE_int_" + dir;
  TString this_hist_name_end = suffix + "end_" + dir;
  TString this_hist_name_inc = suffix + "inc_" + dir;
  cout << "[Draw_Eslice_distribution_inc_true] Hist int : " << this_hist_name_int << endl;
  cout << "[Draw_Eslice_distribution_inc_true] Hist inc : " << this_hist_name_inc << endl;
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
  if((TH1D*)gDirectory -> Get(this_hist_name_inc + "_1")){
    maphist[this_hist_name_inc] = (TH1D*)gDirectory -> Get(this_hist_name_inc + "_1") -> Clone();
    maphist[this_hist_name_inc] -> Add((TH1D*)gDirectory -> Get(this_hist_name_inc + "_2"));
    Rebin_with_overflow(this_hist_name_inc, N_KE_bins, KE_binning);
  }
  else maphist[this_hist_name_inc] = nullptr;

  if(maphist[this_hist_name_int] == nullptr || maphist[this_hist_name_inc] == nullptr) return;

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
  Draw_Eslice_distribution_inc_true("_" + file_prefix + "_" + P_str + "GeV" + file_suffix, dir, interaction, "KE [MeV]", P_str, 0., 1000., false);
}

void Draw_Xsec(){

  setTDRStyle();

  TString filename = "PionXsec_1.0GeV_QE_xsec.root";
  TString file_suffix = "_QE_xsec.root";
  //Run_Draw_Slice_distribution("PionXsec", file_suffix, "QE_", "1.0", "beam_window_Preweight_piOnly");
  //Run_Draw_Slice_distribution("PionXsec", file_suffix, "", "1.0", "beam_window_Preweight_piOnly");


  // == E slice
  //mc_PionXsec_1.0GeV_beam_study.root
  file_suffix = "_Eslice_test.root";
  file_suffix = "_beam_study.root";
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
