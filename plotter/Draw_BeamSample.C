#include "canvas_margin.h"
#include "mylib.h"

void Draw_MC_shape_comparison(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, double ymax){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  double mc_max = -1.;
  double muon_max = -1.;
  const int N_stages = 4;
  //TString stages[N_stages] = {"AfterTarget", "COLL1", "BPROF1", "BPROF2", "BPROF3", "TRIG1", "BPROFEXT", "BPROF4", "TRIG2", "NP04front", "NP04FieldCage"};
  //TString stages[N_stages] = {"AfterTarget", "COLL1", "BPROF3", "TRIG1", "BPROFEXT", "BPROF4", "TRIG2", "NP04front", "NP04FieldCage"};
  TString stages[N_stages] = {"AfterTarget", "BPROF3", "NP04front", "NP04FieldCage"};


  for(int i = 0; i < N_stages; i++){
    TString this_hist_name_all = histname + "_" + stages[i];
    TString this_hist_name_pion = histname + "_" + stages[i] + "_PID211";
    TString this_hist_name_muon = histname + "_" + stages[i] + "_PID13";

    maphist[this_hist_name_all] = (TH1D*)gDirectory -> Get(this_hist_name_all) -> Clone();
    maphist[this_hist_name_pion] = (TH1D*)gDirectory -> Get(this_hist_name_pion) -> Clone();
    maphist[this_hist_name_muon] = (TH1D*)gDirectory -> Get(this_hist_name_muon) -> Clone();

    maphist[this_hist_name_all] -> Rebin(rebin);
    maphist[this_hist_name_pion] -> Rebin(rebin);
    maphist[this_hist_name_muon] -> Rebin(rebin);

    //maphist[this_hist_name] -> Scale(1. / maphist[this_hist_name] -> Integral());
    double this_max = maphist[this_hist_name_all] -> GetMaximum();
    double this_muon_max = maphist[this_hist_name_muon]-> GetMaximum();
    if(mc_max < this_max) mc_max = this_max;
    if(muon_max < this_muon_max) muon_max = this_muon_max;
  }

  TString nameofhistogram = histname + "Draw_MC_shape_comparison" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",600,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetLogy(); // == logy
  
  // == Draw all
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  //template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> GetYaxis() -> SetRangeUser(1., mc_max * 100.); // == logy
  //template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.20, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> SetNColumns(3);

  //Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  Int_t colour_array[] = {0, 632, 800, 416, 600};
  for(int i = 0; i < N_stages; i++){
    TString this_hist_name = histname + "_" + stages[i];
    TString this_legend_str = stages[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i + 1]);
    maphist[this_hist_name] -> SetLineStyle(1);
    maphist[this_hist_name] -> SetLineWidth(3);
    maphist[this_hist_name] -> Draw("histsame");
    maplegend[legend]->AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();
  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/BeamSample_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  // == Draw pion
  template_h -> Draw();
  for(int i = 0; i < N_stages; i++){
    TString this_hist_name = histname + "_" + stages[i] + "_PID211";
    maphist[this_hist_name] -> SetLineColor(colour_array[i + 1]);
    maphist[this_hist_name] -> SetLineStyle(1);
    maphist[this_hist_name] -> SetLineWidth(3);
    maphist[this_hist_name] -> Draw("histsame");
  }
  maplegend[legend] -> Draw("same");
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  gPad->RedrawAxis();
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/BeamSample_" + histname + "_pion_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  // == Draw muon
  
  template_h -> Draw();
  for(int i = 0; i < N_stages; i++){
    TString this_hist_name = histname + "_" + stages[i] + "_PID13";
    maphist[this_hist_name] -> SetLineColor(colour_array[i + 1]);
    maphist[this_hist_name] -> SetLineStyle(1);
    maphist[this_hist_name] -> SetLineWidth(3);
    maphist[this_hist_name] -> Draw("histsame");
  }
  template_h -> GetYaxis() -> SetRangeUser(1.0, muon_max * 100.);
  maplegend[legend] -> Draw("same");
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  gPad->RedrawAxis();
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/BeamSample_" + histname + "_muon_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  f_mc -> Close();

}

void Draw_pion_muon_comparison(TString filename, TString histname, TString TitleX, TString beam_P, double xmin, double xmax, double rebin, double ymax){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TH1D* pion_hist = (TH1D*)gDirectory -> Get(histname + "_PID211");
  TH1D* muon_hist = (TH1D*)gDirectory ->Get(histname + "_PID13");
  pion_hist -> Rebin(rebin);
  muon_hist -> Rebin(rebin);
  double mc_max = pion_hist -> GetMaximum();

  TString nameofhistogram = histname + "Draw_MC_shape_comparison" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",600,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetLogy(); // == logy

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0.1, mc_max * 10.); // == logy
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.60, 0.70, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);

  pion_hist -> SetLineColor(kBlue);
  muon_hist -> SetLineColor(kRed);

  pion_hist -> SetLineWidth(3);
  muon_hist -> SetLineWidth(3);

  pion_hist -> Draw("histsame");
  muon_hist -> Draw("histsame");

  cout << "N pion : " << pion_hist -> Integral() << ", N muon : " << muon_hist -> Integral() << ", ratio : " << muon_hist -> Integral()/pion_hist -> Integral() << endl;

  maplegend[legend] -> AddEntry(pion_hist, "#pi^{+}", "l");
  maplegend[legend] -> AddEntry(muon_hist, "#mu^{+}", "l");
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();
  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_histname;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_histname.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_histname.SetTextSize(0.050);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_histname.DrawLatex(0.55, 0.6, histname);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/BeamSample_pion_vs_muon_" + histname + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);
}

void Draw_VirtualDetector(TString filename, TString histname, TString stage, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TString hist_pion_str = "VirtualDetector_" + histname + "_" + stage + "_211";
  TString hist_muon_str = "VirtualDetector_" + histname+ "_" +stage +"_13";

  TH1D* pion_hist = (TH1D*)gDirectory -> Get(hist_pion_str) -> Clone();
  TH1D* muon_hist = (TH1D*)gDirectory ->Get(hist_muon_str) -> Clone();
  pion_hist -> Rebin(rebin);
  muon_hist -> Rebin(rebin);
  double pion_max = pion_hist -> GetMaximum();
  double muon_max = muon_hist -> GetMaximum();
  double mc_max = max(pion_max, muon_max);

  double N_pion = pion_hist -> Integral();
  double N_muon = muon_hist -> Integral();
  double ratio_pion = N_pion / (N_pion + N_muon);
  double ratio_muon = N_muon / (N_pion + N_muon);
  TString ratio_pion_str = Form("%.3f", ratio_pion);
  TString ratio_muon_str = Form("%.3f", ratio_muon); 
  cout << "[" << stage << "] pion : muon = " << ratio_pion << " : " << ratio_muon << endl;
  
  TString nameofhistogram = histname + "Draw_VirtualDetector" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",600,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  //mapcanvas[canvas] -> SetLogy(); // == logy 

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0.1, mc_max * 10.); // == logy
  template_h -> GetYaxis() -> SetRangeUser(0.1, 100000);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.60, 0.80, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);

  pion_hist -> SetLineColor(kBlue);
  muon_hist -> SetLineColor(kRed);

  pion_hist -> SetLineWidth(3);
  muon_hist -> SetLineWidth(3);

  pion_hist -> Draw("histsame");
  muon_hist -> Draw("histsame");

  maplegend[legend] -> AddEntry(pion_hist, "#pi^{+}", "l");
  maplegend[legend] -> AddEntry(muon_hist, "#mu^{+}", "l");
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();
  mapcanvas[canvas] -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_histname, latex_ratio;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_histname.SetNDC();
  latex_ratio.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_histname.SetTextSize(0.050);
  latex_ratio.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_histname.DrawLatex(0.55, 0.75, stage);
  latex_ratio.DrawLatex(0.55, 0.70, "#pi^{+} : #mu^{+} = " + ratio_pion_str + " : " + ratio_muon_str);
  TString pdfname;  
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/VirtualDetector_pion_vs_muon_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);
  mapcanvas[canvas]  -> Close();
}

void Study_true_beam(TString filename, TString histname, TString stage, TString TitleX, TString beam_P, double xmin, double xmax, double rebin){
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  TString hist_pion_str = "VirtualDetector_" + histname + "_" + stage + "_211";
  TString hist_muon_str = "VirtualDetector_" + histname+ "_" +stage +"_13";

  TH1D* pion_hist = (TH1D*)gDirectory -> Get(hist_pion_str) -> Clone();
  TH1D* muon_hist = (TH1D*)gDirectory ->Get(hist_muon_str) -> Clone();
  pion_hist -> Rebin(rebin);
  muon_hist -> Rebin(rebin);
  double pion_max = pion_hist -> GetMaximum();
  double muon_max = muon_hist -> GetMaximum();
  double mc_max = max(pion_max, muon_max);

  double N_pion = pion_hist -> Integral();
  double N_muon = muon_hist -> Integral();
  double ratio_pion = N_pion / (N_pion + N_muon);
  double ratio_muon = N_muon / (N_pion + N_muon);
  TString ratio_pion_str = Form("%.3f", ratio_pion);
  TString ratio_muon_str = Form("%.3f", ratio_muon);
  cout << "[" << stage << "] pion : muon = " << ratio_pion << " : " << ratio_muon << endl;

  TString nameofhistogram = histname + "Study_true_beam" + beam_P;
  TString canvas = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",600,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleSize(0.05);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h -> GetYaxis() -> SetRangeUser(0., 2000.);
  template_h -> Draw();

  maplegend[legend] = new TLegend(0.60, 0.80, 0.90, 0.90);
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);

  pion_hist -> SetLineColor(kGreen);
  muon_hist -> SetLineColor(kBlue);

  pion_hist -> SetLineWidth(3);
  muon_hist -> SetLineWidth(3);

  pion_hist -> Draw("histsame");
  muon_hist -> Draw("histsame");

  maplegend[legend] -> AddEntry(pion_hist, "#pi^{+}", "l");
  maplegend[legend] -> AddEntry(muon_hist, "#mu^{+}", "l");
  maplegend[legend] -> Draw("same");

  TF1 *gaus_pion = new TF1("gaus_pion", "gaus", 800., 1200.);
  pion_hist -> Fit(gaus_pion, "R", "", 800., 1200.);
  gaus_pion -> SetLineColor(kSpring+3);
  gaus_pion -> SetLineStyle(7);
  gaus_pion -> Draw("lsame");

  TF1 *gaus_muon_right_tail =  new TF1("gaus_muon_right_tail", "gaus", 950., 1200.);
  muon_hist -> Fit(gaus_muon_right_tail, "R", "", 950., 1200.);
  gaus_muon_right_tail -> SetLineColor(kCyan);
  gaus_muon_right_tail -> SetLineStyle(7);
  //gaus_muon_right_tail -> Draw("lsame");

  TF1 *gaus_subtrac = new TF1("gaus_subtrac", "gaus", 800., 1200.);


  // == MC generation of decay-in-flight muon P
  int N_trial = 10000000.;
  double P_muon_CoM = (mass_pion * mass_pion - mass_muon * mass_muon) / (2.0 * mass_pion);
  double E_muon_CoM = sqrt(P_muon_CoM * P_muon_CoM + mass_muon * mass_muon);
  int N_bin = (xmax - xmin) / (10. * rebin);
  TH1D* expected_DiF_muon_P = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D* h_muon_fit_right_tail = new TH1D("", "", N_bin, xmin, xmax);
  TH1D* h_muon_fit_right_tail_basis = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D* h_pion_fit_basis = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TF1 *gaus_muon_right_tail_ext = new TF1("gaus_muon_right_tail_ext", "gaus", 0., 1200.);
  gaus_muon_right_tail_ext -> SetParameter(0, gaus_muon_right_tail -> GetParameter(0));
  gaus_muon_right_tail_ext -> SetParameter(1, gaus_muon_right_tail -> GetParameter(1));
  gaus_muon_right_tail_ext -> SetParameter(2, gaus_muon_right_tail -> GetParameter(2));  
  TH1D* h_cos_theta_CoM = new TH1D("cos_theta_CoM", "cos_theta_CoM", 200., -2., 2.);
  TH1D* h_theta_lab = new TH1D("theta_lab", "theta_lab", 4000., 0., 4.0);
  int N_not_decay = 0.;
  TH1D * expected_DiF_muon_P_2p5 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_2p4 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_2p2 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_1p8 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_1p2 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_0p5 = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_1p8_right = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_1p2_right = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * expected_DiF_muon_P_0p5_right = new TH1D("", "", 10000. / rebin, 0., 100000.);
  TH1D * sum_basis = new TH1D("", "", 10000. / rebin, 0., 100000.);
  for(int i = 0; i < N_trial; i++){
    double this_pion_P = gaus_pion -> GetRandom(800., 1200.);
    
    // == CoM
    double cos_ran = (2.0 * gRandom->Rndm() - 1.0);
    double sin_ran = sqrt(1. - cos_ran * cos_ran);
    double Pz_muon_CoM = P_muon_CoM * cos_ran;
    double Px_muon_CoM = P_muon_CoM * sin_ran;
    h_cos_theta_CoM -> Fill(cos_ran);

    // == Lorentz boost and decay using random number
    double gamma = sqrt(this_pion_P * this_pion_P + mass_pion * mass_pion ) / mass_pion;
    double ctau = 7.827; // = [m]
    double gammactau = gamma * ctau;
    double distance_BPROF3to4 = 14.5; // = [m];
    double P_decay = 1. - exp(-1. * distance_BPROF3to4 / gammactau);
    //cout << "P_decay : " << P_decay << endl;
    bool decay = (gRandom->Rndm() < P_decay);

    double gamma_beta = this_pion_P / mass_pion;
    double Pz_muon_lab = gamma * Pz_muon_CoM + gamma_beta * E_muon_CoM;
    double P_muon_lab = sqrt(Pz_muon_lab * Pz_muon_lab + Px_muon_CoM * Px_muon_CoM);
    double E_muon_lab = gamma * E_muon_CoM + gamma_beta * Pz_muon_CoM;;
    double cos_theta_lab = Pz_muon_lab / P_muon_lab; 
    double theta_lab = acos(cos_theta_lab) * 180. / TMath::Pi();
    //h_theta_lab -> Fill(theta_lab);

    // == For double checking
    double E_pion = sqrt(this_pion_P * this_pion_P + mass_pion * mass_pion);
    double gm = gamma * mass_pion;
    
    if(decay){
      h_theta_lab -> Fill(theta_lab);
      expected_DiF_muon_P -> Fill(P_muon_lab);
      if(theta_lab < 2.5) expected_DiF_muon_P_2p5 -> Fill(P_muon_lab);
      if(theta_lab < 2.4) expected_DiF_muon_P_2p4 -> Fill(P_muon_lab);
      if(theta_lab < 2.2) expected_DiF_muon_P_2p2 -> Fill(P_muon_lab);
      if(theta_lab < 1.8) expected_DiF_muon_P_1p8 -> Fill(P_muon_lab);
      if(theta_lab < 1.2) expected_DiF_muon_P_1p2 -> Fill(P_muon_lab);
      if(theta_lab < 0.5) expected_DiF_muon_P_0p5 -> Fill(P_muon_lab);
      if(P_muon_lab > 750.){
	if(theta_lab < 1.8) expected_DiF_muon_P_1p8_right -> Fill(P_muon_lab);
	if(theta_lab < 1.2) expected_DiF_muon_P_1p2_right -> Fill(P_muon_lab);
	if(theta_lab < 0.5) expected_DiF_muon_P_0p5_right -> Fill(P_muon_lab);
      }

    }
    else N_not_decay++;
    /*
    cout << "P_muon_CoM : " << P_muon_CoM << ", this_pion_P : " << this_pion_P << ", gbm : " << gamma_beta * mass_pion
	 << ", E_pion : " << E_pion << ", gm : " << gm << ", E_muon_lab : " << E_muon_lab << ", Pz_muon_lab : " <<Pz_muon_lab << ", P_muon_lab : " << P_muon_lab << endl;
    */
 
    h_muon_fit_right_tail -> Fill(gaus_muon_right_tail_ext -> GetRandom(0., 1200.));
    h_muon_fit_right_tail_basis -> Fill(gaus_muon_right_tail_ext -> GetRandom(0., 1200.));
    h_pion_fit_basis -> Fill(gaus_pion -> GetRandom(0., 2000.));
  }
  cout << "Scale for basis : " << gaus_pion -> Integral(800., 1200.) / (N_trial + 0.) << endl;
  expected_DiF_muon_P_2p5 -> Scale(0.1 * gaus_pion -> Integral(800., 1200.) / (N_trial + 0.));
  expected_DiF_muon_P_2p4 -> Scale(0.1 * gaus_pion -> Integral(800., 1200.) / (N_trial + 0.));
  expected_DiF_muon_P_2p2 -> Scale(0.1 * gaus_pion -> Integral(800., 1200.) / (N_trial + 0.));
  expected_DiF_muon_P_1p8 -> Scale(0.1 * gaus_pion -> Integral(800., 1200.) / (N_trial + 0.));
  expected_DiF_muon_P_0p5 -> Scale(0.1 * gaus_pion -> Integral(800., 1200.) / (N_trial + 0.));

  h_muon_fit_right_tail_basis -> Scale( 600.  / h_muon_fit_right_tail_basis -> GetMaximum() );
  h_pion_fit_basis -> Scale(0.1 * gaus_pion -> GetParameter(0) / h_pion_fit_basis -> GetMaximum());
 
  double max_MC_muon = expected_DiF_muon_P -> GetMaximum();
  //expected_DiF_muon_P -> Scale(160. / max_MC_muon);
  double bin_size = 100000. / (10000. / rebin);
  double scale = gaus_pion -> Integral(800., 1200.) / (N_not_decay + 0.);
  cout << "scale : " << scale << ", bin_size: " << bin_size << endl;
  cout << "expected_DiF_muon_P -> Integral() : " << expected_DiF_muon_P -> Integral() << endl;
  expected_DiF_muon_P -> Scale( scale / bin_size);
  cout << "gaus_pion -> Integral(800., 1200.) : " << gaus_pion -> Integral(800., 1200.) << ", expected_DiF_muon_P -> Integral() : " << expected_DiF_muon_P -> Integral() << endl;
  expected_DiF_muon_P -> SetLineColor(kRed);
  expected_DiF_muon_P -> SetLineStyle(7);
  expected_DiF_muon_P -> Draw("histsame");

  TH1D* h_central_muon_exp = (TH1D*)muon_hist -> Clone();
  h_central_muon_exp -> Add(expected_DiF_muon_P, -1);
  h_central_muon_exp -> SetLineColor(kYellow);
  h_central_muon_exp -> Draw("histsame");

  h_central_muon_exp -> Fit(gaus_subtrac, "R", "", 800., 1200.);
  gaus_subtrac -> SetLineColor(kCyan);
  gaus_subtrac -> SetLineStyle(7);
  gaus_subtrac -> Draw("lsame");

  gPad->RedrawAxis();
  mapcanvas[canvas] -> cd();

  TString fit_pion_str = Form("(#mu, #sigma) = (%.2f, %.2f)", gaus_pion -> GetParameter(1), gaus_pion -> GetParameter(2));
  TString fit_subtrac_str = Form("(#mu, #sigma) = (%.2f, %.2f)", gaus_subtrac -> GetParameter(1), gaus_subtrac -> GetParameter(2));

  TLegend *l = new TLegend(0.2 , 0.6, 0.5 , 0.9); 
  l -> AddEntry(gaus_pion, fit_pion_str, "l");
  l -> AddEntry(gaus_subtrac, fit_subtrac_str, "l");
  l -> Draw("same");

  TLatex latex_ArgoNeuT, latex_data_POT, latex_histname, latex_ratio;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_histname.SetNDC();
  latex_ratio.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_histname.SetTextSize(0.050);
  latex_ratio.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.61, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  latex_histname.DrawLatex(0.55, 0.75, stage);
  latex_ratio.DrawLatex(0.55, 0.70, "#pi^{+} : #mu^{+} = " + ratio_pion_str + " : " + ratio_muon_str);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/Fitting_virtual_pion_vs_muon_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);
  mapcanvas[canvas]  -> Close();
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  c->cd();
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  TH1D *template_h2 = new TH1D("", "", 3000, -3., 3.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h2 -> SetStats(0);
  template_h2 -> GetXaxis() -> SetTitle("cos#theta_{C.o.M}");
  template_h2 -> GetXaxis() -> SetTitleSize(0.05);
  template_h2 -> GetXaxis() -> SetLabelSize(0.035);
  template_h2 -> GetYaxis() -> SetTitle("Events");
  template_h2 -> GetYaxis() -> SetLabelSize(0.035);
  template_h2 -> Draw(); 

  //h_cos_theta_CoM -> Scale(1./ h_cos_theta_CoM -> Integral());
  mc_max = h_cos_theta_CoM -> GetMaximum();
  template_h2 -> GetXaxis() -> SetRangeUser(-1.2, 1.2);
  template_h2 -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h2 -> Draw();
  h_cos_theta_CoM -> Draw("histsame");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/MC_cos_theta_CoM_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  //h_theta_lab -> Scale(1. / h_theta_lab -> Integral () );
  mc_max = h_theta_lab -> GetMaximum();
  template_h2 -> GetXaxis() -> SetTitle("#theta_{Lab.} (degree)");
  template_h2 -> GetXaxis() -> SetRangeUser(0., 4.0);
  template_h2 -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h2 -> Draw();
  h_theta_lab -> Draw("histsame");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/MC_theta_lab_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);
  
  mc_max = expected_DiF_muon_P -> GetMaximum();
  TH1D *template_h3 = new TH1D("", "", 1500, 0., 1500.);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h3 -> SetStats(0);
  template_h3-> GetXaxis() -> SetTitle("P_{Lab, #mu} (MeV)");
  template_h3 -> GetXaxis() -> SetRangeUser(0., 1500.);
  template_h3 -> GetYaxis() -> SetRangeUser(0., mc_max * 1.5);
  template_h3 -> Draw();
  expected_DiF_muon_P -> SetLineColor(kBlack);
  expected_DiF_muon_P -> SetLineStyle(1);
  expected_DiF_muon_P -> Draw("histsame");
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/MC_P_muon_lab_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);

  template_h3 -> Draw();
  muon_hist -> Draw("histsame");
 
  h_pion_fit_basis -> SetLineColor(kRed);
  h_pion_fit_basis -> SetLineStyle(7);
  h_pion_fit_basis -> Scale(0.05 * muon_hist -> Integral() / h_pion_fit_basis -> Integral());
  h_pion_fit_basis -> Draw("histsame");
  expected_DiF_muon_P_2p5 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_2p5 -> Integral());
  expected_DiF_muon_P_2p4 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_2p4 -> Integral());
  expected_DiF_muon_P_2p2 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_2p2 -> Integral());
  expected_DiF_muon_P_1p8 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_1p8 -> Integral());
  expected_DiF_muon_P_1p2 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_1p2 -> Integral());
  expected_DiF_muon_P_0p5 -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_0p5 -> Integral());
  expected_DiF_muon_P_1p8_right -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_1p8_right -> Integral());
  expected_DiF_muon_P_1p2_right -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_1p2_right -> Integral());
  expected_DiF_muon_P_0p5_right -> Scale(0.05 * pion_hist -> Integral() / expected_DiF_muon_P_0p5_right -> Integral());
  expected_DiF_muon_P_2p5 -> SetLineColor(kRed);
  expected_DiF_muon_P_2p4 -> SetLineColor(kRed);
  expected_DiF_muon_P_2p2 -> SetLineColor(kRed);
  expected_DiF_muon_P_1p8 -> SetLineColor(kRed);
  expected_DiF_muon_P_0p5 -> SetLineColor(kRed);
  expected_DiF_muon_P_2p5 -> SetLineStyle(7);
  expected_DiF_muon_P_2p4 -> SetLineStyle(7);
  expected_DiF_muon_P_2p2 -> SetLineStyle(7);
  expected_DiF_muon_P_1p8 -> SetLineStyle(7);
  expected_DiF_muon_P_0p5 -> SetLineStyle(7);
  expected_DiF_muon_P_2p5 -> Draw("histsame");
  expected_DiF_muon_P_2p4 -> Draw("histsame");
  expected_DiF_muon_P_2p2 -> Draw("histsame");
  expected_DiF_muon_P_1p8 -> Draw("histsame");
  expected_DiF_muon_P_0p5 -> Draw("histsame");

  sum_basis -> Add(h_pion_fit_basis);
  sum_basis -> Add(expected_DiF_muon_P_2p5);
  sum_basis -> Add(expected_DiF_muon_P_2p4);
  sum_basis -> Add(expected_DiF_muon_P_2p2);
  sum_basis -> Add(expected_DiF_muon_P_1p8);
  sum_basis -> Add(expected_DiF_muon_P_0p5);
  sum_basis -> Add(expected_DiF_muon_P_1p8_right);
  sum_basis -> Add(expected_DiF_muon_P_1p2_right);
  sum_basis -> Add(expected_DiF_muon_P_0p5_right);

  sum_basis -> SetLineColor(kOrange);
  sum_basis -> Draw("histsame");

  pdfname = WORKING_DIR + "/output/plots/BeamStudy/Fit_P_muon_lab_" + histname + "_" + stage + "_" + beam_P + "GeV.pdf";
  c -> SaveAs(pdfname);
  c -> Close();

  TFile *f_out = new TFile(WORKING_DIR + "/output/root/Template_" + stage + ".root", "RECREATE");
  f_out -> cd();
  muon_hist -> SetName("data");
  expected_DiF_muon_P_2p5 -> SetName("Theta_2p5_basis");
  expected_DiF_muon_P_2p4 -> SetName("Theta_2p4_basis");
  expected_DiF_muon_P_2p2 -> SetName("Theta_2p2_basis");
  expected_DiF_muon_P_1p8 -> SetName("Theta_1p8_basis");
  expected_DiF_muon_P_1p2 -> SetName("Theta_1p2_basis");
  expected_DiF_muon_P_0p5 -> SetName("Theta_0p5_basis");
  expected_DiF_muon_P_1p8_right -> SetName("Theta_1p8_basis_right");
  expected_DiF_muon_P_1p2_right -> SetName("Theta_1p2_basis_right");
  expected_DiF_muon_P_0p5_right -> SetName("Theta_0p5_basis_right");
  h_pion_fit_basis -> SetName("Pion_fit");
  muon_hist -> Write();
  expected_DiF_muon_P_2p5 -> Write();
  expected_DiF_muon_P_2p4 -> Write();
  expected_DiF_muon_P_2p2 -> Write();
  expected_DiF_muon_P_1p8 -> Write();
  expected_DiF_muon_P_1p2 -> Write();
  expected_DiF_muon_P_0p5 -> Write();
  expected_DiF_muon_P_1p8_right -> Write();
  expected_DiF_muon_P_1p2_right -> Write();
  expected_DiF_muon_P_0p5_right -> Write();
  h_pion_fit_basis -> Write();

  f_out -> Close();
}

void Run_Draw_VirtualDetector(TString filename, TString beam_P){
  const int N_stages = 7;
  TString stages[N_stages] = {"AfterTarget", "BPROF1", "BPROF3", "BPROF4", "COLL1", "TRIG1", "TRIG2"};
  for(int i = 0; i < N_stages; i++){
    if(i == 6) Study_true_beam(filename, "P", stages[i], "P (MeV/c)", beam_P, 0., 1800., 2.);
    //Draw_VirtualDetector(filename, "P", stages[i], "P (MeV/c)", beam_P, 0., 1800., 1.);
  }

}

void Draw_BeamSample(){

  setTDRStyle();
  TString file_suffix = "_0.5M.root";

  /*
  Draw_MC_shape_comparison("_BeamSampleAna_0.5GeV" + file_suffix, "P", "P (MeV/c)", "0.5", 0., 1000., 10., 0.4);
  Draw_MC_shape_comparison("_BeamSampleAna_1.0GeV" + file_suffix, "P", "P (MeV/c)", "1.0", 0., 1500., 10., 0.4);
  Draw_MC_shape_comparison("_BeamSampleAna_2.0GeV" + file_suffix, "P", "P (MeV/c)", "2.0", 0., 2500., 10., 0.4);
  Draw_pion_muon_comparison("_BeamSampleAna_1.0GeV" + file_suffix, "P_AfterTarget", "P (MeV/c)", "1.0", 600., 2000., 20., 0.4);
  */

  Run_Draw_VirtualDetector("_BeamSampleAna_1.0GeV" + file_suffix, "1.0");
}

