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


void Draw_BeamSample(){

  setTDRStyle();
  TString file_suffix = ".root";

  Draw_MC_shape_comparison("_BeamSampleAna_0.5GeV" + file_suffix, "P", "P (MeV/c)", "0.5", 0., 1000., 10., 0.4);
  Draw_MC_shape_comparison("_BeamSampleAna_1.0GeV" + file_suffix, "P", "P (MeV/c)", "1.0", 0., 1500., 10., 0.4);
  Draw_MC_shape_comparison("_BeamSampleAna_2.0GeV" + file_suffix, "P", "P (MeV/c)", "2.0", 0., 2500., 10., 0.4);

}

