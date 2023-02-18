#include "canvas_margin.h"
#include "mylib.h"

void Draw_FakeData_IsQE(TString filename, TString dir, TString histname, TString title_X, TString beam_P, double xmin, double xmax, double rebin_x){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";

  // == Call MC Histograms
  TFile *f_mc = new TFile(root_file_path + filename);
  gDirectory -> Cd(dir);
  int N_type = N_pi_type;
  //histname = "hdaughter_reco_" + histname;
  for(int i = 0; i < N_type; i++){
    TString this_histname = Form(histname + "_" + dir + "_%d", i);
    TString this_histname_IsQE = Form(histname + "_" + dir + "_%d_IsQE", i);

    if((TH1D*)gDirectory -> Get(this_histname)){
      maphist[this_histname] = (TH1D*)gDirectory -> Get(this_histname) -> Clone();
      maphist[this_histname] -> Rebin(rebin_x);
    }
    else maphist[this_histname] = nullptr;

    if((TH1D*)gDirectory -> Get(this_histname_IsQE)){
      maphist[this_histname_IsQE] = (TH1D*)gDirectory -> Get(this_histname_IsQE) -> Clone();
      maphist[this_histname_IsQE] -> Rebin(rebin_x);
    }
    else maphist[this_histname_IsQE] = nullptr;
  }

  TH1D *h_fakedata = (TH1D*)maphist[histname + "_" + dir + "_0"] -> Clone();
  h_fakedata -> Add(maphist[histname +"_" + dir + "_0_IsQE"]);

  TH1D *h_IsQE = (TH1D*)maphist[histname +"_" + dir + "_1_IsQE"] -> Clone();
  for(int i = 2; i < N_type; i++){
    TString this_histname_IsQE = Form(histname + "_" + dir + "_%d_IsQE", i);
    if(maphist[this_histname_IsQE] != nullptr){
      h_IsQE -> Add(maphist[this_histname_IsQE]);
    }
  }
  TH1D *h_MC_Sum = (TH1D*)h_IsQE -> Clone();

  // == Draw
  TCanvas *c = new TCanvas("", "", 800, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TPad *pad_1 = new TPad("", "", 0, 0.25, 1, 1);
  pad_1 -> SetTopMargin( 0.07 );
  pad_1 -> SetBottomMargin( 0.05 );
  pad_1 -> SetLeftMargin( 0.15 );
  pad_1 -> SetRightMargin( 0.03 );
  pad_1 -> Draw();
  pad_1 -> cd();

  double y_max = h_fakedata -> GetMaximum();
  TH1D *pad1_template = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(3);
  pad1_template -> SetStats(0);
  pad1_template -> GetXaxis() -> SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitleSize(0);
  pad1_template -> GetYaxis() -> SetLabelSize(0.05);
  pad1_template -> GetYaxis() -> SetTitleSize(0.07);
  pad1_template -> GetYaxis() -> SetTitleOffset(1.02);
  pad1_template -> GetYaxis() -> SetTitle("Events");
  pad1_template -> GetYaxis() -> SetRangeUser(0., y_max * 1.8);
  pad1_template -> Draw("hist");

  TLegend *l_pad1 = new TLegend(0.20, 0.70, 0.90, 0.90);
  l_pad1 -> SetFillColor(kWhite);
  l_pad1 -> SetLineColor(kWhite);
  l_pad1 -> SetBorderSize(1);
  l_pad1 -> SetFillStyle(4000);
  l_pad1 -> SetShadowColor(0);
  l_pad1 -> SetEntrySeparation(0.3);
  l_pad1 -> SetNColumns(3);
  
  Int_t colour_array[] = {0, 632, 800, 867, 600, 416, 901, 432, 400, 920};
  THStack * h_mc_stack = new THStack("", "");
  h_mc_stack -> Add(h_IsQE);
  h_IsQE -> SetFillColor(kMagenta);  
  TString IsQE_N_event = Form("%.1f", h_IsQE -> Integral());
  
  for(int i = 1; i < N_type; i++){
    TString this_hist_name = Form(histname + "_" + dir + "_%d", i);
    if(maphist[this_hist_name] != nullptr){
      TString this_N_event = Form("%.1f", maphist[this_hist_name] -> Integral());
      TString this_legend_str = "";
      this_legend_str = pi_type_str[i] + " " + this_N_event;
      maphist[this_hist_name] -> SetLineColor(colour_array[i]);
      maphist[this_hist_name] -> SetFillColor(colour_array[i]);
      h_mc_stack -> Add(maphist[this_hist_name]);
      h_MC_Sum -> Add(maphist[this_hist_name]);
      l_pad1 -> AddEntry(maphist[this_hist_name], this_legend_str, "f");
    }
  }
  l_pad1 -> AddEntry(h_IsQE, "#Delta E_{QE}^{true} < 200 MeV " + IsQE_N_event, "f"); 

  h_fakedata -> SetLineColor(kBlack);
  h_fakedata -> SetMarkerColor(kBlack);
  h_fakedata -> SetMarkerStyle(20);
  h_fakedata -> SetMarkerSize(1);
  h_fakedata -> SetLineWidth(2);
  TString data_N_event = Form("%.1f", h_fakedata -> Integral());
  l_pad1 -> AddEntry(h_fakedata, "Fake Data " + data_N_event, "lp");
  h_mc_stack -> Draw("histsame");
  h_fakedata -> Draw("epsame");
  l_pad1 -> Draw("same");
  gPad->RedrawAxis();

  // == Bottom pad
  c -> cd();
  TPad *pad_2 = new TPad("", "", 0, 0, 1, 0.25);
  pad_2 -> SetTopMargin( 0.05 );
  pad_2 -> SetBottomMargin( 0.4 );
  pad_2 -> SetLeftMargin( 0.15 );
  pad_2 -> SetRightMargin( 0.03 );
  pad_2 -> Draw();
  pad_2 -> cd();
  TH1D * pad2_template = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  pad2_template -> Draw("hist");
  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(title_X);
  pad2_template -> GetXaxis() -> SetTitleSize(0.15);
  pad2_template -> GetXaxis() -> SetLabelSize(0.125);
  pad2_template -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.4);
  pad2_template -> GetYaxis() -> SetLabelSize(0.09);
  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.0, 3.0);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");

  TH1D * data_mc_ratio = (TH1D*)h_fakedata -> Clone();
  data_mc_ratio -> Divide(h_MC_Sum);
  data_mc_ratio -> Draw("psame");

  TLegend *l_pad2 = new TLegend(0.2, 0.85, 0.6, 0.95);
  l_pad2 -> SetBorderSize(0);
  l_pad2 -> SetNColumns(3);
  l_pad2 -> AddEntry(data_mc_ratio, "Fake./MC.", "lp");

  TLine *pad2_line = new TLine(xmin, 1, xmax, 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

  l_pad2 -> Draw("same");

  gPad->RedrawAxis();

  // == Latex
  c -> cd();
  TLatex latex_ProtoDUNE, latex_data;
  latex_ProtoDUNE.SetNDC();
  latex_data.SetNDC();
  latex_ProtoDUNE.SetTextSize(0.035);
  latex_data.SetTextSize(0.035);
  latex_ProtoDUNE.DrawLatex(0.15, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data.DrawLatex(0.63, 0.96, "Run 1, " + beam_P + " GeV/c Beam");
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/PionXsec/MC_vs_FakeData/MC_vs_FakeData_" + histname + "_" + beam_P + "GeV_" + dir + ".pdf";
  c -> SaveAs(pdfname);

  f_mc -> Close();
  c -> Close();
}

void Run_Draw_FakeData_IsQE(TString filename, TString dir, TString histname, TString title_X, TString beam_P, double xmin, double xmax, double rebin_x){
  Draw_FakeData_IsQE(filename, dir, histname, title_X, beam_P, xmin, xmax, rebin_x);
}

void Draw_FakeData(){

  cout << "============================" << endl;
  cout << "=========== Start ==========" << endl;
  cout << "============================" << endl;

  setTDRStyle();
  TString file_suffix = ".root";
  // mc_PionXsec_1.0GeV_IsQE.root
  // hdaughter_reco_likelihood_fit_pion_EQEmE_NC_40_noweight_0
  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "likelihood_fit_pion_EQEmE_NC_40", "#Delta E_{QE} [MeV]", "1.0", -3000., 3000., 40.);
  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "likelihood_fit_pion_EQEmE_NC_4", "#Delta E_{QE} [MeV]", "1.0", -3000., 3000., 40.);

  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "gaus_fit_pion_EQEmE_NC_40", "#Delta E_{QE} [MeV]", "1.0", -3000., 3000., 40.);
  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "gaus_fit_pion_EQEmE_NC_4", "#Delta E_{QE} [MeV]", "1.0", -3000., 3000., 40.);

  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "stop_pion_EQEmE_NC_40", "#Delta E_{QE} [MeV]", "1.0", -1000., 2000., 40.);
  Run_Draw_FakeData_IsQE("mc_PionXsec_1.0GeV_IsQE.root", "noweight", "stop_pion_EQEmE_NC_4", "#Delta E_{QE} [MeV]", "1.0", -1000., 2000., 40.);
}
