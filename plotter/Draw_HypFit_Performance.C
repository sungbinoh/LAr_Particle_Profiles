#include "canvas_margin.h"
#include "mylib.h"

void Draw_1D_Truth_Performance(TString filename, TString dir, TString particle, TString method, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){

  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "range") method_for_latex = "CSDA(range)";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Truth_Performance] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_mc = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "KE_Res_1D_" +particle + "_" + method+ "_" + KE_range + "_";
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;


  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Truth_Performance] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      cout << "N_part : " << N_part << ", this_NHits_low_string : " <<this_NHits_low_string << ", this_NHits_high_string : " << this_NHits_high_string << endl;
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  
  // == Draw all
  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/KE_Res_1D_" + particle + "_" + method+ "_" + KE_range + ".pdf";
  c -> SaveAs(pdfname);

  f_mc -> Close();
  delete l;
  c -> Close();
}

void Draw_1D_Fit_Efficiency(TString filename, TString dir, TString particle, TString method, TString TitleX, double xmin, double xmax, double rebin){
    
    TString particle_for_latex = "";
    TString method_for_latex = "";
    if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
    if(particle.Contains("proton")) particle_for_latex = "p^{+}";
    if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
    if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

    cout << "[Draw_1D_Truth_Performance] Start" << endl;
    TString input_file_dir = getenv("LArProf_WD");
    TString root_file_path =input_file_dir + "/output/root/";
    TFile *f_mc = new TFile(root_file_path + "mc" + filename);
    gDirectory -> Cd(dir);

    TString range_ratio_str[6] = {"LessThan0p1", "0p1to0p3", "0p3to0p5", "0p5to0p7", "0p7to0p9", "BiggerThan0p9"};
    TString range_ratio_str_latex[6] = {"< 0.1", "0.1 - 0.3", "0.3 - 0.5", "0.5 - 0.7", "0.7 - 0.9", "> 0.9"};
    int color_array[6] = {632, 800, 400, 416, 600, 880};

    vector<TString> valid_range_ratio_vtr;
    valid_range_ratio_vtr.clear();
    vector<TString> valid_range_ratio_latex_vtr;
    valid_range_ratio_latex_vtr.clear();
    vector<int> valid_range_ratio_color_vtr;
    valid_range_ratio_color_vtr.clear();

    TString range_hist_name_prefix = "KE_eff_" + particle + "_range_";
    TString method_hist_name_prefix = "KE_eff_" + particle + "_" + method + "_";
    for(int i = 0; i < 6; i++){
      TString this_range_hist_name = range_hist_name_prefix + range_ratio_str[i] + "_" + dir;
      TString this_method_hist_name = method_hist_name_prefix + range_ratio_str[i] + "_" + dir;
      if((TH1D*)gDirectory -> Get(this_range_hist_name)){
	maphist[this_range_hist_name] = (TH1D*)gDirectory -> Get(this_range_hist_name) -> Clone();
	maphist[this_range_hist_name] -> Rebin(rebin);
      }
      else maphist[this_range_hist_name] = nullptr;

      if((TH1D*)gDirectory -> Get(this_method_hist_name)){
        maphist[this_method_hist_name] = (TH1D*)gDirectory -> Get(this_method_hist_name) -> Clone();
        maphist[this_method_hist_name] -> Rebin(rebin);
        valid_range_ratio_vtr.push_back(range_ratio_str[i]);
	valid_range_ratio_latex_vtr.push_back(range_ratio_str_latex[i]);
	valid_range_ratio_color_vtr.push_back(color_array[i]);
      }
      else maphist[this_method_hist_name] = nullptr;
    }

    if(valid_range_ratio_vtr.size() == 0) return;

    // == Make overall efficiency
    TH1D *h_range_all = (TH1D*)maphist[range_hist_name_prefix + valid_range_ratio_vtr.at(0) + "_" + dir] -> Clone();
    for(unsigned int i = 1; i < valid_range_ratio_vtr.size(); i++){
      h_range_all -> Add(maphist[range_hist_name_prefix + valid_range_ratio_vtr.at(i) + "_" + dir]);
    }
    TH1D *h_method_all = (TH1D*)maphist[method_hist_name_prefix + valid_range_ratio_vtr.at(0) + "_" + dir] -> Clone();
    for(unsigned int i = 1; i < valid_range_ratio_vtr.size(); i++){
      h_method_all -> Add(maphist[method_hist_name_prefix+ valid_range_ratio_vtr.at(i) + "_" + dir]);
    }
    h_method_all -> Divide(h_range_all);

    // == Make efficiency histograms
    for(unsigned int i = 0; i < valid_range_ratio_vtr.size(); i++){
      TString this_range_hist_name = range_hist_name_prefix + valid_range_ratio_vtr.at(i) + "_" + dir;
      TString this_method_hist_name = method_hist_name_prefix + valid_range_ratio_vtr.at(i) + "_" + dir;
      maphist[this_method_hist_name] -> Divide(maphist[this_range_hist_name]);
    }
    

    TCanvas *c = new TCanvas("","", 600, 800);
    canvas_margin(c);
    gStyle -> SetOptStat(1111);

    // == Draw all
    TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
    gStyle->SetOptTitle(0);
    gStyle->SetLineWidth(2);
    template_h -> SetStats(0);
    template_h -> GetXaxis() -> SetTitle(TitleX);
    template_h -> GetXaxis() -> SetTitleOffset(1.3);
    template_h -> GetXaxis() -> SetTitleSize(0.035);
    template_h -> GetXaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetTitle("Efficiency");
    template_h -> GetYaxis() -> SetLabelSize(0.035);
    template_h -> GetYaxis() -> SetRangeUser(0., 1.5);
    template_h -> Draw();

    TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
    l -> SetFillColor(kWhite);
    l -> SetLineColor(kWhite);
    l -> SetBorderSize(1);
    l -> SetFillStyle(1001);
    l -> SetShadowColor(0);
    l -> SetEntrySeparation(0.3);
    l -> SetNColumns(2);

    h_method_all -> SetLineColor(kBlack);
    h_method_all -> SetLineWidth(3);
    h_method_all -> Draw("histsame");
    l -> AddEntry(h_method_all, "Overall", "l");
    
    for(unsigned int i = 0; i < valid_range_ratio_vtr.size(); i++){
      TString this_method_hist_name = method_hist_name_prefix + valid_range_ratio_vtr.at(i) + "_" + dir;
      maphist[this_method_hist_name] -> SetLineColor(valid_range_ratio_color_vtr.at(i));
      maphist[this_method_hist_name] -> SetLineWidth(3);
      maphist[this_method_hist_name] -> SetLineStyle(7);
      maphist[this_method_hist_name] -> Draw("histsame");
      l -> AddEntry(maphist[this_method_hist_name], "L_{Reco.} / L_{CSDA} : " + valid_range_ratio_latex_vtr.at(i), "l");
    }

    h_method_all -> Draw("histsame");
    l -> Draw("same");

    gPad->RedrawAxis();
    c -> cd();
    TLatex latex_ArgoNeuT, latex_data_POT;
    latex_ArgoNeuT.SetNDC();
    latex_data_POT.SetNDC();
    latex_ArgoNeuT.SetTextSize(0.035);
    latex_data_POT.SetTextSize(0.035);
    latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
    latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
    TString pdfname;
    TString WORKING_DIR = getenv("LArProf_WD");
    pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Efficiency_" + particle + "_" + method+ ".pdf";
    c -> SaveAs(pdfname);

    f_mc -> Close();
    delete l;
    c -> Close();

}



void Draw_1D_Reco_Performance(TString filename, TString dir, TString particle, TString method, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){

  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Reco_Performance] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_data = new TFile(root_file_path + "data" + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "Reco_KE_Res_1D_" +particle + "_" + method+ "_" + KE_range + "_";
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;

  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Reco_Performance] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Reco_KE_Res_1D_" + particle + "_" + method+ "_" + KE_range + ".pdf";
  c -> SaveAs(pdfname);

  f_data -> Close();
  delete l;
  c -> Close();



}

void Draw_1D_Beam_Performance(TString data_or_mc, TString filename, TString dir, TString particle, TString method, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){

  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "range") method_for_latex = "CSDA(range)";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Beam_Performance] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_data = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "Beam_Reco_KE_Res_1D_" +particle + "_" + method+ "_" + KE_range + "_";
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;

  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Beam_Performance] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Beam_Reco_" + data_or_mc + "_KE_Res_1D_" + particle + "_" + method+ "_" + KE_range + ".pdf";
  c -> SaveAs(pdfname);

  f_data -> Close();
  delete l;
  c -> Close();

}

void Draw_1D_Beam_Performance_MC_PID(TString data_or_mc, TString filename, TString dir, TString particle, TString method, TString PID, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){
  
  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "range") method_for_latex = "CSDA(range)";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Beam_Performance_MC_PID] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_data = new TFile(root_file_path + data_or_mc + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "Beam_Reco_KE_Res_1D_" + particle + "_" + PID + "_" + method + "_" + KE_range + "_";
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    //cout << "[Draw_1D_Beam_Performance_MC_PID] this_hist_name : " << this_hist_name << endl;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;

  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Beam_Performance_MC_PID] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Beam_Reco_" + data_or_mc + "_KE_Res_1D_" + particle + "_" + PID + "_" + method + "_" + KE_range + ".pdf";
  c -> SaveAs(pdfname);

  f_data -> Close();
  delete l;
  c -> Close();

}

void Draw_1D_Beam_True_Performance(TString filename, TString dir, TString particle, TString method, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){

  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "range") method_for_latex = "CSDA(range)";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Beam_Performance_MC_PID] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_data = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "Beam_true_KE_Res_1D_" + particle + "_" + method + "_" + KE_range + "_";
  if(filename.Contains("Micheless")) hist_name_prefix = "Micheless_" + hist_name_prefix;
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;

  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Beam_Performance_MC_PID] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Beam_true_KE_Res_1D_" + particle + "_" + method + "_" + KE_range + ".pdf";
  if(filename.Contains("Micheless")) pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Micheless_Beam_true_KE_Res_1D_" + particle + "_" + method + "_" + KE_range + ".pdf";

  c -> SaveAs(pdfname);

  f_data -> Close();
  delete l;
  c -> Close();
}

void Draw_1D_Beam_True_Performance_PID(TString filename, TString dir, TString particle, TString method, TString PID, TString KE_range, TString TitleX, double xmin, double xmax, double rebin, TString KE_ranges_legend){
  
  TString particle_for_latex = "";
  TString method_for_latex = "";
  if(particle.Contains("pion")) particle_for_latex = "#pi^{+}";
  if(particle.Contains("proton")) particle_for_latex = "p^{+}";
  if(method == "range") method_for_latex = "CSDA(range)";
  if(method == "gaussian") method_for_latex = "Fitted (Gaussian)";
  if(method == "likelihood") method_for_latex = "Fitted (Likelihood)";

  cout << "[Draw_1D_Beam_Performance_MC_PID] Start" << endl;
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_data = new TFile(root_file_path + "mc" + filename);
  gDirectory -> Cd(dir);

  TString NHits_str[31] = {""};
  NHits_str[0] = "NHits15to30";
  for(int i = 1; i < 31; i++){
    TString this_NHits_str = Form("NHits%dto%d", i * 30, (i + 1) * 30);
    NHits_str[i] = this_NHits_str;
  }

  vector<TString> valid_NHits_str;
  valid_NHits_str.clear();

  TString hist_name_prefix = "Beam_true_KE_Res_1D_" + particle + "_" + PID + "_" + method + "_" + KE_range + "_";
  if(filename.Contains("Micheless")) hist_name_prefix = "Micheless_" + hist_name_prefix;
  for(int i = 0; i < 31; i++){
    TString this_hist_name = hist_name_prefix + NHits_str[i] + "_" + dir;
    if((TH1D*)gDirectory -> Get(this_hist_name)){
      maphist[this_hist_name] = (TH1D*)gDirectory -> Get(this_hist_name) -> Clone();
      maphist[this_hist_name] -> Rebin(rebin);
      valid_NHits_str.push_back(NHits_str[i]);
    }
    else maphist[this_hist_name] = nullptr;
  }

  if(valid_NHits_str.size() == 0) return;

  TH1D *h_all = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(0) + "_" + dir] -> Clone();
  for(unsigned int i = 1; i < valid_NHits_str.size(); i++){
    h_all -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
  }
  double y_max = h_all -> GetMaximum();

  cout << "[Draw_1D_Beam_Performance_MC_PID] Merging NHits for valid_NHits_str.size() : " << valid_NHits_str.size() << endl;

  int N_per_group = 1 + valid_NHits_str.size() / 4;
  int N_groups = 0;
  TString NHits_low_strings[4] = {""};
  TString NHits_high_strings[4] = {""};
  for(unsigned int i = 0; i < valid_NHits_str.size(); i++){
    int this_group = i / N_per_group;
    TString this_group_str = Form("%d", this_group);
    if(i % N_per_group == 0){
      maphist[particle + method + KE_range + this_group_str] = (TH1D*)maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir] -> Clone();
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      int N_part = tx -> GetEntries();
      TString this_NHits_low_string =  ((TObjString *)(tx->At(1)))->String();
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      this_NHits_low_string.Remove(0, 1);
      NHits_low_strings[this_group] = this_NHits_low_string;
      NHits_high_strings[this_group] = this_NHits_high_string;
      N_groups++;
    }
    else{
      maphist[particle + method + KE_range + this_group_str] -> Add(maphist[hist_name_prefix + valid_NHits_str.at(i) + "_" + dir]);
      TObjArray *tx = valid_NHits_str.at(i).Tokenize("to");
      TString this_NHits_high_string =  ((TObjString *)(tx->At(2)))->String();
      NHits_high_strings[this_group] = this_NHits_high_string;
    }
  }

  TCanvas *c = new TCanvas("","", 600, 800);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);

  TH1D *template_h = new TH1D("", "", 1, xmin, xmax);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> SetStats(0);
  template_h -> GetXaxis() -> SetTitle(TitleX);
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetXaxis() -> SetTitleSize(0.035);
  template_h -> GetXaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetLabelSize(0.035);
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  TLegend *l = new TLegend(0.20, 0.70, 0.90, 0.90);
  l -> SetFillColor(kWhite);
  l -> SetLineColor(kWhite);
  l -> SetBorderSize(1);
  l -> SetFillStyle(1001);
  l -> SetShadowColor(0);
  l -> SetEntrySeparation(0.3);
  l -> SetNColumns(3);

  h_all -> SetLineColor(kBlack);
  h_all -> SetLineWidth(3);
  h_all -> Draw("histsame");

  Int_t colour_array[] = {632, 800, 416, 600};
  for(unsigned int i = 0; i < N_groups; i++){
    TString this_NHits_4_str = Form("%d", i);
    TString this_hist_name = particle + method + KE_range + this_NHits_4_str;
    TString this_legend_str = "N_{Hits} " + NHits_low_strings[i] + " - " + NHits_high_strings[i];
    maphist[this_hist_name] -> SetLineColor(colour_array[i]);
    maphist[this_hist_name] -> SetLineStyle(7);
    maphist[this_hist_name] -> SetLineWidth(2);
    maphist[this_hist_name] -> Draw("histsame");
    l -> AddEntry(maphist[this_hist_name], this_legend_str, "l");
  }

  l -> AddEntry(h_all, "All", "l");
  l -> Draw("same");

  gPad->RedrawAxis();
  c -> cd();
  TLatex latex_ArgoNeuT, latex_data_POT, latex_KE_range;
  latex_ArgoNeuT.SetNDC();
  latex_data_POT.SetNDC();
  latex_KE_range.SetNDC();
  latex_ArgoNeuT.SetTextSize(0.035);
  latex_data_POT.SetTextSize(0.035);
  latex_KE_range.SetTextSize(0.035);
  latex_ArgoNeuT.DrawLatex(0.57, 0.96, "#font[62]{ProtoDUNE-SP} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_data_POT.DrawLatex(0.158, 0.96, particle_for_latex + ", " + method_for_latex);
  latex_KE_range.DrawLatex(0.20, 0.90, KE_ranges_legend);
  TString pdfname;
  TString WORKING_DIR = getenv("LArProf_WD");
  pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Beam_true_KE_Res_1D_" + particle + "_" + PID + "_" + method + "_" + KE_range + ".pdf";
  if(filename.Contains("Micheless")) pdfname = WORKING_DIR + "/output/plots/HypFit_Performance/Micheless_Beam_true_KE_Res_1D_" + particle + "_" + PID + "_" + method + "_" + KE_range + ".pdf";

  c -> SaveAs(pdfname);

  f_data -> Close();
  delete l;
  c -> Close();
}


void Run_Draw_1D_Truth_Performance(TString filename, TString dir){
  TString KE_ranges[11] = {"KE0to50", "KE50to100", "KE100to200", "KE200to300", "KE300to400", "KE400to500", "KE500to600", "KE600to700", "KE700to800", "KE800to900", "KE900to1000"};
  TString KE_ranges_legend[11] = {"KE 0 - 50 (MeV)", "KE 50 - 100 (MeV)", "KE 100 - 200 (MeV)", "KE 200 - 300 (MeV)", "KE 300 - 400 (MeV)", "KE 400 - 500 (MeV)", "KE 500 - 600 (MeV)",
                                  "KE 600 - 700 (MeV)", "KE 700 - 800 (MeV)", "KE 800 - 900 (MeV)", "KE 900 - 1000 (MeV)"};

  TString interactions[3] = {"", "_stopped", "_interacted"};
  TString methods[3] = {"range", "gaussian", "likelihood"};

  for(int i = 0; i < 11; i ++){
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 3; k++){
	if(filename.Contains("pion")) Draw_1D_Truth_Performance(filename, dir, "pion" + interactions[j], methods[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{True}}{KE_{True}}", -2., 2., 5., KE_ranges_legend[i]);
        if(filename.Contains("proton")) Draw_1D_Truth_Performance(filename, dir, "proton" + interactions[j], methods[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{True}}{KE_{True}}", -2., 2., 5., KE_ranges_legend[i]);
      }
    }
  }

  for(int j = 1; j < 3; j++){
    if(filename.Contains("pion")) Draw_1D_Fit_Efficiency(filename, dir, "pion", methods[j], "KE [MeV]", 0., 180., 10.);
    if(filename.Contains("proton"))Draw_1D_Fit_Efficiency(filename, dir, "proton", methods[j], "KE [MeV]", 0., 1800., 10.);
  }
}

void Run_Draw_1D_Reco_Performance(TString filename, TString dir){
  TString KE_ranges[11] = {"KE0to50", "KE50to100", "KE100to200", "KE200to300", "KE300to400", "KE400to500", "KE500to600", "KE600to700", "KE700to800", "KE800to900", "KE900to1000"};
  TString KE_ranges_legend[11] = {"KE 0 - 50 (MeV)", "KE 50 - 100 (MeV)", "KE 100 - 200 (MeV)", "KE 200 - 300 (MeV)", "KE 300 - 400 (MeV)", "KE 400 - 500 (MeV)", "KE 500 - 600 (MeV)",
                                  "KE 600 - 700 (MeV)", "KE 700 - 800 (MeV)", "KE 800 - 900 (MeV)", "KE 900 - 1000 (MeV)"};

  TString methods[2] = {"gaussian", "likelihood"};
  
  for(int i = 0; i < 11; i ++){
    for(int j = 0; j < 2; j++){
      Draw_1D_Reco_Performance(filename, dir, "pion", methods[j], KE_ranges[i], "#frac{KE_{Fit} - KE_{Range}}{KE_{Range}}", -2., 2., 10., KE_ranges_legend[i]);
    }
  }
}

void Run_Draw_1D_Beam_Performance(TString data_or_mc, TString filename, TString dir){
  const int N_KE = 17;
  TString KE_ranges[N_KE] = {"KE0to50", "KE50to100", "KE100to200", "KE200to300", "KE300to400", "KE400to500", "KE500to600", "KE600to700", "KE700to800", "KE800to900", "KE900to1000",
			   "KE1000to1100", "KE1100to1200", "KE1200to1300", "KE1300to1400", "KE1400to1500", "KE1500to1600"};
  TString KE_ranges_legend[N_KE] = {"KE 0 - 50 (MeV)", "KE 50 - 100 (MeV)", "KE 100 - 200 (MeV)", "KE 200 - 300 (MeV)", "KE 300 - 400 (MeV)", "KE 400 - 500 (MeV)", "KE 500 - 600 (MeV)",
				    "KE 600 - 700 (MeV)", "KE 700 - 800 (MeV)", "KE 800 - 900 (MeV)", "KE 900 - 1000 (MeV)",
				    "KE 1000 - 1100 (MeV)", "KE 1100 - 1200 (MeV)", "KE 1200 - 1300 (MeV)", "KE 1300 - 1400 (MeV)", "KE 1400 - 1500 (MeV)", "KE 1500 - 1600 (MeV)"};

  TString methods[3] = {"range", "gaussian", "likelihood"};
  TString PID[2] = {"PID13", "PID211"};
  for(int i = 0; i < N_KE; i ++){
    for(int j = 0; j < 3; j++){
      Draw_1D_Beam_Performance(data_or_mc, filename, dir, "pion", methods[j], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
      Draw_1D_Beam_Performance(data_or_mc, filename, dir, "proton", methods[j], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
      if(data_or_mc == "mc"){
	for(int k = 0; k < 2; k++){
	  Draw_1D_Beam_Performance_MC_PID(data_or_mc, filename, dir, "proton", methods[j], PID[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
	  Draw_1D_Beam_Performance_MC_PID(data_or_mc, filename, dir, "pion", methods[j], PID[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
	}
      }
    }
  }
}

void Run_Draw_1D_Beam_True_Performance(TString filename, TString dir){
  const int N_KE = 17;
  TString KE_ranges[N_KE] = {"KE0to50", "KE50to100", "KE100to200", "KE200to300", "KE300to400", "KE400to500", "KE500to600", "KE600to700", "KE700to800", "KE800to900", "KE900to1000",
			     "KE1000to1100", "KE1100to1200", "KE1200to1300", "KE1300to1400", "KE1400to1500", "KE1500to1600"};
  TString KE_ranges_legend[N_KE] = {"KE 0 - 50 (MeV)", "KE 50 - 100 (MeV)", "KE 100 - 200 (MeV)", "KE 200 - 300 (MeV)", "KE 300 - 400 (MeV)", "KE 400 - 500 (MeV)", "KE 500 - 600 (MeV)",
                                    "KE 600 - 700 (MeV)", "KE 700 - 800 (MeV)", "KE 800 - 900 (MeV)", "KE 900 - 1000 (MeV)",
                                    "KE 1000 - 1100 (MeV)", "KE 1100 - 1200 (MeV)", "KE 1200 - 1300 (MeV)", "KE 1300 - 1400 (MeV)", "KE 1400 - 1500 (MeV)", "KE 1500 - 1600 (MeV)"};

  TString methods[3] = {"range", "gaussian", "likelihood"};
  TString PID[2] = {"PID13", "PID211"};
  for(int i = 0; i < N_KE; i ++){
    for(int j = 0; j < 3; j++){
      Draw_1D_Beam_True_Performance(filename, dir, "pion", methods[j], KE_ranges[i], "#frac{KE_{Meas.} - KE_{True, ff}}{KE_{True, ff}}", -2., 2., 5., KE_ranges_legend[i]);
      Draw_1D_Beam_True_Performance(filename, dir, "proton", methods[j], KE_ranges[i], "#frac{KE_{Meas.} - KE_{True, ff}}{KE_{True, ff}}", -2., 2., 5., KE_ranges_legend[i]);
      for(int k = 0; k < 2; k++){
	Draw_1D_Beam_True_Performance_PID(filename, dir, "pion", methods[j], PID[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{True, ff}}{KE_{True, ff}}", -2., 2., 5., KE_ranges_legend[i]);
      }
      //Draw_1D_Beam_True_Performance(filename, dir, "proton", methods[j], PID[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
      //Draw_1D_Beam_True_Performance(filename, dir, "pion", methods[j], PID[k], KE_ranges[i], "#frac{KE_{Meas.} - KE_{Beam Inst.}}{KE_{Beam Inst.}}", -2., 2., 5., KE_ranges_legend[i]);
    }
  }
}

void Draw_HypFit_Performance(){

  setTDRStyle();
  TString file_suffix = ".root";
  //Run_Draw_1D_Truth_Performance("_HypFit_1.0GeV" + file_suffix, "beam_window_noweight");
  //Run_Draw_1D_Truth_Performance("_HypFit_merged" + file_suffix, "beam_window_noweight");
  //Run_Draw_1D_Truth_Performance("_HypFit_1.0GeV_daughter_pion" + file_suffix, "beam_window_noweight");
  //Run_Draw_1D_Truth_Performance("_HypFit_1.0GeV_daughter_proton" + file_suffix, "beam_window_noweight");
  Run_Draw_1D_Truth_Performance("_HypFit_daughter_proton" + file_suffix, "beam_window_noweight");


  //Run_Draw_1D_Reco_Performance("_HypFit_1.0GeV" + file_suffix, "beam_window_noweight");

  //Run_Draw_1D_Beam_Performance("data", "_HypFit_0.5GeV_beam_pion" + file_suffix, "beam_window_noweight");
  //Run_Draw_1D_Beam_Performance("mc", "_HypFit_0.5GeV_beam_pion" + file_suffix, "beam_window_noweight");

  //Run_Draw_1D_Beam_True_Performance("_HypFit_0.5GeV_beam_pion" + file_suffix, "beam_window_noweight");
  //Run_Draw_1D_Beam_True_Performance("_HypFit_0.5GeV_beam_pion_Micheless" + file_suffix, "beam_window_noweight");


  //Run_Draw_1D_Beam_Performance("_HypFit_merged_beam_proton" + file_suffix, "beam_window_noweight");
}

