#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"

TF1 * KE_remove_beam_plug(double p[], int PID, double KE_min, double KE_max, TString key = ""){

  TF1 * empty = new TF1("", "pol2", KE_min, KE_max);
  empty -> SetParameter(0, 0.);
  empty -> SetParameter(1, 0.);
  empty -> SetParameter(2, 0.);

  double p_beam_plug[2];
  if(PID == 211){
    p_beam_plug[0] = 9.8 + 0.642037;
    //p_beam_plug[0] = 9.8;
    p_beam_plug[1] = 0.0002149;
  }
  else if(PID == 2212){
    p_beam_plug[0] = 22.06 + 4.15558;
    //p_beam_plug[0] = 22.06;
    p_beam_plug[1] = -0.01351;
  }
  else return empty;

  if(key == "proton_data") p_beam_plug[0] = p_beam_plug[0] + 5.98156;
  //p_beam_plug[0] = 0.;
  //p_beam_plug[1] = 0.;

  double this_p[4];
  this_p[0] = p[0] - p_beam_plug[0];
  this_p[1] = p[1] - p_beam_plug[1];
  this_p[2] = p[2];
  this_p[3] = p[3];

  TF1 * out = new TF1("", "pol2", KE_min, KE_max);
  out -> SetParameters(this_p);

  return out;
}

TF1 * KE_to_P(double p[], int PID, double KE_min, double KE_max, TString key = ""){

  TF1 * empty = new TF1("", "pol2", KE_min, KE_max);
  empty -> SetParameter(0, 0.);
  empty -> SetParameter(1, 0.);
  empty -> SetParameter(2, 0.);

  double mass = mass_pion;
  double p_beam_plug[2];
  if(PID == 211){
    p_beam_plug[0] = 9.8 + 0.642037;
    //p_beam_plug[0] = 9.8;
    p_beam_plug[1] = 0.0002149;
  }
  else if(PID == 2212){
    mass = mass_proton;
    p_beam_plug[0] = 22.06 + 4.15558;
    //p_beam_plug[0] = 22.06;
    p_beam_plug[1] = -0.01351;
  }
  else return empty;

  if(key == "proton_data") p_beam_plug[0] = p_beam_plug[0] + 5.98156;
  //p_beam_plug[0] = 0.;
  //p_beam_plug[1] = 0.;

  double this_p[4];
  this_p[0] = p[0] - p_beam_plug[0];
  this_p[1] = p[1] - p_beam_plug[1];
  this_p[2] = p[2];
  this_p[3] = p[3];

  double P_min = sqrt(pow(KE_min + mass, 2.) - pow(mass, 2.));
  double P_max = sqrt(pow(KE_max + mass, 2.) - pow(mass, 2.));

  TF1 * out = new TF1("", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_min, P_max);
  out -> SetParameters(this_p);

  return out;
}

void Test_MC_Truth(){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);

  double x_low = 200.; 
  double x_high = 1500.;

  TH1D* template_h = new TH1D("", "", 1., x_low, x_high);
  template_h -> SetStats(0);
   gStyle->SetOptTitle(0);   
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle("KE_{Beam Inst.} or P_{Beam Inst.} [MeV]");
  template_h -> GetYaxis() -> SetTitle("#Delta KE or #Delta P [MeV]");
  template_h -> GetYaxis() -> SetRangeUser(-40., 120.);
  template_h -> Draw();

  double parameter_pion_central[4] = {222.3, -0.6348, 0.0004500, 139.57};
  double parameter_pion_up[4] = {289.8, -0.8049, 0.0005592, 139.57};
  double parameter_pion_down[4] = {178.3, -0.5222, 0.0003757, 139.57};

  double parameter_proton_central[4] = {39.78, -0.2396, 0.0004498, 938.272};
  double parameter_proton_up[4] = {60.52, -0.3404, 0.0005857, 938.272};
  double parameter_proton_down[4] = {17.92, -0.1334, 0.0003075, 938.272};

  double KE_pion_low = 700.;
  double KE_pion_high = 1100.;
  TF1 * f_pion_KE = KE_remove_beam_plug(parameter_pion_central, 211, KE_pion_low, KE_pion_high);
  TF1 * f_pion_KE_up = KE_remove_beam_plug(parameter_pion_up, 211, KE_pion_low, KE_pion_high);
  TF1 * f_pion_KE_down = KE_remove_beam_plug(parameter_pion_down, 211, KE_pion_low, KE_pion_high);
  f_pion_KE -> SetLineColor(kRed);
  f_pion_KE_up -> SetLineColor(kRed);
  f_pion_KE_down -> SetLineColor(kRed);
  f_pion_KE_up -> SetLineStyle(5);
  f_pion_KE_down -> SetLineStyle(5);
  f_pion_KE -> Draw("lsame");
  f_pion_KE_up -> Draw("lsame");
  f_pion_KE_down -> Draw("lsame");

  double KE_proton_low = 300.;
  double KE_proton_high = 600.;
  TF1 * f_proton_KE = KE_remove_beam_plug(parameter_proton_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_proton_KE_up = KE_remove_beam_plug(parameter_proton_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_proton_KE_down = KE_remove_beam_plug(parameter_proton_down, 2212, KE_proton_low, KE_proton_high);
  f_proton_KE -> SetLineColor(kBlue);
  f_proton_KE_up -> SetLineColor(kBlue);
  f_proton_KE_down -> SetLineColor(kBlue);
  f_proton_KE_up -> SetLineStyle(5);
  f_proton_KE_down -> SetLineStyle(5);
  f_proton_KE -> Draw("lsame");
  f_proton_KE_up -> Draw("lsame");
  f_proton_KE_down -> Draw("lsame");

  cout << "[Test_MC_Truth] f_pion(870.) : " << f_pion_KE -> Eval(870.) << ", f_proton(433.) : " << f_proton_KE -> Eval(433.) << endl;

  TF1 * f_pion_P = KE_to_P(parameter_pion_central, 211, KE_pion_low, KE_pion_high);
  TF1 * f_pion_P_up = KE_to_P(parameter_pion_up, 211, KE_pion_low, KE_pion_high);
  TF1 * f_pion_P_down = KE_to_P(parameter_pion_down, 211, KE_pion_low, KE_pion_high);
  f_pion_P -> SetLineColor(kOrange);
  f_pion_P_up -> SetLineColor(kOrange);
  f_pion_P_down -> SetLineColor(kOrange);
  f_pion_P_up -> SetLineStyle(5);
  f_pion_P_down -> SetLineStyle(5);
  f_pion_P -> Draw("lsame");
  f_pion_P_up -> Draw("lsame");
  f_pion_P_down -> Draw("lsame"); 

  TF1 * f_proton_P = KE_to_P(parameter_proton_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_proton_P_up = KE_to_P(parameter_proton_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_proton_P_down = KE_to_P(parameter_proton_down, 2212, KE_proton_low, KE_proton_high);
  f_proton_P -> SetLineColor(kCyan);
  f_proton_P_up -> SetLineColor(kCyan);
  f_proton_P_down -> SetLineColor(kCyan);
  f_proton_P_up -> SetLineStyle(5);
  f_proton_P_down -> SetLineStyle(5);
  f_proton_P -> Draw("lsame");
  f_proton_P_up -> Draw("lsame");
  f_proton_P_down -> Draw("lsame");

  /*
  for(int i = 0; i < 100; i++){
    double P = P_proton_low + 5. * i;
    cout << "f_proton_P -> Eval(" << P << ") : " << f_proton_P -> Eval(P) << endl;
  }

  for(int i = 0; i < 100; i++){
    double P = P_pion_low + 5. * i;
    cout << "f_pion_P -> Eval(" << P << ") : " << f_pion_P -> Eval(P) << endl;
  }
  */

  TLegend *l = new TLegend(0.18, 0.77, 0.92, 0.92);
  l -> SetNColumns(2);
  l -> AddEntry(f_proton_KE, "p #DeltaKE(KE)", "l");
  l -> AddEntry(f_pion_KE, "#pi^{+} #DeltaKE(KE)", "l");
  l -> AddEntry(f_proton_P,  "p #DeltaP(P)", "l");
  l -> AddEntry(f_pion_P,  "#pi^{+} #DeltaP(P)", "l");
  l -> Draw("same");

  TString pdfname = "";
  pdfname = "./output/plots/BeamStudy/Upstream_Eloss/convert/Spectrometer_Bias_mc_true.pdf";
  c -> SaveAs(pdfname);

  c -> Close();

}


void Overlap_Proton_Results(){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);

  double x_low = 200.;
  double x_high = 1500.;

  TH1D* template_h = new TH1D("", "", 1., x_low, x_high);
  template_h -> SetStats(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle("KE_{Beam Inst.} or P_{Beam Inst.} [MeV]");
  template_h -> GetYaxis() -> SetTitle("#Delta KE or #Delta P [MeV]");
  template_h -> GetYaxis() -> SetRangeUser(-40., 120.);
  template_h -> Draw();

  double KE_proton_low = 300.;
  double KE_proton_high = 600.;

  double p_elas_true_central[4] = {39.78, -0.2396, 0.0004498, 938.272};
  double p_elas_true_up[4] = {60.52, -0.3404, 0.0005857, 938.272};
  double p_elas_true_down[4] = {17.92, -0.1334, 0.0003075, 938.272};

  double p_all_fitted_central[4] = {37.57, -0.2144, 0.0004282, 938.272};
  double p_all_fitted_up[4] = {61.69, -0.3305, 0.0005820, 938.272};
  double p_all_fitted_down[4] = {12.66, -0.09451, 0.0002699, 938.272};

  double p_data_fitted_central[4] = {47.35, -0.1748, 0.0003067, 938.272};
  double p_data_fitted_up[4] = {65.99, -0.2688, 0.0004365, 938.272};
  double p_data_fitted_down[4] = {28.27, -0.07865, 0.0001744, 938.272};

  TF1 * f_KE_elas_true_central = KE_remove_beam_plug(p_elas_true_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_KE_elas_true_up = KE_remove_beam_plug(p_elas_true_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_KE_elas_true_down = KE_remove_beam_plug(p_elas_true_down, 2212, KE_proton_low, KE_proton_high);
  f_KE_elas_true_central -> SetLineColor(kRed);
  f_KE_elas_true_up -> SetLineColor(kRed);
  f_KE_elas_true_down -> SetLineColor(kRed);
  f_KE_elas_true_up -> SetLineStyle(5);
  f_KE_elas_true_down -> SetLineStyle(5);
  f_KE_elas_true_central -> Draw("lsame");
  f_KE_elas_true_up -> Draw("lsame");
  f_KE_elas_true_down -> Draw("lsame");  

  TF1 * f_KE_all_fitted_central = KE_remove_beam_plug(p_all_fitted_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_KE_all_fitted_up = KE_remove_beam_plug(p_all_fitted_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_KE_all_fitted_down = KE_remove_beam_plug(p_all_fitted_down, 2212, KE_proton_low, KE_proton_high);
  f_KE_all_fitted_central -> SetLineColor(kGreen);
  f_KE_all_fitted_up -> SetLineColor(kGreen);
  f_KE_all_fitted_down -> SetLineColor(kGreen);
  f_KE_all_fitted_up -> SetLineStyle(5);
  f_KE_all_fitted_down -> SetLineStyle(5);
  f_KE_all_fitted_central -> Draw("lsame");
  f_KE_all_fitted_up -> Draw("lsame");
  f_KE_all_fitted_down -> Draw("lsame");

  TF1 * f_KE_data_fitted_central = KE_remove_beam_plug(p_data_fitted_central, 2212, KE_proton_low, KE_proton_high, "proton_data");
  TF1 * f_KE_data_fitted_up = KE_remove_beam_plug(p_data_fitted_up, 2212, KE_proton_low, KE_proton_high, "proton_data");
  TF1 * f_KE_data_fitted_down = KE_remove_beam_plug(p_data_fitted_down, 2212, KE_proton_low, KE_proton_high, "proton_data");
  f_KE_data_fitted_central -> SetLineColor(kBlack);
  f_KE_data_fitted_up -> SetLineColor(kBlack);
  f_KE_data_fitted_down -> SetLineColor(kBlack);
  f_KE_data_fitted_up -> SetLineStyle(5);
  f_KE_data_fitted_down -> SetLineStyle(5);
  f_KE_data_fitted_central -> Draw("lsame");
  f_KE_data_fitted_up -> Draw("lsame");
  f_KE_data_fitted_down -> Draw("lsame");
  
  TF1 * f_p_elas_true_central = KE_to_P(p_elas_true_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_p_elas_true_up = KE_to_P(p_elas_true_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_p_elas_true_down = KE_to_P(p_elas_true_down, 2212, KE_proton_low, KE_proton_high);
  f_p_elas_true_central -> SetLineColor(kRed);
  f_p_elas_true_up -> SetLineColor(kRed);
  f_p_elas_true_down -> SetLineColor(kRed);
  f_p_elas_true_up -> SetLineStyle(5);
  f_p_elas_true_down -> SetLineStyle(5);
  f_p_elas_true_central -> Draw("lsame");
  f_p_elas_true_up -> Draw("lsame");
  f_p_elas_true_down -> Draw("lsame");

  TF1 * f_p_all_fitted_central = KE_to_P(p_all_fitted_central, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_p_all_fitted_up = KE_to_P(p_all_fitted_up, 2212, KE_proton_low, KE_proton_high);
  TF1 * f_p_all_fitted_down = KE_to_P(p_all_fitted_down, 2212, KE_proton_low, KE_proton_high);
  f_p_all_fitted_central -> SetLineColor(kGreen);
  f_p_all_fitted_up -> SetLineColor(kGreen);
  f_p_all_fitted_down -> SetLineColor(kGreen);
  f_p_all_fitted_up -> SetLineStyle(5);
  f_p_all_fitted_down -> SetLineStyle(5);
  f_p_all_fitted_central -> Draw("lsame");
  f_p_all_fitted_up -> Draw("lsame");
  f_p_all_fitted_down -> Draw("lsame");

  TF1 * f_p_data_fitted_central = KE_to_P(p_data_fitted_central, 2212, KE_proton_low, KE_proton_high, "proton_data");
  TF1 * f_p_data_fitted_up = KE_to_P(p_data_fitted_up, 2212, KE_proton_low, KE_proton_high, "proton_data");
  TF1 * f_p_data_fitted_down = KE_to_P(p_data_fitted_down, 2212, KE_proton_low, KE_proton_high, "proton_data");
  f_p_data_fitted_central -> SetLineColor(kBlack);
  f_p_data_fitted_up -> SetLineColor(kBlack);
  f_p_data_fitted_down -> SetLineColor(kBlack);
  f_p_data_fitted_up -> SetLineStyle(5);
  f_p_data_fitted_down -> SetLineStyle(5);
  f_p_data_fitted_central -> Draw("lsame");
  f_p_data_fitted_up -> Draw("lsame");
  f_p_data_fitted_down -> Draw("lsame");

  double x_KE_all_fitted_zero = f_KE_all_fitted_central -> GetX(0., 300., 600.);
  cout << "[Overlap_Proton_Results] f_KE_all_fitted_central -> GetX(0., 300., 600.) : " << f_KE_all_fitted_central -> GetX(0., 300., 600.) << endl;
  cout << "[Overlap_Proton_Results] f_KE_data_fitted_central -> Eval(x_KE_all_fitted_zero) : " << f_KE_data_fitted_central -> Eval(x_KE_all_fitted_zero) << endl;

  TLegend *l = new TLegend(0.18, 0.77, 0.92, 0.92);
  l -> SetNColumns(2);
  l -> AddEntry(f_p_elas_true_central, "p with KE_{ff}^{true} (Elas. Only MC)", "l");
  l -> AddEntry(f_p_all_fitted_central, "p with KE_{ff}^{fitted} (All MC)", "l");
  l -> AddEntry(f_p_data_fitted_central,  "p with KE_{ff}^{fitted} (Data)", "l");
  l -> Draw("same");

  TString pdfname = "";
  pdfname = "./output/plots/BeamStudy/Upstream_Eloss/convert/Spectrometer_Bias_protons.pdf";
  c -> SaveAs(pdfname);

  c -> Close();
}

void Compare_Muon_Pion_in_P(){

  TCanvas *c = new TCanvas("", "", 800, 600);
  canvas_margin(c);
  gStyle -> SetOptStat(1111);
  gStyle->SetOptFit(0);

  double x_low = 700.;
  double x_high = 1300.;

  TH1D* template_h = new TH1D("", "", 1., x_low, x_high);
  template_h -> SetStats(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  template_h -> GetXaxis() -> SetTitle("P_{Beam Inst.} [MeV]");
  template_h -> GetYaxis() -> SetTitle("#Delta P [MeV]");
  template_h -> GetYaxis() -> SetRangeUser(-40., 100.);
  template_h -> Draw();  

  double muon_P_true_central[3] = {132.8, -0.3533, 0.0002326};
  double muon_P_true_up[3] = {394.6, -0.9012, 0.0005285};
  double muon_P_true_down[3] = {-76.58, 0.08846, -0.000009925};

  double muon_P_range_central[3] = {279.8, -0.6707, 0.0004078};
  double muon_P_range_up[3] = {397.1, -0.9096, 0.0005322};
  double muon_P_range_down[3] = {151.8, -0.4091, 0.0002713};

  double pion_P_true_central[3] = {268.7, -0.6600, 0.0004022};
  double pion_P_true_up[3] = {331.7, -0.7969, 0.0004777};
  double pion_P_true_down[3] = {223.2, -0.5612, 0.0003475};

  TF1 * f_muon_P_true_central = new TF1("", "pol2", 800., 1200.);
  TF1 * f_muon_P_true_up = new TF1("", "pol2", 800., 1200.);
  TF1 * f_muon_P_true_down = new TF1("", "pol2", 800., 1200.);
  f_muon_P_true_central -> SetParameters(muon_P_true_central);
  f_muon_P_true_up -> SetParameters(muon_P_true_up);
  f_muon_P_true_down -> SetParameters(muon_P_true_down);
  f_muon_P_true_central -> SetLineColor(kRed);
  f_muon_P_true_up -> SetLineColor(kRed);
  f_muon_P_true_down -> SetLineColor(kRed);
  f_muon_P_true_up -> SetLineStyle(5);
  f_muon_P_true_down -> SetLineStyle(5);
  f_muon_P_true_central -> Draw("lsame");
  f_muon_P_true_up -> Draw("lsame");
  f_muon_P_true_down -> Draw("lsame");

  TF1 * f_muon_P_range_central = new TF1("", "pol2", 800., 1200.);
  TF1 * f_muon_P_range_up = new TF1("", "pol2", 800., 1200.);
  TF1 * f_muon_P_range_down = new TF1("", "pol2", 800., 1200.);
  f_muon_P_range_central -> SetParameters(muon_P_range_central);
  f_muon_P_range_up -> SetParameters(muon_P_range_up);
  f_muon_P_range_down -> SetParameters(muon_P_range_down);
  f_muon_P_range_central -> SetLineColor(kGreen);
  f_muon_P_range_up -> SetLineColor(kGreen);
  f_muon_P_range_down -> SetLineColor(kGreen);
  f_muon_P_range_up -> SetLineStyle(5);
  f_muon_P_range_down -> SetLineStyle(5);
  f_muon_P_range_central -> Draw("lsame");
  f_muon_P_range_up -> Draw("lsame");
  f_muon_P_range_down -> Draw("lsame");

  TF1 * f_pion_P_true_central = new TF1("", "pol2", 800., 1200.);
  TF1 * f_pion_P_true_up = new TF1("", "pol2", 800., 1200.);
  TF1 * f_pion_P_true_down = new TF1("", "pol2", 800., 1200.);
  f_pion_P_true_central -> SetParameters(pion_P_true_central);
  f_pion_P_true_up -> SetParameters(pion_P_true_up);
  f_pion_P_true_down -> SetParameters(pion_P_true_down);
  f_pion_P_true_central -> SetLineColor(kBlack);
  f_pion_P_true_up -> SetLineColor(kBlack);
  f_pion_P_true_down -> SetLineColor(kBlack);
  f_pion_P_true_up -> SetLineStyle(5);
  f_pion_P_true_down -> SetLineStyle(5);
  f_pion_P_true_central -> Draw("lsame");
  f_pion_P_true_up -> Draw("lsame");
  f_pion_P_true_down -> Draw("lsame");
  
  TLegend *l = new TLegend(0.18, 0.77, 0.92, 0.92);
  l -> SetNColumns(3);
  l -> AddEntry(f_muon_P_true_central, "#mu with P_{ff}^{true}", "l");
  l -> AddEntry(f_muon_P_range_central, "#mu with P_{ff}^{range}", "l");
  l -> AddEntry(f_pion_P_true_central,  "#pi with P_{ff}^{true}", "l");
  l -> Draw("same");

  TString pdfname = "";
  pdfname = "./output/plots/BeamStudy/Upstream_Eloss/convert/Comparison_muon_pion_in_P.pdf";
  c -> SaveAs(pdfname);

  c -> Close();  


}

void Convert_Spectrometer_Bias(){
  setTDRStyle();

  Test_MC_Truth();
  Overlap_Proton_Results(); 
  Compare_Muon_Pion_in_P();
}
