#include "canvas_margin.h"
#include "mylib.h"
#include "LanGausFit.h"

TF1 * KE_to_P(double p[], int PID, double KE_min, double KE_max){

  TF1 * empty = new TF1("", "pol2", xmin, xmax);
  out -> SetParameter(0, 0.);
  out -> SetParameter(1, 0.);
  out -> SetParameter(2, 0.);

  double mass = mass_pion;
  double p_beam_plug[2];
  if(PID == 211){
    p_beam_plug[0] = 9.8;
    p_beam_plug[1] = 0.0002149;
  }
  else if(PID == 2212){
    mass = mass_proton;
    p_beam_plug[0] = 22.06;
    p_beam_plug[1] = -0.01351;
  }
  else return out;

  p[0] = p[0] - p_beam_plug[0];
  p[1] = p[1] - p_beam_plug[1];

  double P_min = sqrt(pow(KE_min + mass, 2.) - pow(mass, 2.));
  double P_max = sqrt(pow(KE_max + mass, 2.) - pow(mass, 2.));

  TF1 * out = new TF1("", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_min, P_max);
  out -> SetParameters(p);

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
  double parameter_pion_beamplug[2] = {9.8, 0.0002149};

  double parameter_proton_central[4] = {39.78, -0.2396, 0.0004498, 938.272};
  double parameter_proton_up[4] = {60.52, -0.3404, 0.0005857, 938.272};
  double parameter_proton_down[4] = {17.92, -0.1334, 0.0003075, 938.272};
  double parameter_proton_beamplug[2] = {22.06, -0.01351};

  for(int i = 0; i < 2; i++){
    parameter_pion_central[i] = parameter_pion_central[i] - parameter_pion_beamplug[i];
    parameter_pion_up[i] = parameter_pion_up[i] - parameter_pion_beamplug[i];
    parameter_pion_down[i] = parameter_pion_down[i] - parameter_pion_beamplug[i];

    parameter_proton_central[i] =parameter_proton_central[i] - parameter_proton_beamplug[i];
    parameter_proton_up[i] = parameter_proton_up[i] - parameter_proton_beamplug[i];
    parameter_proton_down[i] = parameter_proton_down[i] - parameter_proton_beamplug[i];
  }


  double KE_pion_low = 700.;
  double KE_pion_high = 1100.;
  TF1 * f_pion_KE = new TF1("f_pion_KE", "pol2", KE_pion_low, KE_pion_high);
  TF1 * f_pion_KE_up = new TF1("f_pion_KE_up", "pol2", KE_pion_low, KE_pion_high);
  TF1 * f_pion_KE_down = new TF1("f_pion_KE_down", "pol2", KE_pion_low, KE_pion_high);
  f_pion_KE -> SetParameters(parameter_pion_central);
  f_pion_KE_up -> SetParameters(parameter_pion_up);
  f_pion_KE_down -> SetParameters(parameter_pion_down);
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
  TF1 * f_proton_KE = new TF1("f_proton_KE", "pol2", KE_proton_low, KE_proton_high);
  TF1 * f_proton_KE_up = new TF1("f_proton_KE_up", "pol2", KE_proton_low, KE_proton_high);
  TF1 * f_proton_KE_down = new TF1("f_proton_KE_down", "pol2", KE_proton_low, KE_proton_high);
  f_proton_KE -> SetParameters(parameter_proton_central);
  f_proton_KE_up -> SetParameters(parameter_proton_up);
  f_proton_KE_down -> SetParameters(parameter_proton_down);
  f_proton_KE -> SetLineColor(kBlue);
  f_proton_KE_up -> SetLineColor(kBlue);
  f_proton_KE_down -> SetLineColor(kBlue);
  f_proton_KE_up -> SetLineStyle(5);
  f_proton_KE_down -> SetLineStyle(5);
  f_proton_KE -> Draw("lsame");
  f_proton_KE_up -> Draw("lsame");
  f_proton_KE_down -> Draw("lsame");

  double P_pion_low = sqrt(pow(KE_pion_low + mass_pion, 2.) - pow(mass_pion, 2.));
  double P_pion_high = sqrt(pow(KE_pion_high + mass_pion, 2.) - pow(mass_pion, 2.));
  double P_proton_low = sqrt(pow(KE_proton_low + mass_proton, 2.) - pow(mass_proton, 2.));
  double P_proton_high = sqrt(pow(KE_proton_high + mass_proton, 2.) - pow(mass_proton, 2.));


  /*
  TF1 * f_pion_P_in_KE = new TF1("f_pion_P_in_KE", "sqrt( pow([0] + (1. + [1]) * x + [2] * x * x + [3], 2.)  - [3] * [3] ) - sqrt( pow(x + [3] , 2.) - [3] * [3] )", KE_pion_low, KE_pion_high);
  f_pion_P_in_KE -> SetParameter(0, 153.8);
  f_pion_P_in_KE -> SetParameter(1, -0.5034);
  f_pion_P_in_KE -> SetParameter(2, 0.0003763);
  f_pion_P_in_KE -> SetParameter(3, mass_pion);
  f_pion_P_in_KE -> SetLineColor(kBlue);
  f_pion_P_in_KE -> Draw("lsame");
  */

  //TF1 * f_pion_P = new TF1("f_pion_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x - [3] * [3]) - [3]) + [2] * (sqrt(x*x - [3] * [3]) - [3]) * (sqrt(x*x - [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x - [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_pion_low, P_pion_high);
  TF1 * f_pion_P = new TF1("f_pion_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_pion_low, P_pion_high);
  TF1 * f_pion_P_up = new TF1("f_pion_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_pion_low, P_pion_high);
  TF1 * f_pion_P_down = new TF1("f_pion_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_pion_low, P_pion_high);
  f_pion_P -> SetParameters(parameter_pion_central);
  f_pion_P_up -> SetParameters(parameter_pion_up);
  f_pion_P_down -> SetParameters(parameter_pion_down);
  f_pion_P -> SetLineColor(kOrange);
  f_pion_P_up -> SetLineColor(kOrange);
  f_pion_P_down -> SetLineColor(kOrange);
  f_pion_P_up -> SetLineStyle(5);
  f_pion_P_down -> SetLineStyle(5);
  //f_pion_P -> Draw("lsame");
  //f_pion_P_up -> Draw("lsame");
  //f_pion_P_down -> Draw("lsame");

  /*
  TF1 * f_proton_P_in_KE = new TF1("f_proton_P_in_KE", "sqrt( pow([0] + (1. + [1]) * x + [2] * x * x + [3], 2.)  - [3] * [3] ) - sqrt( pow(x + [3] , 2.) - [3] * [3] )", KE_proton_low, KE_proton_high);
  f_proton_P_in_KE -> SetParameter(0,  29.27);
  f_proton_P_in_KE -> SetParameter(1,  -0.28189);
  f_proton_P_in_KE -> SetParameter(2,  0.0005164);
  f_proton_P_in_KE -> SetParameter(3,  mass_proton);
  f_proton_P_in_KE-> SetLineColor(kBlue);
  f_proton_P_in_KE -> Draw("lsame");
  */

  TF1 * f_proton_P = new TF1("f_proton_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_proton_low, P_proton_high);
  TF1 * f_proton_P_up = new TF1("f_proton_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_proton_low, P_proton_high);
  TF1 * f_proton_P_down = new TF1("f_proton_P", "sqrt( pow([0] + (1. + [1]) * (sqrt(x*x + [3] * [3]) - [3]) + [2] * (sqrt(x*x + [3] * [3]) - [3]) * (sqrt(x*x + [3] * [3]) - [3]) + [3], 2.)  - [3] * [3] ) - sqrt( pow((sqrt(x*x + [3] * [3]) - [3]) + [3] , 2.) - [3] * [3] )", P_proton_low, P_proton_high);
  f_proton_P -> SetParameters(parameter_proton_central);
  f_proton_P_up -> SetParameters(parameter_proton_up);
  f_proton_P_down -> SetParameters(parameter_proton_down);
  f_proton_P -> SetLineColor(kCyan);
  f_proton_P_up -> SetLineColor(kCyan);
  f_proton_P_down -> SetLineColor(kCyan);
  f_proton_P_up -> SetLineStyle(5);
  f_proton_P_down -> SetLineStyle(5);
  //f_proton_P -> Draw("lsame");
  //f_proton_P_up -> Draw("lsame");
  //f_proton_P_down -> Draw("lsame");

  for(int i = 0; i < 100; i++){
    double P = P_proton_low + 5. * i;
    cout << "f_proton_P -> Eval(" << P << ") : " << f_proton_P -> Eval(P) << endl;
  }

  for(int i = 0; i < 100; i++){
    double P = P_pion_low + 5. * i;
    cout << "f_pion_P -> Eval(" << P << ") : " << f_pion_P -> Eval(P) << endl;
  }
  
  TLegend *l = new TLegend(0.18, 0.77, 0.92, 0.92);
  l -> SetNColumns(2);
  l -> AddEntry(f_proton_KE, "p #DeltaKE(KE)", "l");
  l -> AddEntry(f_pion_KE, "#pi^{+} #DeltaKE(KE)", "l");
  //l -> AddEntry(f_proton_P,  "p #DeltaP(P)", "l");
  //l -> AddEntry(f_pion_P,  "#pi^{+} #DeltaP(P)", "l");
  l -> Draw("same");

  TString pdfname = "";
  pdfname = "./output/plots/BeamStudy/Upstream_Eloss/convert/Spectrometer_Bias_mc_true.pdf";
  c -> SaveAs(pdfname);

  c -> Close();

}

void Convert_Spectrometer_Bias(){
  setTDRStyle();

  Test_MC_Truth();
}
