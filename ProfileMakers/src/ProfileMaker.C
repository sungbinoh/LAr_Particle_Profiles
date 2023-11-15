#include "ProfileMaker.h"
#include "TF1.h"

const double pitch = 0.65; // [cm]
void ProfileMaker::Execute(){
  cout << "[ProfileMaker::Execute] Start" << endl;

  // == Test functions
  //KE_to_ResLength_BB(200, mass_muon);
  cout << ResLength_to_KE_BB(KE_to_ResLength_BB(200, mass_muon), mass_muon) << endl;
  cout << ResLength_to_KE_BB(KE_to_ResLength_BB(20, mass_pion), mass_pion) << endl;
  cout << ResLength_to_KE_BB(3.5, mass_pion) << endl;
  double c = KE_to_Momentum(ResLength_to_KE_BB(38, mass_pion), mass_pion);
  cout << "dEdx" << dEdx.dEdx_Bethe_Bloch(400, mass_pion) << endl;
  double b = KE_to_Momentum(771.932, mass_pion);
  KE_to_ResLength_BB(300., mass_pion);
  cout << ResLength_to_KE_BB(20., mass_pion) << endl;
  double a = KE_to_Momentum(500., mass_pion);
  a = KE_to_Momentum(400., mass_pion);
  a = KE_to_Momentum(1000., mass_pion);

  // Wmax
  cout << "[Get_Wmax] dEdx.Get_Wmax(1000., mass_pion) : " << dEdx.Get_Wmax(1000., mass_pion) << endl;
  cout << "[Get_Wmax] dEdx.Get_Wmax(500., mass_pion) : " << dEdx.Get_Wmax(500., mass_pion) << endl;
  cout << "[Get_Wmax] dEdx.Get_Wmax(100., mass_pion) : " << dEdx.Get_Wmax(100., mass_pion) << endl;
  cout << "[Get_Wmax] dEdx.Get_Wmax(50., mass_pion) : " << dEdx.Get_Wmax(50., mass_pion) << endl;

  // == Range to KE
  
  // == KE to range
  cout << "[KE_to_ResLength_BB] proton KE 40 MeV (P " << KE_to_Momentum(40, 938.272) << " MeV) : " << KE_to_ResLength_BB(40., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 50 MeV (P " << KE_to_Momentum(50, 938.272) << " MeV) : " << KE_to_ResLength_BB(50., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 60 MeV (P " << KE_to_Momentum(60, 938.272) << " MeV) : " << KE_to_ResLength_BB(60., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 70 MeV (P " << KE_to_Momentum(70, 938.272) << " MeV) : " << KE_to_ResLength_BB(70., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 80 MeV (P " << KE_to_Momentum(80, 938.272) << " MeV) : " << KE_to_ResLength_BB(80., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 90 MeV (P " << KE_to_Momentum(90, 938.272) << " MeV) : " << KE_to_ResLength_BB(90., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 100 MeV (P " << KE_to_Momentum(100, 938.272) << " MeV) : " << KE_to_ResLength_BB(100., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 150 MeV (P " << KE_to_Momentum(150, 938.272) << " MeV) : " << KE_to_ResLength_BB(150., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 200 MeV (P " << KE_to_Momentum(200, 938.272) << " MeV) : " << KE_to_ResLength_BB(200., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 250 MeV (P " << KE_to_Momentum(250, 938.272) << " MeV) : " << KE_to_ResLength_BB(250., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 300 MeV (P " << KE_to_Momentum(300, 938.272) << " MeV) : " << KE_to_ResLength_BB(300., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 500 MeV (P " << KE_to_Momentum(500, 938.272) << " MeV) : " << KE_to_ResLength_BB(500., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 700 MeV (P " << KE_to_Momentum(700, 938.272) << " MeV) : " << KE_to_ResLength_BB(700., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 1000 MeV (P " << KE_to_Momentum(1000, 938.272) << " MeV) : " << KE_to_ResLength_BB(1000., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 1500 MeV (P " << KE_to_Momentum(1500, 938.272) << " MeV) : " << KE_to_ResLength_BB(1500., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 2000 MeV (P " << KE_to_Momentum(2000, 938.272) << " MeV) : " << KE_to_ResLength_BB(2000., 938.272) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] proton KE 2300 MeV (P " << KE_to_Momentum(2300, 938.272) << " MeV) : " << KE_to_ResLength_BB(2300., 938.272) << " cm" << endl;


  cout << "[KE_to_ResLength_BB] pion KE 50 MeV (P " << KE_to_Momentum(50, mass_pion) << " MeV) : " << KE_to_ResLength_BB(50., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 100 MeV (P " << KE_to_Momentum(100, mass_pion) << " MeV) : " << KE_to_ResLength_BB(100., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 150 MeV (P " << KE_to_Momentum(150, mass_pion) << " MeV) : " << KE_to_ResLength_BB(150., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 200 MeV (P " << KE_to_Momentum(200, mass_pion) << " MeV) : " << KE_to_ResLength_BB(200., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 250 MeV (P " << KE_to_Momentum(250, mass_pion) << " MeV) : " << KE_to_ResLength_BB(250., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 300 MeV (P " << KE_to_Momentum(300, mass_pion) << " MeV) : " << KE_to_ResLength_BB(300., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 400 MeV (P " << KE_to_Momentum(400, mass_pion) << " MeV) : " << KE_to_ResLength_BB(400., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 500 MeV (P " << KE_to_Momentum(500, mass_pion) << " MeV) : " << KE_to_ResLength_BB(500., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 700 MeV (P " << KE_to_Momentum(700, mass_pion) << " MeV) : " << KE_to_ResLength_BB(700., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 1000 MeV (P " << KE_to_Momentum(1000, mass_pion) << " MeV) : " << KE_to_ResLength_BB(1000., mass_pion) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] pion KE 1500 MeV (P " << KE_to_Momentum(1500, mass_pion) << " MeV) : " << KE_to_ResLength_BB(1500., mass_pion) << " cm" << endl;

  cout << "[KE_to_ResLength_BB] muon KE 100 MeV (P " << KE_to_Momentum(100, mass_muon) << " MeV) : " << KE_to_ResLength_BB(100., mass_muon) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] muon KE 150 MeV (P " << KE_to_Momentum(150, mass_muon) << " MeV) : " << KE_to_ResLength_BB(150., mass_muon) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] muon KE 250 MeV (P " << KE_to_Momentum(250, mass_muon) << " MeV) : " << KE_to_ResLength_BB(250., mass_muon) << " cm" << endl;
  cout << "[KE_to_ResLength_BB] muon KE 450 MeV (P " << KE_to_Momentum(450, mass_muon) << " MeV) : " << KE_to_ResLength_BB(450., mass_muon) << " cm" << endl;


  cout << "[P_to_ResLength_BB] muon Momentum 1000 MeV/c (KE " << Momentum_to_KE(1000., mass_muon) << " MeV/c) : " << KE_to_ResLength_BB(Momentum_to_KE(1000., mass_muon), mass_muon) << " cm" << endl;
  cout << "[P_to_ResLength_BB] pion Momentum 1000 MeV/c (KE " << Momentum_to_KE(1000., mass_pion) << " MeV/c) : " << KE_to_ResLength_BB(Momentum_to_KE(1000., mass_pion), mass_pion) << " cm" << endl;
  cout << "[P_to_ResLength_BB] proton Momentum 6000 MeV/c (KE " << Momentum_to_KE(6000., mass_proton) << " MeV/c) : " << KE_to_ResLength_BB(Momentum_to_KE(6000., mass_pion), mass_proton) << " cm" << endl;


  // == KE to dE/dx BB
  cout << "<dE/dx> pion KE = 20 MeV (P = " << KE_to_Momentum(20, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(20, mass_pion) << "\t" << dEdx.dEdx_Landau_Vavilov(20, 0.65, mass_pion) << endl;
  cout << "<dE/dx> pion KE = 200 MeV (P = " << KE_to_Momentum(200, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(200, mass_pion) << "\t" << dEdx.dEdx_Landau_Vavilov(200, 0.65, mass_pion) << endl;
  cout << "<dE/dx> pion KE = 500 MeV (P = " << KE_to_Momentum(500, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(500, mass_pion) << "\t" << dEdx.dEdx_Landau_Vavilov(500, 0.65, mass_pion) << endl;
  cout << "<dE/dx> pion KE = 947 MeV (P = " << KE_to_Momentum(947, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(947, mass_pion) << "\t" << dEdx.dEdx_Landau_Vavilov(947, 0.65, mass_pion) << endl;
  cout << "<dE/dx> proton KE = 20 MeV (P = " << KE_to_Momentum(20, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(20, mass_proton) << "\t" << dEdx.dEdx_Landau_Vavilov(20, 0.65, mass_proton) << endl;
  cout << "<dE/dx> proton KE = 200 MeV (P = " << KE_to_Momentum(200, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(200, mass_proton) << "\t" << dEdx.dEdx_Landau_Vavilov(200, 0.65, mass_proton) << endl;
  cout << "<dE/dx> proton KE = 500 MeV (P = " << KE_to_Momentum(500, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(500, mass_proton) << "\t" << dEdx.dEdx_Landau_Vavilov(500, 0.65, mass_proton) << endl;
  /*
  KE_to_ResLength_LV(200, mass_muon);
  KE_to_ResLength_LV(300, mass_muon);
  KE_to_ResLength_LV(400, mass_muon);
  KE_to_ResLength_LV(800, mass_muon);
  */
  /*
  KE_to_ResLength_Urban(200, mass_muon);
  KE_to_ResLength_Urban(300, mass_muon);
  KE_to_ResLength_Urban(400, mass_muon);
  KE_to_ResLength_Urban(800, mass_muon);
  */
  outfile -> cd();

  const int N_KEs = 500;
  vector<double> Sum_dEdx;
  vector<double> KEs;
  Sum_dEdx.resize(N_KEs);
  KEs.resize(N_KEs);
  for(int i_KE = 0; i_KE < N_KEs; i_KE++){
    double this_KE = 50. + 2. * (i_KE + 0.);
    double this_Sum_dEdx = Sum_dEdx_length_KE0(this_KE, mass_muon, 40.);
    KEs[i_KE] = this_KE;
    Sum_dEdx[i_KE] = this_Sum_dEdx;
    //cout << "this_KE : " << this_KE << ", this_Sum_dEdx : " << this_Sum_dEdx << endl;
  }

  TGraph *gr_sum_dedx = new TGraph(N_KEs, &KEs[0], &Sum_dEdx[0]);
  gr_sum_dedx -> SetName("Sum_dEdx_50cm");
  gr_sum_dedx -> Write();

  /*
  Validate_Landau(200, mass_muon);
  Validate_Landau(50, mass_muon);
  Validate_Landau(10, mass_muon);

  Mimic_GEANT(200, mass_muon);
  Mimic_GEANT(1000, mass_muon);
  */

  ///////////////////////////
  // == Produce profiles
  ///////////////////////////
  //Produce_Profile_all();
  //Produce_Profile("pion", mass_pion);
  //Produce_Profile("proton", mass_proton);
  //Produce_Profile("muon", mass_muon);

  Produce_KE_vs_dEdx("pion", mass_pion);
  //Produce_KE_vs_dEdx("proton", mass_proton);
  //Produce_kappa("pion_0p65cm", mass_pion, 0.65);
  //Produce_kappa("proton_0p65cm", mass_proton, 0.65);
  //Produce_KE_vs_Range("pion", mass_pion);
 
  ///////////////////////////
  // == Produce PDFs
  ///////////////////////////
  Produce_dEdx_PDF("dEdx_PDF_pion_100MeV_0p1cm", mass_pion, 100., 0.1, 0., 5., 2.0);
  Produce_dEdx_PDF("dEdx_PDF_pion_80MeV_0p1cm", mass_pion, 80., 0.1, 0., 10., 2.0);
  Produce_dEdx_PDF("dEdx_PDF_pion_1MeV_0p1cm", mass_pion, 1., 0.1, 50., 90., 2.0);

  /*
  Produce_dEdx_PDF("dEdx_PDF_proton_50MeV_0p65cm", mass_proton, 50, 0.65, 0., 20., 1.5);
  Produce_dEdx_PDF("dEdx_PDF_proton_80MeV_0p65cm", mass_proton, 80, 0.65, 0., 20., 1.5);
  Produce_dEdx_PDF("dEdx_PDF_proton_100MeV_0p65cm", mass_proton, 100, 0.65, 0., 20., 1.5);

  Produce_dEdx_PDF("dEdx_PDF_pion_947MeV_0p65cm", mass_pion, 947, 0.65, 1., 5., 2.0);
  Produce_dEdx_PDF("dEdx_PDF_pion_500MeV_0p65cm", mass_pion, 500, 0.65, 1., 5., 2.0);
  Produce_dEdx_PDF("dEdx_PDF_proton_500MeV_0p65cm", mass_proton, 500, 0.65, 1., 5., 1.5);
  Produce_dEdx_PDF("dEdx_PDF_pion_200MeV_0p65cm", mass_pion, 200, 0.65, 1., 5., 2.0);
  Produce_dEdx_PDF("dEdx_PDF_proton_200MeV_0p65cm", mass_proton, 200, 0.65, 2., 8., 1.5);
  Produce_dEdx_PDF("dEdx_PDF_pion_20MeV_0p65cm", mass_pion, 20, 0.65, 4., 8., 1.5);
  Produce_dEdx_PDF("dEdx_PDF_proton_20MeV_0p65cm", mass_proton, 20, 0.65, 23., 28., 1.5);

  ///////////////////////////
  // == Produce Likelihoods
  ///////////////////////////
  Produce_dEdx_likelihood("dEdx_likelihood_pion_1p65MeVcm_0p65cm", mass_pion, 1.65, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_pion_1p80MeVcm_0p65cm", mass_pion, 1.80, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_pion_5p75MeVcm_0p65cm", mass_pion, 5.75, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_pion_5p88MeVcm_0p65cm", mass_pion, 5.88, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_pion_5p99MeVcm_0p65cm", mass_pion, 5.99, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_proton_2p49MeVcm_0p65cm", mass_proton, 2.49, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_proton_4p40MeVcm_0p65cm", mass_proton, 4.40, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_proton_25p2MeVcm_0p65cm", mass_proton, 25.2, 0.65, 10., 30., 2.0);

  Produce_dEdx_likelihood("dEdx_likelihood_pion_2p15MeVcm_0p65cm", mass_proton, 2.15, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_proton_2p79MeVcm_0p65cm", mass_proton, 2.79, 0.65, 1., 2000., 2.0);
  Produce_dEdx_likelihood("dEdx_likelihood_proton_25p2MeVcm_0p65cm", mass_proton, 25.2, 0.65, 1., 50., 2.0); 

  Produce_Range_from_Momentum_Gaussian("Range_from_Momentum_mu1007MeV_sigma68p17_pion", mass_pion, 1007., 68.17);
  */
}

/*
void ProfileMaker::Produce_Profile_all(){



}
*/
void ProfileMaker::Produce_KE_vs_dEdx(TString name, double mass){
  
  vector<double> KE_vec;
  vector<double> dEdx_BB_vec;
  double KE_min = 5;
  double KE_max = 2000.;
  int N_step = 10000;
  double step = (KE_max - KE_min) / (N_step + 0.);
  for(int i = 0; i < N_step; i++){
    double this_KE = KE_min + step * (i + 0.);
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(this_KE, mass);
    KE_vec.push_back(this_KE);
    dEdx_BB_vec.push_back(this_dEdx);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);

  TH1D *template_h = new TH1D("", "", 1., 0., 2000.);
  template_h ->GetXaxis() -> SetTitle("KE [MeV]");
  template_h ->GetYaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetYaxis() -> SetRangeUser(1.5, 3.);
  template_h -> Draw();

  TGraph *dEdx_range = new TGraph(N_step, &KE_vec[0], &dEdx_BB_vec[0]);
  dEdx_range -> SetName(name + "_range_vs_dEdx");
  dEdx_range -> SetLineColor(kRed);
  dEdx_range -> SetLineWidth(3);
  dEdx_range -> Draw("same");

  c -> SaveAs("./output/" + name + "_KE_vs_dEdx.pdf");
  c -> Close();
}

void ProfileMaker::Produce_Profile(TString name, double mass){
  // == Range vs dE/dx
  vector<double> KE_vec;
  vector<double> KE_to_mom_vec;
  vector<double> KE_to_range_BB_vec;
  vector<double> KE_to_dEdx_BB_vec;
  for(int i = 0; i < 1000; i++){
    double this_KE = 10.+ 2. * (i + 0.);
    double this_mom = KE_to_Momentum(this_KE, mass);
    double this_range = KE_to_ResLength_BB(this_KE, mass);
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(this_KE, mass);
    KE_vec.push_back(this_KE);
    KE_to_mom_vec.push_back(this_mom);
    KE_to_range_BB_vec.push_back(this_range);
    KE_to_dEdx_BB_vec.push_back(this_dEdx);
  }
  
  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  TH1D *template_h = new TH1D("", "", 1., 0.1, 1000.);
  template_h ->GetXaxis() -> SetTitle("Range [cm]");
  template_h ->GetYaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h -> GetYaxis() -> SetRangeUser(0., 10.);
  template_h -> Draw();
  TGraph *dEdx_range = new TGraph(1000, &KE_to_range_BB_vec[0], &KE_to_dEdx_BB_vec[0]);
  dEdx_range -> SetName(name + "_range_vs_dEdx");
  dEdx_range -> SetLineColor(kRed);
  dEdx_range -> SetLineWidth(3);
  dEdx_range -> Draw("same");

  c -> SaveAs("./output/" + name + "_range_vs_dEdx.pdf");
  c -> Close();
  dEdx_range -> Write();

  KE_vec.clear();
  KE_to_mom_vec.clear();
  KE_to_range_BB_vec.clear();
  KE_to_dEdx_BB_vec.clear();
}



void ProfileMaker::Produce_Range_from_Momentum_Gaussian(TString name, double mass, double mean, double sigma){
  TF1 *this_gaus = new TF1("this_gaus", "gaus", mean * 0.5, mean * 1.5);
  this_gaus -> SetParameter(0, 10.);
  this_gaus -> SetParameter(1, mean);
  this_gaus -> SetParameter(2, sigma);
  TH1D * h_P = new TH1D("", "", 300., 0., 1500.);
  TH1D * h_range_muon =new TH1D("", "", 160., 0., 800.);
  TH1D * h_range_pion = new TH1D("", "", 160., 0., 800.);
  TH1D * h_range_proton = new TH1D("", "", 80., 0., 800.);
  int N_trial = 100000.;
  for(int i = 0; i < N_trial; i++){
    if(i % 1000 == 0) cout << i << " / " << N_trial << endl;
    double this_P = this_gaus -> GetRandom();
    h_P -> Fill(this_P);

    double this_muon_range = KE_to_ResLength_BB(Momentum_to_KE(this_P, mass_muon), mass_muon);
    double this_pion_range = KE_to_ResLength_BB(Momentum_to_KE(this_P, mass_pion), mass_pion);
    double this_proton_range = KE_to_ResLength_BB(Momentum_to_KE(this_P, mass_proton), mass_proton);
    h_range_muon -> Fill(this_muon_range);
    h_range_pion -> Fill(this_pion_range);
    h_range_proton -> Fill(this_proton_range);
    //cout << "this_P ; " << this_P << ", KE_pion : " << Momentum_to_KE(this_P, mass_pion) << ", this_pion_range : " << this_pion_range << ", KE_proton : " << Momentum_to_KE(this_P, mass_proton) << ", this_proton_range : " << this_proton_range << endl;
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);

  // == Draw h_P
  double y_max = h_P -> GetMaximum();
  y_max = std::max(y_max, h_range_muon -> GetMaximum());
  y_max = std::max(y_max, h_range_pion -> GetMaximum());
  y_max = std::max(y_max, h_range_proton -> GetMaximum());
  TH1D *template_h = new TH1D("", "", 1., 0., 1400.);
  template_h -> GetXaxis() -> SetTitle("Momentum [MeV/c] || Range [cm]");
  template_h -> GetXaxis() -> SetTitleOffset(1.3);
  template_h -> GetYaxis() -> SetTitle("Events");
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.5);
  template_h -> Draw();

  h_P -> SetLineColor(kRed);
  h_P -> SetLineWidth(2);
  h_P -> Draw("histsame");
  h_range_muon -> SetLineColor(kCyan);
  h_range_muon -> SetLineWidth(2);
  h_range_muon -> Draw("histsame");
  h_range_pion -> SetLineColor(kBlue);
  h_range_pion -> SetLineWidth(2);
  h_range_pion -> Draw("histsame");
  h_range_proton -> SetLineColor(kGreen);
  h_range_proton -> SetLineWidth(2);
  h_range_proton -> Draw("histsame");

  c -> SaveAs("./output/plots/BeamStudy/" + name + ".pdf");
  
  c -> Close();
}

void ProfileMaker::Produce_kappa(TString name, double mass, double width){
  vector<double> KE_vec;
  vector<double> kappa_vec;
  for(int i = 0; i < 1000; i++){
    double this_KE = 1.+ 2. * (i + 0.);
    double this_xi = dEdx.Get_Landau_xi(this_KE, width, mass);
    double this_Wmax = dEdx.Get_Wmax(this_KE, mass);
    double this_kappa = this_xi / this_Wmax;
    KE_vec.push_back(this_KE);
    kappa_vec.push_back(this_kappa);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  c -> SetLogx();
  c -> SetLogy();
  TH1D *template_h = new TH1D("", "", 1., 0.1, 2000.);
  template_h ->GetXaxis() -> SetTitle("KE [MeV]");
  template_h ->GetXaxis() -> SetTitleOffset(1.3);
  template_h ->GetYaxis() -> SetTitle("#kappa ");
  template_h -> GetYaxis() -> SetRangeUser(0.0001, 1000.);
  template_h -> Draw();
  TGraph *KE_kappa = new TGraph(1000, &KE_vec[0], &kappa_vec[0]);
  KE_kappa -> SetName(name + "_KE_vs_Kappa");
  KE_kappa -> SetLineColor(kRed);
  KE_kappa -> SetLineWidth(3);
  KE_kappa -> Draw("same");
  TLine *line1 = new TLine(0.1, 0.01, 2000., 0.01);
  line1 -> SetLineColor(kBlue);
  line1 -> SetLineStyle(7);
  line1 -> Draw("same");
  TLine *line2 = new TLine(0.1, 10., 2000., 10.);
  line2 -> SetLineColor(kBlue);
  line2 -> SetLineStyle(7);
  line2 -> Draw("same");

  c -> SaveAs("./output/" + name + "_KE_vs_Kappa.pdf");

  KE_kappa -> Write();

  KE_vec.clear();
  kappa_vec.clear();
}

void ProfileMaker::Produce_dEdx_PDF(TString name, double mass, double KE, double width, double xmin, double xmax, double ymax){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  
  double this_xi = dEdx.Get_Landau_xi(KE, width, mass);
  double this_Wmax = dEdx.Get_Wmax(KE, mass);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = dEdx.dEdx_Bethe_Bloch(KE, mass);

  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, width};
  cout << "[ProfileMaker::Produce_dEdx_PDF] par[0] : " << par[0] << ", par[1] : " << par[1] << endl;

  TF1 *this_PDF = dEdx.dEdx_PDF(par);
  this_PDF -> SetNpx(1000);
  cout << "[ProfileMaker::Produce_dEdx_PDF] this_PDF -> Eval(2.1) : " << this_PDF -> Eval(2.1) << endl;
  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  //c -> SetLogx();
  TH1D *template_h = new TH1D("", "", 1., xmin, xmax);
  template_h ->GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h ->GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> Draw();

  double this_MPV =  dEdx.dEdx_Landau_Vavilov(KE, width, mass);
  TLine *l_mean = new TLine(this_dEdx_BB, 0., this_dEdx_BB, ymax);
  l_mean -> SetLineColor(kGreen);
  l_mean -> SetLineStyle(7);
  l_mean -> SetLineWidth(3);
  l_mean -> Draw("same");
  TLine *l_MPV = new TLine(this_MPV, 0., this_MPV, ymax);
  l_MPV -> SetLineColor(kBlue);
  l_MPV -> SetLineStyle(6);
  l_MPV -> Draw("same");
  
  this_PDF -> Draw("lsame");
  cout << "[ProfileMaker::Produce_dEdx_PDF] Integral : " << this_PDF -> Integral(0., 100.) << endl;
  cout << "[ProfileMaker::Produce_dEdx_PDF] this_PDF.Mean() : " << this_PDF -> Mean(0., 100.) << ", this_dEdx_BB : " << this_dEdx_BB << endl;


  TLegend *l = new TLegend(0.6, 0.6, 0.85, 0.85);
  l -> AddEntry(l_mean, "<dE/dx>", "l");
  l -> AddEntry(l_MPV, "(dE/dx)_{MPV}", "l");
  l -> SetLineColor(kWhite);
  l -> Draw("same");

  TLatex latex_kappa, latex_info;
  latex_kappa.SetNDC();
  latex_info.SetNDC();
  latex_kappa.SetTextSize(0.040);
  latex_info.SetTextSize(0.040);
  TString PDF_region_str = "";
  if(this_kappa < 0.01) PDF_region_str = "Landau";
  else if(this_kappa > 10.) PDF_region_str = "Gaussain";
  else PDF_region_str = "Vavilov";
  TString kappa_str = Form("%.3f", this_kappa);
  TString width_str = Form("%.2f", width);
  latex_kappa.DrawLatex(0.10, 0.92, "#kappa = " + kappa_str + " (" + PDF_region_str + ")");
  latex_info.DrawLatex(0.52, 0.92, "LAr, #rho = 1.39 g/cm^{3}, #deltax = " + width_str);

  c -> SaveAs("./output/" + name + "_dEdx_PDF.pdf");

}

void ProfileMaker::Produce_dEdx_likelihood(TString name, double mass, double dEdx_value, double width, double xmin, double xmax, double ymax){

  vector<double> KE_vec;
  vector<double> likelihood_vec;

  double maxium_l = -1.;
  double maximum_l_KE = 1.;
  for(int i = 0; i < 20000; i++){
    double this_KE = 1.+ 0.1 * (i + 0.);
  
    double gamma = (this_KE/mass)+1.0;
    double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  
    double this_xi = dEdx.Get_Landau_xi(this_KE, width, mass);
    double this_Wmax = dEdx.Get_Wmax(this_KE, mass);
    double this_kappa = this_xi / this_Wmax;
    double this_dEdx_BB = dEdx.dEdx_Bethe_Bloch(this_KE, mass);
    double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, width};
    TF1 *this_PDF = dEdx.dEdx_PDF(par);
    double this_likelihood = this_PDF -> Eval(dEdx_value);

    KE_vec.push_back(this_KE);
    likelihood_vec.push_back(this_likelihood);
    if(maxium_l < this_likelihood){
      maxium_l = this_likelihood;
      maximum_l_KE = this_KE;
    }
  }
  
 
  TCanvas *c = new TCanvas("", "", 800, 600);
  c -> SetLogx();
  gStyle->SetOptStat(0);

  double y_max = *max_element(likelihood_vec.begin(), likelihood_vec.end());
  TH1D *template_h = new TH1D("", "", 1., xmin, xmax);
  template_h ->GetXaxis() -> SetTitle("KE [MeV]");
  template_h ->GetYaxis() -> SetTitle("Likelihood");
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.3);
  template_h -> Draw();

  TGraph *L_gr = new TGraph(20000, &KE_vec[0], &likelihood_vec[0]);
  L_gr -> SetName(name + "_KE_vs_Kappa");
  L_gr -> SetLineColor(kRed);
  L_gr -> SetLineWidth(3);
  L_gr -> Draw("same");

  TString max_KE_str = Form("%.1f", maximum_l_KE);
  TLine *l_max = new TLine(maximum_l_KE, 0., maximum_l_KE, y_max);
  l_max -> SetLineColor(kBlue);
  l_max -> SetLineStyle(7);
  l_max -> Draw("same");
  
  TLegend *l = new TLegend(0.30, 0.75, 0.90, 0.90);
  l -> AddEntry(l_max, "Maximum Likelihood KE : " + max_KE_str + " [MeV]", "l");
  //l -> SetLineColor(kWhite);
  l -> Draw("same");

  c -> SaveAs("./output/" + name + "_KE_Likelihood.pdf");
 
  KE_vec.clear();
  likelihood_vec.clear();
  c -> Close();
}

void ProfileMaker::Produce_KE_vs_Range(TString name, double mass){

  cout << "[ProfileMaker::Produce_KE_vs_Range] Start for " << name << endl;

  double KE_low = 10.;
  double KE_high = 2000.;
  double KE_step = 10.;
  int N_steps = (KE_high - KE_low) / KE_step;
  vector<double> KE_vec, range_vec;
  for(int i = 0; i < N_steps; i++){
    double this_KE = KE_low +(i + 0.) * KE_step;
    double this_range = KE_to_ResLength_BB(this_KE, mass);
    KE_vec.push_back(this_KE);
    range_vec.push_back(this_range);
  }

  TCanvas *c = new TCanvas("", "", 800, 600);
  //c -> SetLogx();
  gStyle->SetOptStat(0);

  double y_max = *max_element(range_vec.begin(), range_vec.end());
  TH1D *template_h = new TH1D("", "", 1., 0., KE_high);
  template_h -> GetXaxis() -> SetTitle("KE [MeV]");
  template_h -> GetYaxis() -> SetTitle("Range [cm]");
  template_h -> GetYaxis() -> SetRangeUser(0., y_max * 1.3);
  template_h -> Draw();

  TGraph *range_gr = new TGraph(N_steps, &KE_vec[0], &range_vec[0]);
  range_gr -> SetName(name + "_KE_vs_Range");
  range_gr -> SetLineColor(kRed);
  range_gr -> SetLineWidth(3);
  range_gr -> Draw("same");

  c -> SaveAs("./output/" + name + "_KE_vs_Range.pdf");

  KE_vec.clear();
  range_vec.clear();
  c -> Close();

}



double ProfileMaker::KE_to_Momentum(double KE, double mass){
  double mom = sqrt(pow(KE, 2) + 2.0 * KE * mass);
  //cout << "[ProfileMaker::KE_to_Momentum] KE\t" << KE << "\tmass\t" << mass << "\tmom\t" << mom << endl;

  return mom;
}

double ProfileMaker::Momentum_to_KE(double P, double mass){
  double KE = sqrt(pow(P, 2) + pow(mass, 2)) - mass;
  return KE;
}

double ProfileMaker::KE_to_ResLength_BB(double KE, double mass){
  // == Res length using Bethe-Bloch (BB)
  double res_length_BB = 0.;
  double KE_BB = KE;
  double step = 0.01; // [cm]
  while(KE_BB > 0.01){
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(KE_BB, mass);
    KE_BB = KE_BB - this_dEdx * step;
    res_length_BB += step;
    
  }
  //cout << "[ProfileMaker::KE_to_ResLength_BB] KE\t" << KE << " [MeV]\t" << "Range\t" << res_length_BB * dEdx.rho << " [g/cm2]\tres_length_BB\t" << res_length_BB << " [cm]" << endl;

  return res_length_BB;
}

double ProfileMaker::ResLength_to_KE_BB(double ResLength, double mass){
  // == KE to ResLength using Bethe-Bloch (BB)
  double KE_BB = 0.1; // [MeV], starting with non-zero energy
  double this_length = 0.;
  double step = 0.01; // [cm]
  bool first = true;
  while(this_length < ResLength){
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(KE_BB, mass);
    KE_BB += this_dEdx * step;
    this_length += step;
    //cout << this_length << ", " << this_dEdx * step << ", " << KE_BB << endl;
  }

  cout << "[ProfileMaker::ResLength_to_KE_BB] ResLength\t" << ResLength  << "\t[cm]\tKE_BB\t" << KE_BB << " [MeV]" << endl;

  return KE_BB;
}

void ProfileMaker::Mimic_GEANT(double KE, double mass){
  double this_xsec = dEdx.Macroscopic_Xsec(KE, mass, 0.4);
  double this_lambda = pitch * this_xsec;
  double CoverI = dEdx.dEdx_Bethe_Bloch(KE, mass) * pitch/ dEdx.I;
  double this_Wmax = dEdx.Get_Wmax(KE, mass);
  //this_Wmax = 0.04;
  int this_N = gRandom->Poisson(this_lambda);

  double dE_BB = dEdx.dEdx_Bethe_Bloch(KE, mass) * pitch;
  cout << "[ProfileMaker::Mimic_GEANT] KE\t" << KE << "\tmass\t" << mass << "\tthis_xsec\t" << this_xsec << "\tthis_lambda\t" << this_lambda
       << "\tthis_N\t" << this_N << "\tdE_BB\t" << dE_BB << "\tdE/dx_BB\t" << dE_BB / pitch << endl;
  for(int i = 0; i < 10000; i++){
    double dE_Urban = 0.;
    for(int j = 0; j < this_N; j++){
      double u_j = gRandom->Uniform(1.);
      double dE_j = dEdx.I / (1 - u_j * ((this_Wmax - dEdx.I)/this_Wmax));
      dE_Urban += dE_j;
    }
    FillHist(Form("Urban_KE%.f", KE), dE_Urban, 1., 100., 0., 10.);
  }
}

void ProfileMaker::KE_to_ResLength_Urban(double KE, double mass){
  // == Res length using the Urban model
  int N_trial = 2000;
  for(int i_Ur = 0; i_Ur < N_trial; i_Ur++){
    double res_length_Ur = 0.;
    double KE_Ur = KE;
    
    while(KE_Ur > 0.1){
      double this_xsec = dEdx.Macroscopic_Xsec(KE_Ur, mass, 0.4);
      double this_lambda = pitch * this_xsec;
      double this_Wmax = dEdx.Get_Wmax(KE_Ur, mass);
      //this_Wmax = 0.04;
      int this_N = gRandom->Poisson(this_lambda);
      double dE_Urban = 0.;
      for(int j = 0; j < this_N; j++){
	double u_j = gRandom->Uniform(1.);
	double dE_j = dEdx.I / (1 - u_j * ((this_Wmax - dEdx.I)/this_Wmax));
	dE_Urban += dE_j;
      }
      KE_Ur = KE_Ur - dE_Urban;
      res_length_Ur += pitch;
      if(KE_Ur < 1.0) break;
    }
    double this_range = res_length_Ur * dEdx.rho;
    TString histname = Form("Range_%.f", KE);
    FillHist(histname, this_range, 1., 5000, 0., 1000.);

    if((i_Ur % (N_trial / 20)) == 0) cout << "[ProfileMaker::KE_to_ResLength_Urban]" << i_Ur << "/" << N_trial << endl;
  }
  

}

void ProfileMaker::KE_to_ResLength_LV(double KE, double mass){
  // == Res length using Landau-Vavilov (LV)
  int N_trial = 2000;
  for(int i_LV = 0; i_LV < N_trial; i_LV++){
    double res_length_LV = 0.;
    double KE_LV = KE;
    
    while(KE_LV > 0.1){
      double this_MPV = dEdx.dEdx_Landau_Vavilov(KE_LV, pitch, mass);
      double this_FWHM = 4 * dEdx.Get_Landau_xi(KE_LV, pitch, mass);
      double mpv_shift = -0.22278298;
      double corrected_mpv = this_MPV - this_FWHM * mpv_shift;
      
      //cout << "[ProfileMaker::KE_to_ResLength] MPV : " << corrected_mpv << ", this_FWHM : " << this_FWHM << endl; 
      TF1 *this_landau = new TF1("landauDistr", "TMath::Landau(x,[0],[1],1)", 0, corrected_mpv + this_FWHM * 50.);
      this_landau -> SetParameter(0, corrected_mpv);
      this_landau -> SetParameter(1, this_FWHM);
      
      double this_dEdx = this_landau->GetRandom();
      delete this_landau;
      //cout << "[ProfileMaker::KE_to_ResLength] KE_LV : " << KE_LV << ", this_dEdx : " << this_dEdx << endl;
      //KE_LV = KE_LV - this_dEdx * pitch;
      KE_LV = KE_LV - this_dEdx;
      res_length_LV += pitch;
      if(KE_LV < 1.0) break;
    }
    double this_range = res_length_LV * dEdx.rho;
    TString histname = Form("Range_%.f", KE);
    FillHist(histname, this_range, 1., 5000, 0., 1000.);
    if((i_LV % (N_trial / 20)) == 0) cout << "[ProfileMaker::KE_to_ResLength LV] KE\t" << KE << " MeV\t" << "Range\t" <<  res_length_LV * dEdx.rho << " [g/cm2]\t" << i_LV << "/" << N_trial << endl;
  }
}

double ProfileMaker::Sum_dEdx_length_KE0(double KE, double mass, double length){
  
  double KE_BB = KE;
  double this_length = 0.;
  double sum_dEdx = 0.;
  while(this_length < length && KE_BB > 0.){
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(KE_BB, mass);
    //cout << "sum_dEdx : " << sum_dEdx << endl;
    sum_dEdx = sum_dEdx + this_dEdx * pitch;
    KE_BB = KE_BB - this_dEdx * pitch;
    this_length += pitch;
  }

  return sum_dEdx;
}

void ProfileMaker::Validate_Landau(double KE, double mass){

  double this_MPV = dEdx.dEdx_Landau_Vavilov(KE, pitch, mass);
  double this_FWHM = 4 * dEdx.Get_Landau_xi(KE, pitch, mass);
  double mpv_shift = -0.22278298;
  double corrected_mpv = this_MPV - this_FWHM * mpv_shift;
  TF1 *landauDistr = new TF1("landauDistr", "TMath::Landau(x,[0],[1],1)", 0, corrected_mpv + this_FWHM * 50.);
  TString graph_name = Form("Landau_%.f", KE);
  landauDistr -> SetName(graph_name);
  landauDistr -> SetParameter(0, corrected_mpv);
  landauDistr -> SetParameter(1, this_FWHM);
  landauDistr -> Write();
  //delete landauDistr;  
  int N_trial = 10000;
  for(int i_LV = 0; i_LV < N_trial; i_LV++){
    double this_dEdx = landauDistr -> GetRandom();
    FillHist("h_" + graph_name, this_dEdx, 1., 100, 0., this_FWHM * 50.);
  }
  cout << "[ProfileMaker::Validate_Landau] KE\t" << KE << "\t" << "dE/dx BB\t" << dEdx.dEdx_Bethe_Bloch(KE, mass) << ", MPV\t" << this_MPV <<"\tFWHM\t" << this_FWHM << endl;
}
/*
void ProfileMaker::Sqrt_s_in_Momentum(TString particle_str, double mass){

  

}
*/
ProfileMaker::ProfileMaker(){

}
ProfileMaker::~ProfileMaker(){

}
