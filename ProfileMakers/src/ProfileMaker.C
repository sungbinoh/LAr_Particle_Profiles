#include "ProfileMaker.h"
#include "TF1.h"

const double pitch = 0.65; // [cm]
void ProfileMaker::Execute(){
  cout << "[ProfileMaker::Execute] Start" << endl;

  double mass_muon = 105.658; // [MeV]
  double mass_pion = 139.57; // [MeV]
  double mass_proton = 938.272; // [MeV]

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

  // == KE to dE/dx BB
  cout << "<dE/dx> pion KE = 20 MeV (P = " << KE_to_Momentum(20, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(20, mass_pion) << endl;
  cout << "<dE/dx> pion KE = 200 MeV (P = " << KE_to_Momentum(200, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(200, mass_pion) << endl;
  cout << "<dE/dx> pion KE = 500 MeV (P = " << KE_to_Momentum(500, mass_pion) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(500, mass_pion) << endl;
  cout << "<dE/dx> proton KE = 20 MeV (P = " << KE_to_Momentum(20, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(20, mass_proton) << endl;
  cout << "<dE/dx> proton KE = 200 MeV (P = " << KE_to_Momentum(200, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(200, mass_proton) << endl;
  cout << "<dE/dx> proton KE = 500 MeV (P = " << KE_to_Momentum(500, mass_proton) << " MeV/c) : " << dEdx.dEdx_Bethe_Bloch(500, mass_proton) << endl;
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

  Validate_Landau(200, mass_muon);
  Validate_Landau(50, mass_muon);
  Validate_Landau(10, mass_muon);

  Mimic_GEANT(200, mass_muon);
  Mimic_GEANT(1000, mass_muon);

  ///////////////////////////
  // == Produce profiles
  ///////////////////////////
  //Produce_Profile("pion", mass_pion);
  //Produce_Profile("proton", mass_proton);
  //Produce_kappa("pion_0p65cm", mass_pion, 0.65);
  //Produce_kappa("proton_0p65cm", mass_proton, 0.65);
  Produce_dEdx_PDF("dEdx_PDF_pion_500MeV_0p65cm", mass_pion, 500, 0.65, 0., 6., 0.22);
  Produce_dEdx_PDF("dEdx_PDF_proton_500MeV_0p65cm", mass_proton, 500, 0.65, 0., 6., 0.25);
  Produce_dEdx_PDF("dEdx_PDF_pion_200MeV_0p65cm", mass_pion, 200, 0.65, 0., 6., 0.22);
  Produce_dEdx_PDF("dEdx_PDF_proton_200MeV_0p65cm", mass_proton, 200, 0.65, 2., 8., 0.4);
  Produce_dEdx_PDF("dEdx_PDF_pion_20MeV_0p65cm", mass_pion, 20, 0.65, 4., 8., 0.5);
  Produce_dEdx_PDF("dEdx_PDF_proton_20MeV_0p65cm", mass_proton, 20, 0.65, 23., 28., 1.2);

}

void ProfileMaker::Produce_Profile(TString name, double mass){

  vector<double> KE_vec;
  vector<double> KE_to_mom_vec;
  vector<double> KE_to_range_BB_vec;
  vector<double> KE_to_dEdx_BB_vec;
  for(int i = 0; i < 1000; i++){
    double this_KE = 3.+ 2. * (i + 0.);
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
  c -> SetLogx();
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

  dEdx_range -> Write();

  KE_vec.clear();
  KE_to_mom_vec.clear();
  KE_to_range_BB_vec.clear();
  KE_to_dEdx_BB_vec.clear();
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

  TCanvas *c = new TCanvas("", "", 800, 600);
  gStyle->SetOptStat(0);
  //c -> SetLogx();
  TH1D *template_h = new TH1D("", "", 1., xmin, xmax);
  template_h ->GetXaxis() -> SetTitle("dE/dx [MeV/cm]");
  template_h ->GetYaxis() -> SetTitle("A.U.");
  template_h -> GetYaxis() -> SetRangeUser(0., ymax);
  template_h -> Draw();

  this_PDF -> Draw("lsame");

  c -> SaveAs("./output/" + name + "_dEdx_PDF.pdf");

}

double ProfileMaker::KE_to_Momentum(double KE, double mass){
  double mom = sqrt(pow(KE, 2) + 2.0 * KE * mass);
  //cout << "[ProfileMaker::KE_to_Momentum] KE\t" << KE << "\tmass\t" << mass << "\tmom\t" << mom << endl;

  return mom;
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

ProfileMaker::ProfileMaker(){

}
ProfileMaker::~ProfileMaker(){

}
