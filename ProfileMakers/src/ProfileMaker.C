#include "ProfileMaker.h"
#include "TF1.h"

const double pitch = 0.65; // [cm]
void ProfileMaker::Execute(){
  cout << "[ProfileMaker::Execute] Start" << endl;

  double mass_muon = 105.658; // [MeV]
  double mass_pion = 139.57; // [MeV]
  double mass_proton = 938.272; // [MeV]

  //KE_to_ResLength_BB(200, mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(200, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(300, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(400, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(800, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(1000, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(2000, mass_muon), mass_muon);
  ResLength_to_KE_BB(KE_to_ResLength_BB(4000, mass_muon), mass_muon);

  ResLength_to_KE_BB(KE_to_ResLength_BB(20, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(50, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(80, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(100, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(200, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(300, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(400, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(800, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(1000, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(2000, mass_pion), mass_pion);
  ResLength_to_KE_BB(KE_to_ResLength_BB(4000, mass_pion), mass_pion);

  KE_to_ResLength_BB(1000., mass_pion);
  KE_to_ResLength_BB(500., mass_pion);
  ResLength_to_KE_BB(20., mass_pion);

  double a = KE_to_Momentum(500., mass_pion);
  a = KE_to_Momentum(400., mass_pion);
  a = KE_to_Momentum(1000., mass_pion);

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

  
}

double ProfileMaker::KE_to_Momentum(double KE, double mass){
  double mom = sqrt(pow(KE, 2) + 2.0 * KE * mass);
  cout << "[ProfileMaker::KE_to_Momentum] KE\t" << KE << "\tmass\t" << mass << "\tmom\t" << mom << endl;

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
  cout << "[ProfileMaker::KE_to_ResLength_BB] KE\t" << KE << " [MeV]\t" << "Range\t" << res_length_BB * dEdx.rho << " [g/cm2]\tres_length_BB\t" << res_length_BB << " [cm]" << endl;

  return res_length_BB;
}

void ProfileMaker::ResLength_to_KE_BB(double ResLength, double mass){
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
