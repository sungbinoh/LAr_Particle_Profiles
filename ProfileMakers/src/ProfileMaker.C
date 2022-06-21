#include "ProfileMaker.h"
#include "TF1.h"

const double pitch = 0.65; // [cm]
void ProfileMaker::Execute(){
  cout << "[ProfileMaker::Execute] Start" << endl;

  double mass_muon = 105.658; // [MeV]
  double mass_pion = 139.57; // [MeV]
  double mass_proton = 938.272; // [MeV]

  KE_to_ResLength_BB(200, mass_muon);
  KE_to_ResLength_BB(300, mass_muon);
  KE_to_ResLength_BB(400, mass_muon);
  KE_to_ResLength_BB(800, mass_muon);

  KE_to_ResLength_LV(200, mass_muon);
  KE_to_ResLength_LV(300, mass_muon);
  KE_to_ResLength_LV(400, mass_muon);
  KE_to_ResLength_LV(800, mass_muon);
  
}

void ProfileMaker::KE_to_ResLength_BB(double KE, double mass){
  // == Res length using Bethe-Bloch (BB)
  double res_length_BB = 0.;
  double KE_BB = KE;
  while(KE_BB > 0.){
    double this_dEdx = dEdx.dEdx_Bethe_Bloch(KE_BB, mass);
    KE_BB = KE_BB - this_dEdx * pitch;
    res_length_BB += pitch;
  }
  cout << "[ProfileMaker::KE_to_ResLength_BB] KE\t" << KE << " MeV\t" << "Range\t" << res_length_BB * dEdx.rho << " [g/cm2]" << endl;

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
      TF1 *this_landau = new TF1("landauDistr", "TMath::Landau(x,[0],[1],1)", 0, 500);
      this_landau -> SetParameter(0, corrected_mpv);
      this_landau -> SetParameter(1, this_FWHM);
      
      double this_dEdx = this_landau->GetRandom();
      delete this_landau;
      //cout << "[ProfileMaker::KE_to_ResLength] KE_LV : " << KE_LV << ", this_dEdx : " << this_dEdx << endl;
      KE_LV = KE_LV - this_dEdx * pitch;
      res_length_LV += pitch;
      if(KE_LV < 1.0) break;
    }
    double this_range = res_length_LV * dEdx.rho;
    TString histname = Form("Range_%.f", KE);
    FillHist(histname, this_range, 1., 5000, 0., 1000.);
    if((i_LV % (N_trial / 20)) == 0) cout << "[ProfileMaker::KE_to_ResLength LV] KE\t" << KE << " MeV\t" << "Range\t" <<  res_length_LV * dEdx.rho << " [g/cm2]\t" << i_LV << "/" << N_trial << endl;
  }
}


ProfileMaker::ProfileMaker(){

}
ProfileMaker::~ProfileMaker(){

}
