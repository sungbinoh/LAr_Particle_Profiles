#ifndef ProfileMaker_h
#define ProfileMaker_h

#include "ProfileMakerCore.h"

class ProfileMaker : public ProfileMakerCore {

public:
  
  void Execute();
  //void Test();
  void Produce_Profile(TString name, double mass);
  void Produce_KE_vs_dEdx(TString name, double mass);
  void Produce_kappa(TString name, double mass, double width);
  void Produce_dEdx_PDF(TString name, double mass, double KE, double width, double xmin, double xmax, double ymax);
  void Produce_dEdx_likelihood(TString name, double mass, double dEdx, double width, double xmin, double xmax, double ymax);
  void Produce_Range_from_Momentum_Gaussian(TString name, double mass, double mean, double sigma);
  double KE_to_Momentum(double KE, double mass);
  double Momentum_to_KE(double P, double mass);
  double KE_to_ResLength_BB(double KE, double mass);
  double ResLength_to_KE_BB(double ResLength, double mass);
  void Mimic_GEANT(double KE, double mass);
  void KE_to_ResLength_Urban(double KE, double mass);
  void KE_to_ResLength_LV(double KE, double mass);
  double Sum_dEdx_length_KE0(double KE, double mass, double length);
  void Validate_Landau(double KE, double mass);
  ProfileMaker();
  virtual ~ProfileMaker();
};

#endif
