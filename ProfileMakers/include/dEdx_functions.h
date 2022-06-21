#ifndef dEdx_functions_h
#define dEdx_functions_h

#include "TMath.h"

class dEdx_functions{

  // == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
  const double rho = 1.39; // [g/cm3], density of LAr   
  const double K = 0.307075; // [MeV cm2 / mol]
  const double A = 39.948; // [g / mol], atomic mass of Ar
  const double I = 188.0e-6; // [MeV], mean excitation energy
  const double Me = 0.511; // [Mev], mass of electron
  // == Parameters for the density correction
  const double density_C = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a = 0.19559;
  const double density_k = 3.0;

 public :
  double Density_Correction(double beta, double gamma);
  double dEdx_Bethe_Bloch(double KE, double mass);
  double Landau_Vavilov(double KE, double dx, double mass);
  double Get_Landau_xi(double KE, double dx, double mass);

  dEdx_functions();
  ~dEdx_functions();
};

#endif
