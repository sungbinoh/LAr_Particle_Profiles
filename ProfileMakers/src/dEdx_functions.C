#include "dEdx_functions.h"

double dEdx_functions::Density_Correction(double beta, double gamma){
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double dEdx_functions::Get_Wmax(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));

  return Wmax;
}

double dEdx_functions::dEdx_Bethe_Bloch(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = Density_Correction(beta, gamma);

  // == dE/dx with the density correction
  double f = rho * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}

double dEdx_functions::Macroscopic_Xsec(double KE, double mass, double r_param = 0.4){
  // == Calculate macroscopic ionization cross section used by the GEANT
  // ==== Nuclear Instruments and Methods in Physics Research A 362 (1995) 416-422
  // ==== Default r_param = 0.4
  double Wmax = Get_Wmax(KE,mass);
  Wmax = 0.04; // == [MeV]
  double dEdx = dEdx_Bethe_Bloch(KE, mass);
  double xsec = dEdx * (Wmax - I) * r_param / (I * Wmax * TMath::Log(Wmax / I));
  double factor = (Wmax - I) / ( Wmax * TMath::Log(Wmax / I));
  double log = TMath::Log(Wmax / I);
  //cout << "[dEdx_functions::Macroscopic_Xsec] Wmax\t" << Wmax << "\tfactor\t" << factor << "\tlog\t" << log << endl;

  return xsec;
}

double dEdx_functions::dEdx_Landau_Vavilov(double KE, double dx, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double delta = Density_Correction(beta, gamma);
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  double a0 = 2.0 * Me * pow(beta * gamma, 2) / I;
  //double this_dpdx = (xi / dx) * (TMath::Log(a0) + + TMath::Log(xi / I) + 0.2 - pow(beta, 2) - delta) ;
  double this_dpdx = (xi) * (TMath::Log(a0) + + TMath::Log(xi / I) + 0.2 - pow(beta, 2) - delta) ; // == MPV of dE, FWHM = 4xi

  return this_dpdx;
}

double dEdx_functions::Get_Landau_xi(double KE, double dx, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  return xi;
}

dEdx_functions::dEdx_functions(){

}

dEdx_functions::~dEdx_functions(){

}
