#include "dEdx_functions.h"

ROOT::Math::VavilovAccurate vav;

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

double dEdx_functions::Landau_dE_cutoff(double KE, double width, double mass){

  double this_cutoff = 100.;
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = Get_Wmax(KE,mass);
  double this_Landau_xi = Get_Landau_xi(KE, width, mass);

  double this_kappa = this_Landau_xi / Wmax;
  if(this_kappa > 0.01) return this_cutoff;

  double lambda_bar = -1. * gamma_prime - beta * beta - TMath::Log(this_Landau_xi / Wmax);
  double lambda_max = 0.51146 + 1.19486 * lambda_bar + (0.465814 + 0.0115374 * lambda_bar) * exp(1.17165 + 0.979242 * lambda_bar);
  //cout << "[dEdx_functions::Landau_dE_cutoff] gamma_prime : " << gamma_prime << ", lambda_bar : " << lambda_bar << ", lambda_max : " << lambda_max << endl;
  
  TF1 *landau_f = new TF1("landau_f", "TMath::Landau(x)", -20., 5000.);
  double repro_lambda_bar = landau_f -> Mean(-20., lambda_max);
  //cout << "[dEdx_functions::Landau_dE_cutoff] repro_lambda_bar : " << repro_lambda_bar << endl;
  //cout << "[dEdx_functions::Landau_dE_cutoff] lambda_bar - repro_lambda_bar : " << lambda_bar - repro_lambda_bar << endl;

  this_cutoff = (lambda_max + gamma_prime + beta*beta + TMath::Log(this_Landau_xi / Wmax)) * (this_Landau_xi / width) + dEdx_Bethe_Bloch(KE, mass);

  return this_cutoff;
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
  //double this_dpdx = (xi / dx) * (TMath::Log(a0) + TMath::Log(xi / I) + 0.2 - pow(beta, 2) - delta) ;
  double this_dpdx = (xi / dx) * (TMath::Log(a0) + TMath::Log(xi / I) + 0.2 - pow(beta, 2) - delta) ; // == MPV of dE, FWHM = 4xi
  
  return this_dpdx;
}

double dEdx_functions::Get_Landau_xi(double KE, double dx, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * dx * 0.5 * K * (18.0 / A) * pow(1. / beta, 2);
  return xi;
}

double dEdx_PDF_setting(double *x, double *par){

  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width}
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =  this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =  TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav =  vav.Pdf(y, par[0], par[1]);
    this_vav =  this_vav / a;
  }
  
  return this_vav;
}

TF1 *dEdx_functions::dEdx_PDF(double *par){

  TF1 *out = new TF1("", dEdx_PDF_setting, -100., 1000., 5);
  //cout << "[dEdx_functions::dEdx_PDF] par[0] : " << par[0] << ", par[1] : " << par[1] << ", par[2] : " << par[2] << endl;
  out -> SetParameter(0, par[0]);
  out -> SetParameter(1, par[1]);
  out -> SetParameter(2, par[2]);
  out -> SetParameter(3, par[3]);
  out -> SetParameter(4, par[4]);
  out -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  return out;
}
/*
double dEdx_functions::Get_Kappa(double KE, double dx, double mass){
  

}
*/
dEdx_functions::dEdx_functions(){

}

dEdx_functions::~dEdx_functions(){

}
