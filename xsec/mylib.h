#ifndef mylib_h
#define mylib_h

ofstream outfile;

using namespace std;
const double alpha = 1 - 0.6827;
// == Debugging Mode
bool debug = true;
bool blind_SR = false;
TString Simulator;
// == Call all needed maps
map<TString, TH1D*> maphist;
map<TString, TH1F*> mapTH1F;
map<TString, TH2D*> maphist2D;
map<TString, TGraph*> map_gr;
map<TString, TGraphAsymmErrors*> map_asym_gr;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, double> map_overflow;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, std::vector<double> > map_bin_vector;
map<TString, std::vector<TString> > map_sample_names;
map<TString, std::vector<double> > map_syst_array;
map<TString, std::vector<double> > map_syst_table;

double mass_muon = 105.658; // [MeV]
double mass_pion = 139.57; // [MeV]
double mass_proton = 938.272; // [MeV]

const double K = 0.307075; // [MeV cm2 / mol]
const double I = 188.0e-6; // [MeV], mean excitation energy
const double Me = 0.511; // [Mev], mass of electron
// == Parameters for the density correction
const double density_C = 5.2146;
const double density_y0 = 0.2;
const double density_y1 = 3.0;
const double density_a = 0.19559;
const double density_k = 3.0;

const double LAr_density = 1.39; // [g/cm3]
const double N_A = 6.02; // [10^23 / mol]
const double M_Ar = 39.948; // [g / mol]
const double xsec_unit = M_Ar / (LAr_density * N_A); // [10^4 mb cm]

const int N_pi_type = 10;
TString pi_type_str[N_pi_type] = {"Fake Data",
				  "#pi^{+}_{Inel.}",
				  "#pi^{+}_{Elas.}",
				  "#mu",
				  "Misid. : cosmic",
				  "Misid. : p",
				  "Misid. : #pi^{+}",
				  "Misid. : #mu",
				  "Misid. : e/#gamma",
				  "Misid. : other"};
const int N_p_type = 9;
TString p_type_str[N_p_type] = {"Fake Data",
				"PInel",
				"PElas",
				"misID:cosmic",
				"misID:p",
				"misID:pi",
				"misID:mu",
				"misID:e/#gamma",
				"misID:other"};

// == Uniform 100 MeV 
//const int N_KE_bins = 15;
//double KE_binning[N_KE_bins] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1500.};

// == Uniiform 50 MeV
const int N_KE_bins = 24;
double KE_binning[N_KE_bins] = {0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.,
				1050., 1100., 1500.};

//const int N_KE_bins = 15;
//double KE_binning[N_KE_bins] = {0., 400., 500., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1100., 1200., 1300.}; 

//const int N_KE_bins = 17;
//double KE_binning[N_KE_bins] = {0., 100., 200., 300., 400., 500., 600., 650., 700., 750., 800., 850., 900., 950., 1000., 1200., 1300.};

const int N_MC = 6;
TString MC_category[N_MC] = {"NC", "NumuCC", "External_NC", "External_NueCC", "External_NumuCC", "NueCC"};

const int N_smear_bit = 8;
TString smear_flags[N_smear_bit] = {"NONE", "P", "Theta", "P_Theta", "Phi", "P_Phi", "Phi_Theta", "All"};

vector<double> vx, vy, vexl, vexh, veyl, veyh;

double dEdx_Bethe_Bloch(double KE, double mass);

TH1D * GetHist(TString hname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;

  return h;

}

void Rebin_with_overflow(TString histname, int N_bin, double binx[]){
  if(debug) cout << "[Rebin_with_overflow] " << histname << endl;
  maphist[histname + "rebin"] = (TH1D*)maphist[histname] -> Rebin(N_bin - 1, histname + "rebin", binx);
  double last_bin = 0.;
  last_bin =  maphist[histname + "rebin"] -> GetBinContent(N_bin - 1) + maphist[histname + "rebin"] -> GetBinContent(N_bin);
  maphist[histname + "rebin"] -> SetBinContent(N_bin - 1, last_bin);
}

TH1D* Make_incident_histogram(TString h_init_name, TString h_end_name){

  cout << "[Make_incident_histogram] Making incident histogram using " << h_init_name << " and " << h_end_name << endl;
  int N_bin = maphist[h_init_name] -> GetNbinsX();
  TH1D *out = (TH1D*) maphist[h_init_name] -> Clone();
  for(int i_bin = 1; i_bin < N_bin + 1; i_bin++){
    double sum_N_end = 0.;
    double sum_N_init = 0.;
    for(int j_bin = 1; j_bin < i_bin + 1; j_bin++){
      //cout << "[sum_N_end] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_end += maphist[h_end_name] -> GetBinContent(j_bin);
    }
    for(int j_bin = 1; j_bin < i_bin + 0; j_bin++){
      //cout << "[sum_N_init] (i, j) : (" << i << ", " << j << ")" << endl;
      sum_N_init += maphist[h_init_name] -> GetBinContent(j_bin);
    }
    double this_content = sum_N_end - sum_N_init;
    double this_err = sqrt(this_content);
    out -> SetBinContent(i_bin, this_content);
    out -> SetBinError(i_bin, this_err);
  }

  return out;
}

TH1D* Make_cross_section_histogram(TString h_inc_name, TString h_int_name){

  cout << "[Make_cross_section_histogram] for " << h_inc_name << " and " << h_int_name << endl;
  int N_bin = maphist[h_inc_name] -> GetNbinsX();
  TH1D * out = (TH1D*)  maphist[h_inc_name] -> Clone();
  for(int i_bin = 1; i_bin < N_bin + 1; i_bin++){
    double this_inc = maphist[h_inc_name] -> GetBinContent(i_bin);
    double this_int = maphist[h_int_name] -> GetBinContent(i_bin);
    double ratio = 1.;
    double ratio_err = 0.;
    if(this_inc > 0. && this_int > 0. && this_inc != this_int){
      ratio = this_inc / (this_inc - this_int);
      double numer_err = sqrt(this_inc);
      double denom_err = sqrt(this_inc + this_int);
      ratio_err = ratio * fabs(1. / numer_err + 1. / denom_err);
    }
    double log_ratio = log(ratio);
    double log_err = log(1. + ratio_err / ratio );
    double this_KE = out -> GetXaxis() -> GetBinCenter(i_bin);
    double this_dE = out -> GetXaxis() -> GetBinWidth(i_bin);
    double this_dEdx = dEdx_Bethe_Bloch(this_KE, mass_pion);
    double this_xsec = 10000. * xsec_unit * log_ratio * this_dEdx / this_dE;
    double this_xsec_err = 10000. * xsec_unit * log_err * this_dEdx / this_dE;
    out -> SetBinContent(i_bin, this_xsec);
    out -> SetBinError(i_bin, this_xsec_err);
    cout << "[Make_cross_section_histogram] " << i_bin << " xsec : " << this_xsec << " +- " << this_xsec_err << ", this_inc : " << this_inc << ", this_int : " << this_int << endl;
  }

  return out;
}

void change_to_pseudo_data(TString current_histname){
  
  int N_bin = maphist[current_histname] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    int current_bin_int =  maphist[current_histname] -> GetBinContent(i);
    double current_error = pow(current_bin_int, 0.5);
    maphist[current_histname] -> SetBinContent(i, current_bin_int);
    maphist[current_histname] -> SetBinError(i, current_error);
  }
}

void Proper_error_data(TString nameofhistogram, int N_bin, double binx[]){
  
  map_asym_gr[nameofhistogram + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + "rebin"));
  
  for(int i = 0; i < N_bin; i++){
    int N = GetHist(nameofhistogram + "rebin") -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYlow(i, (N-L) ); // / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXlow(i, 0);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYhigh(i, (U-N) );// / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[nameofhistogram + "correct_error"] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y; // ( binx[i + 1] - binx[i] );
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPoint(i, current_x, modified_y);

      
      //if(debug) cout << "i : " << i << ", current_x : " << current_x << ", current_y : " << current_y << ", modified_y : " << modified_y << endl;

      veyl.push_back( (N-L)); // / (binx[i + 1] - binx[i]) );
      veyh.push_back( (U-N)); // / (binx[i + 1] - binx[i]));
    }
    else{
      double zerodata_err_low = 0.1 ;// / (binx[i + 1] - binx[i]);
      double zerodata_err_high = 1.8 ;// / (binx[i + 1] - binx[i]);

      veyl.push_back(zerodata_err_low);
      veyh.push_back(zerodata_err_high);
      
      double current_x = GetHist(nameofhistogram + "rebin") -> GetBinCenter(i + 1);
      double current_ex = GetHist(nameofhistogram + "rebin") -> GetBinWidth(i + 1) /2.;

      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYlow(i, zerodata_err_low);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEYhigh(i, zerodata_err_high);
      map_asym_gr[nameofhistogram + "correct_error"] -> SetPointEXhigh(i, 0.);

      vx.push_back(current_x);
      vexl.push_back(current_ex);
      vexh.push_back(current_ex);
    }
  }//end for i on N_bin
    
}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, bool YErrorZero=false){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins();
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=0; i<Nbins; i++){
    x[i] = xaxis->GetBinCenter(i+1);
    y[i] = hist->GetBinContent(i+1);
    xlow[i] = xaxis->GetBinCenter(i+1)-xaxis->GetBinLowEdge(i+1);
    xup[i] = xaxis->GetBinUpEdge(i+1)-xaxis->GetBinCenter(i+1);
    ylow[i] = hist->GetBinError(i+1);
    yup[i] = hist->GetBinError(i+1);
    if(YErrorZero){
      ylow[i] = 0;
      yup[i] = 0;
    }
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left, int n_x_shift, int i_x_shift){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);

    double dx = (xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left))/double(n_x_shift+1);
    x[i-1] = xaxis->GetBinLowEdge(i+n_skip_x_left) + double(i_x_shift+1) * dx;
    xlow[i-1] = double(i_x_shift+1) * dx;
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-x[i-1];
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}


TGraphAsymmErrors* GraphSubtract(TGraphAsymmErrors *a, TGraphAsymmErrors *b, bool Rel){

  //==== do a-b

  int NX = a->GetN();
  TGraphAsymmErrors* gr_out = (TGraphAsymmErrors*)a->Clone();

  for(int i=0; i<NX; i++){

    double a_x, a_y, b_x, b_y;

    a->GetPoint(i, a_x, a_y);
    b->GetPoint(i, b_x, b_y);

    if(Rel==true){
      gr_out->SetPoint( i, a_x, (a_y-b_y)/b_y );
    }
    else{
      gr_out->SetPoint( i, a_x, a_y-b_y );
    }

  }

  return gr_out;

}

void RemoveLargeError(TGraphAsymmErrors *a){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if(y+yerr_high==1.){
      a->SetPointEYhigh( i, yerr_low );
    }

  }

}

void ScaleGraph(TGraphAsymmErrors *a, double c){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    a->SetPoint(i, x, c*y);
    a->SetPointEYlow(i, c*yerr_low);
    a->SetPointEYhigh(i, c*yerr_high);

  }

}



double GetMaximum(TH1D* hist){

  TAxis *xaxis = hist->GetXaxis();

  double maxval(-1.);
  for(int i=1; i<=xaxis->GetNbins(); i++){
    if( hist->GetBinContent(i) + hist->GetBinError(i) > maxval ){
      maxval = hist->GetBinContent(i) + hist->GetBinError(i);
    }
  }

  return maxval;

}

double GetMaximum(TGraphAsymmErrors *a){

  int NX = a->GetN();

  double maxval(-9999.);
  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if( y+yerr_high > maxval ){
      maxval = y+yerr_high;
    }

  }

  return maxval;

}

double GetYieldSystematics(TH1D *hist){

  int n_syst = hist->GetXaxis()->GetNbins();
  int n_source = (n_syst-1)/2;

  //==== Bin index
  //==== i=1 : central
  //==== i=2 : _MuonEn_up
  //==== i=3 : _MuonEn_down
  //==== -> n_syst = 3
  //==== -> n_source = (n_syst-1)/2 = (3-1)/2 = 1

  double y_central = hist->GetBinContent(1);

  double sum_syst = 0.;
  for(int i=1; i<=n_source; i++){
    double yield_up = hist->GetBinContent(i*2);
    double yield_down = hist->GetBinContent(i*2+1);

    double syst_up = fabs(yield_up-y_central);
    double syst_down = fabs(yield_down-y_central);

    sum_syst += 0.5*(syst_up*syst_up+syst_down*syst_down);

  }

  return sqrt(sum_syst);

}

TDirectory *MakeTemporaryDirectory(){

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}

void AddPhantomZero(double a, TString align, int digit_int, int digit_frac){

  if(align=="r"){
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    for(int i=0; i<digit_int-(number_maxdigit+1); i++) cout << "\\phantom{0}";
    cout << std::fixed<<std::setprecision(digit_frac) << a;
  }

  else if(align=="l"){
    int target_total_digit = digit_int+digit_frac;
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    cout << std::fixed<<std::setprecision(digit_frac) << a;
    for(int i=0; i<target_total_digit-(number_maxdigit+1)-digit_frac; i++) cout << "\\phantom{0}";
  }

}

// === dE/dx
double Density_Correction(double beta, double gamma){
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

double Get_Wmax(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));

  return Wmax;
}

double dEdx_Bethe_Bloch(double KE, double mass){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * Me * pow(beta * gamma, 2)) / (1.0 + 2.0 * Me * (gamma / mass) + pow((Me / mass),2));
  double delta = Density_Correction(beta, gamma);

  // == dE/dx with the density correction
  double f = LAr_density * K * (18.0 / M_Ar) * pow(1. / beta, 2);
  double a0 = 0.5 * TMath::Log(2.0 * Me * pow(beta * gamma, 2) * Wmax / (I * I));
  double this_dEdx = f * ( a0 - pow(beta, 2) - delta / 2.0); // [MeV/cm]

  return this_dEdx;
}

#endif
