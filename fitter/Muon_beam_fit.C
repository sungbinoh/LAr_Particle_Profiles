#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;
 
TH1 *makeTH1();
TTree *makeTTree();

void Fit_test(){

  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_in = new TFile(root_file_path + "Template_TRIG2.root");
  TH1D * data = (TH1D*)gDirectory->Get("data");
  TH1D * pion_fit = (TH1D*)gDirectory->Get("Pion_fit");
  TH1D * expected_DiF_muon_P_2p5 = (TH1D*)gDirectory->Get("Theta_2p5_basis");
  TH1D * expected_DiF_muon_P_2p4 = (TH1D*)gDirectory->Get("Theta_2p4_basis");
  TH1D * expected_DiF_muon_P_2p2 = (TH1D*)gDirectory->Get("Theta_2p2_basis");
  TH1D * expected_DiF_muon_P_1p8 = (TH1D*)gDirectory->Get("Theta_1p8_basis");
  TH1D * expected_DiF_muon_P_1p2 = (TH1D*)gDirectory->Get("Theta_1p2_basis");
  TH1D * expected_DiF_muon_P_0p5 = (TH1D*)gDirectory->Get("Theta_0p5_basis");
  TH1D * expected_DiF_muon_P_1p8_right = (TH1D*)gDirectory->Get("Theta_1p8_basis_right");
  TH1D * expected_DiF_muon_P_1p2_right = (TH1D*)gDirectory->Get("Theta_1p2_basis_right");
  TH1D * expected_DiF_muon_P_0p5_right = (TH1D*)gDirectory->Get("Theta_0p5_basis_right");
  TH1D * hist_test = (TH1D*) expected_DiF_muon_P_1p2 -> Clone();
  hist_test -> Add(expected_DiF_muon_P_1p8_right);
  hist_test -> Add(expected_DiF_muon_P_2p4);

  RooRealVar x("x", "x", 500., 1200.);
  RooDataHist temp_data("data", "data", x, Import(*data));
  RooDataHist temp_pion_fit("temp_pion_fit", "temp_pion_fit", x, Import(*pion_fit));
  RooDataHist temp_P_2p5("P_2p5", "P_2p5", x, Import(*expected_DiF_muon_P_2p5));
  RooDataHist temp_P_2p4("P_2p4", "P_2p4", x, Import(*expected_DiF_muon_P_2p4));
  RooDataHist temp_P_2p2("P_2p2", "P_2p2", x, Import(*expected_DiF_muon_P_2p2));
  RooDataHist temp_P_1p8("P_1p8", "P_1p8", x, Import(*expected_DiF_muon_P_1p8));
  RooDataHist temp_P_1p2("P_1p2", "P_1p2", x, Import(*expected_DiF_muon_P_1p2));
  RooDataHist temp_P_0p5("P_0p5", "P_0p5", x, Import(*expected_DiF_muon_P_0p5));
  RooDataHist temp_P_1p8_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));
  RooDataHist temp_P_1p2_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));
  RooDataHist temp_P_0p5_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));
  RooDataHist temp_test("temp_test", "temp_test", x, Import(*hist_test));
  
  RooHistPdf pdf_data("pdf_data", "pdf_data", x, temp_data, 2);
  RooHistPdf pdf_pion_fit("pdf_pion_fit", "pdf_pion_fit", x, temp_pion_fit, 2);
  RooHistPdf pdf_P_2p5("pdf_P_2p5", "pdf_P_2p5", x, temp_P_2p5, 2);
  RooHistPdf pdf_P_2p4("pdf_P_2p4", "pdf_P_2p4", x, temp_P_2p4, 2);
  RooHistPdf pdf_P_2p2("pdf_P_2p2", "pdf_P_2p2", x, temp_P_2p2, 2);
  RooHistPdf pdf_P_1p8("pdf_P_1p8", "pdf_P_1p8", x, temp_P_1p8, 2);
  RooHistPdf pdf_P_1p2("pdf_P_1p2", "pdf_P_1p2", x, temp_P_1p2, 2);
  RooHistPdf pdf_P_0p5("pdf_P_0p5", "pdf_P_0p5", x, temp_P_0p5, 2);
  RooHistPdf pdf_P_1p8_right("pdf_P_1p8_right", "pdf_P_1p8_right", x, temp_P_1p8_right, 2);
  RooHistPdf pdf_P_1p2_right("pdf_P_1p2_right", "pdf_P_1p2_right", x, temp_P_1p2_right, 2);
  RooHistPdf pdf_P_0p5_right("pdf_P_0p5_right", "pdf_P_0p5_right", x, temp_P_0p5_right, 2);

  double nominal = 0.05;
  double max = 1.;

  RooRealVar frac_2p5("frac_2p5", "frac_2p5", 0.5, 0., max);
  RooRealVar frac_2p4("frac_2p4", "frac_2p4", nominal, 0., max);
  RooRealVar frac_2p2("frac_2p2", "frac_2p2", nominal, 0., max);
  RooRealVar frac_1p8("frac_1p8", "frac_1p8", nominal, 0., max);
  RooRealVar frac_1p2("frac_1p2", "frac_1p2", nominal, 0., max);
  RooRealVar frac_0p5("frac_0p5", "frac_0p5", nominal, 0., max);
  RooRealVar frac_1p8_right("frac_1p8_right", "frac_1p8_right", 0.5, 0., max);
  RooRealVar frac_1p2_right("frac_1p2_right", "frac_1p2_right", 0.5, 0., max);
  RooRealVar frac_0p5_right("frac_0p5_right", "frac_0p5_right", 0.5, 0., max);

  RooPlot *xframe = x.frame(Title("Example of composite pdf=(sig1+sig2)+bkg"));

  //RooAddPdf model_test("model_test", "model_test", RooArgList(pdf_P_2p5, pdf_P_2p4), RooArgList(frac_2p5), true);
  //RooAddPdf model_test("model_test", "model_test", RooArgList(pdf_P_2p5, pdf_P_2p4, pdf_P_2p4, pdf_P_1p2), RooArgList(frac_2p5, frac_2p4, frac_2p2), true);
  //RooAddPdf model_test("model", "sum_all", RooArgList(pdf_P_2p5, pdf_P_2p4, pdf_P_2p2, pdf_P_1p8, pdf_P_1p2, pdf_P_0p5, pdf_P_1p8_right, pdf_P_1p2_right, pdf_P_0p5_right, pdf_pion_fit),
  //		       RooArgList(frac_2p5, frac_2p4, frac_2p2, frac_1p8, frac_1p2, frac_0p5, frac_1p8_right, frac_1p2_right), true);
  RooAddPdf model_test("model", "sum_all", RooArgList(pdf_P_2p5, pdf_P_2p4, pdf_P_2p2, pdf_P_1p8, pdf_P_1p2, pdf_P_0p5, pdf_P_1p8_right, pdf_P_1p2_right,  pdf_pion_fit),
                       RooArgList(frac_2p5, frac_2p4, frac_2p2, frac_1p8, frac_1p2, frac_0p5, frac_1p8_right, frac_1p2_right), true);
  //RooDataHist *data_test = model_test.generate(x, 100000);
  model_test.fitTo(temp_data);
  //model_test.fitTo(temp_P_2p2);
  //temp_P_2p2.plotOn(xframe);
  //data_test->plotOn(xframe);
  temp_data.plotOn(xframe, LineColor(kBlack), LineStyle(kSolid));
  model_test.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  model_test.plotOn(xframe, Components(pdf_P_2p5), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_2p4), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_2p2), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_1p8), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_1p2), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_0p5), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_1p8_right), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_1p2_right), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_P_0p5_right), LineStyle(kDashed));
  model_test.plotOn(xframe, Components(pdf_pion_fit), LineStyle(kDashed));

  TCanvas *c = new TCanvas("rf706_histpdf", "rf706_histpdf", 800, 400);
  xframe->Draw();
  TString WORKING_DIR = getenv("LArProf_WD");
  TString pdfname;
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/Template_fit_Muon_1GeV.pdf";
  c -> SaveAs(pdfname);


  f_in -> Close();
}

void Muon_beam_fit(){

  Fit_test();
  /*
  TString input_file_dir = getenv("LArProf_WD");
  TString root_file_path =input_file_dir + "/output/root/";
  TFile *f_in = new TFile(root_file_path + "Template_TRIG2.root");
  TH1D * data = (TH1D*)gDirectory->Get("data");
  TH1D * expected_DiF_muon_P_2p5 = (TH1D*)gDirectory->Get("Theta_2p5_basis");
  TH1D * expected_DiF_muon_P_2p4 = (TH1D*)gDirectory->Get("Theta_2p4_basis");
  TH1D * expected_DiF_muon_P_2p2 = (TH1D*)gDirectory->Get("Theta_2p2_basis");
  TH1D * expected_DiF_muon_P_1p8 = (TH1D*)gDirectory->Get("Theta_1p8_basis");
  TH1D * expected_DiF_muon_P_1p2 = (TH1D*)gDirectory->Get("Theta_1p2_basis");
  TH1D * expected_DiF_muon_P_0p5 = (TH1D*)gDirectory->Get("Theta_0p5_basis");
  TH1D * expected_DiF_muon_P_1p8_right = (TH1D*)gDirectory->Get("Theta_1p8_basis_right");
  TH1D * expected_DiF_muon_P_1p2_right = (TH1D*)gDirectory->Get("Theta_1p2_basis_right");
  TH1D * expected_DiF_muon_P_0p5_right = (TH1D*)gDirectory->Get("Theta_0p5_basis_right");
  
  RooRealVar x("x", "x", 0., 1500.);
  RooDataHist temp_data("data", "data", x, Import(*data));
  RooDataHist temp_P_2p5("P_2p5", "P_2p5", x, Import(*expected_DiF_muon_P_2p5));
  RooDataHist temp_P_2p4("P_2p4", "P_2p4", x, Import(*expected_DiF_muon_P_2p4));
  RooDataHist temp_P_2p2("P_2p2", "P_2p2", x, Import(*expected_DiF_muon_P_2p2));
  RooDataHist temp_P_1p8("P_1p8", "P_1p8", x, Import(*expected_DiF_muon_P_1p8));
  RooDataHist temp_P_1p2("P_1p2", "P_1p2", x, Import(*expected_DiF_muon_P_1p2));
  RooDataHist temp_P_0p5("P_0p5", "P_0p5", x, Import(*expected_DiF_muon_P_0p5));
  RooDataHist temp_P_1p8_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));
  RooDataHist temp_P_1p2_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));
  RooDataHist temp_P_0p5_right("P_1p8_right", "P_1p8_right", x, Import(*expected_DiF_muon_P_1p8_right));

  RooHistPdf pdf_data("pdf_data", "pdf_data", x, temp_data, 2);
  RooHistPdf pdf_P_2p5("pdf_P_2p5", "pdf_P_2p5", x, temp_P_2p5, 2);
  RooHistPdf pdf_P_2p4("pdf_P_2p4", "pdf_P_2p4", x, temp_P_2p4, 2);
  RooHistPdf pdf_P_2p2("pdf_P_2p2", "pdf_P_2p2", x, temp_P_2p2, 2);
  RooHistPdf pdf_P_1p8("pdf_P_1p8", "pdf_P_1p8", x, temp_P_1p8, 2);
  RooHistPdf pdf_P_1p2("pdf_P_1p2", "pdf_P_1p2", x, temp_P_1p2, 2);
  RooHistPdf pdf_P_0p5("pdf_P_0p5", "pdf_P_0p5", x, temp_P_0p5, 2);
  RooHistPdf pdf_P_1p8_right("pdf_P_1p8_right", "pdf_P_1p8_right", x, temp_P_1p8_right, 2);
  RooHistPdf pdf_P_1p2_right("pdf_P_1p2_right", "pdf_P_1p2_right", x, temp_P_1p2_right, 2);
  RooHistPdf pdf_P_0p5_right("pdf_P_0p5_right", "pdf_P_0p5_right", x, temp_P_0p5_right, 2);

  double nominal = 1.;
  double max = 1.;
  
  RooRealVar frac_2p5("frac_2p5", "frac_2p5", nominal, 0., max);
  RooRealVar frac_2p4("frac_2p4", "frac_2p4", nominal, 0., max);
  RooRealVar frac_2p2("frac_2p2", "frac_2p2", nominal, 0., max);
  RooRealVar frac_1p8("frac_1p8", "frac_1p8", nominal, 0., max);
  RooRealVar frac_1p2("frac_1p2", "frac_1p2", nominal, 0., max);
  RooRealVar frac_0p5("frac_0p5", "frac_0p5", nominal, 0., max);
  RooRealVar frac_1p8_right("frac_1p8_right", "frac_1p8_right", nominal, 0., max);
  RooRealVar frac_1p2_right("frac_1p2_right", "frac_1p2_right", nominal, 0., max);
  RooRealVar frac_0p5_right("frac_0p5_right", "frac_0p5_right", nominal, 0., max);

  RooPlot *xframe = x.frame(Title("Example of composite pdf=(sig1+sig2)+bkg"));
  pdf_P_2p5.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_2p4.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_2p2.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_1p8.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_1p2.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_0p5.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_1p8_right.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_1p2_right.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  pdf_P_0p5_right.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));

  temp_P_2p5.plotOn(xframe, LineColor(kBlue), LineStyle(kSolid));

  RooAddPdf model("model", "sum_all", RooArgList(pdf_P_2p5, pdf_P_2p4, pdf_P_2p2, pdf_P_1p8, pdf_P_1p2, pdf_P_0p5, pdf_P_1p8_right, pdf_P_1p2_right, pdf_P_0p5_right),
		  RooArgList(frac_2p5, frac_2p4, frac_2p2, frac_1p8, frac_1p2, frac_0p5, frac_1p8_right, frac_1p2_right), true);

  
  model.fitTo(temp_data, Extended());
  
  //model.plotOn(xframe, LineColor(kRed), LineStyle(kDashed));
  //model.plotOn(xframe, Components(RooArgSet(frac_2p5, frac_2p4, frac_2p2, frac_1p8, frac_1p2, frac_0p5, frac_1p8_right, frac_1p2_right)), LineColor(kRed), LineStyle(kDashed));
  model.Print("t");
  //temp_data.plotOn(xframe, LineColor(kBlack), LineStyle(kSolid));
  //pdf_P_2p5.plotOn(xframe, LineColor(kGreen), LineStyle(kSolid));
  

  RooAddPdf model_test("model_test", "model_test", RooArgList(pdf_P_2p5, pdf_P_2p4), RooArgList(frac_2p5), true);
  RooDataSet *data_test = model_test.generate(x, 1000);
  model_test.fitTo(*data_test);



  TCanvas *c = new TCanvas("rf706_histpdf", "rf706_histpdf", 800, 400);
  xframe->Draw();
  TString WORKING_DIR = getenv("LArProf_WD");
  TString pdfname;
  pdfname = WORKING_DIR + "/output/plots/BeamStudy/Template_fit_Muon_1GeV.pdf";
  c -> SaveAs(pdfname);


  f_in -> Close();
  */
}
