#ifndef LANGAUSFIT_H
#define LANGAUSFIT_H

#include "TH1.h"
#include "TF1.h"

TF1 *langaufit(TH1 *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status);

TF1 *runlangaufit(TH1 *his, int plane);

#endif
