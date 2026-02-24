#ifndef MULTGAUSTEST_H
#define MULTGAUSTEST_H


double Gaus_fit(double *dim, double *par);
double decay_parent(double *dim, double *par);

Double_t decay_background(Double_t *dim, Double_t *par);
Double_t decay_nbeta(Double_t *dim, Double_t *par);


//44SGated
void Decayp(TH1 *hist);

//44S
void Decay_exp(TH1 *hist);

//46Cl
void Decayb(TH1 *hist);

void MultGausFit(TH1 *hist, double lower, double upper);

#endif