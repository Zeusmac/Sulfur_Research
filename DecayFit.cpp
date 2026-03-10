#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>

#include "DecayNetwork.h"

DecayNetwork net;

double decay_model(double *x, double *par)
{

    double t = x[0];

    double lambda_parent  = par[0];
    double lambda_daughter = par[1];
    double lambda_bn = par[2];
    double lambda_gd = par[3];

    double N0 = par[4];
    double background = par[5];

    net.setLambda("Parent",lambda_parent);
    net.setLambda("Daughter",lambda_daughter);
    net.setLambda("BnDaughter",lambda_bn);
    net.setLambda("Granddaughter",lambda_gd);

    net.setInitial("Parent",N0);

    net.solve(t);

    double total =
        net.getPopulation("Parent")
      + net.getPopulation("Daughter")
      + net.getPopulation("BnDaughter")
      + net.getPopulation("Granddaughter")
      + background;

    return total;
}

int main()
{

    //---------------------------------
    // build decay network
    //---------------------------------

    net.addIsotope("Parent");
    net.addIsotope("Daughter");
    net.addIsotope("BnDaughter");
    net.addIsotope("Granddaughter");

    net.addDecay("Parent","Daughter",0.85);
    net.addDecay("Parent","BnDaughter",0.15);
    net.addDecay("Daughter","Granddaughter",1.0);

    //---------------------------------
    // load histogram
    //---------------------------------

    TFile file("Mar_5_26/decay_curve44S.root");

    TH1D *h = (TH1D*)file.Get("decay_curve44S");

    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    //---------------------------------
    // create fit function
    //---------------------------------

    TF1 *fit = new TF1("fit",decay_model,xmin,xmax,6);

    fit->SetParNames(
        "lambda_parent",
        "lambda_daughter",
        "lambda_bn",
        "lambda_gd",
        "N0",
        "background"
    );

    //---------------------------------
    // initial guesses
    //---------------------------------

    fit->SetParameters(
        0.5,
        0.3,
        0.4,
        0.2,
        1000,
        5
    );

    //---------------------------------
    // parameter bounds
    //---------------------------------

    fit->SetParLimits(0,0,10);
    fit->SetParLimits(1,0,10);
    fit->SetParLimits(2,0,10);
    fit->SetParLimits(3,0,10);

    fit->SetParLimits(4,0,1e7);
    fit->SetParLimits(5,0,1000);

    //---------------------------------
    // perform fit
    //---------------------------------

    h->Fit(fit,"R");

    //---------------------------------
    // draw result
    //---------------------------------

    TCanvas *c1 = new TCanvas("c1","Decay Fit",800,600);

    h->Draw();
    fit->Draw("same");

    c1->SaveAs("decay_fit.png");

}
