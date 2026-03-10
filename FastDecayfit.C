#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include "FastDecay.h"

// create global network
FastDecayNetwork net(500,20.0);

// branching fraction for beta-n
double Pn = 0.15;

double decay_model(double *x, double *par)
{
    double t = x[0];

    double lambda_p = par[0];
    double lambda_d = par[1];
    double lambda_bn = par[2];
    double lambda_gd = par[3];

    double N0 = par[4];
    double background = par[5];

    // set parameters
    net.setLambda("Parent",lambda_p);
    net.setLambda("Daughter",lambda_d);
    net.setLambda("BnDaughter",lambda_bn);
    net.setLambda("Granddaughter",lambda_gd);

    net.setInitial("Parent",N0);
    net.setInitial("Daughter",0);
    net.setInitial("BnDaughter",0);
    net.setInitial("Granddaughter",0);

    // build solution once
    net.buildSolution();

    double total =
        net.getPopulation("Parent",t)
      + net.getPopulation("Daughter",t)
      + net.getPopulation("BnDaughter",t)
      + net.getPopulation("Granddaughter",t)
      + background;

    return total;
}

void fast_decay_fit()
{
    // build decay network
    net.addIsotope("Parent");
    net.addIsotope("Daughter");
    net.addIsotope("BnDaughter");
    net.addIsotope("Granddaughter");

    net.addDecay("Parent","Daughter",1-Pn);
    net.addDecay("Parent","BnDaughter",Pn);
    net.addDecay("Daughter","Granddaughter",1.0);

    // load histogram
    TFile f("Mar_5_26/decay_curve44S.root");
    TH1D *h = (TH1D*)f.Get("decay_curve44S");

    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    // TF1 fit
    TF1 *fit = new TF1("fit",decay_model,xmin,xmax,6);
    fit->SetParNames("lambda_p","lambda_d","lambda_bn","lambda_gd","N0","background");

    // initial guesses
    fit->SetParameters(0.5,0.3,0.2,0.1,1000,1);

    // bounds
    fit->SetParLimits(0,0,10);
    fit->SetParLimits(1,0,10);
    fit->SetParLimits(2,0,10);
    fit->SetParLimits(3,0,10);
    fit->SetParLimits(4,0,1e6);
    fit->SetParLimits(5,0,1000);

    // perform fit
    h->Fit(fit,"R");

    // draw
    TCanvas *c = new TCanvas();
    h->Draw();
    fit->Draw("same");
    c->SaveAs("fast_decay_fit.png");
}
