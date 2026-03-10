#include <iostream>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMinuit.h"

using namespace std;

static double tmax = 1.0;
static int nSteps = 200;

struct Populations
{
    double Np;
    double Nd;
    double Nbn;
    double Ngd;
};

Populations RK4Solve(double t,
                     double lambda_p,
                     double lambda_d,
                     double lambda_bn,
                     double lambda_gd,
                     double Pn,
                     double N0)
{

    double dt = t / nSteps;

    Populations pop;
    pop.Np = N0;
    pop.Nd = 0;
    pop.Nbn = 0;
    pop.Ngd = 0;

    for(int i=0;i<nSteps;i++)
    {

        double dNp = -lambda_p * pop.Np;

        double dNd =
            (1-Pn)*lambda_p*pop.Np
            - lambda_d*pop.Nd;

        double dNbn =
            Pn*lambda_p*pop.Np
            - lambda_bn*pop.Nbn;

        double dNgd =
            lambda_bn*pop.Nbn
            - lambda_gd*pop.Ngd;

        pop.Np += dNp * dt;
        pop.Nd += dNd * dt;
        pop.Nbn += dNbn * dt;
        pop.Ngd += dNgd * dt;
    }

    return pop;
}

double DecayModel(double *x, double *par)
{

    double t = x[0];

    double lambda_p  = par[0];
    double lambda_d  = par[1];
    double lambda_bn = par[2];
    double lambda_gd = par[3];

    double N0 = par[4];
    double bg = par[5];

    double Pn = par[6];

    double eff_p  = par[7];
    double eff_d  = par[8];
    double eff_bn = par[9];
    double eff_gd = par[10];

    Populations pop =
        RK4Solve(t,
                 lambda_p,
                 lambda_d,
                 lambda_bn,
                 lambda_gd,
                 Pn,
                 N0);

    double counts =
        eff_p  * lambda_p  * pop.Np +
        eff_d  * lambda_d  * pop.Nd +
        eff_bn * lambda_bn * pop.Nbn +
        eff_gd * lambda_gd * pop.Ngd +
        bg;

    return counts;
}

int main()
{

    TFile *f = new TFile("Mar_5_26/decay_curve44S.root");
    TH1D *h = (TH1D*)f->Get("decay_curve44S");
    h -> Rebin(10);

    tmax = h->GetXaxis()->GetXmax();

    TF1 *fit = new TF1("fit",DecayModel,0,tmax,11);

    fit->SetParNames(
        "lambda_p",
        "lambda_d",
        "lambda_bn",
        "lambda_gd",
        "N0",
        "bg",
        "Pn",
        "eff_p",
        "eff_d",
        "eff_bn",
        "eff_gd"
    );

    fit->SetParameters(
        0.01,//lam p
        0.005,// lam d
        0.002, // lam bn
        0.001, // lam gd
        220000,//N0
        220000,//bckgnd
        0.2, // pn
        0.5,
        0.5,
        0.5,
        0.5
    );

    fit->SetParLimits(0,1e-5,1);
    fit->SetParLimits(1,1e-5,1);
    fit->SetParLimits(2,1e-5,1);
    fit->SetParLimits(3,1e-5,1);

    fit->SetParLimits(4,0,1e7);
    fit->SetParLimits(5,0,1e6);

    fit->SetParLimits(6,0,1);

    fit->SetParLimits(7,0,1);
    fit->SetParLimits(8,0,1);
    fit->SetParLimits(9,0,1);
    fit->SetParLimits(10,0,1);

    h->Fit("fit","R");

    TCanvas *c = new TCanvas();
    h->Draw();
    fit->Draw("same");

    c->SaveAs("numerical_decay_fit.png");

    cout << "\nFit complete\n";

}
