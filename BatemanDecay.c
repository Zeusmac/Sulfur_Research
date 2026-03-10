#include <cstring>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

class DecayNetwork {

private:

    int maxN;
    int n;

    char **names;

    double *lambda;
    double *N0;
    double **branch;

public:

    DecayNetwork(int maxTerms)
    {
        maxN = maxTerms;
        n = 0;

        names = new char*[maxN];
        lambda = new double[maxN];
        N0 = new double[maxN];

        branch = new double*[maxN];

        for(int i=0;i<maxN;i++)
        {
            names[i] = new char[64];
            branch[i] = new double[maxN];

            lambda[i] = 0;
            N0[i] = 0;

            for(int j=0;j<maxN;j++)
                branch[i][j] = 0;
        }
    }

    ~DecayNetwork()
    {
        for(int i=0;i<maxN;i++)
        {
            delete[] names[i];
            delete[] branch[i];
        }

        delete[] names;
        delete[] branch;
        delete[] lambda;
        delete[] N0;
    }

    int findIndex(const char *name)
    {
        for(int i=0;i<n;i++)
        {
            if(strcmp(names[i],name)==0)
                return i;
        }
        return -1;
    }

    int addNucleus(const char *name)
    {
        int idx = findIndex(name);

        if(idx >= 0)
            return idx;

        strcpy(names[n],name);
        n++;

        return n-1;
    }

    void setLambda(const char *name,double val)
    {
        int i = findIndex(name);
        if(i>=0) lambda[i] = val;
    }

    void setInitial(const char *name,double val)
    {
        int i = findIndex(name);
        if(i>=0) N0[i] = val;
    }

    void addDecay(const char *parent,const char *daughter,double br)
    {
        int p = addNucleus(parent);
        int d = addNucleus(daughter);

        branch[p][d] = br;
    }

private:

    void derivatives(double *N,double *dNdt)
    {
        for(int i=0;i<n;i++)
        {
            dNdt[i] = -lambda[i]*N[i];

            for(int j=0;j<n;j++)
                dNdt[i] += branch[j][i]*lambda[j]*N[j];
        }
    }

public:

    void solve(double t,double *result)
    {
        int steps = 300;
        double dt = t/steps;

        double *N = new double[n];
        double *k1 = new double[n];
        double *k2 = new double[n];
        double *k3 = new double[n];
        double *k4 = new double[n];
        double *temp = new double[n];

        for(int i=0;i<n;i++)
            N[i] = N0[i];

        for(int s=0;s<steps;s++)
        {

            derivatives(N,k1);

            for(int i=0;i<n;i++)
                temp[i] = N[i] + 0.5*dt*k1[i];

            derivatives(temp,k2);

            for(int i=0;i<n;i++)
                temp[i] = N[i] + 0.5*dt*k2[i];

            derivatives(temp,k3);

            for(int i=0;i<n;i++)
                temp[i] = N[i] + dt*k3[i];

            derivatives(temp,k4);

            for(int i=0;i<n;i++)
                N[i] += dt/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
        }

        for(int i=0;i<n;i++)
            result[i] = N[i];

        delete[] N;
        delete[] k1;
        delete[] k2;
        delete[] k3;
        delete[] k4;
        delete[] temp;
    }

};

DecayNetwork net;
double decay_model(double *x, double *par)
{
    double t = x[0];

    // parameters
    double lambda_p  = par[0];
    double lambda_d  = par[1];
    double lambda_bn = par[2];
    double lambda_gd = par[3];

    double N0 = par[4];
    double background = par[5];

    // update network parameters
    net.setLambda("Parent",lambda_p);
    net.setLambda("Daughter",lambda_d);
    net.setLambda("BnDaughter",lambda_bn);
    net.setLambda("Granddaughter",lambda_gd);

    net.setInitial("Parent",N0);

    // solve Bateman system
    net.solve(t);

    double counts =
          net.getPopulation("Parent")
        + net.getPopulation("Daughter")
        + net.getPopulation("BnDaughter")
        + net.getPopulation("Granddaughter")
        + background;

    return counts;
}


int main()
{
    // -------------------------
    // build decay network
    // -------------------------

    net.addIsotope("Parent");
    net.addIsotope("Parent");
    net.addIsotope("Daughter");
    net.addIsotope("BnDaughter");
    net.addIsotope("Granddaughter");

    net.addDecay("Parent","Daughter",0.8);
    net.addDecay("Parent","BnDaughter",0.2);
    net.addDecay("Daughter","Granddaughter",1.0);

    // -------------------------
    // load histogram
    // -------------------------

    TFile f("Mar_5_26/decay_curve44S.root");
    TH1D *h = (TH1D*)f.Get("decay_curve44S");

    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    // -------------------------
    // define fit function
    // -------------------------

    TF1 *fit = new TF1("fit",decay_model,xmin,xmax,6);

    fit->SetParNames(
        "lambda_p",
        "lambda_d",
        "lambda_bn",
        "lambda_gd",
        "N0",
        "background"
    );

    // initial guesses
    fit->SetParameters(
        0.5,
        0.3,
        0.2,
        0.1,
        1000,
        10
    );

    // bounds
    fit->SetParLimits(0,0,10);
    fit->SetParLimits(1,0,10);
    fit->SetParLimits(2,0,10);
    fit->SetParLimits(3,0,10);
    fit->SetParLimits(4,0,1e6);
    fit->SetParLimits(5,0,1000);

    // -------------------------
    // perform fit
    // -------------------------

    h->Fit(fit,"R");

    // -------------------------
    // draw result
    // -------------------------

    TCanvas *c = new TCanvas();
    h->Draw();
    fit->Draw("same");

    c->SaveAs("fit.png");

    return 0;
}
