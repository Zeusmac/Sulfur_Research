#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
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
	    // k1
	    double k1p  = -lambda_p * pop.Np;

	    double k1d  =
		(1-Pn)*lambda_p*pop.Np
		- lambda_d*pop.Nd;

	    double k1bn =
		Pn*lambda_p*pop.Np
		- lambda_bn*pop.Nbn;

	    double k1gd =
		lambda_bn*pop.Nbn
		- lambda_gd*pop.Ngd;


	    // k2 (half step)
	    double Np2  = pop.Np  + 0.5*dt*k1p;
	    double Nd2  = pop.Nd  + 0.5*dt*k1d;
	    double Nbn2 = pop.Nbn + 0.5*dt*k1bn;
	    double Ngd2 = pop.Ngd + 0.5*dt*k1gd;

	    double k2p  = -lambda_p * Np2;

	    double k2d  =
		(1-Pn)*lambda_p*Np2
		- lambda_d*Nd2;

	    double k2bn =
		Pn*lambda_p*Np2
		- lambda_bn*Nbn2;

	    double k2gd =
		lambda_bn*Nbn2
		- lambda_gd*Ngd2;


	    // k3 (half step again)
	    double Np3  = pop.Np  + 0.5*dt*k2p;
	    double Nd3  = pop.Nd  + 0.5*dt*k2d;
	    double Nbn3 = pop.Nbn + 0.5*dt*k2bn;
	    double Ngd3 = pop.Ngd + 0.5*dt*k2gd;

	    double k3p  = -lambda_p * Np3;

	    double k3d  =
		(1-Pn)*lambda_p*Np3
		- lambda_d*Nd3;

	    double k3bn =
		Pn*lambda_p*Np3
		- lambda_bn*Nbn3;

	    double k3gd =
		lambda_bn*Nbn3
		- lambda_gd*Ngd3;


	    // k4 (full step)
	    double Np4  = pop.Np  + dt*k3p;
	    double Nd4  = pop.Nd  + dt*k3d;
	    double Nbn4 = pop.Nbn + dt*k3bn;
	    double Ngd4 = pop.Ngd + dt*k3gd;

	    double k4p  = -lambda_p * Np4;

	    double k4d  =
		(1-Pn)*lambda_p*Np4
		- lambda_d*Nd4;

	    double k4bn =
		Pn*lambda_p*Np4
		- lambda_bn*Nbn4;

	    double k4gd =
		lambda_bn*Nbn4
		- lambda_gd*Ngd4;


	    // final update
	    pop.Np  += dt/6.0 * (k1p  + 2*k2p  + 2*k3p  + k4p);
	    pop.Nd  += dt/6.0 * (k1d  + 2*k2d  + 2*k3d  + k4d);
	    pop.Nbn += dt/6.0 * (k1bn + 2*k2bn + 2*k3bn + k4bn);
	    pop.Ngd += dt/6.0 * (k1gd + 2*k2gd + 2*k3gd + k4gd);
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

void GetPars(vector<string> &Parnames, vector<double> &pars, vector<pair<double,double>> &bounds, char *filename)
{

	string line;
	stringstream ss;
	fstream file(filename, ios::in);
	if(file.is_open()) 
	{
		while(getline(file, line))
		{
		char f = line[0];
        	if (line.size() == 0) continue;
        	if (f == '#') continue;
        	if (f == ' ') continue;

		ss.clear();
                ss.str(line);

		string Names;
		double par, low, high;
		ss >> Names >> par >> low >> high;
		pair<double, double> b = {low, high};
		Parnames.push_back(Names);
		pars.push_back(par);
		bounds.push_back(b);
		} 	
	}
	else
	{
		cerr << "cant open decayPars.txt"<< endl;
	}
	file.close();
}

int main(int argc, char* argv[])
{
    vector<double> pars;
    vector<string> Parnames;
    vector<pair<double, double>> bounds;
    GetPars(Parnames,pars, bounds, argv[2]);
		
    int bin = stoi(argv[1]);
    TFile* f = new TFile(argv[3]);
    TH1D* h = (TH1D*)f->Get("decay_curve44S");

    h -> Rebin(bin);
    printf( "%i ms per bin\n", bin);
                                    
    tmax = h->GetXaxis()->GetXmax();

    TF1 *fit = new TF1("fit",DecayModel, 0,tmax,11);
    
    for(int i = 0; i < 11; i++)
    {
	fit -> SetParName(i, Parnames[i].c_str());
	fit -> SetParameter(i, pars[i]);
	fit -> SetParLimits(i, bounds[i].first, bounds[i].second);
	printf("%s		Inital: %f	bounds: [%f	, 	%f]\n", Parnames[i].c_str(), pars[i], bounds[i].first, bounds[i].second);
    }
 
    h->Fit("fit","R","",1 ,tmax);

    Double_t decayhl = TMath::Log(2)/(fit ->GetParameter(0));
    Double_t decay_error = (TMath::Log(2)/TMath::Power((fit -> GetParameter(0)),2))*(fit -> GetParError(0));
    printf("Decay: %7fms +_ %7fms \n", decayhl, decay_error);

    Double_t decay_daughter = TMath::Log(2)/(fit -> GetParameter(1));
    Double_t decay_daughter_error = (TMath::Log(2)/TMath::Power((fit ->GetParameter(1)),2))*(fit -> GetParError(1));
    printf("Decay_daughter: %7fms +_ %7fms \n", decay_daughter, decay_daughter_error);
 
    Double_t decay_betan = TMath::Log(2)/(fit -> GetParameter(2));
    Double_t decay_betan_error = (TMath::Log(2)/TMath::Power((fit ->GetParameter(2)),2))*(fit -> GetParError(2));
    printf("Decay_Nbeta: %7fms +_ %7fms \n", decay_betan, decay_betan_error);
    
    Double_t decay_gdaughter = TMath::Log(2)/(fit -> GetParameter(3));
    Double_t decay_gdaughter_error = (TMath::Log(2)/TMath::Power((fit ->GetParameter(3)),2))*(fit -> GetParError(3));
    printf("Decay_gdaughter: %7fms +_ %7fms \n", decay_gdaughter, decay_gdaughter_error);

    Double_t chi2 = fit -> GetChisquare();
    double_t NDF = fit -> GetNDF();
    printf("Chi^2/NDF: %7f\n", chi2/NDF);

    TCanvas *c = new TCanvas("decayfit","Fit of 44S");
    c->cd();
    h -> Draw("E");
    fit -> Draw("same");
    c->SaveAs("numerical_decay_fit.root");
    
    c->SaveAs("numerical_decay_fit.png");
    
    cout << "\nFit complete\n";

}
