#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

// ---------------------------
// Time-window normalization
// ---------------------------
double Tcorr = 1.0;
double WindowNorm(double lambda){ return 1.0 / (1.0 - exp(-lambda*Tcorr)); }

// ---------------------------
// Read parameters from file
// ---------------------------
void GetPars(vector<string> &names, vector<double> &pars, vector<pair<double,double>> &bounds, string filename){
    ifstream file(filename);
    if(!file.is_open()){ cerr<<"Cannot open parameter file: "<<filename<<endl; exit(1);}
    string line;
    while(getline(file,line)){
        if(line.empty() || line[0]=='#' || line[0]==' ') continue;
        stringstream ss(line);
        string name; double X0, lo, hi;
        ss >> name >> X0 >> lo >> hi;
        names.push_back(name);
        pars.push_back(X0);
        bounds.push_back({lo,hi});
    }
    file.close();
}

// ---------------------------
// Full Bateman decay model
// Parameters:
// 0-4: lambda_p, lambda_d, lambda_bn, lambda_gd, lambda_2n_d
// 5: lambda_bn_gd (new for beta-n granddaughter)
// 6: N0
// 7: bg
// 8-13: efficiencies eff_p, eff_d, eff_bn, eff_gd, eff_2n_d, eff_bn_gd
// 14-15: branching ratios Pn, P2n
// ---------------------------
double TotalModelFull(double *x, double *par){
    double t = x[0];

    double lambda_p      = par[0];
    double lambda_d      = par[1];
    double lambda_bn     = par[2];
    double lambda_gd     = par[3];
    double lambda_2n_d   = par[4];
    double lambda_bn_gd  = par[5];
    double N0            = par[6];
    double bg            = par[7];

    double eff_p         = par[8];
    double eff_d         = par[9];
    double eff_bn        = par[10];
    double eff_gd        = par[11];
    double eff_2n_d      = par[12];
    double eff_bn_gd     = par[13];

    double Pn            = par[14];
    double P2n           = par[15];

    // physical limits
    if(Pn<0) Pn=0; if(Pn>1) Pn=1;
    if(P2n<0) P2n=0; if(P2n>1) P2n=1;
    double branch_d = 1.0 - Pn - P2n;
    if(branch_d<0) branch_d=0;

    double rate = 0.0;

    // -------- Parent --------
    double Np = N0 * exp(-lambda_p*t);
    rate += WindowNorm(lambda_p) * eff_p * lambda_p * Np;

    // -------- Daughter --------
    double Nd = branch_d * N0 * (lambda_p/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
    rate += WindowNorm(lambda_d) * eff_d * lambda_d * Nd;

    // -------- Granddaughter --------
    double denom1 = (lambda_d-lambda_p)*(lambda_gd-lambda_p);
    double denom2 = (lambda_p-lambda_d)*(lambda_gd-lambda_d);
    double denom3 = (lambda_p-lambda_gd)*(lambda_d-lambda_gd);
    if(fabs(denom1)<1e-20) denom1=1e-20;
    if(fabs(denom2)<1e-20) denom2=1e-20;
    if(fabs(denom3)<1e-20) denom3=1e-20;
    double Ngd = branch_d * N0 * lambda_p * lambda_d * (
        exp(-lambda_p*t)/denom1 +
        exp(-lambda_d*t)/denom2 +
        exp(-lambda_gd*t)/denom3
    );
    rate += WindowNorm(lambda_gd) * eff_gd * lambda_gd * Ngd;

    // -------- Beta-n daughter --------
    double Nbn = Pn * N0 * (lambda_p/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
    rate += WindowNorm(lambda_bn) * eff_bn * lambda_bn * Nbn;

    // -------- Beta-n granddaughter (separate) --------
    double denom1_bn = (lambda_bn-lambda_p)*(lambda_bn_gd-lambda_p);
    double denom2_bn = (lambda_p-lambda_bn)*(lambda_bn_gd-lambda_bn);
    double denom3_bn = (lambda_p-lambda_bn_gd)*(lambda_bn-lambda_bn_gd);
    if(fabs(denom1_bn)<1e-20) denom1_bn=1e-20;
    if(fabs(denom2_bn)<1e-20) denom2_bn=1e-20;
    if(fabs(denom3_bn)<1e-20) denom3_bn=1e-20;
    double Nbgd = Pn * N0 * lambda_p * lambda_bn * (
        exp(-lambda_p*t)/denom1_bn +
        exp(-lambda_bn*t)/denom2_bn +
        exp(-lambda_bn_gd*t)/denom3_bn
    );
    rate += WindowNorm(lambda_bn_gd) * eff_bn_gd * lambda_bn_gd * Nbgd;

    // -------- Beta-2n daughter --------
    double N2n = P2n * N0 * (lambda_p/(lambda_2n_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_2n_d*t));
    rate += WindowNorm(lambda_2n_d) * eff_2n_d * lambda_2n_d * N2n;

    // -------- Background --------
    rate += bg;

    return rate;
}

// ---------------------------
int main(int argc,char* argv[]){
    if(argc<5){ cout<<"Usage: ./bateman_fit <rebin> <rootfile> <histname> <paramfile>\n"; return 1;}
    int rebin = stoi(argv[1]);
    string filename = argv[2];
    string histname = argv[3];
    string paramfile = argv[4];

    vector<string> names;
    vector<double> pars;
    vector<pair<double,double>> bounds;
    GetPars(names,pars,bounds,paramfile);
    int npars = pars.size();

    TFile *f = TFile::Open(filename.c_str());
    if(!f || f->IsZombie()){ cerr<<"Error opening ROOT file\n"; return 1;}
    TH1D* h=(TH1D*)f->Get(histname.c_str());
    if(!h){ cerr<<"Histogram not found\n"; return 1;}
    h->Rebin(rebin);

    Tcorr = h->GetXaxis()->GetXmax();
    double xmin = 1; //h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    TF1 *fit = new TF1("fit",TotalModelFull,xmin,xmax,npars);
    for(int i=0;i<npars;i++){
        fit->SetParName(i,names[i].c_str());
        fit->SetParameter(i,pars[i]);
        fit->SetParLimits(i,bounds[i].first,bounds[i].second);
    }

    h->Fit("fit","R");

    cout<<"\nHalf-lives (ms):"<<endl;
    for(int i=0;i<5;i++){ // first 5 parameters are decay constants
        double lambda = fit->GetParameter(i);
        double hl = TMath::Log(2)/lambda;
        double hl_err = TMath::Log(2)/pow(lambda,2)*fit->GetParError(i);
        cout<<fit->GetParName(i)<<": "<<hl<<" ± "<<hl_err<<endl;
    }

    cout<<"Chi2/NDF = "<<fit->GetChisquare()<<" / "<<fit->GetNDF()
        <<" = "<<fit->GetChisquare()/fit->GetNDF()<<endl;

    int nbins = h->GetNbinsX();
    TH1D *hfit_total = new TH1D("fit_total","Total Fit",nbins,xmin,xmax);
    for(int i=1;i<=nbins;i++){
        double x=h->GetBinCenter(i);
        hfit_total->SetBinContent(i,fit->Eval(x));
    }

    TCanvas *c = new TCanvas("c","Bateman Fit",1200,800);
    h->SetLineColor(kBlack); h->Draw("E");
    fit -> Draw("same");
    //hfit_total->SetLineColor(kOrange+2); hfit_total->Draw("same");
    TLegend *leg=new TLegend(0.1,0.1,0.2,0.2);
    leg->AddEntry(h,"Data","l");
    leg->AddEntry(hfit_total,"Total Fit","l");
    leg->Draw();
    c->SaveAs("bateman_fit.png");

    TFile *fout = new TFile("bateman_fit_output.root","RECREATE");
    h->Write();
    hfit_total->Write();
    c -> Write();
    fout->Close();

    cout<<"Fit and histogram saved."<<endl;
    return 0;
}
