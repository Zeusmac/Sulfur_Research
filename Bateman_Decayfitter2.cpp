#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"

using namespace std;

double Tcorr = 1.0;
double WindowNorm(double lambda){ return 1.0; }

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

string FormatNumber(double x){
    std::ostringstream ss;
    double ax = fabs(x);
    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss << std::scientific << std::setprecision(6) << x;
    else
        ss << std::fixed << std::setprecision(6) << x;
    return ss.str();
}

// Draw a TPaveText box with chi2 and half-lives
void DrawFitStatsBox(TCanvas* c,TF1* fit, double xmin, double xmax){
    TPaveText* box = new TPaveText(0.65,0.65,0.95,0.95,"NDC");
    box->SetFillColor(0);
    box->SetBorderSize(1);
    box->SetTextAlign(12);
    box->AddText(("Chi2/NDF: "+FormatNumber(fit->GetChisquare())+" / "+to_string(fit->GetNDF())+
                 " = "+FormatNumber(fit->GetChisquare()/fit->GetNDF())).c_str());

    for(int i=0;i<5;i++){ // first 5 parameters are decay constants
        double lambda = fit->GetParameter(i);
        double hl = TMath::Log(2)/lambda;
        double hl_err = TMath::Log(2)/pow(lambda,2)*fit->GetParError(i);
        string text = fit->GetParName(i) + string(": ") + FormatNumber(hl)+" ± "+FormatNumber(hl_err);
        box->AddText(text.c_str());
    }
    box->Draw();
}

// ---------------------------
// Full Bateman decay model
// ---------------------------
double TotalModelFull(double *x, double *par){
    double t = x[0];
    if(t<0) return par[7]; // background

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

    // -------- Beta-n granddaughter --------
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
    if(argc<8){ cout<<"Usage: ./bateman_fit <rebin> <low> <high> <rootfile> <histname> <paramfile> <result filename>\n"; return 1;}
    int rebin = stoi(argv[1]);
    string filename = argv[4];
    string histname = argv[5];
    string paramfile = argv[6];

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

    double xmin = stod(argv[2]);
    double xmax = stod(argv[3]);

    TF1 *fit = new TF1("fit",TotalModelFull,xmin,xmax,npars);
    for(int i=0;i<npars;i++){
        fit->SetParName(i,names[i].c_str());
        fit->SetParameter(i,pars[i]);
        fit->SetParLimits(i,bounds[i].first,bounds[i].second);
    }

    TFitResultPtr result =  h->Fit("fit","SR","",xmin,xmax);
    double init[16];
    for(int i=0;i<fit->GetNpar();i++)
	    init[i] = fit->GetParameter(i);

TCanvas *c = new TCanvas("c","Bateman Fit",1200,800);
c->cd();

h->SetTitle("Bateman Fit;Time (ms);Counts");
h->Draw("E");

// total fit
fit->SetLineColor(kRed);
fit->Draw("same");

// -------- Parent component --------
TF1* f_parent = new TF1("parent","[0]*exp(-[1]*x)",xmin,xmax);
f_parent->SetParameters(N0*eff_p*lambda_p,lambda_p);
f_parent->SetLineColor(kBlue);
f_parent->SetLineStyle(2);
f_parent->Draw("same");

// -------- Daughter component --------
TF1* f_daughter = new TF1("daughter",
"[0]*(exp(-[1]*x)-exp(-[2]*x)) + [3]",
xmin,xmax);

f_daughter->SetParameters(scale_d,lambda_p,lambda_d,bg);
f_daughter->SetLineColor(kGreen+2);
f_daughter->SetLineStyle(2);
f_daughter->Draw("same");

// -------- Background --------
TF1* f_bg = new TF1("bg","[0]",xmin,xmax);
f_bg->SetParameter(0,bg);
f_bg->SetLineColor(kMagenta);
f_bg->SetLineStyle(3);
f_bg->Draw("same");


// legend
TLegend *leg = new TLegend(0.15,0.15,0.35,0.35);
leg->SetBorderSize(0);
leg->SetFillStyle(0);

leg->AddEntry(h,"Data","lep");
leg->AddEntry(fit,"Total Fit","l");
leg->AddEntry(f_parent,"Parent","l");
leg->AddEntry(f_daughter,"Daughter","l");
leg->AddEntry(f_bg,"Background","l");

leg->Draw();

// stats box
DrawFitStatsBox(c,fit,xmin,xmax);

c->Update();   // IMPORTANT
    c->SaveAs("bateman_fit.png");

    TFile *fout = new TFile("bateman_fit_output.root","RECREATE");
    h->Write();
    fit->Write();
    f_parent->Write();
    f_daughter->Write();
    f_bg->Write();
    c->Write();
    fout->Close();

    cout<<"Fit and histogram saved."<<endl;
    return 0;
}
