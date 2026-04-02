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
#include "TMath.h"

using namespace std;

// ---------------------------
// Time-window normalization
// ---------------------------
double Tcorr = 1.0;
double WindowNorm(double lambda){ return 1.0;}// /(1.0 - exp(-lambda*Tcorr)); }

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
string FormatNumber(double x)
{
    std::ostringstream ss;

    double ax = fabs(x);

    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss << std::scientific << std::setprecision(6) << x;
    else
        ss << std::fixed << std::setprecision(6) << x;

    return ss.str();
}


void PrintFitResults(TH1* hist,
                       TF1* fit,
                       TFitResultPtr r,
                       const double* initial_pars,
                       double ms_per_bin,
                       const char* filename)
{

    std::ofstream file(filename);

    if(!file.is_open())
    {
        std::cout<<"Could not open output file\n";
        return;
    }

    int npar = fit->GetNpar();

    double xmin = fit->GetXmin();
    double xmax = fit->GetXmax();

    file<<"\n================ FIT REPORT ================\n\n";

    file<<"Histogram           : "<<hist->GetName()<<"\n";
    file<<"Fit Function        : "<<fit->GetName()<<"\n";

    file<<"\nFit Range\n";
    file<<"-------------------------------------------\n";
    file<<"xmin = "<<FormatNumber(xmin)<<"\n";
    file<<"xmax = "<<FormatNumber(xmax)<<"\n";

    file<<"\nBin Width\n";
    file<<"-------------------------------------------\n";
    file<<"ms per bin = "<<FormatNumber(ms_per_bin)<<"\n";

    file<<"\nFit Statistics\n";
    file<<"-------------------------------------------\n";

    file<<"Chi2   = "<<FormatNumber(r->Chi2())<<"\n";
    file<<"NDF    = "<<r->Ndf()<<"\n";
    file<<"EDM    = "<<FormatNumber(r->Edm())<<"\n";
    file<<"NCalls = "<<r->NCalls()<<"\n";

    file<<"\nParameter Table\n";
    file<<"----------------------------------------------------------------------------------\n";

    file<<std::left
        <<std::setw(18)<<"Parameter"
        <<std::setw(16)<<"Initial"
        <<std::setw(16)<<"Value"
        <<std::setw(16)<<"Error"
        <<std::setw(16)<<"LowerBound"
        <<std::setw(16)<<"UpperBound"
        <<"Status\n";

    file<<"----------------------------------------------------------------------------------\n";

    for(int i=0;i<npar;i++)
    {

        const char* name = fit->GetParName(i);

        double val = r->Parameter(i);
        double err = r->ParError(i);

        double lo,hi;
        fit->GetParLimits(i,lo,hi);

        bool fixed = r->IsParameterFixed(i);

        file<<std::setw(18)<<name;
        file<<std::setw(16)<<FormatNumber(initial_pars[i]);
        file<<std::setw(16)<<FormatNumber(val);

        if(fixed)
            file<<std::setw(16)<<"---";
        else
            file<<std::setw(16)<<FormatNumber(err);

        file<<std::setw(16)<<FormatNumber(lo);
        file<<std::setw(16)<<FormatNumber(hi);

        if(fixed)
            file<<"fixed";
        else
            file<<"limited";

        file<<"\n";
    }

    file<<"----------------------------------------------------------------------------------\n";

    file<<"\nHalf-lives (ms)\n";
    file<<"-------------------------------------------\n";

    for(int i=0;i<npar;i++)
    {

        std::string name = fit->GetParName(i);

        if(name.find("lambda") != std::string::npos)
        {

            double lambda = r->Parameter(i);
            double err = r->ParError(i);

            double half = log(2.0)/lambda;

            double half_err = 0;
            if(lambda != 0)
                half_err = (err/lambda)*half;

            file<<std::setw(18)<<name
                <<FormatNumber(half)
                <<" ± "
                <<FormatNumber(half_err)
                <<"\n";
        }
    }

    file<<"\n-------------------------------------------\n";

    file<<"Chi2/NDF = "
        <<FormatNumber(r->Chi2())
        <<" / "
        <<r->Ndf()
        <<" = "
        <<FormatNumber(r->Chi2()/r->Ndf())
        <<"\n";

    file<<"\n===========================================\n";

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
     

    if(t < 0) return bg;
   // if(t == 0) return 0;

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
// Component functions with background
// ---------------------------
double ParentComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double N0       = par[6];
    double eff_p    = par[8];
    double bg       = par[7];
    double Np = N0*exp(-lambda_p*t);
    return eff_p*lambda_p*Np + bg;
}

double DaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double lambda_d = par[1];
    double N0       = par[6];
    double eff_d    = par[9];
    double Pn      = par[14];
    double P2n     = par[15];
    double branch_d = 1.0 - Pn - P2n;
    if(branch_d<0) branch_d=0;
    double Nd = branch_d * N0 * (lambda_p/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
    return eff_d*lambda_d*Nd + par[7];
}

// Similarly define for other components: Granddaughter, Beta-n daughter, Beta-n granddaughter, Beta-2n daughter
// For brevity, I’ll do Granddaughter:
double GranddaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double lambda_d = par[1];
    double lambda_gd = par[3];
    double N0       = par[6];
    double eff_gd   = par[11];
    double Pn      = par[14];
    double P2n     = par[15];
    double branch_d = 1.0 - Pn - P2n;
    if(branch_d<0) branch_d=0;

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
    return eff_gd*lambda_gd*Ngd + par[7];
}
// Beta-n daughter + background
double BetaNDaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p   = par[0];
    double lambda_bn  = par[2];
    double N0         = par[6];
    double eff_bn     = par[10];
    double Pn         = par[14];
    double bg         = par[7];

    double Nbn = Pn * N0 * (lambda_p/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
    return eff_bn * lambda_bn * Nbn + bg;
}

// Beta-n granddaughter + background
double BetaNGranddaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p      = par[0];
    double lambda_bn     = par[2];
    double lambda_bn_gd  = par[5];
    double N0            = par[6];
    double eff_bn_gd     = par[13];
    double Pn            = par[14];
    double bg            = par[7];

    double denom1 = (lambda_bn-lambda_p)*(lambda_bn_gd-lambda_p);
    double denom2 = (lambda_p-lambda_bn)*(lambda_bn_gd-lambda_bn);
    double denom3 = (lambda_p-lambda_bn_gd)*(lambda_bn-lambda_bn_gd);
    if(fabs(denom1)<1e-20) denom1=1e-20;
    if(fabs(denom2)<1e-20) denom2=1e-20;
    if(fabs(denom3)<1e-20) denom3=1e-20;

    double Nbgd = Pn * N0 * lambda_p * lambda_bn * (
        exp(-lambda_p*t)/denom1 +
        exp(-lambda_bn*t)/denom2 +
        exp(-lambda_bn_gd*t)/denom3
    );
    return eff_bn_gd * lambda_bn_gd * Nbgd + bg;
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

    Tcorr = h->GetXaxis()->GetXmax();
    double xmin = stoi(argv[2]); // h->GetXaxis()->GetXmin();
    double xmax =stoi(argv[3]);//h->GetXaxis()->GetXmax();

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
    PrintFitResults(h, fit, result, init, rebin, argv[7]);	

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
    c -> cd();
    h->SetLineColor(kBlack); h->Draw("E");
    fit -> Draw("same");
	    // Draw individual components
	TF1* f_parent = new TF1("Parent", ParentComponent, xmin, xmax, npars);
	TF1* f_daughter = new TF1("Daughter", DaughterComponent, xmin, xmax, npars);
	TF1* f_granddaughter = new TF1("Granddaughter", GranddaughterComponent, xmin, xmax, npars);
	TF1* f_bn_daughter = new TF1("BetaN_Daughter", BetaNDaughterComponent, xmin, xmax, npars);
	TF1* f_bn_granddaughter = new TF1("BetaN_Granddaughter", BetaNGranddaughterComponent, xmin, xmax, npars);

	// Copy parameters from fit
	for(int i=0;i<npars;i++){
	    f_parent->SetParameter(i, fit->GetParameter(i));
	    f_daughter->SetParameter(i, fit->GetParameter(i));
	    f_granddaughter->SetParameter(i, fit->GetParameter(i));
	    f_bn_daughter->SetParameter(i, fit->GetParameter(i));
            f_bn_granddaughter->SetParameter(i, fit->GetParameter(i));
	}

	// Set line colors
	f_parent->SetLineColor(kRed);
	f_daughter->SetLineColor(kBlue);
	f_granddaughter->SetLineColor(kGreen);
	f_bn_daughter->SetLineColor(kBlue);
	f_bn_granddaughter->SetLineColor(kGreen);

	// Draw everything
	h->Draw("E");
	fit->Draw("same");               // total fit
	f_parent->Draw("same");          // parent
	f_daughter->Draw("same");        // daughter
	f_granddaughter->Draw("same");   // granddaughter
	f_bn_daughter->Draw("same");           // beta-n daughter
	f_bn_granddaughter->Draw("same");      // beta-n granddaughter
	// Legend
	TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
	leg->AddEntry(h,"Data","lep");
	leg->AddEntry(fit,"Total Fit","l");
	leg->AddEntry(f_parent,"Parent + bg","l");
	leg->AddEntry(f_daughter,"Daughter + bg","l");
	leg->AddEntry(f_granddaughter,"Granddaughter + bg","l");
leg->AddEntry(f_bn_daughter,"Beta-n Daughter + bg","l");
leg->AddEntry(f_bn_granddaughter,"Beta-n Granddaughter + bg","l");
	leg->Draw();
    c->SaveAs("bateman_fit.png");

    TFile *fout = new TFile("bateman_fit_output.root","RECREATE");
    h->Write();
    //hfit_total->Write();
    c -> Write();
    fout->Close();

    cout<<"Fit and histogram saved."<<endl;
    return 0;
}
