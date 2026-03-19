#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

using namespace std;

class BatemanActivity {

private:

    vector<double> lambda;

public:

    BatemanActivity(const vector<double>& lambdas){
        lambda = lambdas;
    }

    // coefficient for Bateman term
    double coeff(int k,int j){

        double prod_lambda = 1.0;

        for(int i=0;i<k;i++)
            prod_lambda *= lambda[i];

        double denom = 1.0;

        for(int m=0;m<=k;m++){
            if(m==j) continue;
            denom *= (lambda[m] - lambda[j]);
        }

        return prod_lambda/denom;
    }

    // population Nk(t)
    double population(int k,double t,double N0){

        double sum = 0.0;

        for(int j=0;j<=k;j++)
            sum += coeff(k,j)*exp(-lambda[j]*t);

        return N0*sum;
    }

    // activity Ai(t)
    double activity(int k,double t,double N0){

        double Nk = population(k,t,N0);

        return lambda[k]*Nk;
    }

};
double decay_model(double *x,double *par){

    double t = x[0];

    double lambda_p  = par[0];
    double lambda_d  = par[1];
    double lambda_gd = par[2];
    double lambda_bn = par[3];
    double lambda_2n = par[4];

    double N0 = par[5];

    double Pn  = par[6];
    double P2n = par[7];

    double eff_p  = par[8];
    double eff_d  = par[9];
    double eff_gd = par[10];
    double eff_bn = par[11];
    double eff_2n = par[12];

    double bg = par[13];

    double Pbeta = 1.0 - Pn - P2n;

    // main beta chain
    vector<double> main_chain = {lambda_p,lambda_d,lambda_gd};
    BatemanActivity mainSolver(main_chain);

    double A_p  = eff_p  * mainSolver.activity(0,t,N0);
    double A_d  = eff_d  * Pbeta * mainSolver.activity(1,t,N0);
    double A_gd = eff_gd * Pbeta * mainSolver.activity(2,t,N0);

    // beta-n chain
    vector<double> bn_chain = {lambda_p,lambda_bn};
    BatemanActivity bnSolver(bn_chain);

    double A_bn = eff_bn * Pn * bnSolver.activity(1,t,N0);

    // beta-2n chain
    vector<double> b2n_chain = {lambda_p,lambda_2n};
    BatemanActivity b2nSolver(b2n_chain);

    double A_2n = eff_2n * P2n * b2nSolver.activity(1,t,N0);

    return A_p + A_d + A_gd + A_bn + A_2n + bg;
}
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

    //Tmax = h->GetXaxis()->GetXmax();
    double xmin = 1; //h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    TF1 *fit = new TF1("fit",decay_model,xmin,xmax,npars);
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
  
