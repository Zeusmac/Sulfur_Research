#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include "TApplication.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
using namespace std;

//////////////////////////////////////////////////////////////
// Utility
//////////////////////////////////////////////////////////////

string FormatNumber(double x){
    std::ostringstream ss;
    double ax = fabs(x);
    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss << scientific << setprecision(6) << x;
    else
        ss << fixed << setprecision(6) << x;
    return ss.str();
}
struct FitConfig {
    string filename;
    string histname;
    double xmin, xmax;
    int rebin;
    string outfile;

    vector<string> names;
    vector<double> init;
    vector<pair<double,double>> bounds;

    // NEW: channel switches
    int use_beta = 1;
    int use_bn   = 1;
    int use_2n   = 1;
};
bool GetPars(FitConfig &cfg, const string &filename)
{
    ifstream file(filename);
    if(!file.is_open()){
        cerr<<"Cannot open config file\n";
        return false;
    }

    string key;

    while(file >> key)
    {
        if(key=="FILE") file >> cfg.filename;
        else if(key=="HIST") file >> cfg.histname;
        else if(key=="RANGE") file >> cfg.xmin >> cfg.xmax;
        else if(key=="REBIN") file >> cfg.rebin;
	else if(key=="USE_BETA") file >> cfg.use_beta;
	else if(key=="USE_BN")   file >> cfg.use_bn;
	else if(key=="USE_2N")   file >> cfg.use_2n;
        else if(key=="OUTPUT") file >> cfg.outfile;

        else if(key[0]=='#'){
            getline(file,key);
        }
        else
        {
            // parameter line
            string name = key;
            double init, lo, hi;
            file >> init >> lo >> hi;

            cfg.names.push_back(name);
            cfg.init.push_back(init);
            cfg.bounds.push_back({lo,hi});
        }
    }

    return true;
}
void PrintFitResultsAppend(TF1* fit, TFitResultPtr r, ofstream &file,
                          const string &label,
                          const vector<double> &init,
                          const vector<pair<double,double>> &bounds,
                          int rebin, double xmin, double xmax)
{
    file << "\n============================================================\n";
    file << "                     " << label << "\n";
    file << "============================================================\n\n";

    // -------------------------
    // Fit summary
    // -------------------------
    file << "Rebin        : " << rebin << "\n";
    file << "Fit range    : [" << xmin << ", " << xmax << "]\n";
    file << "Chi2 / NDF   : " << FormatNumber(r->Chi2())
         << " / " << r->Ndf()
         << "  =  " << FormatNumber(r->Chi2()/r->Ndf()) << "\n";
    file << "EDM          : " << FormatNumber(r->Edm()) << "\n";
    file << "NCalls       : " << r->NCalls() << "\n";
    file << "Probability  : " << FormatNumber(r->Prob()) << "\n\n";

    // -------------------------
    // Fit status
    // -------------------------
    int fitStatus = r->Status();
    int covStatus = r->CovMatrixStatus();

    file << "Fit Status   : " << fitStatus;
    if(fitStatus == 0) file << "  (OK)";
    else               file << "  (PROBLEM)";
    file << "\n";

    file << "Cov Matrix   : " << covStatus;
    if(covStatus == 3) file << "  (FULL ACCURATE)";
    else if(covStatus == 2) file << "  (FORCED POS-DEF)";
    else if(covStatus == 1) file << "  (APPROXIMATE)";
    else file << "  (NOT VALID)";
    file << "\n\n";

    if(fitStatus != 0 || covStatus < 3){
        file << "⚠ WARNING: Fit may be unreliable\n\n";
    }

    // -------------------------
    // Parameter table
    // -------------------------
    file << left
         << setw(16) << "Parameter" << "\t"
         << setw(14) << "Value" << "\t"
         << setw(24) << "Error (- / +)" << "\t"
         << setw(14) << "Initial" << "\t"
         << setw(20) << "Bounds" << "\t"
         << "Status\n";

    file << string(90,'-') << "\n";

    for(int i=0;i<fit->GetNpar();i++){
        bool fixed = r->IsParameterFixed(i);

        stringstream bnd;
        bnd << "[" << bounds[i].first << ", " << bounds[i].second << "]";

        string errStr;

        if(fixed){
            errStr = "---";
        }
        else{
            double el = r->LowerError(i);
            double eh = r->UpperError(i);

            if(el == 0 && eh == 0){
                double e = r->ParError(i);
                errStr = "± " + FormatNumber(e);
            }
            else{
                stringstream ss;
                ss << "-" << FormatNumber(fabs(el))
                   << "  +" << FormatNumber(eh);
                errStr = ss.str();
            }
        }

        file << setw(16) << fit->GetParName(i) << "\t"
             << setw(14) << FormatNumber(r->Parameter(i)) << "\t"
             << setw(24) << errStr << "\t"
             << setw(14) << FormatNumber(init[i]) << "\t"
             << setw(20) << bnd.str() << "\t"
             << (fixed ? "fixed" : "free") << "\n";
    }

    // -------------------------
    // Half-lives
    // -------------------------
    file << "\n------------------ Half-lives (ms) ------------------\n";

    for(int i=0;i<fit->GetNpar();i++){
        string name = fit->GetParName(i);

        if(name.find("lambda") != string::npos){
            double lam = r->Parameter(i);

            double el = r->LowerError(i);
            double eh = r->UpperError(i);

            double t = TMath::Log(2)/lam;

            string errStr;

            if(el == 0 && eh == 0){
                double e = r->ParError(i);
                double terr = (e/lam)*t;
                errStr = "± " + FormatNumber(terr);
            }
            else{
                double t_low  = log(2)/(lam + eh);
                double t_high = log(2)/(lam + el);

                double err_low  = t - t_low;
                double err_high = t_high - t;

                stringstream ss;
                ss << "-" << FormatNumber(fabs(err_low))
                   << "  +" << FormatNumber(err_high);

                errStr = ss.str();
            }

            file << setw(16) << name
                 << setw(14) << FormatNumber(t)
                 << errStr << "\n";
        }
    }

    // -------------------------
    // Correlation matrix
    // -------------------------
    int n = fit->GetNpar();

    file << "\n---------------- Correlation Matrix ----------------\n\n";

    file << setw(12) << "";
    for(int j=0;j<n;j++)
        file << setw(12) << fit->GetParName(j);
    file << "\n";

    for(int i=0;i<n;i++){
        file << setw(12) << fit->GetParName(i);
        for(int j=0;j<n;j++){
            file << setw(12) << FormatNumber(r->Correlation(i,j));
        }
        file << "\n";
    }

    file << "\n============================================================\n";
}
void DrawContoursAdvanced(TFitResultPtr result, const string& outroot="contours.root") {
    TFile* fout = new TFile(outroot.c_str(), "RECREATE");

    int npar = result->NPar();  // <-- correct
    for (int i = 0; i < npar; ++i) {
        for (int j = i+1; j < npar; ++j) {
            unsigned int npoints = 1000;
            vector<double> x(npoints), y(npoints);
  	    if(result->IsParameterFixed(i) || result->IsParameterFixed(j)) continue;
	    if(!result->CovMatrixStatus()){
 	    std::cout << "Covariance matrix not valid, contours cannot be computed." << std::endl;
   	    return;
	    }
            bool ok = result->Contour(i, j, npoints, x.data(), y.data(), 0.683);
            if (!ok){    std::cout << "Contour failed for pair " 
              << result->GetParameterName(i) << " vs " 
              << result->GetParameterName(j) << std::endl;
    		continue;
		};

            TGraph* gr = new TGraph(npoints, x.data(), y.data());
            gr->SetTitle(Form("Contour %s vs %s; %s; %s",
                              result->GetParameterName(i).c_str(),
                              result->GetParameterName(j).c_str(),
                              result->GetParameterName(i).c_str(),
                              result->GetParameterName(j).c_str()));
            gr->SetLineColor(kRed);
            gr->SetLineWidth(2);
            gr->Write(Form("Contour_%s_vs_%s", 
                           result->GetParameterName(i).c_str(),
                           result->GetParameterName(j).c_str()));
        }
    }

    fout->Close();
    cout << "All contours saved to " << outroot << endl;
}
FitConfig *gCfg = nullptr;
//////////////////////////////////////////////////////////////
// ORIGINAL FULL BATEMAN MODEL
//////////////////////////////////////////////////////////////

Double_t TotalModelFull(double *x, double *par)
{
    double t = x[0];
    if(t < 0) return par[7];
    // ---------------- CHANNEL SWITCHES ----------------
    int use_beta = gCfg->use_beta;
    int use_bn   = gCfg->use_bn;
    int use_2n   = gCfg->use_2n;
    // ---------------- PARAMETERS ----------------
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


    // ---------------- SAFETY ----------------
    if(Pn < 0) Pn = 0; if(Pn > 1) Pn = 1;
    if(P2n < 0) P2n = 0; if(P2n > 1) P2n = 1;

    double raw_beta = 1.0 - Pn - P2n;
    if(raw_beta < 0) raw_beta = 0;

    // ---------------- NORMALIZATION ----------------
    double total_branch = 0;

    if(use_beta) total_branch += raw_beta;
    if(use_bn)   total_branch += Pn;
    if(use_2n)   total_branch += P2n;

    if(total_branch <= 0) return 1e30;

    double branch_d = use_beta ? raw_beta / total_branch : 0;
    double branch_bn = use_bn ? Pn / total_branch : 0;
    double branch_2n = use_2n ? P2n / total_branch : 0;

    Double_t rate = 0.0;

    // ---------------- PARENT ----------------
    double Np = N0 * exp(-lambda_p*t);
    rate += eff_p * lambda_p * Np;

    // ---------------- DAUGHTER ----------------
    if(use_beta){
        double denom = (lambda_d - lambda_p);
        if(fabs(denom) < 1e-12) denom = 1e-12;

        double Nd = branch_d * N0 * (lambda_p/denom) *
                    (exp(-lambda_p*t)-exp(-lambda_d*t));

        rate += eff_d * lambda_d * Nd;
    }

    // ---------------- BETA-n ----------------
    if(use_bn){
        double denom = (lambda_bn - lambda_p);
        if(fabs(denom) < 1e-12) denom = 1e-12;

        double Nbn = branch_bn * N0 * (lambda_p/denom) *
                     (exp(-lambda_p*t)-exp(-lambda_bn*t));

        rate += eff_bn * lambda_bn * Nbn;
    }

    // ---------------- BETA-2n ----------------
    if(use_2n){
        double denom = (lambda_2n_d - lambda_p);
        if(fabs(denom) < 1e-12) denom = 1e-12;

        double N2n = branch_2n * N0 * (lambda_p/denom) *
                     (exp(-lambda_p*t)-exp(-lambda_2n_d*t));

        rate += eff_2n_d * lambda_2n_d * N2n;
    }

    // ---------------- BACKGROUND ----------------
    rate += bg;

    // ---------------- SAFETY: prevent negative ----------------
    if(rate < 0) return 1e30;

    return rate;
}
//////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////

int main(int argc,char* argv[])
{
	if(argc < 2){
		std::cout<<"Usage: ./fit config.txt\n";
		return 1;
	}
	FitConfig cfg;
	gCfg = &cfg;
	if(!GetPars(cfg, argv[1])) return 1;
        TApplication app("app", &argc, argv);
	TFile *f = TFile::Open(cfg.filename.c_str());
	TH1D* h = (TH1D*)f->Get(cfg.histname.c_str());

	h->Rebin(cfg.rebin);

	TF1 *fit = new TF1("fit",TotalModelFull,cfg.xmin,cfg.xmax,cfg.init.size());

	for(int i=0;i<cfg.init.size();i++){
		fit->SetParName(i,cfg.names[i].c_str());
		fit->SetParameter(i,cfg.init[i]);
		fit->SetParLimits(i,cfg.bounds[i].first,cfg.bounds[i].second);
	}

	TFitResultPtr result = h->Fit("fit","S R M","",cfg.xmin,cfg.xmax);

	// Save results
	ofstream out(cfg.outfile);
	PrintFitResultsAppend(fit,result,out,"Fit",
		      cfg.init,cfg.bounds,
		      cfg.rebin,cfg.xmin,cfg.xmax);
	out.close();
	int p1=0,p2=1;
	double maxCorr=0;

	for(int i=0;i<fit->GetNpar();i++){
		for(int j=i+1;j<fit->GetNpar();j++){
			double c = fabs(result->Correlation(i,j));
			if(c > maxCorr){
			    maxCorr = c;
			    p1=i; p2=j;
			}
		}
	}
	DrawContoursAdvanced(result, "bateman_fit_contours.root");
	TCanvas *c = new TCanvas("c","Fit",800,600);
	h->Draw("E");
	fit->Draw("same");
	c -> Update();
	app.Run();
	
	c->SaveAs("fit.root");

	cout<<"Done. Results saved.\n";

	return 0;
}
