#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TApplication.h"
using namespace std;
//////////////////////////////////////////////////////////////
// Utility
//////////////////////////////////////////////////////////////
double chi2(TH1D* hist, TF1* eq){
     double chi2_sum = 0;
     for (int i = 1; i <= hist -> GetNbinsX(); ++i) {  // loop over bins 1..NbinsX
         double x = hist->GetBinCenter(i);
         double c = hist->GetBinContent(i);
         if (c <= 0)  continue;  // skip empty bins
         if (x < eq->GetXmin() || x > eq->GetXmax()) continue;
 
         double model = eq->Eval(x);
         double err2 = hist->GetBinError(i);
         if (err2 <= 0) continue;  // avoid division by zero
         chi2_sum += (model - c) * (model - c) / c;
     }
	return chi2_sum;
}
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
void PrintFitResultsAppend(TH1D* h,TF1* fit, TFitResultPtr r, ofstream &file,
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
    file << "\n--------------------My chi2/eff_d/eff_bn -----------------------\n";
    double a = fit-> GetParameter(6);
    file << "eff_d: "<< (a) << endl;
    file << "eff_bn: "<< (1-a) << endl;
    file << chi2(h,fit) << endl;

    // -------------------------
    // Correlation matrix
    // -------------------------


//    int n = fit->GetNpar();	
	//file << r->Print("V") << endl;
 //file << "\n---------------- Correlation Matrix ----------------\n\n";

 //   file << setw(12) << "";
 //   for(int j=0;j<n;j++)
 //       file << setw(12) << fit->GetParName(j);
 //   file << "\n";

 //   for(int i=0;i<n;i++){
      //  file << setw(12) << fit->GetParName(i);
    //    for(int j=0;j<n;j++){
    //        file << setw(12) << FormatNumber(r->Correlation(i,j));
  //      }
//        file << "\n";
//    }

    file << "\n============================================================\n";
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
    double lambda_b   	 = par[3]; 
    double Nb		 = par[4];
    double N0            = par[5];
    double eff_d         = par[6];
    double eff_bn        =  (1 - par[6]);
    

    // physical limits
    double rate = 0;

    // -------- Parent --------
    double Np = N0 * exp(-lambda_p*t);
    rate +=  Np;

    // -------- Daughter --------
    double Nd = N0 * (1/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
    rate +=  eff_d * lambda_d * Nd;

    // -------- Beta-n daughter --------
    double Nbn =  N0 * (1/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
    rate += eff_bn * lambda_bn * Nbn;

   // --------- expo bckgrnd----------
    double Nbg = Nb * exp(-lambda_b*t);
    rate += Nbg;

    if(t < 0) return Nb*exp(lambda_b*t);
    return rate;
}
// ---------------------------
// Component functions with background
// ---------------------------
double ParentComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double N0       = par[5];
    double bg       = par[6];
    double Np = N0 * exp(-lambda_p*t);
    return Np;
}

double DaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double lambda_d = par[1];
    double N0       = par[5];
    double eff_d         = par[6];
    double Nd =  N0 * (1/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
    return eff_d* lambda_d * Nd;
}

// Similarly define for other components: Granddaughter, Beta-n daughter, Beta-n granddaughter, Beta-2n daughter
// Beta-n daughter + background
double BetaNDaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p   = par[0];
    double lambda_bn  = par[2];
    double N0         = par[5];
    double eff_bn        =  (1 - par[6]);
    double Nbn =  N0 * (1/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
    return eff_bn * lambda_bn * Nbn;
}

double Expobg(double *x, double *par){
    double t = x[0];
    double lambda_b = par[3];
    double Nb       = par[4];
    double Nbg = Nb*exp(-lambda_b*t);
    return Nbg;
}
// ---------------------------
int main(int argc,char* argv[])
{
	if(argc < 2){
		std::cout<<"Usage: ./fit config.txt <FitOptions> <setlogy>\n";
		std::cout<< "./fit -h for options description" << endl;
		return 1;
	}
	if(std::string(argv[1]) == "-h"){
		std::cout<< "fit options: 0 Default Migrad with HESSE\n 1 Default Migrad with MINOS\n 2 Improved Migrad with HESSE\n 3 Improved Migrad with MINOS\n To SetLogY == 1" << endl;
		return 1;
	}  
	FitConfig cfg;
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
	TString fitoptions = "S R";
	int option = 0;
	if(!(argv[2] == NULL)) option = stoi(argv[2]); 	
	if(option == 1) fitoptions + " E" ;   
	if(option == 2) fitoptions + " M";  
	if(option == 3) fitoptions + " M E";  
	TFitResultPtr result = h->Fit("fit",fitoptions,"",cfg.xmin,cfg.xmax);
 
	// Save results
	ofstream out(cfg.outfile + ".txt");
	PrintFitResultsAppend(h,fit,result,out,"Fit",
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
	TCanvas *c = new TCanvas("c","Fit",800,600);
	    // Draw individual components
	TF1* f_parent = new TF1("Parent", ParentComponent, -1000, cfg.xmax, cfg.init.size());
	TF1* f_daughter = new TF1("Daughter", DaughterComponent, -1000, cfg.xmax, cfg.init.size());
	TF1* f_bn_daughter = new TF1("BetaN_Daughter", BetaNDaughterComponent, -1000, cfg.xmax, cfg.init.size());
	TF1* f_expbg = new TF1("expbg", Expobg, -1000, cfg.xmax, cfg.init.size());
	
	// Copy parameters from fit
	for(int i=0;i<fit->GetNpar();i++){
	    f_parent->SetParameter(i, fit->GetParameter(i));
	    f_daughter->SetParameter(i, fit->GetParameter(i));
	    f_bn_daughter->SetParameter(i, fit->GetParameter(i));
	    f_expbg->SetParameter(i, fit->GetParameter(i));
	}
	// Set line colors
	c -> cd();
	f_parent->SetLineColor(kOrange);
	f_daughter->SetLineColor(kYellow + 1);
	f_expbg->SetLineColor(kGreen);
	f_bn_daughter->SetLineColor(kBlue - 4);
	// Draw everything
	int logy = 0;
	if(!(argv[3] == NULL)) logy = stoi(argv[3]); 
	if(logy == 1) c -> SetLogy();
	h->Draw("E");
	//int MaxY = h -> GetBinContent(1);
//	h -> GetYaxis()->SetRangeUser(0,MaxY + 4000);
	fit->Draw("same");               // total fit
	f_parent->Draw("same");          // parent
	f_daughter->Draw("same");        // daughter
	f_expbg->Draw("same");   // granddaughter
	f_bn_daughter->Draw("same");           // beta-n daughter
	// Legend
	TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
	leg->AddEntry(h,"Data","lep");
	leg->AddEntry(fit,"Total Fit","l");
	leg->AddEntry(f_parent,"Parent","l");
	leg->AddEntry(f_daughter,"Daughter","l");
	leg->AddEntry(f_bn_daughter,"Beta-n Daughter","l");
	leg->AddEntry(f_expbg,"Expo Backgrnd","l");
	//leg->AddEntry(f_bck,"Bck","l");
	leg->Draw();
	c -> Update();
	// Save current stdout
	FILE* old_stdout = stdout;

	// Redirect stdout to file
	FILE* file_stdout = freopen((cfg.outfile + "Matrix.txt").c_str(), "w", stdout);

	// Call verbose print
	result->Print("V");

	// Flush and restore stdout
	fflush(stdout);
	stdout = old_stdout;
	TFile* outf = new TFile((cfg.outfile + ".root").c_str(), "RECREATE");
	if (!f || f->IsZombie()) {
   	     std::cerr << "Error opening file\n";
	}
	c -> Write();
	h -> Write();
	fit ->Write();
	f_parent-> Write();
	f_daughter -> Write();
	f_expbg -> Write();
	f_bn_daughter -> Write();
	outf -> Close();
	app.Run();
	return 0;
}
