#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>

#include "TBrowser.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
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
TGraph* ChiVsPar(TH1D* hist, TF1* funk, TFitResultPtr r,
                 int par_number, double stepsize)
{
    TGraph *gr = new TGraph();

    double p0 = funk->GetParameter(par_number);
    double chi2_min = chi2(hist, funk);
    double NDF = r -> Ndf();
    gr->AddPoint(p0, chi2_min/NDF); // Δχ² = 0 at minimum
    cout << chi2_min << endl;
    TF1 *f_up = (TF1*)funk->Clone();
    TF1 *f_dn = (TF1*)funk->Clone();
    double p1 = p0;
    double p2 = p0;	 
    bool go_up = true, go_dn = true;

    while (go_up || go_dn) {

        // --- upward scan ---
        if (go_up) {
            p1 += stepsize;
            f_up->SetParameter(par_number, p1);

            double chi = chi2(hist, f_up);
            double dchi = fabs(chi - chi2_min);
//	    cout << chi << " " << dchi << endl;
            gr->AddPoint(p1, chi/NDF);

            if (dchi >= 2) go_up = false;
        }

        // --- downward scan ---
        if (go_dn) {
            p2 -= stepsize;
            f_dn->SetParameter(par_number, p2);

            double chi = chi2(hist, f_dn);
            double dchi = fabs(chi - chi2_min);
            gr->AddPoint(p2, chi/NDF);

            if (dchi >= 2) go_dn = false;
        }
    }
    std::cout << "Name: " << gr->GetName() << std::endl;
    std::cout << "N points: " << gr->GetN() << std::endl;
    gr -> Sort();
    return gr;
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
class bfitter {
    public:
	    string filename;
	    string histname;
	    double xmin, xmax;
	    int rebin;
	    string outfile;

	    vector<string> names;
	    vector<double> init;
	    vector<pair<double,double>> bounds;

	    // NEW: channel switches
	    int use_expo_bg = 0;
	    int use_flat_bg = 0; 
	    int use_beta = 1;
	    int use_bn   = 1;
	    int use_2n   = 1;

	bool GetPars(const string &Config_filename)
	{
	    ifstream file(Config_filename);
	    if(!file.is_open()){
		cerr<<"Cannot open config file\n";
		return false;
	    }

	    string key;
		while(file >> key)
		{
		    if(key[0]=='#'){
			file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		    }

		    if(key=="FILE"){
			file >> filename;
			continue;
		    }
		    else if(key=="HIST"){
			file >> histname;
			continue;
		    }
		    else if(key=="RANGE"){
			file >> xmin >> xmax;
			continue;
		    }
		    else if(key=="REBIN"){
			file >> rebin;
			continue;
		    }
		    else if(key=="OUTPUT"){
			file >> outfile;
			continue;
		    }
		    else if(key=="USE_BN_CHAN"){
			file >> use_bn;
			continue;
		    }
		    else if(key=="USE_EXPO_BG"){
			file >> use_expo_bg;
			continue;
		    }
		    else if(key=="USE_FLAT_BG"){
			file >> use_flat_bg;
			continue;
		    }

		    // parameter line ONLY if none matched
		    string par_name = key;
		    double inital, lo, hi;

		    if(!(file >> inital >> lo >> hi)){
			cerr << "Error reading parameter line for: " << par_name << endl;
			break;
		    }

		    names.push_back(par_name);
		    init.push_back(inital);
		    bounds.push_back({lo,hi});
		}
	    return true;
	}

	double TotalModelFull(double *x, double *par){
	    double t = x[0];

	    double lambda_p      = par[0];
	    double lambda_d      = par[1];
	    double lambda_bn     = par[2];
	    double N0            = par[3];
	    double eff_d         = par[4];
	    double eff_bn        =  (1 - par[4]);
	    double bg 		= par[5];
	    double lambda_b 	= par[6]; 
	    double Nb		 = par[7];

	    // physical limits
	    double rate = 0;	
	    if(use_flat_bg == 1) rate += bg;

	    // -------- Parent --------
	    double Np = N0 * exp(-lambda_p*t);
	    rate +=  Np;

	    // -------- Daughter --------
	    double Nd = N0 * (1/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
	    rate +=  eff_d * lambda_d * Nd;

	    // -------- Beta-n daughter --------
	    double Nbn =  N0 * (1/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
	    if(use_expo_bg == 1) rate += eff_bn * lambda_bn * Nbn;

	   // --------- expo bckgrnd----------
	    double Nbg = Nb * exp(-lambda_b*t);

	    if(use_expo_bg == 1) rate += Nbg;

	    if(t < 0 && (use_expo_bg == 1)) return Nb*exp(lambda_b*t);
	    return rate;
	}
};

// ---------------------------
// Component functions with background
// ---------------------------
double ParentComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double N0       = par[3];
    double bg 	    = par[5];
    double Np = N0 * exp(-lambda_p*t) + bg;
    return Np;
}

double DaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p = par[0];
    double lambda_d = par[1];
    double N0       = par[3];
   double bg       = par[5];
    double eff_d         = par[4];
    double Nd =  N0 * (1/(lambda_d-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_d*t));
    return eff_d* lambda_d * Nd + bg;
}

// Similarly define for other components: Granddaughter, Beta-n daughter, Beta-n granddaughter, Beta-2n daughter
// Beta-n daughter + background
double BetaNDaughterComponent(double *x, double *par){
    double t = x[0];
    double lambda_p   = par[0];
    double lambda_bn  = par[2];
    double N0         = par[3];
    double bg       = par[5];
    double eff_bn        =  (1 - par[4]);
    double Nbn =  N0 * (1/(lambda_bn-lambda_p)) * (exp(-lambda_p*t)-exp(-lambda_bn*t));
    return eff_bn * lambda_bn * Nbn + bg;
}

double Expobg(double *x, double *par){
    double t = x[0];
    double lambda_b = par[6];
    double Nb       = par[7];
    double bg       = par[5];
    double Nbg = Nb*exp(-lambda_b*t);
    return Nbg + bg;
}

void PrintFitResultsAppend(TH1D* h,TF1* fit, TFitResultPtr r, ofstream &file,
                          string filename,TString fitoptions,const string &label,
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
    file << filename << endl; 
    file << fitoptions << endl;
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
    double a = fit-> GetParameter(4);
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
int main(int argc,char* argv[])
{
	if(argc < 2){
		std::cout<<"Usage: ./fit config.txt <FitOptions> <setlogy>\n";
		std::cout<< "./fit -h for options description" << endl;
		return 1;
	}
	if(std::string(argv[1]) == "-h"){
		std::cout<< "fit options: 0\n Default Migrad with HESSE\n 1 Default Migrad with MINOS\n 2 Improved Migrad with HESSE\n 3 Improved Migrad with MINOS\n 4 LogLiklihood fit\n 5  logliklihood Mirgad\n Logliklihood fit with MINOs\n To SetLogY == 1" << endl;
		return 1;
	}  
	bfitter cfg;
	double parent_halflife;
	double parent_halflife_err;
	double daughter_halflife;
	double daughter_halflife_err;
	if(!cfg.GetPars(argv[1])) return 1;
	TFile *f = TFile::Open(cfg.filename.c_str());
	TH1D* h = (TH1D*)f->Get(cfg.histname.c_str());
	h ->Sumw2();
	h->Rebin(cfg.rebin);
	int t = h->GetNbinsX();
	TGraphErrors* gr = new TGraphErrors(t);

	for (int i = 0; i < t; i++) {
    		double x = h ->GetBinCenter(i+1);
    		double y = h ->GetBinContent(i+1);

    		double ey = TMath::Sqrt(y);   // Poisson error
    		double ex = 0;//cfg.rebin;         // usually zero unless you want bin width

    		gr->SetPoint(i, x, y);
    		gr->SetPointError(i, ex, ey);
	}	
	//TGraphErrors *gr = new TGraphErrors(h);
	cout << "N parameters: " << cfg.init.size() << endl;
	TF1 *fit = new TF1("fit",
    		[&](double *x, double *par) {
        		return cfg.TotalModelFull(x, par);
    		},
    			cfg.xmin, cfg.xmax, cfg.init.size());
	
	for(int i=0;i<cfg.init.size();i++){
		fit->SetParName(i,cfg.names[i].c_str());
		fit->SetParameter(i,cfg.init[i]);
		fit->SetParLimits(i,cfg.bounds[i].first,cfg.bounds[i].second);
	}
	TString fitoptions = "S R EX0";
	int option = 0;
	if(!(argv[2] == NULL)) option = stoi(argv[2]); 	
	if(option == 1) fitoptions += " E" ;   
	if(option == 2) fitoptions += " M";  
	if(option == 3) fitoptions += " M E"; 
	if(option == 4) fitoptions += " L"; 
	if(option == 5) fitoptions += " M L"; 
	if(option == 6) fitoptions += " L E"; 
	cout <<" fitoptions :" << fitoptions << endl;
	
	TFitResultPtr result = gr->Fit(fit,fitoptions,"",cfg.xmin,cfg.xmax);
 	
	parent_halflife = TMath::Log(2)/fit -> GetParameter(0);
	daughter_halflife = TMath::Log(2)/fit -> GetParameter(1);
	
	// Save results
	ofstream out(cfg.outfile + ".txt");
	PrintFitResultsAppend(h,fit,result,out,cfg.filename.c_str(),fitoptions,"Fit",
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
	
	TF1* f_bg = new TF1("bg","[0]", cfg.xmax);
	TF1* f_expbg = new TF1("expbg", Expobg, 0, cfg.xmax, cfg.init.size());
	TF1* f_p = new TF1("Parent", ParentComponent, 0, cfg.xmax, cfg.init.size());
	TF1* f_d = new TF1("Daughter", DaughterComponent, 0, cfg.xmax, cfg.init.size());
	TF1* f_bnd = new TF1("BetaN_Daughter", BetaNDaughterComponent, 0, cfg.xmax, cfg.init.size());
	// Copy parameters from fit
	for(int i=0;i<fit->GetNpar();i++){
	    f_p->SetParameter(i, fit->GetParameter(i));
	    f_d->SetParameter(i, fit->GetParameter(i));
	    f_bnd->SetParameter(i, fit->GetParameter(i));
	    if(cfg.use_expo_bg == 1) f_expbg->SetParameter(i, fit->GetParameter(i));
	}
	if(cfg.use_flat_bg == 1) f_bg -> SetParameters(fit -> GetParameter(5));

	// Set line colors
	c -> cd();
	f_p->SetLineColor(kOrange);
	f_d->SetLineColor(kYellow + 1);
	f_expbg->SetLineColor(kGreen);
	f_bnd->SetLineColor(kBlue - 4);

	// Draw everything
	int logy = 0;
	if(!(argv[3] == NULL)) logy = stoi(argv[3]); 
	if(logy == 1) c -> SetLogy();
	gr->Draw("AP");
	//gr->SetMinimum(0.1); // to let log scale work 
	gr->SetMinimum(h->GetMinimum()*.77777775);
	gr->SetMaximum(h->GetMaximum());	
	fit -> Draw("same");
	f_p->Draw("same");          // parent
	f_d->Draw("same");        // daughter
	if(cfg.use_expo_bg == true) f_expbg->Draw("same");   // granddaughter
	if(cfg.use_bn == 1) f_bnd->Draw("same");           // beta-n daughter
	if(cfg.use_flat_bg == true) f_bg -> Draw("same"); 
	// Legend
	TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(h,"Data","lep");
	leg->AddEntry(fit,"Total Fit","l");
	leg->AddEntry(f_p,"Parent","l");
	leg->AddEntry(f_d,"Daughter","l");
	if(cfg.use_bn == 1) leg->AddEntry(f_bnd,"Beta-n Daughter","l");
	if(cfg.use_expo_bg == 1) leg->AddEntry(f_expbg,"Expo Backgrnd","l");
	if(cfg.use_flat_bg == 1) leg->AddEntry(f_bg,"Bck","l");
	leg->AddEntry(f_p, Form("Decay: %.4fms", parent_halflife), "");
	leg->AddEntry(f_d, Form("Decay: %.4fms", daughter_halflife), "");
    	leg->AddEntry((TObject*)0, Form("#chi^{2}: %.4f",result -> Chi2()), "");

	leg->Draw();
	c -> Update();
	gr -> SetTitle("44S");
	gr->GetYaxis()->SetTitle(Form("Counts| %i ms per bin",cfg.rebin));
	gr->GetXaxis()->SetTitle("ms");
	gr->GetXaxis()->CenterTitle(true);  	
	gr->GetYaxis()->CenterTitle(true);  	
//	TGraph * grChi = ChiVsPar(h,fit,result,0,.001);// making a graph of chi2 vs lambda_p
//	grChi->SetName("ChiVsPar");
//	grChi->SetTitle("Delta Chi^2 vs Parameter;Parameter;#Delta#chi^{2}");

	// Save current stdout
	FILE* old_stdout = stdout;

	// Redirect stdout to file
	FILE* file_stdout = freopen((cfg.outfile + "Matrix.txt").c_str(), "w", stdout);

	// Call verbose print
	result->Print("V");

	// Flush and restore stdout
	fflush(stdout);
	stdout = old_stdout;
	
	TApplication app("app", &argc, argv);
	TFile* outf = new TFile((cfg.outfile + ".root").c_str(), "RECREATE");
	if (!f || f->IsZombie()) {
   	     std::cerr << "Error opening file\n";
	}
//	grChi -> Write("ChivsPar");
	c -> Write();
	gr -> Write();
	fit ->Write();
	f_p-> Write();
	f_d -> Write();
        if(cfg.use_expo_bg == 1)f_expbg -> Write();
	f_bnd -> Write();
	if(cfg.use_expo_bg == 1)f_bg -> Write();
	TBrowser * n = new TBrowser();
	app.Run();
	outf -> Close();
	return 0;
}
