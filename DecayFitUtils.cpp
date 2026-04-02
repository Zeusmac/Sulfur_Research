#include "DecayFitUtils.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "TMath.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFitResult.h"
using namespace std;

// ==============================
// Utility
// ==============================
string FormatNumber(double x){
    ostringstream ss;
    double ax = fabs(x);
    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss << scientific << setprecision(6) << x;
    else
        ss << fixed << setprecision(6) << x;
    return ss.str();
}

double SafeExp(double x){
    if(x < -700) return 0;
    if(x > 700) return exp(700);
    return exp(x);
}

double LambdaToHalfLife(double lambda){
    if(lambda <= 0) return 0;
    return log(2)/lambda;
}

// ==============================
// Config
// ==============================
bool ReadConfig(const string& fname, FitConfig &cfg)
{
    ifstream file(fname);
    if(!file.is_open()) return false;

    string key;
    while(file >> key)
    {
        if(key=="FILE") file >> cfg.filename;
        else if(key=="HIST") file >> cfg.histname;
        else if(key=="FIT_RANGE") file >> cfg.xmin >> cfg.xmax;
        else if(key=="REBIN") file >> cfg.rebin;

        else if(key=="GATE1") file >> cfg.gate1_low >> cfg.gate1_high;
        else if(key=="GATE2") file >> cfg.gate2_low >> cfg.gate2_high;

        else if(key[0]=='#') file.ignore(1000,'\n');

        else {
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

// ==============================
// TF1 setup
// ==============================
void SetupTF1(TF1* f,
              const vector<string>& names,
              const vector<double>& init,
              const vector<pair<double,double>>& bounds)
{
    for(size_t i=0;i<init.size();i++){
        f->SetParName(i,names[i].c_str());
        f->SetParameter(i,init[i]);
        f->SetParLimits(i,bounds[i].first,bounds[i].second);
    }
}

void CopyParameters(TF1* src, TF1* dest)
{
    for(int i=0;i<src->GetNpar();i++)
        dest->SetParameter(i, src->GetParameter(i));
}

// ==============================
// Histogram helpers
// ==============================
TH1D* MakeProjectionX(TH2D* h2, const string& name,
                      int ylow, int yhigh, int rebin)
{
    TH1D* h = h2->ProjectionX(name.c_str(), ylow, yhigh);
    h->Sumw2();
    h->Rebin(rebin);
    return h;
}

TH1D* SubtractHistograms(TH1D* h1, TH1D* h2, const string& name)
{
    TH1D* h = (TH1D*)h1->Clone(name.c_str());
    h->Add(h2,-1);
    return h;
}

// ==============================
// Fit helper
// ==============================
TFitResultPtr FitHistogram(TH1D* h, TF1* f,
                           double xmin, double xmax,
                           const string& options)
{
    return h->Fit(f, options.c_str(), "", xmin, xmax);
}

// ==============================
// Correlation
// ==============================
pair<int,int> FindMaxCorrelation(TF1* fit, TFitResultPtr r)
{
    int p1=0,p2=1;
    double maxCorr=0;

    for(int i=0;i<fit->GetNpar();i++){
        for(int j=i+1;j<fit->GetNpar();j++){
            double c = fabs(r->Correlation(i,j));
            if(c > maxCorr){
                maxCorr = c;
                p1=i; p2=j;
            }
        }
    }
    return {p1,p2};
}

// ==============================
// Results printing
// ==============================
void PrintFitResultsAppend(TF1* fit, TFitResultPtr r, std::ofstream &file,
                          const std::string &label,
                          const std::vector<double> &init,
                          const std::vector<std::pair<double,double>> &bounds)
{
    file << "\n============================================================\n";
    file << "                     " << label << "\n";
    file << "============================================================\n\n";

    // =============================
    // Fit summary
    // =============================
    file << "Chi2 / NDF   : " << FormatNumber(r->Chi2())
         << " / " << r->Ndf()
         << "  =  " << FormatNumber(r->Chi2()/r->Ndf()) << "\n";
    file << "EDM          : " << FormatNumber(r->Edm()) << "\n";
    file << "NCalls       : " << r->NCalls() << "\n";
    file << "Probability  : " << FormatNumber(r->Prob()) << "\n\n";

    int fitStatus = r->Status();
    int covStatus = r->CovMatrixStatus();

    file << "Fit Status   : " << fitStatus
         << ((fitStatus==0) ? " (OK)" : " (PROBLEM)") << "\n";

    file << "Cov Matrix   : " << covStatus
         << ((covStatus==3) ? " (FULL ACCURATE)" :
             (covStatus==2) ? " (FORCED POS-DEF)" :
             (covStatus==1) ? " (APPROXIMATE)" :
                              " (NOT VALID)") << "\n\n";

    if(fitStatus != 0 || covStatus < 3){
        file << "WARNING: Fit may be unreliable\n\n";
    }

    // =============================
    // Parameter table
    // =============================
    file << std::left
         << std::setw(16) << "Parameter" << "\t"
         << std::setw(14) << "Value" << "\t"
         << std::setw(24) << "Error (- / +)" << "\t"
         << std::setw(14) << "Initial" << "\t"
         << std::setw(20) << "Bounds" << "\t"
         << "Status\n";

    file << std::string(90,'-') << "\n";

    for(int i=0;i<fit->GetNpar();i++){
        bool fixed = r->IsParameterFixed(i);

        std::stringstream bnd;
        bnd << "[" << bounds[i].first << ", " << bounds[i].second << "]";

        std::string errStr;

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
                std::stringstream ss;
                ss << "-" << FormatNumber(fabs(el))
                   << "  +" << FormatNumber(eh);
                errStr = ss.str();
            }
        }

        file << std::setw(16) << fit->GetParName(i) << "\t"
             << std::setw(14) << FormatNumber(r->Parameter(i)) << "\t"
             << std::setw(24) << errStr << "\t"
             << std::setw(14) << FormatNumber(init[i]) << "\t"
             << std::setw(20) << bnd.str() << "\t"
             << (fixed ? "fixed" : "free") << "\n";
    }

    // =============================
    // Half-lives (ASYMMETRIC)
    // =============================
    file << "\n------------------ Half-lives ------------------\n";

    for(int i=0;i<fit->GetNpar();i++){
        std::string name = fit->GetParName(i);

        if(name.find("lambda") != std::string::npos){
            double lam = r->Parameter(i);

            if(lam <= 0) continue;

            double el = r->LowerError(i);
            double eh = r->UpperError(i);

            double t = log(2)/lam;

            std::string errStr;

            if(el == 0 && eh == 0){
                double e = r->ParError(i);
                double terr = (e/lam)*t;
                errStr = "± " + FormatNumber(terr);
            }
            else{
                // Propagate asymmetric errors properly
                double t_low  = log(2)/(lam + eh);  // upper lambda → lower t
                double t_high = log(2)/(lam + el);  // lower lambda → upper t

                double err_low  = t - t_low;
                double err_high = t_high - t;

                std::stringstream ss;
                ss << "-" << FormatNumber(fabs(err_low))
                   << "  +" << FormatNumber(err_high);

                errStr = ss.str();
            }

            file << std::setw(16) << name
                 << std::setw(14) << FormatNumber(t)
                 << errStr << "\n";
        }
    }

    // =============================
    // Covariance matrix (NEW)
    // =============================
    int n = fit->GetNpar();

    file << "\n---------------- Covariance Matrix ----------------\n\n";

    file << std::setw(12) << "";
    for(int j=0;j<n;j++)
        file << std::setw(12) << fit->GetParName(j);
    file << "\n";

    for(int i=0;i<n;i++){
        file << std::setw(12) << fit->GetParName(i);
        for(int j=0;j<n;j++){
            file << std::setw(12)
                 << FormatNumber(r->CovMatrix(i,j));
        }
        file << "\n";
    }

    // =============================
    // Correlation matrix
    // =============================
    file << "\n---------------- Correlation Matrix ----------------\n\n";

    file << std::setw(12) << "";
    for(int j=0;j<n;j++)
        file << std::setw(12) << fit->GetParName(j);
    file << "\n";

    for(int i=0;i<n;i++){
        file << std::setw(12) << fit->GetParName(i);
        for(int j=0;j<n;j++){
            file << std::setw(12)
                 << FormatNumber(r->Correlation(i,j));
        }
        file << "\n";
    }

    file << "\n============================================================\n";
}
// ==============================
// Stats box (GENERALIZED)
// ==============================
void DrawFitStatsBox(TF1* fit,
                     const vector<int>& lambda_indices,
                     double x1, double y1)
{
    TPaveText *pave = new TPaveText(x1,y1,x1+0.25,y1+0.25,"NDC");
    pave->SetFillColor(0);

    double chi2_ndf = fit->GetChisquare()/fit->GetNDF();
    pave->AddText(Form("#chi^{2}/NDF = %.2f",chi2_ndf));

    for(int idx : lambda_indices){
        double lam = fit->GetParameter(idx);
        double lam_err = fit->GetParError(idx);

        double t = LambdaToHalfLife(lam);
        double terr = (lam_err/lam)*t;

        pave->AddText(Form("t_{1/2} = %.2f ± %.2f", t, terr));
    }

    pave->Draw();
}

// ==============================
// Component builder
// ==============================
TF1* MakeComponent(const string& name,
                   Double_t (*func)(double*, double*),
                   TF1* fit,
                   int npar,
                   double xmin, double xmax,
                   int color, int style)
{
    TF1* f = new TF1(name.c_str(), func, xmin, xmax, npar);
    CopyParameters(fit,f);
    f->SetLineColor(color);
    f->SetLineStyle(style);
    return f;
}
