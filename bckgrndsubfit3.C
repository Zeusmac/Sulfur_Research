#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>


using namespace std;

//////////////////////////////////////////////////////////////
// Config struct
//////////////////////////////////////////////////////////////

struct FitConfig
{
    string filename;
    string histname;

    double Xmin, Xmax;
    int gate1_low, gate1_high;
    int gate2_low, gate2_high;
    int bin;

    vector<string> names1, names2;
    vector<double> init1, init2;
    vector<pair<double,double>> bounds1, bounds2;

    bool use_parent1 = true;
    bool use_daughter1 = true;
    bool use_parent2 = true;
    bool use_daughter2 = true;
};

//////////////////////////////////////////////////////////////
// Read config
//////////////////////////////////////////////////////////////

bool ReadConfig(const char* fname, FitConfig &cfg)
{
    ifstream file(fname);
    if(!file.is_open()){
        cout<<"Could not open config file: "<<fname<<"\n";
        return false;
    }

    string key;
    while(file >> key)
    {
        if(key=="#") file.ignore(numeric_limits<streamsize>::max(), '\n');

        else if(key=="FILE") file >> cfg.filename;
        else if(key=="HIST") file >> cfg.histname;
        else if(key=="FIT_RANGE") file >> cfg.Xmin >> cfg.Xmax;
        else if(key=="GATE1") file >> cfg.gate1_low >> cfg.gate1_high;
        else if(key=="GATE2") file >> cfg.gate2_low >> cfg.gate2_high;
        else if(key=="BIN") file >> cfg.bin;

        else if(key=="INIT1") {
            cfg.init1.clear();
            double v; for(int i=0;i<5;i++){ file>>v; cfg.init1.push_back(v);}
        }
        else if(key=="INIT2") {
            cfg.init2.clear();
            double v; for(int i=0;i<5;i++){ file>>v; cfg.init2.push_back(v);}
        }

        else if(key=="BOUNDS1") {
            cfg.bounds1.clear();
            double minv,maxv;
            for(int i=0;i<5;i++){ file>>minv>>maxv; cfg.bounds1.push_back({minv,maxv});}
        }
        else if(key=="BOUNDS2") {
            cfg.bounds2.clear();
            double minv,maxv;
            for(int i=0;i<5;i++){ file>>minv>>maxv; cfg.bounds2.push_back({minv,maxv});}
        }

        else if(key=="NAMES1") {
            cfg.names1.clear();
            string s; for(int i=0;i<5;i++){ file>>s; cfg.names1.push_back(s);}
        }
        else if(key=="NAMES2") {
            cfg.names2.clear();
            string s; for(int i=0;i<5;i++){ file>>s; cfg.names2.push_back(s);}
        }

        else if(key=="USE_PARENT1") file >> cfg.use_parent1;
        else if(key=="USE_EXPBCK1") file >> cfg.use_daughter1;
        else if(key=="USE_PARENT2") file >> cfg.use_parent2;
        else if(key=="USE_EXPBCK2") file >> cfg.use_daughter2;
    }

    return true;
}

//////////////////////////////////////////////////////////////
// Filename generator
//////////////////////////////////////////////////////////////

string MakeOutputName(const FitConfig& cfg)
{
    ostringstream name;
    name << "fit_"
         << "g" << cfg.gate1_low << "-" << cfg.gate1_high
         << "_bg" << cfg.gate2_low << "-" << cfg.gate2_high
         << "_bin" << cfg.bin
         << ".txt";
    return name.str();
}

//////////////////////////////////////////////////////////////
// Decay functions
//////////////////////////////////////////////////////////////

double decay_parent(double* x, double* par)
{
    return par[0] * TMath::Exp(-par[1]*x[0]);
}

double decay_model(double* x, double* par, bool use_parent, bool use_daughter)
{ 
    double val = par[4]; // Background
    if(x[0] < 0 ) return par[4];
    if(use_parent) val += decay_parent(x,&par[0]);
    if(use_daughter) val += decay_parent(x,&par[2]);
    return val;
}
double ParentComponent(double* x, double* par)
{
    double t = x[0];
    double lambda_p = par[1];
    double N0       = par[0];
    double bg       = par[4];
    double Np = N0*exp(-lambda_p*t);
    return Np + bg;
}

//////////////////////////////////////////////////////////////
// Formatter
//////////////////////////////////////////////////////////////

string FormatNumber(double x)
{
    ostringstream ss;
    double ax = fabs(x);
    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss<<scientific<<setprecision(6)<<x;
    else
        ss<<fixed<<setprecision(6)<<x;
    return ss.str();
}

//////////////////////////////////////////////////////////////
// Print fit results
//////////////////////////////////////////////////////////////

void PrintFitResultsAppend(TF1* fit, TFitResultPtr r, ofstream &file,
                          const string &label,
                          const vector<double> &init,
                          const vector<pair<double,double>> &bounds)
{
    file << "\n================ " << label << " =================\n\n";

    file << "Chi2 / NDF   : " << FormatNumber(r->Chi2())
         << " / " << r->Ndf()
         << "  =  " << FormatNumber(r->Chi2()/r->Ndf()) << "\n";
    file << "EDM          : " << FormatNumber(r->Edm()) << "\n";
    file << "NCalls       : " << r->NCalls() << "\n";
    file << "Probability  : " << FormatNumber(r->Prob()) << "\n\n";

    int fitStatus = r->Status();
    int covStatus = r->CovMatrixStatus();
    file << "Fit Status   : " << fitStatus << ((fitStatus==0)?" (OK)":" (PROBLEM)") << "\n";
    file << "Cov Matrix   : " << covStatus
         << ((covStatus==3)?" (FULL ACCURATE)":
             (covStatus==2)?" (FORCED POS-DEF)":
             (covStatus==1)?" (APPROXIMATE)":" (NOT VALID)") << "\n\n";

    // Parameters
    file << left << setw(16) << "Parameter" << "\t"
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
        if(fixed) errStr = "---";
        else{
            double el = r->LowerError(i);
            double eh = r->UpperError(i);
            if(el==0 && eh==0){
                double e = r->ParError(i);
                errStr = "± "+FormatNumber(e);
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
             << (fixed?"fixed":"free") << "\n";
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
    // Correlation matrix
    int n = fit->GetNpar();
    file << "\n---- Correlation Matrix ----\n";
    file << setw(12) << "";
    for(int j=0;j<n;j++) file << setw(12) << fit->GetParName(j);
    file << "\n";
    for(int i=0;i<n;i++){
        file << setw(12) << fit->GetParName(i);
        for(int j=0;j<n;j++)
            file << setw(12) << FormatNumber(r->Correlation(i,j));
        file << "\n";
    }
    file << "========================================\n";
}

//////////////////////////////////////////////////////////////
// Stats box
//////////////////////////////////////////////////////////////

void DrawFitStatsBox(TF1* fit, bool use_parent, bool use_daughter, double x1=0.35, double y1=0.65)
{
    double chi2 = fit->GetChisquare();
    double ndf  = fit->GetNDF();
    double chi2_ndf = (ndf>0)? chi2/ndf : 0;

    double half_p=0, half_p_err=0;
    double half_b=0, half_b_err=0;

    if(use_parent){
        double lam = fit->GetParameter(1);
        double lam_err = fit->GetParError(1);
        half_p = TMath::Log(2)/lam;
        half_p_err = (TMath::Log(2)/(lam*lam))*lam_err;
    }
    if(use_daughter){
        double lam = fit->GetParameter(3);
        double lam_err = fit->GetParError(3);
        half_b = TMath::Log(2)/lam;
        half_b_err = (TMath::Log(2)/(lam*lam))*lam_err;
    }

    TPaveText *pave = new TPaveText(x1,y1,x1+0.25,y1+0.25,"NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->SetTextFont(42);
    pave->SetTextSize(0.06);

    pave->AddText(Form("#chi^{2}/NDF = %.2f",chi2_ndf));
    if(use_parent) pave->AddText(Form("t_{1/2}^{(p)} = %.2f #pm %.2f", half_p, half_p_err));
    if(use_daughter) pave->AddText(Form("t_{1/2}^{(b)} = %.2f #pm %.2f", half_b, half_b_err));

    pave->Draw();
}

//////////////////////////////////////////////////////////////
// Main function
//////////////////////////////////////////////////////////////

void bckgrndsubfit3(const char* configfile="bckgrndsub_config.txt")
{
    FitConfig cfg;
    if(!ReadConfig(configfile,cfg)) return;

    // Open file
    TFile* file1 = new TFile(cfg.filename.c_str());
    TH2D* dtime = (TH2D*)file1->Get(cfg.histname.c_str());

    // Projections
    TH1D* gamma   = dtime->ProjectionX("gamma",cfg.gate1_low,cfg.gate1_high);
    TH1D* bckgrnd = dtime->ProjectionX("bckgrnd",cfg.gate2_low,cfg.gate2_high);

    gamma->Sumw2(); bckgrnd->Sumw2();
    gamma->Rebin(cfg.bin); bckgrnd->Rebin(cfg.bin);

    TH1D* sub = (TH1D*)gamma->Clone("subtracted");
    sub->Add(bckgrnd,-1);

    gamma->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);
    bckgrnd->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);
    sub->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);

    // --- Fits ---
    TF1* fit1 = new TF1("fit1",
        [&](double* x, double* p){ return decay_model(x,p,cfg.use_parent1,cfg.use_daughter1); },
        cfg.Xmin,cfg.Xmax,5);

    TF1* fit2 = new TF1("fit2",[&](double* x, double* p){ return decay_model(x,p,cfg.use_parent1,cfg.use_daughter1); },
         cfg.Xmin,cfg.Xmax,5);


    // Parameter names
    for(int i=0;i<5;i++){
        fit1->SetParName(i,cfg.names1[i].c_str());
        fit2->SetParName(i,cfg.names2[i].c_str());

        fit1->SetParameter(i,cfg.init1[i]);
        fit2->SetParameter(i,cfg.init2[i]);

        fit1->SetParLimits(i,cfg.bounds1[i].first,cfg.bounds1[i].second);
        fit2->SetParLimits(i,cfg.bounds2[i].first,cfg.bounds2[i].second);
    }

    TFitResultPtr r1 = gamma->Fit(fit1,"S R L E");
    TFitResultPtr r2 = sub->Fit(fit2,"S R L E");
    	
    // --- Components ---
    TF1* p1 = new TF1("p1", ParentComponent, cfg.Xmin,cfg.Xmax,3);
    TF1* b1 = new TF1("b1", ParentComponent, cfg.Xmin,cfg.Xmax,3);
    TF1* c1 = new TF1("c1","[0]", cfg.Xmin,cfg.Xmax);

    p1->SetParameters(fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(4));
    b1->SetParameters(fit1->GetParameter(2), fit1->GetParameter(3), fit1->GetParameter(4));
    c1->SetParameter(0, fit1->GetParameter(4));

    p1->SetLineColor(kRed);   p1->SetLineStyle(2);
    b1->SetLineColor(kBlue);  b1->SetLineStyle(2);
    c1->SetLineColor(kGreen+2); c1->SetLineStyle(2);

    TF1* p2 = new TF1("p2", ParentComponent, cfg.Xmin,cfg.Xmax,3);
    TF1* b2 = new TF1("b2", ParentComponent, cfg.Xmin,cfg.Xmax,3);
    TF1* c2 = new TF1("c2","[0]", cfg.Xmin,cfg.Xmax);

    p2->SetParameters(fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(4));
    b2->SetParameters(fit2->GetParameter(2), fit2->GetParameter(3), fit2->GetParameter(4));
    c2->SetParameter(0, fit2->GetParameter(4));

    p2->SetLineColor(kRed);   p2->SetLineStyle(2);
    b2->SetLineColor(kBlue);  b2->SetLineStyle(2);
    c2->SetLineColor(kGreen+2); c2->SetLineStyle(2);

    // --- Titles ---
    gamma->SetTitle(Form("Gamma Gate %d-%d",cfg.gate1_low,cfg.gate1_high));
    bckgrnd->SetTitle(Form("Background Gate %d-%d",cfg.gate2_low,cfg.gate2_high));
    sub->SetTitle("Background Subtracted");

    // --- Save results ---
    string outname = MakeOutputName(cfg);
    ofstream out(outname);
    cout<<"Saving results to: "<<outname<<endl;

    PrintFitResultsAppend(fit1,r1,out,"Gamma Fit",cfg.init1,cfg.bounds1);
    PrintFitResultsAppend(fit2,r2,out,"Subtracted Fit",cfg.init2,cfg.bounds2);

    out.close();

    // --- Canvas ---
    TCanvas* c = new TCanvas("c","Fits",800,1200);
    c->Divide(1,3);

    c->cd(1);
    TLegend* leg = new TLegend(0.55,0.6,0.88,0.88);
    leg->AddEntry(gamma,"Data","lep");
    leg->AddEntry(fit1,"Total Fit","l");
    if(cfg.use_parent1) leg->AddEntry(p1,"Parent","l");
    if(cfg.use_daughter1) leg->AddEntry(b1,"Daughter","l");
    leg->AddEntry(c1,"Constant Bkg","l");
    leg->Draw();

    gamma->Draw("E");
    fit1->Draw("same"); 
    if(cfg.use_parent1) p1->Draw("same");
    if(cfg.use_daughter1) b1->Draw("same");

    DrawFitStatsBox(fit1,cfg.use_parent1,cfg.use_daughter1);

    c->cd(2);
    bckgrnd->Draw("E");

    c->cd(3);
    sub->Draw("E");
    //fit2->Draw("same");
    if(cfg.use_parent2) p2->Draw("same");
   // if(cfg.use_daughter2) b2->Draw("same");
    c2->Draw("same");

    DrawFitStatsBox(fit2,cfg.use_parent2,cfg.use_daughter2);

    string rootfile = outname.substr(0,outname.length()-4)+".root";
    TFile * outfile = new TFile(rootfile.c_str(),"Recreate");
    gamma -> Write();
    bckgrnd -> Write();
    sub -> Write();
    c->Write();
    outfile -> Close();
    cout<<"Canvas saved as: "<<rootfile<<endl;
}
