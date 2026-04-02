#include "DecayFitUtils.h"

#include "TApplication.h"
#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"

// ============================================
// Simple decay model
// ============================================
double decay_parent(double* x, double* par)
{
    return par[0] * exp(-par[1]*x[0]);
}

double decay_model(double* x, double* par,
                   bool use_parent, bool use_daughter)
{
    double val = par[4];
    if(x[0] < 0) return val;

    if(use_parent)   val += decay_parent(x,&par[0]);
    if(use_daughter) val += decay_parent(x,&par[2]);

    return val;
}

// ============================================
// MAIN
// ============================================
int main(int argc,char* argv[])
{
    if(argc < 2){
        std::cout<<"Usage: ./bkgfit config.txt\n";
        return 1;
    }

    FitConfig cfg;
    if(!ReadConfig(argv[1], cfg)) return 1;

    TApplication app("app",&argc,argv);

    // Load file
    TFile* file = TFile::Open(cfg.filename.c_str());
    TH2D* h2 = (TH2D*)file->Get(cfg.histname.c_str());

    // Projections (UTIL FUNCTIONS 🚀)
    TH1D* gamma = MakeProjectionX(h2,"gamma",
                                 cfg.gate1_low,cfg.gate1_high,
                                 cfg.rebin);

    TH1D* bg = MakeProjectionX(h2,"bg",
                              cfg.gate2_low,cfg.gate2_high,
                              cfg.rebin);

    TH1D* sub = SubtractHistograms(gamma,bg,"sub");

    gamma->GetXaxis()->SetRangeUser(cfg.xmin,cfg.xmax);
    sub->GetXaxis()->SetRangeUser(cfg.xmin,cfg.xmax);

    // Fit functions
    TF1* fit1 = new TF1("fit1",
        [&](double* x,double* p){
            return decay_model(x,p,cfg.use_parent,cfg.use_daughter);
        },
        cfg.xmin,cfg.xmax,cfg.init.size());

    TF1* fit2 = new TF1("fit2",
        [&](double* x,double* p){
            return decay_model(x,p,cfg.use_parent,cfg.use_daughter);
        },
        cfg.xmin,cfg.xmax,cfg.init.size());

    // Setup params (UTIL 🚀)
    SetupTF1(fit1, cfg.names, cfg.init, cfg.bounds);
    SetupTF1(fit2, cfg.names, cfg.init, cfg.bounds);

    // Fit (UTIL 🚀)
    auto r1 = FitHistogram(gamma, fit1, cfg.xmin, cfg.xmax);
    auto r2 = FitHistogram(sub,   fit2, cfg.xmin, cfg.xmax);

    // Output name
    std::string outname = "bkg_fit.txt";
    std::ofstream out(outname);

    PrintFitResultsAppend(fit1,r1,out,"Gamma Fit",
                          cfg.init,cfg.bounds);

    PrintFitResultsAppend(fit2,r2,out,"Subtracted Fit",
                          cfg.init,cfg.bounds);

    out.close();

    // Components (UTIL 🚀)
    TF1* p1 = MakeComponent("p1", decay_parent, fit1,
                           fit1->GetNpar(),
                           cfg.xmin,cfg.xmax,
                           kRed,2);

    // Canvas
    TCanvas* c = new TCanvas("c","Fits",800,1000);
    c->Divide(1,2);

    c->cd(1);
    gamma->Draw("E");
    fit1->Draw("same");
    p1->Draw("same");

    DrawFitStatsBox(fit1,{1,3});

    c->cd(2);
    sub->Draw("E");
    fit2->Draw("same");

    DrawFitStatsBox(fit2,{1,3});

    c->SaveAs("bkg_fit.root");

    app.Run();
    return 0;
}
