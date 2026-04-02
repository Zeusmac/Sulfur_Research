#include "DecayFitUtils.h"

#include "TApplication.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"

FitConfig *gCfg = nullptr;

// ============================================
// Your full Bateman model (UNCHANGED)
// ============================================
Double_t TotalModelFull(double *x, double *par)
{
    double t = x[0];
    if(t < 0) return par[7];

    int use_beta = gCfg->use_parent;
    int use_bn   = gCfg->use_daughter;

    double lambda_p = par[0];
    double lambda_d = par[1];
    double N0       = par[6];
    double bg       = par[7];

    double rate = N0 * exp(-lambda_p*t);

    if(use_beta){
        double denom = (lambda_d - lambda_p);
        if(fabs(denom)<1e-12) denom=1e-12;

        rate += N0*(lambda_p/denom)*
                (exp(-lambda_p*t)-exp(-lambda_d*t));
    }

    return rate + bg;
}

// ============================================
// MAIN
// ============================================
int main(int argc,char* argv[])
{
    if(argc < 2){
        std::cout<<"Usage: ./fit config.txt\n";
        return 1;
    }

    FitConfig cfg;
    gCfg = &cfg;

    if(!ReadConfig(argv[1], cfg)) return 1;

    TApplication app("app",&argc,argv);

    // Load histogram
    TFile *f = TFile::Open(cfg.filename.c_str());
    TH1D* h = (TH1D*)f->Get(cfg.histname.c_str());

    h->Rebin(cfg.rebin);
    h->GetXaxis()->SetRangeUser(cfg.xmin,cfg.xmax);

    // Fit function
    TF1 *fit = new TF1("fit",TotalModelFull,
                       cfg.xmin,cfg.xmax,cfg.init.size());

    SetupTF1(fit, cfg.names, cfg.init, cfg.bounds);

    // Fit
    auto result = FitHistogram(h, fit, cfg.xmin, cfg.xmax);

    // Save results
    std::ofstream out("fit_results.txt");
    PrintFitResultsAppend(fit,result,out,"Bateman Fit",
                          cfg.init,cfg.bounds);
    out.close();

    // Components (example)
    TF1* parent = MakeComponent("parent", TotalModelFull,
                               fit, fit->GetNpar(),
                               cfg.xmin,cfg.xmax,
                               kRed,2);

    // Canvas
    TCanvas *c = new TCanvas("c","Fit",800,600);

    h->Draw("E");
    fit->Draw("same");
    parent->Draw("same");

    DrawFitStatsBox(fit,{0,1}); // lambda indices

    c->SaveAs("bateman_fit.root");

    app.Run();
    return 0;
}
