#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"

//////////////////////////////////////////////////////////////
// Decay functions
//////////////////////////////////////////////////////////////

double decay_parent(double *dim, double *par)
{
    return par[0] * TMath::Exp(-par[1]*dim[0]);
}

double decay_model(double* dim, double* par)
{
    double x = dim[0];
    if(x < 0) return par[4];

    return decay_parent(dim,&par[0]) +
           decay_parent(dim,&par[2]) +
           par[4];
}

//////////////////////////////////////////////////////////////
// Stats box
//////////////////////////////////////////////////////////////

void DrawFitStatsBox(TF1* fit, double x1=0.35, double y1=0.65)
{
    double chi2 = fit->GetChisquare();
    double ndf  = fit->GetNDF();
    double chi2_ndf = (ndf > 0) ? chi2/ndf : 0;

    double lambda_p = fit->GetParameter(1);
    double lambda_p_err = fit->GetParError(1);

    double half_p = TMath::Log(2)/lambda_p;
    double half_p_err = (TMath::Log(2)/(lambda_p*lambda_p))*lambda_p_err;

    double lambda_b = fit->GetParameter(3);
    double lambda_b_err = fit->GetParError(3);

    double half_b = TMath::Log(2)/lambda_b;
    double half_b_err = (TMath::Log(2)/(lambda_b*lambda_b))*lambda_b_err;

    TPaveText *pave = new TPaveText(x1, y1, x1+0.25, y1+0.25, "NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->SetTextFont(42);
    pave->SetTextSize(0.03);

    pave->AddText(Form("#chi^{2}/NDF = %.1f / %d = %.2f", chi2, (int)ndf, chi2_ndf));
    pave->AddText(Form("t_{1/2}^{(p)} = %.2f #pm %.2f ms", half_p, half_p_err));
    pave->AddText(Form("t_{1/2}^{(b)} = %.2f #pm %.2f ms", half_b, half_b_err));

    pave->Draw();
}

//////////////////////////////////////////////////////////////
// Main macro
//////////////////////////////////////////////////////////////

void bckgrndsubfit2()
{
    int bin = 10;
    double Xmin = 1;
    double Xmax = 1000;

    TFile *file1 = new TFile("HistRootFiles/e21062_44S.root");

    if(!file1 || file1->IsZombie()){
        std::cout<<"File failed to open"<<std::endl;
        return;
    }

    TH2D *dtime = (TH2D*)file1->Get("dtime");
    if(!dtime){
        std::cout<<"Histogram not found"<<std::endl;
        return;
    }

    //////////////////////////////////////////////////////////
    // Projections
    //////////////////////////////////////////////////////////

    TH1D* gamma   = dtime->ProjectionX("gamma",327,330);
    TH1D* bckgrnd = dtime->ProjectionX("bckgrnd",333,336);

    gamma->Rebin(bin);
    bckgrnd->Rebin(bin);

    gamma->SetDirectory(0);
    bckgrnd->SetDirectory(0);

    //////////////////////////////////////////////////////////
    // Background subtraction
    //////////////////////////////////////////////////////////

    TH1D *sub = (TH1D*)gamma->Clone("subtracted");
    sub->Add(bckgrnd,-1);

    //////////////////////////////////////////////////////////
    // Set range
    //////////////////////////////////////////////////////////

    gamma->GetXaxis()->SetRangeUser(Xmin,Xmax);
    bckgrnd->GetXaxis()->SetRangeUser(Xmin,Xmax);
    sub->GetXaxis()->SetRangeUser(Xmin,Xmax);

    //////////////////////////////////////////////////////////
    // Fit functions
    //////////////////////////////////////////////////////////

    TF1 *parent = new TF1("parent", decay_model, Xmin, Xmax, 5);
    parent->SetParameters(1000,0.001,1000,0.0001,200);
    parent->SetParName(0,"Amplitude");
    parent->SetParName(1,"lambda_p");
    parent->SetParName(2,"bckgrnd source");
    parent->SetParName(3,"lambda_b");
    parent->SetParName(4,"Background");

    parent->SetParameters(1000,0.001,1000,.0001, 200);
    parent->SetParLimits(0,0,13000);
    parent->SetParLimits(1,.00001,.006);
    parent->SetParLimits(2,0,10000);
    parent->SetParLimits(3,.00001,1);
    parent->SetParLimits(4,0,20000);

    gamma->Fit(parent,"SR","",Xmin,Xmax);

    TF1 *parentsub = new TF1("parentsub", decay_model, Xmin, Xmax, 5);
    parentsub->SetParameters(1000,0.001,1000,0.0001,200);
    parentsub->SetParName(0,"Amplitude");
    parentsub->SetParName(1,"lambda_p");
    parentsub->SetParName(2,"bckgrnd source");
    parentsub->SetParName(3,"lambda_b");
    parentsub->SetParName(4,"Background");

    parentsub->SetParameters(1000,0.001,1000,.0001, 200);
    parentsub->SetParLimits(0,0,13000);
    parentsub->SetParLimits(1,.00001,.006);
    parentsub->SetParLimits(2,0,10000);
    parentsub->SetParLimits(3,.00001,1);
    parentsub->SetParLimits(4,0,20000);

    //////////////////////////////////////////////////////////
    // Fit
    //////////////////////////////////////////////////////////

    sub->Fit(parentsub,"SR","",Xmin,Xmax);

    //////////////////////////////////////////////////////////
    // Components (gamma)
    //////////////////////////////////////////////////////////

    //TF1* comp_parent = new TF1("comp_parent",
      //  [](double *x, double *par){ return par[0]*TMath::Exp(-par[1]*x[0]); },
       // Xmin, Xmax, 2);

    TF1* comp_b = new TF1("comp_b",
        [](double *x, double *par){ return par[0]*TMath::Exp(-par[1]*x[0]) + par[2]; },
        Xmin, Xmax, 3);

    TF1* comp_const = new TF1("comp_const",
        [](double *x, double *par){ return par[0]; },
        Xmin, Xmax, 1);

    //comp_parent->SetParameters(parent->GetParameter(0), parent->GetParameter(1));
    comp_b->SetParameters(parent->GetParameter(2), parent->GetParameter(3), parent->GetParameter(4));
    comp_const->SetParameter(0, parent->GetParameter(4));

    //comp_parent->SetLineColor(kBlue);
    //comp_parent->SetLineStyle(2);

    comp_b->SetLineColor(kGreen+2);
    comp_b->SetLineStyle(2);

    comp_const->SetLineColor(kMagenta);
    comp_const->SetLineStyle(3);

    //////////////////////////////////////////////////////////
    // Components (subtracted)
    //////////////////////////////////////////////////////////
	TF1* comp_b_sub = new TF1(
	    "comp_b_sub",
	    [](double *x, double *par){ return par[0]*TMath::Exp(-par[1]*x[0])+par[2]; },
	    Xmin, Xmax, 3);
	comp_b_sub->SetParameters(parentsub->GetParameter(2), parentsub->GetParameter(3),parentsub->GetParameter(4));
	comp_b_sub->SetLineColor(kGreen+2);
	comp_b_sub->SetLineStyle(2);

	TF1* comp_const_sub = new TF1(
	    "comp_const_sub",
	    [](double *x, double *par){ return par[0]; },
	    Xmin, Xmax, 1
	);
	comp_const_sub->SetParameter(0, parentsub->GetParameter(4));
	comp_const_sub->SetLineColor(kMagenta);
	comp_const_sub->SetLineStyle(3);

    //////////////////////////////////////////////////////////
    // Canvas
    //////////////////////////////////////////////////////////

    TCanvas *c = new TCanvas("c","Fits",800,1200);
    c->Divide(1,3);

    //////////////////////////////////////////////////////////
    // Pad 1: gamma
    //////////////////////////////////////////////////////////

    c->cd(1);
    //gPad->SetLogy();
	gamma->SetTitle("Gamma Gate 327-330");
	gamma->GetXaxis()->SetTitle("Time (ms)");
	gamma->GetYaxis()->SetTitle("Counts");
    gamma->Draw("E");
    parent->Draw("same");

    //comp_parent->Draw("same");
    comp_b->Draw("same");
    comp_const->Draw("same");

    DrawFitStatsBox(parent);

    TLegend* leg = new TLegend(0.15,0.15,0.40,0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.025);

    leg->AddEntry(gamma,"Data","lep");
    leg->AddEntry(parent,"Total Fit","l");
    //leg->AddEntry(comp_parent,"Parent","l");
    leg->AddEntry(comp_b,"expo background","l");
    leg->AddEntry(comp_const,"Background","l");
    leg->Draw();

    //////////////////////////////////////////////////////////
    // Pad 2: background
    //////////////////////////////////////////////////////////

    c->cd(2);
  //  gPad->SetLogy();
    bckgrnd->SetTitle("Background Gate 333-336");
    bckgrnd->GetXaxis()->SetTitle("Time (ms)");
    bckgrnd->GetYaxis()->SetTitle("Counts");
    bckgrnd->SetLineColor(kRed);
    bckgrnd->Draw("E");

    //////////////////////////////////////////////////////////
    // Pad 3: subtracted
    //////////////////////////////////////////////////////////

    c->cd(3);
//    gPad->SetLogy();
    sub->SetTitle("Background Subtracted");
    sub->GetXaxis()->SetTitle("Time (ms)");
    sub->GetYaxis()->SetTitle("Counts");
    sub->Draw("E");
    parentsub->Draw("same");

//    comp_parent_sub->Draw("same");
    comp_b_sub->Draw("same");
    comp_const_sub->Draw("same");

    DrawFitStatsBox(parentsub);

    //////////////////////////////////////////////////////////
    // Save
    //////////////////////////////////////////////////////////

    c->SaveAs("background_subtraction_fit.root");
}
