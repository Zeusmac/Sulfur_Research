#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

double decay_parent(double *dim, double *par)
{
  /* function for fitting decay curves 
  a is amplitude of source,
  b is tau of the source,
  bkgrd is backround,
  this first term is the exponetial decay 
  */
  double a = par[0];
  double b = par[1];
  double x = dim[0];

  return  a *  TMath::Exp(-(b) * (x));
}
double decay_model(double* dim, double* par)
{
  double x = dim[0];
  double bckgrnd = par[4];
  if( x < 0 ) return bckgrnd;
  return decay_parent(dim,&par[0]) + decay_parent(dim,&par[2]) + bckgrnd;
}
std::string FormatNumber(double x)
{
    std::ostringstream ss;

    double ax = fabs(x);

    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss << std::scientific << std::setprecision(6) << x;
    else
        ss << std::fixed << std::setprecision(6) << x;

    return ss.str();
}


void SaveFullFitReport(TH1* hist,
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
void bckgrndsubfit()
{
    int bin = 10;
    Double_t Xmin = 0;	

    TFile *file1 = new TFile("HistRootFiles/e21062_44S.root");

    if(!file1 || file1->IsZombie()){
        std::cout<<"File failed to open"<<std::endl;
        return;
    }

    TH2D *dtime = (TH2D*)file1->Get("dtimeN");

    if(!dtime){
        std::cout<<"Histogram decay_curve not found"<<std::endl;
        file1->ls();
        return;
    }

     //dtime -> Rebin(10);
    // projections
    TH1D* gamma   = dtime->ProjectionX("gamma",327,330);
    TH1D* bckgrnd   = dtime->ProjectionX("bckgrnd",333,336);

    bckgrnd -> Rebin(bin);
    gamma -> Rebin(bin);		
   gamma->SetDirectory(0);
   bckgrnd->SetDirectory(0);

    // background subtraction
    TH1D *sub = (TH1D*)gamma->Clone("subtracted");
    //sub -> Rebin(bin);	
    sub->Add(bckgrnd,-1);
    sub -> Draw("E");	
    // decay fit function
    double Xmax = 2000;//sub->GetXaxis()->GetXmax();

    TF1 *parent = new TF1(
        "parent",
        decay_model,
        Xmin,
        Xmax,
	5);
    parent->SetParName(0,"Amplitude");
    parent->SetParName(1,"lambda_p");
    parent->SetParName(2,"bckgrnd source");
    parent->SetParName(3,"lambda_b");
    parent->SetParName(4,"Background");

    parent->SetParameters(1000,0.001,1000,.0001,200);
    parent->SetParLimits(0,0,13000);
    parent->SetParLimits(1,.00001,.006);
    parent->SetParLimits(2,0,10000);
    parent->SetParLimits(3,.00001,1);
    parent->SetParLimits(4,0,20000);
    
    TF1 *parentsub = new TF1(
        "parentsub",
        decay_model,
        Xmin,
        Xmax,
	5);
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

    TFitResultPtr gresult = gamma->Fit(parent,"SR","", Xmin, Xmax);

    double init[16];
    for(int i=0;i<parent->GetNpar();i++)
	    init[i] = parent->GetParameter(i);

    SaveFullFitReport(gamma, parent, gresult, init, bin, "329gamma_fit_report.txt");
    
    TFitResultPtr subresult = sub->Fit(parentsub,"SR","",Xmin, Xmax);
    
    double initsub[16];
    for(int i=0;i<parentsub->GetNpar();i++)
	    initsub[i] = parentsub->GetParameter(i);

    SaveFullFitReport(sub, parentsub, subresult, init, bin, "329bckgrndSub_fit_report.txt");
    // half-life calculation
    double lambda = parent->GetParameter(1);
    double lambda_err = parent->GetParError(1);

    double half = TMath::Log(2)/lambda;
    double half_err = (TMath::Log(2)/(lambda*lambda))*lambda_err;

    printf("Half-life = %.4f ms ± %.4f ms\n",half,half_err);

    std::cout<<"Chi2/NDF = "
             <<parent->GetChisquare()<<" / "
             <<parent->GetNDF()<<" = "
             <<parent->GetChisquare()/parent->GetNDF()
             <<std::endl;

    // canvas with three plots
    TCanvas *c = new TCanvas("c","Background subtraction",600,1200);
    c->Divide(1,3);
    	
    c->cd(1);
    gamma->SetTitle("Gamma Gate 327-330");
    //gamma->GetXaxis() -> SetRange(1000,2000);
    gamma->Draw("E");
    parent->Draw("same");

    c->cd(2);
    bckgrnd->SetTitle("Background Gate 333-336");
    //bckgrnd->GetXaxis() -> SetRange(1000,2000);
    bckgrnd->SetLineColor(kRed);
    bckgrnd->Draw("E");

    c->cd(3);
    sub->SetTitle("Background Subtracted");
    //sub->GetXaxis() -> SetRange(1000,2000);
    sub->Draw("E");
    parentsub->Draw("same");

    c->SaveAs("background_subtraction_fit.root");
}
