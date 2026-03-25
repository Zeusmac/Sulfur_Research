#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>


//////////////////////////////////////////////////////////////
// Config struct
//////////////////////////////////////////////////////////////

struct FitConfig
{
    std::string filename;
    std::string histname;

    double Xmin, Xmax;

    int gate1_low, gate1_high;
    int gate2_low, gate2_high;

    int bin;

    double init1[5];
    double init2[5];

    double bounds1[10];
    double bounds2[10];
};

//////////////////////////////////////////////////////////////
// Read config
//////////////////////////////////////////////////////////////

bool ReadConfig(const char* fname, FitConfig &cfg)
{
    std::ifstream file(fname);
    if(!file.is_open()){
        std::cout<<"Could not open config file\n";
        return false;
    }

    std::string key;

    while(file >> key)
    {   
	if(key=="#") continue;
        if(key=="FILE") file >> cfg.filename;
        else if(key=="HIST") file >> cfg.histname;
        else if(key=="FIT_RANGE") file >> cfg.Xmin >> cfg.Xmax;
        else if(key=="GATE1") file >> cfg.gate1_low >> cfg.gate1_high;
        else if(key=="GATE2") file >> cfg.gate2_low >> cfg.gate2_high;
        else if(key=="BIN") file >> cfg.bin;

        else if(key=="INIT1")
            for(int i=0;i<5;i++) file >> cfg.init1[i];

        else if(key=="INIT2")
            for(int i=0;i<5;i++) file >> cfg.init2[i];

        else if(key=="BOUNDS1")
            for(int i=0;i<10;i++) file >> cfg.bounds1[i];

        else if(key=="BOUNDS2")
            for(int i=0;i<10;i++) file >> cfg.bounds2[i];
    }

    return true;
}

//////////////////////////////////////////////////////////////
// Filename generator
//////////////////////////////////////////////////////////////

std::string MakeOutputName(const FitConfig& cfg)
{
    std::ostringstream name;

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
// Formatter
//////////////////////////////////////////////////////////////

std::string FormatNumber(double x)
{
    std::ostringstream ss;
    double ax = fabs(x);

    if(ax > 1e4 || (ax < 1e-3 && ax > 0))
        ss<<std::scientific<<std::setprecision(6)<<x;
    else
        ss<<std::fixed<<std::setprecision(6)<<x;

    return ss.str();
}

//////////////////////////////////////////////////////////////
// Print results (enhanced)
//////////////////////////////////////////////////////////////
void PrintFitResultsAppend(TF1* fit, TFitResultPtr r, std::ofstream &file, const std::string &label,
                           const double init[5], const double bounds[10], int bin,
                           int gate1_low, int gate1_high, int gate2_low, int gate2_high,
                           double fitLow, double fitHigh)
{
    file<<"\n========== "<<label<<" ==========\n\n";

    // General info
    file<<"Bin number: "<<bin<<"\n";
    file<<"Gate1 range: ["<<gate1_low<<","<<gate1_high<<"]\n";
    file<<"Gate2 range: ["<<gate2_low<<","<<gate2_high<<"]\n";
    file<<"Fit range: ["<<fitLow<<","<<fitHigh<<"]\n\n";

    file<<"Chi2   = "<<FormatNumber(r->Chi2())<<"\n";
    file<<"NDF    = "<<r->Ndf()<<"\n";
    file<<"EDM    = "<<FormatNumber(r->Edm())<<"\n";
    file<<"NCalls = "<<r->NCalls()<<"\n\n";

    file<<std::left
        <<std::setw(18)<<"Parameter"
        <<std::setw(18)<<"Value"
        <<std::setw(18)<<"Error"
        <<std::setw(18)<<"Initial"
        <<std::setw(18)<<"Bounds"
        <<"Status\n";

    for(int i=0;i<fit->GetNpar();i++)
    {
        bool fixed = r->IsParameterFixed(i);
        std::stringstream bnd;
        bnd<<"["<<bounds[2*i]<<","<<bounds[2*i+1]<<"]";

        file<<std::setw(18)<<fit->GetParName(i)
            <<std::setw(18)<<FormatNumber(r->Parameter(i))
            <<std::setw(18)<<(fixed ? "---" : FormatNumber(r->ParError(i)))
            <<std::setw(18)<<FormatNumber(init[i])
            <<std::setw(18)<<bnd.str()
            <<(fixed ? "fixed" : "limited")<<"\n";
    }

    // Half-lives
    file<<"\nHalf-lives (ms)\n";
    for(int i=0;i<fit->GetNpar();i++)
    {
        std::string name = fit->GetParName(i);
        if(name.find("lambda") != std::string::npos)
        {
            double lam = r->Parameter(i);
            double err = r->ParError(i);
            double t = log(2)/lam;
            double terr = (err/lam)*t;
            file<<name<<": "<<FormatNumber(t)<<" ± "<<FormatNumber(terr)<<"\n";
        }
    }

    file<<"Chi2/NDF = "<<FormatNumber(r->Chi2()/r->Ndf())<<"\n\n";

    // Correlation matrix
    int n = fit->GetNpar();
    file<<"Correlation Matrix:\n";
    file<<std::setw(12)<<"";
    for(int j=0;j<n;j++) file<<std::setw(12)<<fit->GetParName(j);
    file<<"\n";

    for(int i=0;i<n;i++)
    {
        file<<std::setw(12)<<fit->GetParName(i);
        for(int j=0;j<n;j++)
        {
            file<<std::setw(12)<<FormatNumber(r->Correlation(i,j));
        }
        file<<"\n";
    }
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
    pave->SetTextSize(0.06);

    pave->AddText(Form("#chi^{2}/NDF = %.2f", chi2_ndf));
    pave->AddText(Form("t_{1/2}^{(p)} = %.2f #pm %.2f", half_p, half_p_err));
    pave->AddText(Form("t_{1/2}^{(b)} = %.2f #pm %.2f", half_b, half_b_err));

    pave->Draw();
}

//////////////////////////////////////////////////////////////
// MAIN
//////////////////////////////////////////////////////////////

void bckgrndfit3(const char* configfile="bckgrndsub_congfig.txt")
{
    FitConfig cfg;

    if(!ReadConfig(configfile,cfg)) return;

    TFile *file1 = new TFile(cfg.filename.c_str());
    TH2D *dtime = (TH2D*)file1->Get(cfg.histname.c_str());

    //////////////////////////////////////////////////////////
    // Projections
    //////////////////////////////////////////////////////////

    TH1D* gamma   = dtime->ProjectionX("gamma",cfg.gate1_low,cfg.gate1_high);
    TH1D* bckgrnd = dtime->ProjectionX("bckgrnd",cfg.gate2_low,cfg.gate2_high);
  
   bckgrnd->Sumw2();	
   gamma -> Sumw2();
    gamma->Rebin(cfg.bin);
    bckgrnd->Rebin(cfg.bin);
    TH1D *sub = (TH1D*)gamma->Clone("subtracted");
    
    sub->Add(bckgrnd,-1);
    bckgrnd->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);
    gamma->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);
    sub->GetXaxis()->SetRangeUser(cfg.Xmin,cfg.Xmax);

    //////////////////////////////////////////////////////////
    // Fits
    //////////////////////////////////////////////////////////

    TF1 *fit1 = new TF1("fit1", decay_model, cfg.Xmin, cfg.Xmax, 5);
    TF1 *fit2 = new TF1("fit2", "[0]*exp(-[1]*x) + [2]", cfg.Xmin, cfg.Xmax, 3);

    fit1->SetParNames("Amplitude","lambda_p","bkg src","lambda_b","Background");
    fit2->SetParNames("Amplitude","lambda_p","bkg src","lambda_b","Background");

    for(int i=0;i<5;i++)
    {
        fit1->SetParameter(i,cfg.init1[i]);
        fit2->SetParameter(i,cfg.init2[i]);

        fit1->SetParLimits(i,cfg.bounds1[2*i],cfg.bounds1[2*i+1]);
        fit2->SetParLimits(i,cfg.bounds2[2*i],cfg.bounds2[2*i+1]);
    }

    TFitResultPtr r2 = bckgrnd->Fit(fit2,"SR");
     

	    // -------- Components for fit1 (gamma) --------
	TF1 *p1 = new TF1("p1", decay_parent, cfg.Xmin, cfg.Xmax, 2);
	TF1 *b1 = new TF1("b1", decay_parent, cfg.Xmin, cfg.Xmax, 2);
	TF1 *c1 = new TF1("c1", "[0]",  cfg.Xmin, cfg.Xmax);

	p1->SetParameters(fit1->GetParameter(0),fit1->GetParameter(1));
	b1->SetParameters(fit1->GetParameter(2), fit1->GetParameter(3));
	c1->SetParameters(fit1->GetParameter(4));

	// styles
	p1->SetLineColor(kRed);
	p1->SetLineStyle(2);

	b1->SetLineColor(kBlue);
	b1->SetLineStyle(2);

	c1->SetLineColor(kGreen+2);
	c1->SetLineStyle(2);


	// -------- Components for fit2 (subtracted) --------
	TF1 *c2 = new TF1("c2","[0]",  cfg.Xmin, cfg.Xmax);

	c2->SetParameters(fit2->GetParameter(4));

	c2->SetLineColor(kGreen+2);
	c2->SetLineStyle(2);    

    //////////////////////////////////////////////////////////
    // Titles
    //////////////////////////////////////////////////////////

    gamma->SetTitle(Form("Gamma Gate %d-%d",cfg.gate1_low,cfg.gate1_high));
    bckgrnd->SetTitle(Form("Background Gate %d-%d",cfg.gate2_low,cfg.gate2_high));
    sub->SetTitle("Background Subtracted");
	// -----------------------------
	// Save results
	// -----------------------------
	std::string outname = MakeOutputName(cfg);
	std::ofstream out(outname);

	std::cout<<"Saving results to: "<<outname<<std::endl;

	// Use enhanced function
	PrintFitResultsAppend(fit1, r1, out, "Gamma Fit",
			      cfg.init1, cfg.bounds1, cfg.bin,
			      cfg.gate1_low, cfg.gate1_high,
			      cfg.gate2_low, cfg.gate2_high,
			      cfg.Xmin, cfg.Xmax);

	PrintFitResultsAppend(fit2, r2, out, "Subtracted Fit",
			      cfg.init2, cfg.bounds2, cfg.bin,
			      cfg.gate1_low, cfg.gate1_high,
			      cfg.gate2_low, cfg.gate2_high,
			      cfg.Xmin, cfg.Xmax);

	out.close();
    //////////////////////////////////////////////////////////
    // Canvas
    //////////////////////////////////////////////////////////

    TCanvas *c = new TCanvas("c","Fits",800,1200);
        c->Divide(1,2);

	c->cd(1);
	TLegend *leg = new TLegend(0.55,0.6,0.88,0.88);
	leg->AddEntry(gamma,"Data","lep");
	leg->AddEntry(fit1,"Total Fit","l");
	leg->AddEntry(p1,"Parent","l");
	leg->AddEntry(b1,"Bkg Source","l");
	leg->AddEntry(c1,"Constant Bkg","l");
	leg->Draw();
	gamma->Draw("E");
	fit1->Draw("same");

	p1->Draw("same");
	b1->Draw("same");
	c1->Draw("same");

	DrawFitStatsBox(fit1);
	c->cd(2);
	bckgrnd->Draw("E");

	DrawFitStatsBox(fit2);
    std::string rootfile = outname.substr(0,outname.length() - 4) + ".root" ;
    std::string pngfile = outname.substr(0,outname.length() - 4) + ".png" ;
    c->SaveAs(rootfile.c_str());
    c->SaveAs(pngfile.c_str());
    
}
