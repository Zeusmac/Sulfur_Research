double decay_parent(double* x, double* par)
{
    return par[0] * TMath::Exp(-par[1]*x[0]);
}

double decay_model(double* x, double* par, bool use_parent, bool use_daughter)
{ 
    double val = par[4]; // Background
    if(x[0] < 0 ) return par[4];
    val += decay_parent(x,&par[0]);
    val += decay_parent(x,&par[2]);
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

void FindChi2Error()
{
	TFile* f = TFile::Open("saved_histograms.root");
	if (!f || f->IsZombie()) return;
	TH1D* hist      = (TH1D*) f->Get("gamma");
	TH1D* bckgrnd    = (TH1D*) f->Get("bckgrnd");
	TH1D* subtracted = (TH1D*) f->Get("subtracted");

	double MaxX = 4000; //hist -> GetXaxis() -> GetXmax();
	double MinX = hist -> GetXaxis() -> GetXmin();
	TF1 * decay = (TF1*) hist -> GetFunction("fit1");
  Int_t n_params = decay->GetNpar();
  for (int i = 0; i < n_params; ++i) 
  {
    Double_t par_min, par_max;
    Double_t par = decay->GetParameter(i);
    std::cout << "Parameter " << i <<  decay->GetParName(i) <<" :" << par  << std::endl;
  }

	double chi2_sum = 0;
	for (int i = MinX; i <= 4980; ++i) 
	{
		double x = hist->GetBinCenter(i);
		Double_t c = hist->GetBinContent(i);
		if (c > 0)  
		{
			Double_t model =  decay -> Eval(x);
			double err2 = hist->GetBinError(i);
			chi2_sum += (((model - c)*(model - c))/(err2));
		}
	}
	cout << chi2_sum << endl;
}
