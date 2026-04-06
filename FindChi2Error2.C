double chi2(TH1D* hist, TF1* eq){
     double chi2_sum = 0;
     for (int i = 1; i <= hist -> GetNbinsX(); ++i) {  // loop over bins 1..NbinsX
         double x = hist->GetBinCenter(i);
         double c = hist->GetBinContent(i);
         if (c <= 0) c = 1.0;  // skip empty bins
         if (x < eq->GetXmin() || x > eq->GetXmax()) continue;
 
         double model = eq->Eval(x);
         double err2 = hist->GetBinError(i);
         if (err2 <= 0) continue;  // avoid division by zero
         chi2_sum += (model - c) * (model - c) / c;
     }
	return chi2_sum;
}

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
  double bkgrd = par[2];
  double x = dim[0];

  if(x < 0) return bkgrd; 
  return  a * TMath::Exp(-(b) * (x)) + bkgrd;
}

void FindChi2Error2()
{
    int lamp = 0.005999;
    TFile* f = TFile::Open("saved_histograms.root");
    if (!f || f->IsZombie()) return;

    TH1D* hist      = (TH1D*) f->Get("gamma");
    TH1D* hist2      = (TH1D*) f->Get("bckgrnd");
    TH1D* hist3      = (TH1D*) f->Get("subtracted");
    if (!hist || hist->IsZombie()){
	cout << "hist no existy" << endl; 
	return;
	}
    TF1* decay      = hist3->GetFunction("fit2");
    if (!decay) {
        std::cout << "No fit function found!\n";
        return;
    }
    TF1 * eq = new TF1("fit",decay_parent,-1000,4000, 3);
    int j = 0;
    int n_params = decay->GetNpar();
    for (int i = 0; i < n_params; ++i) {
        double par = decay->GetParameter(i);
        std::cout << "Parameter " << i << " (" << decay->GetParName(i) << ") = " << par << std::endl;
	if(i == 2 || i == 3) continue;
	eq -> SetParameter(j, par);
	j++;
    }
    eq -> SetParameter(1,lamp);
    double chi2_sum = 0;
    for (int i = 1; i <= hist -> GetNbinsX(); ++i) {  // loop over bins 1..NbinsX
        double x = hist->GetBinCenter(i);
        double c = hist->GetBinContent(i);
	double c2 = hist2->GetBinContent(i);
	double c3 = hist3-> GetBinContent(i);
        if (c <= 0) c = 1.0;
	if (c2 <= 0) c2 = 1.0;  // skip empty bins
        if (x < eq->GetXmin() || x > eq->GetXmax()) continue;

        double model = eq->Eval(x);
        double err2 = c2 + c;
        if (err2 <= 0) continue;  // avoid division by zero
        chi2_sum += (model - c3) * (model - c3) / err2;
    }

    std::cout << "Chi2 = " << chi2_sum/165 << std::endl;
}
