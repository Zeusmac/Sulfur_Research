void FindChi2Error2()
{
    TFile* f = TFile::Open("saved_histograms.root");
    if (!f || f->IsZombie()) return;

    TH1D* hist      = (TH1D*) f->Get("sub");
    if (!f || hist->IsZombie()){
	cout << "hist no existy" << endl; 
	return;
	}
    TF1* decay      = hist->GetFunction("fit2");
    if (!decay) {
        std::cout << "No fit function found!\n";
        return;
    }
	TF1 * eq = new TF1("fit", "[0] * exp(-[1]*x) + [2]",0,4000);
	int j = 0;
    int n_params = decay->GetNpar();
    for (int i = 0; i < n_params; ++i) {
        double par = decay->GetParameter(i);
        std::cout << "Parameter " << i << " (" << decay->GetParName(i) << ") = " << par << std::endl;
	if(i == 2 || i == 3) continue;
	eq -> SetParameter(j, par);
	j++;
    }

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

    std::cout << "Chi2 = " << chi2_sum << std::endl;
}
