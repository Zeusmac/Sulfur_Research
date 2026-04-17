#include <sstream>
#include <fstream>


double chi2(TH1D* hist, TF1* eq){
     double chi2_sum = 0;
     for (int i = 1; i <= hist -> GetNbinsX(); ++i) {  // loop over bins 1..NbinsX
         double x = hist->GetBinCenter(i);
         double c = hist->GetBinContent(i);
         if (c <= 0)  continue;  // skip empty bins
         if (x < eq->GetXmin() || x > eq->GetXmax()) continue;
 
         double model = eq->Eval(x);
         double err2 = hist->GetBinError(i);
         if (err2 <= 0) continue;  // avoid division by zero
         chi2_sum += (model - c) * (model - c) / c;
     }
	return chi2_sum;
}


void FitGammaPeaksROOT() {
    // -------------------------------
    // Load histogram
    // -------------------------------
    TFile *file = new TFile("gamma_spec44S_main.root");   // change this
    TH1D *h = (TH1D*)file->Get("gamma_spec");        // change this
    if (!h) {
        std::cout << "Histogram not found!" << std::endl;
        return;
    }
    TGraphErrors * gr = new TGraphErrors(h); 
     h-> Sumw2(kTRUE);
//    h -> Rebin(2);
    TFitResultPtr result[50];	
    // -------------------------------
    // Peak search using TSpectrum
    // -------------------------------
    TSpectrum *s = new TSpectrum(50); // max peaks
    Int_t nPeaks = s->Search(h, 2, "", 0.05);

    std::cout << "Found " << nPeaks << " peaks\n";

    Double_t *xPeaks = s->GetPositionX();

    // -------------------------------
    // Draw histogram
    // -------------------------------
    TCanvas *c = new TCanvas("c", "Gamma Fits", 1000, 700);
    gr->Draw("AP");

    std::vector<TF1*> fits;

    // -------------------------------
    // Fit each peak
    // -------------------------------
    for (int i = 0; i < nPeaks; i++) {

        double peakX = xPeaks[i];
	if (h->GetBinContent(h->FindBin(peakX)) < 50) continue;
	std::sort(xPeaks, xPeaks + nPeaks);
        // Define fit window
	double sigma_guess = 2.0; // or from calibration

	double xmin = peakX - 2.3*sigma_guess;
	double xmax = peakX + 2.3*sigma_guess;
        // Create function: gaus + linear
        TF1 *f = new TF1(Form("fit_%d", i), "gaus(0) + pol1(3)", xmin, xmax);

        // Initial guesses
        int bin = h->FindBin(peakX);
        double height = h->GetBinContent(bin);
        f->SetParName(0, "a");
        f->SetParName(1, "x0");
        f->SetParName(2, "sigma");
        f->SetParName(3, "slope");
        f->SetParName(4, "intercept");
	
        f->SetParameters(
            height,     // amplitude
            peakX,      // mean
            2.0,        // sigma
            0,          // slope
            h->GetBinContent(bin) * 0.000001  // intercept
        );

        // Optional parameter limits (helps stability)
        f->SetParLimits(0, 0, 1e9);             // amplitude > 0
        f->SetParLimits(2, 0.2, 10);            // sigma reasonable
	f->SetParLimits(1, xmin,xmax);


        // Fit only in range
 	result[i] = gr->Fit(f, "R S M L+");

        fits.push_back(f);
    }
	
    // -------------------------------
    // Print results
    // -------------------------------
    ofstream out("gamma_peakfits.txt");
	 if(!out.is_open()){
         cerr<<"Cannot open config file\n";
         return;
     }	
 
    TF1* bg = new TF1("bkg", "pol1");
	
    for (size_t i = 0; i < fits.size(); i++) {
        TF1 *f = fits[i];

	bg -> SetParameters(f->GetParameter(3), f->GetParameter(4));
	double a = f -> GetParameter(0);
	double aerr = f-> GetParError(0);
	double sig = f -> GetParameter(2);
	double sigerr = f -> GetParError(2);
	double peak_count =  a * sig * TMath::Sqrt(2 * TMath::Pi());
	double peak_counterr = peak_count * (aerr/a + sig/sigerr);
	double x0 = f->GetParameter(1);
	double bck_counts = bg -> Integral(x0 - 2.5*2,x0 + 2.5*2);	
	double NDF = result[i] -> Ndf();
	double p_value = result[i] -> Prob();
	 
        out << "\nPeak " << i << ":\n";
	out << " p value =\t" << p_value << std::endl;
	out << " Amplitude =\t " << a << " ± " << aerr << std::endl; 
        out << " Mean   =\t" << x0
                  << " ± " << f->GetParError(1) << "\n";
        out << " Sigma  =\t" << sig
                  << " ± " << sigerr << "\n";
	out << " bg_slope  =\t" << f -> GetParameter(3) << " ± " << f -> GetParError(3) << std::endl;
	out << " bg_intercept  =\t" << f -> GetParameter(4) << " ± " <<f -> GetParError(4) << std::endl;
        out << "  Chi2/NDF =\t" << result[i]->Chi2() / NDF << "\n";
	out << " mChi2/NDF = \t" << chi2(h,f)/NDF << "\n";
	out << " FWHM =\t" << 2.355*sig 
			  << " ± " << 2.355*sigerr << std::endl;
	out << " Peak Counts =\t" << peak_count << " ± " << peak_counterr << std::endl; 
	out << " Peak Counts/total_count =\t" << peak_count/(h -> GetEntries()) << std::endl; 
	out << " signal/sqrt(noise) =\t" << (peak_count/TMath::Sqrt(bck_counts)) << "sigma" << std::endl; 
    }	
    out.close();
    c->Update();
}
