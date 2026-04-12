#include <sstream>
#include <fstream>


void FitGammaPeaksROOT() {
    // -------------------------------
    // Load histogram
    // -------------------------------
    TFile *file = new TFile("gamma_spec_main_1000ms.root");   // change this
    TH1D *h = (TH1D*)file->Get("Gamma_spec_main_1000ms");        // change this
    if (!h) {
        std::cout << "Histogram not found!" << std::endl;
        return;
    }

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
    h->Draw();

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

	double xmin = peakX - 2.5*sigma_guess;
	double xmax = peakX + 2.5*sigma_guess;
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
            h->GetBinContent(bin) * 0.0000001  // intercept
        );

        // Optional parameter limits (helps stability)
        f->SetParLimits(0, 0, 1e9);             // amplitude > 0
        f->SetParLimits(2, 0.2, 10);            // sigma reasonable
	f->SetParLimits(1, xmin,xmax);


        // Fit only in range
        h->Fit(f, "R L+");

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
	double bck_counts = bg -> Integral(f->GetParameter(1) - 2.5*2,f->GetParameter(1) + 2.5*2);	
	 
        out << "\nPeak " << i << ":\n";
	out << " Amplitude  = " << a << " ± " << aerr << std::endl; 
        out << "  Mean   = " << f->GetParameter(1)
                  << " ± " << f->GetParError(1) << "\n";
        out << "  Sigma  = " << sig
                  << " ± " << sigerr << "\n";
	out << " bg_slope  = " << f -> GetParameter(3) << " ± " << f -> GetParError(3) << std::endl;
	out << " bg_intercept  = " << f -> GetParameter(4) << f -> GetParError(4) << std::endl;
        out << "  Chi2/NDF = "
                  << f->GetChisquare() / f->GetNDF() << "\n";
	out << " FWHM	= " << 2.355*sig 
			  << " ± " << 2.355*sigerr << std::endl;
	out << " Peak Counts = " << peak_count << " ± " << peak_counterr << std::endl; 
	out << " Peak Counts/total_count = " << peak_count/(h -> GetEntries()) << " ± " << peak_counterr << std::endl; 
	out << " signal/sqrt(noise) = " << peak_count/(bck_counts) << std::endl; 
    }	
    out.close();
    c->Update();
}
