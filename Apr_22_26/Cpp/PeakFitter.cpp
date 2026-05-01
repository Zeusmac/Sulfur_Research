#include "PeakFitter.h"

#include <algorithm>
#include <fstream>
#include <cmath>

#include "TSpectrum.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TDirectory.h"

using namespace std;

PeakFitter::PeakFitter(GammaDB& db,
                       PeakTracker& tracker,
                       TDirectory* sigmaDir,
                       TDirectory* peakDir)
    : db(db), tracker(tracker),
      sigmaDir(sigmaDir), peakDir(peakDir)
{}

// -----------------------------
bool PeakFitter::IsTimeHistogram(const string& name) {
    return (name.find("_") != string::npos &&
            isdigit(name[name.find("_") + 1]));
}

// -----------------------------
pair<double,double> PeakFitter::ExtractTimeRange(const string& name) {
    double t1 = 0, t2 = 0;
    sscanf(name.c_str(), "%*[^0-9]%lf_%lf", &t1, &t2);
    return {t1, t2};
}

// -----------------------------
void PeakFitter::FitHistogram(TH1* h, TFile* fout) {

    if (!h) return;

    TLatex latex;
    latex.SetTextSize(0.025);

    string hname = h->GetName();
     // Canvas
    TCanvas *c = new TCanvas(Form("c_%s", hname.c_str()), hname.c_str(), 1000, 700);

      // TEXT OUTPUT
    string txtname = "../Gamma_fits/" + hname + "_fit.txt";
    ofstream out(txtname);
     if (!out.is_open()) {
        cerr << "Cannot open output file\n";
        return;
    }

    TH1* h_work = (TH1*)h->Clone(Form("%s_work", hname.c_str()));
    h_work -> Sumw2();
    TSpectrum s(1000);

    map<string,int> isotope_counts;

    TH1 *bkg = s.Background(h_work, 14);
    int nPeaks = s.Search(h_work, 2, "", 0.02);
    h_work->Add(bkg, -1);

    double *xPeaks = s.GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);
    
    out << "====================================\n";
    out << "Histogram: " << hname << "\n";
    out << "Found " << nPeaks << " peaks\n";

    TF1 bg("bkg", "pol1");    

    vector<double> energies, sigmas, sigma_errs;

    for (int i = 0; i < nPeaks; i++) {

        double x = xPeaks[i];
        int bin = h_work->FindBin(x);

        if (h_work->GetBinContent(bin) < 50) continue;

        TF1 *f = new TF1(Form("f_%d", i), "gaus(0)+pol1(3)", x-6, x+6);

        f->SetParameters(h_work->GetBinContent(bin), x, 2.0, 0, 1);

        TFitResultPtr r = h_work->Fit(f, "R S Q L +");

        double a = f->GetParameter(0);
        double aerr = f->GetParError(0);
        double sig = f->GetParameter(2);
        double sigerr = f->GetParError(2);
        double x0 = f->GetParameter(1);

        double peak_count = a * sig * sqrt(2*TMath::Pi());
        double peak_counterr = peak_count * sqrt(
            pow(aerr/a,2) + pow(sigerr/sig,2));

        double bck_counts = bg.Integral(x0 - 2.5*sig, x0 + 2.5*sig);
        double SNR = peak_count / sqrt(bck_counts);
        double chi2ndf = r->Chi2() / r->Ndf();
        double p_value = r->Prob();

        //skips bad peaks
        if (h->GetBinContent(bin) < 50) continue;
        if (sig >= 500) continue;
        if (SNR < 10.0) continue;
        // -----------------------------
        // Matching
        // -----------------------------
        auto matches = db.Match(x0, 5);
        // -----------------------------
        // Output text
        // -----------------------------
        out << "\nPeak " << i << ":\n";
        out << " p value =\t" << p_value << endl;

        string label = Form("%.1f", x0);

	    out << "Fit Status   : " << r->Status() << "\n";
	    out << "EDM          : " << r->Edm() << "\n";
        out << " Amplitude =\t " << a << " ± " << aerr << endl;
        out << " Mean   =\t" << x0 << " ± " << f->GetParError(1) << endl;
	    out << " Possible matches:\n";

	    if (matches.empty()) 
        {
    		out << "   No isotope match within tolerance\n";
		} else 
        {
    	    for (auto &m : matches) {

        	    out << "   Isotope: " << m.isotope
            	    << " | Line: " << m.energy << " keV"
            	    << " | ΔE = " << m.deltaE << " keV\n";
    	}

	        // best match (top of sorted list)
	        out << " Best match: "
		    << matches[0].isotope
		    << " (" << matches[0].energy << " keV)\n";
	    }
        out << " Sigma  =\t" << sig << " ± " << sigerr << endl;
        out << " bg_slope  =\t" << f->GetParameter(3) << endl;
        out << " bg_intercept  =\t" << f->GetParameter(4) << endl;
        out << " Chi2/NDF =\t" << chi2ndf << endl;
        out << " FWHM =\t" << 2.355*sig << endl;
        out << " Peak Counts =\t" << peak_count << " ± " << peak_counterr << endl;
        out << " signal/sqrt(noise) =\t" << SNR << endl;
 	
        // -----------------------------
        // Draw label
        // -----------------------------
        double y = f->Eval(x0);
        latex.DrawLatex(x0, y * 1.1, label.c_str());

        energies.push_back(x0);
        sigmas.push_back(sig);
        sigma_errs.push_back(sigerr);

        double time = -1;
        if (IsTimeHistogram(hname)) {
            auto [t1,t2] = ExtractTimeRange(hname);
            time = 0.5*(t1+t2);
        }

        if (time >= 0)
            tracker.Add(x0, time, peak_count, 0);
    }

    // -----------------------------
    // Isotope confirmation
    // -----------------------------
    out << "\n=== Isotope confirmation ===\n";

    for (auto &p : isotope_counts) {
        if (p.second >= 2) {
            out << "CONFIRMED: " << p.first
                << " (" << p.second << " lines)\n";
        } else {
            out << "Possible: " << p.first << endl;
        }
    }
    out.close();

    // -------- Sigma vs Energy plot --------
    TGraphErrors* siggr = new TGraphErrors(energies.size());

    for (size_t i = 0; i < energies.size(); i++) {
        siggr->SetPoint(i, energies[i], sigmas[i]);
        siggr->SetPointError(i, 0, sigma_errs[i]);
    }

    siggr->SetName(Form("sigma_vs_energy_%s", hname.c_str()));
    siggr->SetTitle("Sigma vs Energy;Energy (keV);Sigma (keV)");
    siggr->SetMarkerStyle(20);

     //sigma vs energy dir
    this->sigmaDir->cd();
    siggr->Write();
    fout -> cd();

    this->peakDir->cd();
    //---------- Peak count vs time -----------
    for (auto &kv : tracker.GetTimes()) {

        double energy = kv.first;
        auto &t = kv.second;
        auto &c = tracker.GetCounts().at(energy);
        auto &e = tracker.GetErrors().at(energy);

        int n = t.size();
        if (n < 2) continue;

        TGraphErrors *peakgr = new TGraphErrors(n);

        for (int i = 0; i < n; i++) {
            peakgr->SetPoint(i, t[i], c[i]);
            peakgr->SetPointError(i, 0, e[i]);
        }

        peakgr->SetName(Form("peak_%d_keV", (int)energy));
        peakgr->Write();
    }

    //Writing to out file
    fout->cd();
    c->Write();   // canvas with labels + fits
    h->Write();   // histogram itself

    c->Update();
}
