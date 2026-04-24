#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TSpectrum.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;

// -----------------------------
// Structs
// -----------------------------
struct GammaLine {
    string isotope;
    double energy;
};

struct PeakInfo {
    double x0;
    double sigma;
    double counts;
};

// -----------------------------
vector<GammaLine> db;
vector<PeakInfo> peaks;

// -----------------------------
// Load DB
// -----------------------------
vector<GammaLine> LoadGammaDB(const string& filename) {
    vector<GammaLine> out;
    ifstream in(filename);

    string iso;
    double energy;

    while (in >> iso >> energy) {
        out.push_back({iso, energy});
    }
    return out;
}

// -----------------------------
// Match peaks
// -----------------------------
vector<GammaLine> MatchPeak(double x0, double sigma) {

    vector<GammaLine> matches;
    double tol = 2.0 * sigma;

    for (auto &line : db) {
        if (fabs(line.energy - x0) < tol) {
            matches.push_back(line);
        }
    }
    return matches;
}

// -----------------------------
// Interactive click
// -----------------------------
void HandleClick(Int_t event, Int_t px, Int_t py, TObject*) {

    if (event != kButton1Down) return;

    double x = gPad->AbsPixeltoX(px);
    x = gPad->PadtoX(x);

    double best_dist = 1e9;
    PeakInfo best;

    for (auto &p : peaks) {
        double d = fabs(p.x0 - x);
        if (d < best_dist) {
            best_dist = d;
            best = p;
        }
    }

    if (best_dist < 5) {
        cout << "\nClicked peak:\n";
        cout << " Energy = " << best.x0 << "\n";
        cout << " Counts = " << best.counts << "\n";

        auto matches = MatchPeak(best.x0, best.sigma);

        cout << " Possible isotopes:\n";
        for (auto &m : matches) {
            cout << "   " << m.isotope << " (" << m.energy << " keV)\n";
        }
    }
}

// -----------------------------
// MAIN
// -----------------------------
void FitGammaPeaksROOT() {

    db = LoadGammaDB("gamma_db.txt");

    TFile *file = new TFile("gamma_spec44S_main.root");
    TH1D *h = (TH1D*)file->Get("gamma_spec");

    if (!h) {
        cout << "Histogram not found!\n";
        return;
    }

    h->Sumw2();

    TGraphErrors *gr = new TGraphErrors(h);

    TSpectrum *s = new TSpectrum(50);
    int nPeaks = s->Search(h, 2, "", 0.05);

    cout << "Found " << nPeaks << " peaks\n";

    double *xPeaks = s->GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);

    TCanvas *c = new TCanvas("c", "Gamma Fits", 1000, 700);
    gr->Draw("AP");

    TLatex latex;
    latex.SetTextSize(0.025);

    vector<TF1*> fits;
    TFitResultPtr results[100];

    ofstream out("gamma_peakfits.txt");
    if (!out.is_open()) {
        cerr << "Cannot open output file\n";
        return;
    }

    TF1* bg = new TF1("bkg", "pol1");

    map<string,int> isotope_counts;

    // -----------------------------
    // Fit loop
    // -----------------------------
    for (int i = 0; i < nPeaks; i++) {

        double peakX = xPeaks[i];
        if (h->GetBinContent(h->FindBin(peakX)) < 50) continue;

        double xmin = peakX - 5;
        double xmax = peakX + 5;

        TF1 *f = new TF1(Form("fit_%d", i), "gaus(0)+pol1(3)", xmin, xmax);

        double height = h->GetBinContent(h->FindBin(peakX));

        f->SetParameters(height, peakX, 2.0, 0, 1);
        f->SetParLimits(0, 0, 1e9);
        f->SetParLimits(2, 0.2, 10);
        f->SetParLimits(1, xmin, xmax);

        results[i] = gr->Fit(f, "R S Q");

        fits.push_back(f);

        double a = f->GetParameter(0);
        double aerr = f->GetParError(0);
        double sig = f->GetParameter(2);
        double sigerr = f->GetParError(2);
        double x0 = f->GetParameter(1);

        double peak_count = a * sig * sqrt(2*TMath::Pi());

        // FIXED ERROR PROPAGATION
        double peak_counterr = peak_count * sqrt(
            pow(aerr/a,2) + pow(sigerr/sig,2)
        );

        bg->SetParameters(f->GetParameter(3), f->GetParameter(4));

        double bck_counts = bg->Integral(x0 - 2.5*sig, x0 + 2.5*sig);

        double NDF = results[i]->Ndf();
        double p_value = results[i]->Prob();

        peaks.push_back({x0, sig, peak_count});

        // -----------------------------
        // Matching
        // -----------------------------
        auto matches = MatchPeak(x0, sig);

        string label = Form("%.1f", x0);

        out << "\nPeak " << i << ":\n";
        out << " p value =\t" << p_value << endl;

        for (auto &m : matches) {
            label += " ";
            label += m.isotope;

            out << " Match: " << m.isotope
                << " (" << m.energy << " keV)\n";

            isotope_counts[m.isotope]++;
        }

        out << " Amplitude =\t " << a << " ± " << aerr << endl;
        out << " Mean   =\t" << x0 << " ± " << f->GetParError(1) << endl;
        out << " Sigma  =\t" << sig << " ± " << sigerr << endl;
        out << " bg_slope  =\t" << f->GetParameter(3) << endl;
        out << " bg_intercept  =\t" << f->GetParameter(4) << endl;
        out << " Chi2/NDF =\t" << results[i]->Chi2() / NDF << endl;
        out << " FWHM =\t" << 2.355*sig << endl;
        out << " Peak Counts =\t" << peak_count << " ± " << peak_counterr << endl;
        out << " signal/sqrt(noise) =\t"
            << peak_count / sqrt(bck_counts) << endl;

        // -----------------------------
        // Draw label
        // -----------------------------
        double y = f->Eval(x0);
        latex.DrawLatex(x0, y*1.1, label.c_str());
    }

    // -----------------------------
    // Multi-peak confirmation
    // -----------------------------
    out << "\n=== Isotope confirmation ===\n";
    cout << "\n=== Isotope confirmation ===\n";

    for (auto &p : isotope_counts) {
        if (p.second >= 2) {
            out << "CONFIRMED: " << p.first
                << " (" << p.second << " lines)\n";

            cout << "CONFIRMED: " << p.first
                 << " (" << p.second << " lines)\n";
        }
        else {
            out << "Possible: " << p.first << endl;
        }
    }

    out.close();

    // -----------------------------
    // Enable clicking
    // -----------------------------
    gPad->AddExec("ex", "HandleClick(event,px,py,obj);");

    c->Update();
}
