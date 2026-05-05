#include "ResolutionAnalyzer.h"

#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

ResolutionAnalyzer::ResolutionAnalyzer(PeakTracker& tracker)
    : tracker(tracker)
{}

// ----------------------------------
void ResolutionAnalyzer::BuildFWHMvsTime(TDirectory* sigmaDir) {

    auto& times  = tracker.GetTimes();
    auto& sigmas = tracker.GetSigmas();     
    auto& errors = tracker.GetSigmaErrors();

    for (auto& kv : times) {

        double energy = kv.first;

        auto& t = kv.second;
        auto& s = sigmas.at(energy);
        auto& e = errors.at(energy);

        int n = t.size();
        if (n < 2) continue;

        TGraphErrors* g = new TGraphErrors(n);

        for (int i = 0; i < n; i++) {
            double fwhm = 2.355 * s[i];
            double fwhm_err = 2.355 * e[i];

            g->SetPoint(i, t[i], fwhm);
            g->SetPointError(i, 0, fwhm_err);
        }
        g->Sort();
        g->SetName(Form("FWHM_vs_time_%d_keV", (int)energy));
        g->SetTitle(Form("FWHM vs Time (%.0f keV);Time (ms);FWHM (keV)", energy));
        g->SetMarkerStyle(20);

        // -------- FIT --------
        TF1* f = new TF1(Form("fit_time_%d", (int)energy), "[0] + [1]*x", 0, 1e6);
        f->SetParameters(1, 0);

        g->Fit(f, "Q");  // quiet fit

        // -------- WRITE SAFELY --------
        sigmaDir->WriteTObject(g);
        sigmaDir->WriteTObject(f);
    }
}