#include "PeakFitter.h"

#include <algorithm>
#include <fstream>
#include <cmath>

#include "TSpectrum.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

PeakFitter::PeakFitter(GammaDB& db_, PeakTracker& tracker_)
    : db(db_), tracker(tracker_) {}

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

    string hname = h->GetName();

    TH1* h_work = (TH1*)h->Clone(Form("%s_work", hname.c_str()));

    TSpectrum s(1000);
    TH1 *bkg = s.Background(h_work, 10);
    h_work->Add(bkg, -1);

    int nPeaks = s.Search(h_work, 2, "", 0.02);
    double *xPeaks = s.GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);

    vector<double> energies, sigmas, sigma_errs;

    for (int i = 0; i < nPeaks; i++) {

        double x = xPeaks[i];
        int bin = h_work->FindBin(x);

        if (h_work->GetBinContent(bin) < 50) continue;

        TF1 *f = new TF1(Form("f_%d", i), "gaus(0)+pol1(3)", x-6, x+6);

        f->SetParameters(h_work->GetBinContent(bin), x, 2.0, 0, 1);

        TFitResultPtr r = h_work->Fit(f, "R S Q L 0");

        double sig = f->GetParameter(2);
        double sigerr = f->GetParError(2);
        double x0 = f->GetParameter(1);

        double a = f->GetParameter(0);
        double peak_counts = a * sig * sqrt(2*TMath::Pi());

        if (sig > 6) continue;

        energies.push_back(x0);
        sigmas.push_back(sig);
        sigma_errs.push_back(sigerr);

        double time = -1;
        if (IsTimeHistogram(hname)) {
            auto [t1,t2] = ExtractTimeRange(hname);
            time = 0.5*(t1+t2);
        }

        if (time >= 0)
            tracker.Add(x0, time, peak_counts, 0);
    }

    // -------- Sigma vs Energy plot --------
    TGraphErrors* g = new TGraphErrors(energies.size());

    for (size_t i = 0; i < energies.size(); i++) {
        g->SetPoint(i, energies[i], sigmas[i]);
        g->SetPointError(i, 0, sigma_errs[i]);
    }

    g->SetName(Form("sigma_vs_energy_%s", hname.c_str()));
    g->SetTitle("Sigma vs Energy;Energy (keV);Sigma (keV)");
    g->SetMarkerStyle(20);

    fout->cd();
    g->Write();
}