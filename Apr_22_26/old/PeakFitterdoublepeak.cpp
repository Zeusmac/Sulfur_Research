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
TF1* MakeGaussianModel(double x0, double sigmaGuess, double width, int mode, int i) {

    if (mode == 0) {
        // single Gaussian
        return new TF1(Form("g1_%d", i),
                       "gaus(0)+pol1(3)",
                       x0 - width, x0 + width);
    }

    if (mode == 1) {
        // double Gaussian (broader peaks / pileup)
        return new TF1(Form("g2_%d", i),
                       "gaus(0)+gaus(3)+pol1(6)",
                       x0 - width, x0 + width);
    }

    return nullptr;
}

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
PeakFitter::FitQuality EvaluateFit(TF1* f, TFitResultPtr& r,
                                   double sig, double snr) {

    PeakFitter::FitQuality q;

    q.chi2ndf = r->Chi2() / std::max(1.0, (double)r->Ndf());
    q.status  = r->Status();
    q.sigma   = sig;
    q.snr     = snr;

    return q;
}
vector<vector<double>> ClusterPeaks(double* xPeaks, int nPeaks, double mergeDistance = 3.0) {

    vector<vector<double>> clusters;

    if (nPeaks == 0) return clusters;

    vector<double> current;
    current.push_back(xPeaks[0]);

    for (int i = 1; i < nPeaks; i++) {

        if (fabs(xPeaks[i] - xPeaks[i-1]) < mergeDistance) {
            current.push_back(xPeaks[i]);
        } else {
            clusters.push_back(current);
            current.clear();
            current.push_back(xPeaks[i]);
        }
    }

    clusters.push_back(current);
    return clusters;
}
void PeakFitter::FitHistogram(TH1* h, TFile* fout) {

    if (!h) return;

    string hname = h->GetName();

    TCanvas *c = new TCanvas(Form("c_%s", hname.c_str()),
                             hname.c_str(), 1000, 700);

    h->Draw("HIST");

    TLatex latex;
    latex.SetTextSize(0.025);

    ofstream out("../Gamma_fits/" + hname + "_fit.txt");

    if (!out.is_open()) return;

    TH1* h_work = (TH1*)h->Clone(Form("%s_work", hname.c_str()));
    h_work->Sumw2();

    TSpectrum s(1000);

    TH1 *bkg = s.Background(h_work, 14);
    int nPeaks = s.Search(h_work, 2, "", 0.02);

    h_work->Add(bkg, -1);

    double *xPeaks = s.GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);

    out << "Histogram: " << hname << "\n";
    out << "Raw TSpectrum peaks: " << nPeaks << "\n";

    // -------------------------------
    // STEP 1: FILTER FAKE PEAKS
    // -------------------------------
    vector<double> cleanPeaks;

    for (int i = 0; i < nPeaks; i++) {

        double x = xPeaks[i];
        int bin = h_work->FindBin(x);

        double content = h_work->GetBinContent(bin);

        if (content < 20) continue;            // noise rejection
        if (x < 5 || x > h->GetXaxis()->GetXmax() - 5) continue;

        cleanPeaks.push_back(x);
    }

    sort(cleanPeaks.begin(), cleanPeaks.end());

    // -------------------------------
    // STEP 2: CLUSTER PEAKS
    // -------------------------------
    vector<vector<double>> clusters;

    double clusterWindow = 10.0; // keV (adaptive could be % energy)

    for (double x : cleanPeaks) {

        bool added = false;

        for (auto &cl : clusters) {
            if (fabs(cl.back() - x) < clusterWindow) {
                cl.push_back(x);
                added = true;
                break;
            }
        }

        if (!added) clusters.push_back({x});
    }

    out << "Clusters formed: " << clusters.size() << "\n";

    vector<double> energies, sigmas, sigma_errs;

    TF1 globalBg("bg", "pol1");

    // -------------------------------
    // STEP 3: FIT EACH CLUSTER
    // -------------------------------
    for (size_t i = 0; i < clusters.size(); i++) {

        auto &cl = clusters[i];

        double xmin = cl.front() - 5;
        double xmax = cl.back() + 5;

        int n = cl.size();

        // ---------------------------
        // MODEL SELECTION
        // ---------------------------
        string model;

        if (n == 1) model = "gaus(0)+pol1(3)";
        else if (n == 2) model = "gaus(0)+gaus(3)+pol1(6)";
        else if (n == 3) model = "gaus(0)+gaus(3)+gaus(6)+pol1(9)";
        else  model = "gaus(0)+gaus(3)+gaus(6)+gaus(9)+pol1(12)";

        TF1 *f = new TF1(Form("fit_%s_%zu", hname.c_str(), i),
                         model.c_str(), xmin, xmax);

        // ---------------------------
        // INITIAL PARAMS
        // ---------------------------
        int p = 0;
        for (double x : cl) {
            int bin = h_work->FindBin(x);
            double amp = h_work->GetBinContent(bin);

            f->SetParameter(p, amp);
            f->SetParameter(p + 1, x);
            f->SetParameter(p + 2, 2.0);
            p += 3;
        }

        // background init
        f->SetParameter(p, 0);
        f->SetParameter(p + 1, 0);
        f->SetParameter(p + 2, 0);

        // -------------------------------
        // STEP 4: ITERATIVE FIT REFINEMENT
        // -------------------------------
        TFitResultPtr r;

        double prevChi2 = 1e9;

        for (int iter = 0; iter < 200; iter++) {

            r = h_work->Fit(f, "R Q S");

            double chi2ndf = (r->Ndf() > 0) ? r->Chi2()/r->Ndf() : 999;

            if (fabs(prevChi2 - chi2ndf) < .001) break;

            prevChi2 = chi2ndf;

            // widen if bad fit
            if (chi2ndf > 5) {
                f->SetRange(xmin - .2, xmax + .2);
            }
        }

        double chi2ndf = (r->Ndf() > 0) ? r->Chi2()/r->Ndf() : 0;

        // -------------------------------
        // EXTRACT BEST PEAK (FIRST GAUSS ONLY)
        // -------------------------------
        double a = f->GetParameter(0);
        double mean = f->GetParameter(1);
        double sigma = f->GetParameter(2);
        double sigerr = f->GetParError(2);

        energies.push_back(mean);
        sigmas.push_back(sigma);
        sigma_errs.push_back(sigerr);

        // -------------------------------
        // LOG OUTPUT
        // -------------------------------
        out << "\nCluster " << i << "\n";
        out << "Peaks in cluster: " << n << "\n";
        out << "Model: " << model << "\n";
        out << "Mean: " << mean << "\n";
        out << "Sigma: " << sigma << "\n";
        out << "Chi2/NDF: " << chi2ndf << "\n";

        // -------------------------------
        // DRAW ONLY FINAL FIT
        // -------------------------------
        f->SetLineColor(kRed);
        f->Draw("SAME");

        latex.DrawLatex(mean, f->Eval(mean)*1.05,
                        Form("%.1f", mean));
    }

    out.close();

    // -------------------------------
    // SAVE CANVAS + HIST
    // -------------------------------
    fout->cd();
    c->Write();
    h->Write();
}