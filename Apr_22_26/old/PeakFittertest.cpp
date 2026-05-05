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

    TLatex latex;
    latex.SetTextSize(0.025);

    ofstream out("../Gamma_fits/" + hname + "_fit.txt");
    if (!out.is_open()) return;

    h->Sumw2();
    h->Draw();

    TH1* h_work = (TH1*)h->Clone(Form("%s_work", hname.c_str()));

    // =========================================================
    // BACKGROUND MODEL (improved)
    // =========================================================
    TF1 *bg = new TF1("bg",
        "[0]*exp(-[1]*x) + pol2(2)",
        h->GetXaxis()->GetXmin(),
        h->GetXaxis()->GetXmax()
    );

    bg->SetParameters(h->GetMaximum(), 0.001, 0, 0, 0);

    h_work->Fit(bg, "Q");

    // subtract background for peak finding
    for (int i = 1; i <= h_work->GetNbinsX(); i++) {
        double x = h_work->GetBinCenter(i);
        h_work->SetBinContent(i,
            h_work->GetBinContent(i) - bg->Eval(x));
    }

    // =========================================================
    // TSPECTRUM FIRST PASS
    // =========================================================
    TSpectrum s(1000);
    int nPeaks = s.Search(h_work, 2, "", 0.02);

    double *xPeaks = s.GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);

    out << "Histogram: " << hname << "\n";
    out << "Initial peaks: " << nPeaks << "\n";

    // =========================================================
    // BUILD SIGMA MODEL σ(E)
    // =========================================================
    vector<double> rawE, rawS;

    for (int i = 0; i < nPeaks; i++) {

        double x = xPeaks[i];
        int bin = h->FindBin(x);

        if (h->GetBinContent(bin) < 50) continue;

        TF1 ftmp("tmp", "gaus", x-5, x+5);
        ftmp.SetParameters(h->GetBinContent(bin), x, 2.0);

        h->Fit(&ftmp, "Q R");

        double sig = fabs(ftmp.GetParameter(2));

        if (sig > 0 && sig < 20) {
            rawE.push_back(x);
            rawS.push_back(sig);
        }
    }

    double a = 1.0, b = 0.01;

    if (rawE.size() > 5) {
        TGraph gr(rawE.size());
        for (size_t i = 0; i < rawE.size(); i++)
            gr.SetPoint(i, rawE[i], rawS[i]);

        TF1 fres("fres", "sqrt([0] + [1]*x)", 0, 3000);
        gr.Fit(&fres, "Q");

        a = fres.GetParameter(0);
        b = fres.GetParameter(1);
    }

    auto SigmaE = [&](double E) {
        return sqrt(a + b*E);
    };

    // =========================================================
    // CLUSTER PEAKS
    // =========================================================
    vector<vector<double>> clusters;

    double mergeTol = 2.5;

    for (int i = 0; i < nPeaks; i++) {

        double x = xPeaks[i];
        bool added = false;

        for (auto &cl : clusters) {
            if (fabs(cl.front() - x) < mergeTol * SigmaE(x)) {
                cl.push_back(x);
                added = true;
                break;
            }
        }

        if (!added)
            clusters.push_back({x});
    }

    // =========================================================
    // FIT FUNCTION BUILDER
    // =========================================================
    auto MakeModel = [&](int nPeaks, double xmin, double xmax, int id) {

        if (nPeaks == 1)
            return new TF1(Form("f_%d", id),
                "gaus(0)+pol2(3)", xmin, xmax);

        if (nPeaks == 2)
            return new TF1(Form("f_%d", id),
                "gaus(0)+gaus(3)+pol2(6)", xmin, xmax);

        return new TF1(Form("f_%d", id),
            "gaus(0)+gaus(3)+gaus(6)+pol2(9)", xmin, xmax);
    };

    // =========================================================
    // FINAL FIT LOOP (ADAPTIVE)
    // =========================================================
    for (size_t i = 0; i < clusters.size(); i++) {

        auto &cl = clusters[i];
        sort(cl.begin(), cl.end());

        double E = cl[cl.size()/2];

        // use tracker sigma if available
        double sigma = -1;
        for (auto &kv : tracker.GetSigmas()) {
            if (fabs(kv.first - E) < 3.0)
                sigma = kv.second.back();
        }

        if (sigma < 0)
            sigma = SigmaE(E);

        double xmin = E - 5*sigma;
        double xmax = E + 5*sigma;

        TF1 *f = MakeModel(cl.size(), xmin, xmax, i);

        // --------------------------
        // INITIALIZE PARAMETERS
        // --------------------------
        for (size_t p = 0; p < cl.size(); p++) {

            double xp = cl[p];
            int bin = h->FindBin(xp);

            double height = h->GetBinContent(bin);

            f->SetParameter(p*3 + 0, height);
            f->SetParameter(p*3 + 1, xp);
            f->SetParameter(p*3 + 2, sigma);

            f->SetParLimits(p*3 + 2, 0.5*sigma, 1.5*sigma);
        }

        // --------------------------
        // ITERATIVE FIT LOOP
        // --------------------------
        TFitResultPtr r;

        for (int iter = 0; iter < 3; iter++) {
            r = h->Fit(f, "R Q S");

            double chi2ndf = r->Chi2() /
                std::max(1.0, (double)r->Ndf());

            if (chi2ndf < 2.0) break;

            // widen sigma if bad fit
            for (size_t p = 0; p < cl.size(); p++) {
                double sig = f->GetParameter(p*3+2);
                f->SetParameter(p*3+2, sig * 1.2);
            }
        }

        // =====================================================
        // RESIDUAL PEAK FINDER (SECOND PASS)
        // =====================================================
        TH1* res = (TH1*)h->Clone("res");
        res->Reset();

        for (int b = 1; b <= h->GetNbinsX(); b++) {
            double x = h->GetBinCenter(b);
            res->SetBinContent(b,
                h->GetBinContent(b) - f->Eval(x));
        }

        TSpectrum s2(50);
        int nNew = s2.Search(res, 2, "", 0.05);

        if (nNew > 0 && cl.size() < 3) {
            // try refit with more peaks
            double *xp2 = s2.GetPositionX();

            vector<double> newCl = cl;
            for (int j = 0; j < nNew; j++)
                newCl.push_back(xp2[j]);

            sort(newCl.begin(), newCl.end());

            TF1 *f2 = MakeModel(newCl.size(), xmin, xmax, 100+i);

            r = h->Fit(f2, "R Q S");
            f = f2;
        }

        // =====================================================
        // FINAL RESULTS
        // =====================================================
        double mean = f->GetParameter(1);
        double sig  = fabs(f->GetParameter(2));

        double amp = f->GetParameter(0);
        double counts = amp * sig * sqrt(2*TMath::Pi());

        double chi2ndf = r->Chi2() /
            std::max(1.0, (double)r->Ndf());

        if (counts < 10) continue;

        out << "\nCluster " << i << "\n";
        out << "Mean: " << mean << "\n";
        out << "Sigma: " << sig << "\n";
        out << "Chi2/NDF: " << chi2ndf << "\n";

        auto matches = db.Match(mean, 5);

        for (auto &m : matches)
            out << "Match: " << m.isotope
                << " (" << m.energy << " keV)\n";

        // DRAW ONLY FINAL FIT
        f->SetLineColor(kRed);
        f->Draw("same");

        latex.DrawLatex(mean, f->Eval(mean)*1.1,
                        Form("%.1f", mean));

        // tracker
        tracker.Add(mean, 0, counts, 0, sig, 0);
    }

    // =========================================================
    // SAVE OUTPUT
    // =========================================================
    fout->cd();
    c->Write();
    h->Write();

    c->Update();

    out.close();
}