#include "FitGammaAnalyzer.h"

#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "TSpectrum.h"
#include "TF1.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

// -----------------------------
FitGammaAnalyzer::FitGammaAnalyzer() {}

// -----------------------------
bool FitGammaAnalyzer::LoadGammaDB(const string& filename) {

    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "ERROR: Cannot open file: " << filename << endl;
        return false;
    }

    db.clear();

    string line;
    while (getline(in, line)) {

        if (line.empty() || line[0] == '#') continue;

        string iso;
        double energy;

        stringstream ss(line);
        if (!(ss >> iso >> energy)) continue;

        db.push_back({iso, energy});
    }

    cout << "Loaded DB entries: " << db.size() << endl;
    return true;
}

// -----------------------------
bool FitGammaAnalyzer::IsTimeHistogram(const string& name) {
    return (name.find("_") != string::npos &&
            isdigit(name[name.find("_") + 1]));
}

// -----------------------------
pair<double,double> FitGammaAnalyzer::ExtractTimeRange(const string& name) {

    double t1 = 0, t2 = 0;
    sscanf(name.c_str(), "%*[^0-9]%lf_%lf", &t1, &t2);
    return {t1, t2};
}

// -----------------------------
double FitGammaAnalyzer::FindOrCreatePeakKey(double energy, double tol) {

    for (auto &kv : peak_times) {
        if (fabs(kv.first - energy) < tol) return kv.first;
    }

    peak_times[energy]  = {};
    peak_counts[energy] = {};
    peak_errors[energy] = {};

    return energy;
}

// -----------------------------
vector<MatchResult> FitGammaAnalyzer::MatchPeak(double peakEnergy) {

    vector<MatchResult> matches;
    double tol = 3.0;

    for (auto &line : db) {
        double dE = fabs(line.energy - peakEnergy);
        if (dE <= tol) {
            matches.push_back({line.isotope, line.energy, dE});
        }
    }

    sort(matches.begin(), matches.end(),
         [](const MatchResult &a, const MatchResult &b) {
             return a.deltaE < b.deltaE;
         });

    return matches;
}

// -----------------------------
void FitGammaAnalyzer::FitSingleHistogram(TH1* h, TFile* fout) {

    if (!h) return;

    string hname = h->GetName();

    ofstream out("Gamma_fits/" + hname + "_fit.txt");

    double time = -1;
    if (IsTimeHistogram(hname)) {
        auto [t1, t2] = ExtractTimeRange(hname);
        time = 0.5 * (t1 + t2);
    }

    if (!out.is_open()) return;

    h->Sumw2();

    TCanvas *c = new TCanvas(Form("c_%s", hname.c_str()), "", 1000, 700);
    TLatex latex;
    latex.SetTextSize(0.025);

    TSpectrum s(1000);

    TH1 *bkghist = s.Background(h, 10);
    TH1D * bcksub_h =  h->Add(bkghist, -1);

    int nPeaks = s.Search(h, 2, "", 0.02);
    double *xPeaks = s.GetPositionX();

    sort(xPeaks, xPeaks + nPeaks);

    TF1 bg("bkg", "pol1");
    map<string,int> isotope_counts;

    for (int i = 0; i < nPeaks; i++) {

        double peakX = xPeaks[i];
        int bin = bcksub_h->FindBin(peakX);

        if (h->GetBinContent(bin) < 50) continue;

        double xmin = peakX - 6;
        double xmax = peakX + 6;

        TF1 *f = new TF1(Form("fit_%s_%d", hname.c_str(), i),
                         "gaus(0)+pol1(3)", xmin, xmax);

        double height = bcksub_h->GetBinContent(bin);

        f->SetParameters(height, peakX, 2.0, 0, 1);

        TFitResultPtr result = bcksub_h->Fit(f, "R S Q L 0");

        double a = f->GetParameter(0);
        double sig = f->GetParameter(2);
        double x0 = f->GetParameter(1);

        double peak_count = a * sig * sqrt(2*TMath::Pi());

        bg.SetParameters(f->GetParameter(3), f->GetParameter(4));
        double bck = bg.Integral(x0 - 2.5*sig, x0 + 2.5*sig);

        if (bck <= 0) continue;

        double SNR = peak_count / sqrt(bck);
        if (SNR < 10) continue;

        auto matches = db.Match(x0);

        string label = Form("%.1f", x0);
        for (auto &m : matches) {
            label += " " + m.isotope;
            isotope_counts[m.isotope]++;
        }

        double y = f->Eval(x0);
        latex.DrawLatex(x0, y * 1.1, label.c_str());

        double key = FindOrCreatePeakKey(x0);

        if (time >= 0) {
            peak_times[key].push_back(time);
            peak_counts[key].push_back(peak_count);
            peak_errors[key].push_back(0);
        }
    }

    fout->cd();
    c->Write();
    h->Write();
}

// -----------------------------
void FitGammaAnalyzer::FitAll(const char* rootfile) {

    TFile *file = new TFile(rootfile);
    if (!file || file->IsZombie()) return;

    TFile *fout = new TFile("AllGammaFits.root", "RECREATE");

    TIter next(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            FitSingleHistogram((TH1*)obj, fout);
        }
    }

    fout->Close();

    TFile *fout2 = new TFile("all_peak_time_graphs.root", "RECREATE");

    for (auto &kv : peak_times) {

        double energy = kv.first;
        auto &t = peak_times[energy];
        auto &c = peak_counts[energy];
        auto &e = peak_errors[energy];

        int n = t.size();
        if (n < 2) continue;

        TGraphErrors *g = new TGraphErrors(n);

        for (int i = 0; i < n; i++) {
            g->SetPoint(i, t[i], c[i]);
            g->SetPointError(i, 0, e[i]);
        }

        g->SetName(Form("peak_%d_keV", (int)energy));
        g->Write();
    }

    fout2->Close();

    cout << "Finished analysis\n";
}
