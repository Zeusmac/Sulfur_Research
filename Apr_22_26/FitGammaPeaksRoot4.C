#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "TSpectrum.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

// -----------------------------
struct GammaLine {
    string isotope;
    double energy;
};
struct MatchResult {
    string isotope;
    double energy;
    double deltaE;
};
vector<GammaLine> db;
map<double, vector<double>> peak_times;
map<double, vector<double>> peak_counts;
map<double, vector<double>> peak_errors;
// -----------------------------
bool IsTimeHistogram(string name) {
    return (name.find("_") != string::npos &&
            isdigit(name[name.find("_") + 1]));
}
pair<double,double> ExtractTimeRange(string name) {

    double t1 = 0, t2 = 0;
 
    // expects: Gamma_1001_1101 
    sscanf(name.c_str(), "%*[^0-9]%lf_%lf", &t1, &t2);

    return {t1, t2};
}
double FindOrCreatePeakKey(double energy, double tol = 3.0) {

    for (auto &kv : peak_times) {
        if (fabs(kv.first - energy) < tol) {
            return kv.first;
        }
    }

    // new peak
    peak_times[energy]  = {};
    peak_counts[energy] = {};
    peak_errors[energy] = {};

    return energy;
}
vector<GammaLine> LoadGammaDB(const string& filename) {

    vector<GammaLine> out;
    ifstream in(filename);

    if (!in.is_open()) {
        cerr << "ERROR: Cannot open file: " << filename << endl;
        return out;
    }

    string line;
    while (getline(in, line)) {

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        string iso;
        double energy;

        stringstream ss(line);

        if (!(ss >> iso >> energy)) {
            cerr << "Skipping bad line: " << line << endl;
            continue;
        }

        out.push_back({iso, energy});
        cout << "Loaded: " << iso << " " << energy << endl;
    }

    cout << "DB size = " << out.size() << endl;
    return out;
}
// -----------------------------
vector<MatchResult> MatchPeak(double peakEnergy) {

    vector<MatchResult> matches;

    double tol = 3.0; // keV window

    for (auto &line : db) {

        double dE = fabs(line.energy - peakEnergy);

        if (dE <= tol) {
            matches.push_back({
                line.isotope,
                line.energy,
                dE
            });
        }
    }

    // sort best match first (smallest energy difference)
    sort(matches.begin(), matches.end(),
         [](const MatchResult &a, const MatchResult &b) {
             return a.deltaE < b.deltaE;
         });

    return matches;
}
// -----------------------------
void FitSingleHistogram(TH1* h, TFile* fout) {

    if (!h) return;

    string hname = h->GetName();
    
    // TEXT OUTPUT
    string txtname = "Gamma_fits/" + hname + "_fit.txt";
    ofstream out(txtname);
	double time = 0;

	if (IsTimeHistogram(hname)) {
	    auto [t1, t2] = ExtractTimeRange(hname);
	    time = 0.5 * (t1 + t2);
	} else {
	    time = -1; // special marker for full spectrum
	}

    if (!out.is_open()) {
        cerr << "Cannot open output file\n";
        return;
    }

    h->Sumw2();

    // Canvas
    TCanvas *c = new TCanvas(Form("c_%s", hname.c_str()), hname.c_str(), 1000, 700);
    //h->Draw();

    TLatex latex;
    latex.SetTextSize(0.025);

    TSpectrum s(1000);
    TH1 *bkghist = s.Background(h, 12);
    int nPeaks = s.Search(h, 2, "", 0.02);
    h -> Add(bkghist,-1);
    double *xPeaks = s.GetPositionX();
    sort(xPeaks, xPeaks + nPeaks);

    out << "====================================\n";
    out << "Histogram: " << hname << "\n";
    out << "Found " << nPeaks << " peaks\n";

    TF1 bg("bkg", "pol1");

    map<string,int> isotope_counts;

    // -----------------------------
    // Fit loop
    // -----------------------------
    for (int i = 0; i < nPeaks; i++) {

        double peakX = xPeaks[i];
        int bin = h->FindBin(peakX);

        if (h->GetBinContent(bin) < 50) continue;

	double sigma_guess = 2;
	double xmin = peakX - 3*sigma_guess;
	double xmax = peakX + 3*sigma_guess;
        // UNIQUE function name
        TF1 *f = new TF1(Form("fit_%s_%d", hname.c_str(), i),
                         "gaus(0)+pol1(3)", xmin, xmax);

        double height = h->GetBinContent(bin);

        f->SetParameters(height, peakX, 2.0, 0, 1);
        f->SetParLimits(0, 0, 1e9);
        f->SetParLimits(2, 0.2, 10);
        f->SetParLimits(1, xmin, xmax);

        TFitResultPtr result = h->Fit(f, "R S Q L+");

        double a = f->GetParameter(0);
        double aerr = f->GetParError(0);
        double sig = f->GetParameter(2);
        double sigerr = f->GetParError(2);
        double x0 = f->GetParameter(1);

        double peak_count = a * sig * sqrt(2*TMath::Pi());

        double peak_counterr = peak_count * sqrt(
            pow(aerr/a,2) + pow(sigerr/sig,2)
        );

        bg.SetParameters(f->GetParameter(3), f->GetParameter(4));

        double bck_counts = bg.Integral(x0 - 2.5*sig, x0 + 2.5*sig);

        double SNR = peak_count / sqrt(bck_counts);
        double chi2ndf = result->Chi2() / result->Ndf();
        double p_value = result->Prob();

        // -----------------------------
        // Filtering
        // -----------------------------
        if (h->GetBinContent(bin) < 50) continue;
        if (sig >= 6) continue;
        if (SNR < 10.0) continue;
//	f -> Draw("same");
        // -----------------------------
        // Matching
        // -----------------------------
        auto matches = MatchPeak(x0);
        // -----------------------------
        // Output text
        // -----------------------------
        out << "\nPeak " << i << ":\n";
        out << " p value =\t" << p_value << endl;

        string label = Form("%.1f", x0);

        for (auto &m : matches) {
            out << " Match: " << m.isotope
                << " (" << m.energy << " keV)\n";

            label += " ";
            label += m.isotope;

            isotope_counts[m.isotope]++;
        }
	out << "Fit Status   : " << result->Status() << "\n";
	out << "EDM          : " << result->Edm() << "\n";
        out << " Amplitude =\t " << a << " ± " << aerr << endl;
        out << " Mean   =\t" << x0 << " ± " << f->GetParError(1) << endl;
	out << " Possible matches:\n";

	if (matches.empty()) {
    		out << "   No isotope match within tolerance\n";
		} else {

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
	double key = FindOrCreatePeakKey(x0);
	if (time >= 0) {
	    peak_times[key].push_back(time);
	    peak_counts[key].push_back(peak_count);
	    peak_errors[key].push_back(peak_counterr);
	}	
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

    // -----------------------------
    // Save to ROOT file
    // -----------------------------
    fout->cd();

    c->Write();   // canvas with labels + fits
    h->Write();   // histogram itself

    c->Update();
}

// -----------------------------
// MAIN
// -----------------------------
void FitGammaAll(const char* rootfile) {

    db = LoadGammaDB("Isotope_energys.txt");
	cout << db.size() << endl;
    TFile *file = new TFile(rootfile);

    if (!file || file->IsZombie()) {
        cout << "Error opening file\n";
        return;
    }

    // ONE output ROOT file
    TFile *fout = new TFile("AllGammaFits.root", "RECREATE");

    TIter next(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {

        TObject *obj = key->ReadObj();

        if (obj->InheritsFrom("TH1")) {

            TH1 *h = (TH1*)obj;

            FitSingleHistogram(h, fout);
        }
    }
    
    fout->Close();	
    TFile *fout2 = new TFile("all_peak_time_graphs.root", "RECREATE");

	int graph_id = 0;

	for (auto &kv : peak_times) {

	    double energy = kv.first;

	    auto &t = peak_times[energy];
	    auto &c = peak_counts[energy];
	    auto &e = peak_errors[energy];

	    int n = t.size();

	    if (n < 2) continue;  // skip useless graphs

	    TGraphErrors *g = new TGraphErrors(n);

	    for (int i = 0; i < n; i++) {
		g->SetPoint(i, t[i], c[i]);
		g->SetPointError(i, 0, e[i]);
	    }

	    g->SetName(Form("peak_%d_keV", (int)energy));
	    g->SetTitle(Form("Peak %.1f keV;Time (ms);Counts", energy));
	    g->SetMarkerStyle(20);

	    g->Write();

	    graph_id++;
	}

	fout2->Close();
    cout << "Finished. All fits saved to AllGammaFits.root\n";
}
