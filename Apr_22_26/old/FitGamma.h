#ifndef FITGAMMA_ANALYZER_H
#define FITGAMMA_ANALYZER_H

#include <vector>
#include <map>
#include <string>
#include <utility>

#include "TH1.h"
#include "TFile.h"
#include "TGraphErrors.h"

struct GammaLine {
    std::string isotope;
    double energy;
};

struct MatchResult {
    std::string isotope;
    double energy;
    double deltaE;
};

class FitGammaAnalyzer {
public:
    FitGammaAnalyzer();

    // Load isotope database
    bool LoadGammaDB(const std::string& filename);

    // Main entry
    void FitAll(const char* rootfile);

private:
    // Internal helpers
    bool IsTimeHistogram(const std::string& name);
    std::pair<double,double> ExtractTimeRange(const std::string& name);
    double FindOrCreatePeakKey(double energy, double tol = 3.0);
    std::vector<MatchResult> MatchPeak(double peakEnergy);

    void FitSingleHistogram(TH1* h, TFile* fout);

    // Data (was global before)
    std::vector<GammaLine> db;

    std::map<double, std::vector<double>> peak_times;
    std::map<double, std::vector<double>> peak_counts;
    std::map<double, std::vector<double>> peak_errors;
};

#endif
