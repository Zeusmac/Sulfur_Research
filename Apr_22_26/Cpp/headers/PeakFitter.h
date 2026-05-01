#ifndef PEAKFITTER_H
#define PEAKFITTER_H

#include "TH1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "GammaDB.h"
#include "PeakTracker.h"

class PeakFitter {
public:
        PeakFitter(GammaDB& db,
               PeakTracker& tracker,
               TDirectory* sigmaDir,
               TDirectory* peakDir);

    void FitHistogram(TH1* h, TFile* fout);

private:
    GammaDB& db;
    PeakTracker& tracker;
    TDirectory* sigmaDir;
    TDirectory* peakDir;

    bool IsTimeHistogram(const std::string& name);
    std::pair<double,double> ExtractTimeRange(const std::string& name);
};

#endif
