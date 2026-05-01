#ifndef FITGAMMA_ANALYZER_H
#define FITGAMMA_ANALYZER_H

#include "GammaDB.h"
#include "PeakTracker.h"
#include "PeakFitter.h"
#include "ResolutionAnalyzer.h"

class FitGammaAnalyzer {
public:
    void Run(const char* rootfile);

private:
    GammaDB db;
    PeakTracker tracker;
};

#endif
