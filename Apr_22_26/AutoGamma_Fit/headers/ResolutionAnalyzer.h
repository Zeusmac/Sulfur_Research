#ifndef RESOLUTION_ANALYZER_H
#define RESOLUTION_ANALYZER_H

#include "PeakTracker.h"
#include "TDirectory.h"

class ResolutionAnalyzer {
public:
    ResolutionAnalyzer(PeakTracker& tracker);

    void BuildFWHMvsTime(TDirectory* sigmaDir);

private:
    PeakTracker& tracker;
};

#endif