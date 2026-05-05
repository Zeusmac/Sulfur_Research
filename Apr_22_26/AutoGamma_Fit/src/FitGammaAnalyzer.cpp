#include "FitGammaAnalyzer.h"
#include "TFile.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include <iostream>

using namespace std;

void FitGammaAnalyzer::Run(const char* rootfile) {

    if (!db.Load("../Isotope_energys.txt")) return;

    TFile *file = new TFile(rootfile);
    if (!file || file->IsZombie()) return;

    TFile *fout = new TFile("AllGammaFits.root", "RECREATE");

     //Root Directories
    TDirectory *sigmaDir = fout->mkdir("sigma_vs_energy");
    TDirectory *peakDir = fout->mkdir("peak_count_vs_time");

    PeakFitter fitter(db, tracker, sigmaDir, peakDir);

    TIter next(file->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();

        fout -> cd();

        if (obj->InheritsFrom("TH1")) {
            fitter.FitHistogram((TH1*)obj, fout);
        }
    }

    fout->Close();

    cout << "Analysis complete\n";
}
