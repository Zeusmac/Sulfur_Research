#include <iostream>

#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TString.h"

void Gamma_Spec_5000(const char* filename) {

    TFile *file = new TFile(filename, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file\n";
        return;
    }

    TH2D *h = (TH2D*)file->Get("dtime");
    if (!h) {
        std::cerr << "Error: histogram 'dtime' not found\n";
        file->Close();
        return;
    }

    TFile *out = new TFile("Gamma_spec100mswindow.root", "RECREATE");
    TH1D *py = h -> ProjectionY("gamma_spec_all",1001,6001);
    py -> Write();
    int window = 100; //window of 100ms
    int initial_bin = 1001;

    TH1D *Proj[500];

    for (int i = 0; i < 50; i++) {

        int bin1 = initial_bin + i * window;
        int bin2 = bin1 + window;

        Proj[i] = h->ProjectionY(
            Form("Gamma_%d_%d", bin1, bin2),
            bin1,
            bin2
        );

        Proj[i]->Write();
    }

    out->Close();
    file->Close();
}
