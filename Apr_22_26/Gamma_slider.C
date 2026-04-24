#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TSlider.h"
#include "TExec.h"

TH1D* Proj[500];
TCanvas *c;
TSlider *slider;
int current_index = 0;

void DrawSlice() {
    c->cd();
    Proj[current_index]->Draw("HIST");
    c->Modified();
    c->Update();
}

void Gamma_slider(const char* filename) {

    TFile *f = TFile::Open(filename);
    TH2D *h = (TH2D*)f->Get("dtime");

    int window = 100;
    int initial_bin = 1001;

    for (int i = 0; i < 50; i++) {
        int b1 = initial_bin + i * window;
        int b2 = b1 + window;
        Proj[i] = h->ProjectionY(Form("h_%d", i), b1, b2);
    }
    c = new TCanvas("c", "slider demo", 800, 600);

    slider = new TSlider("slider", "index", 0.1, 0.02, 0.9, 0.06);
    slider->SetRange(0, 1);

    DrawSlice();

    // manual update loop (interactive polling)
    while (true) {
        gSystem->ProcessEvents();

        int new_index = slider->GetMinimum() * 500;
        if (new_index != current_index) {
            current_index = new_index;
            DrawSlice();
        }
    }
}
