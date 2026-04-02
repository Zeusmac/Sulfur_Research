void save_histograms_from_canvas() {
    // --- Open the ROOT file and get the canvas ---
    TFile* f = TFile::Open("Mar_26_26/BestFits/Npeak328fit_bin30.root");
    if (!f || f->IsZombie()) {
        std::cout << "Cannot open file!\n";
        return;
    }

    TCanvas* c = (TCanvas*) f->Get("c");
    if (!c) {
        std::cout << "Canvas not found!\n";
        return;
    }

    // --- Prepare vector for cloned histograms ---
    std::vector<TH1*> hists;

    // --- Loop over pads and primitives ---
    TIter nextPad(c->GetListOfPrimitives());
    while (TObject* padObj = nextPad()) {
        TPad* pad = dynamic_cast<TPad*>(padObj);
        if (!pad) continue;

        TIter nextObj(pad->GetListOfPrimitives());
        while (TObject* obj = nextObj()) {
            TH1* h = dynamic_cast<TH1*>(obj);
            if (h) {
                // Clone histogram to detach from original file
                TH1* clone = (TH1*) h->Clone();
                clone->SetDirectory(0);  // detach from file
                hists.push_back(clone);
                std::cout << "Cloned histogram: " << clone->GetName() << std::endl;
            }
        }
    }

    if (hists.empty()) {
        std::cout << "No histograms found!\n";
        return;
    }

    // --- Save all histograms to a new ROOT file ---
    TFile* outFile = new TFile("saved_histograms.root", "RECREATE");
    for (TH1* h : hists) {
        h->Write();
    }
    outFile->Close();

    std::cout << "Saved " << hists.size() << " histograms to saved_histograms.root" << std::endl;
}
