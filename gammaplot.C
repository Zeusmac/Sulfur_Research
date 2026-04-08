void gammaplot()
{
    // --- Open ROOT files ---
    TFile * f = new TFile("Mar_26_26/Bestfits/fit_g2783-2790_bg2800-2807_bin10.root");
    TFile * f1 = new TFile("Mar_26_26/Bestfits/328peakfit_bin30.root");
    TFile * f2 = new TFile("Mar_26_26/Bestfits/Npeak328fit_bin30.root");

    // --- Get histograms ---
    TCanvas * peak2789 = (TCanvas*) f -> Get("c");
    TCanvas * peak328  = (TCanvas*) f1 -> Get("c");
    TCanvas * peak328N = (TCanvas*) f2 -> Get("c");

    // --- Axis formatting ---

		
    // --- Create canvas and legends ---
    TCanvas* ng = new TCanvas("GammaGated","GammaGated",1200,900);
    ng->Divide(3,1);

    // --- Fill each pad ---
    ng->cd(1);
    peak2789->DrawClonePad();
    ng->cd(2);
    peak328->DrawClonePad();
    ng->cd(3);
    peak328N->DrawClonePad();

    ng->Update();
}
