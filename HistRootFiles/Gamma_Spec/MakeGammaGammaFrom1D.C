void MakeGammaGammaFrom1D(TH1D* h1)
{
    int nbins = h1->GetNbinsX();
    double xmin = h1->GetXaxis()->GetXmin();
    double xmax = h1->GetXaxis()->GetXmax();

    TH2D *GamGam = new TH2D("GamGam", "Gamma-Gamma (Combinatorial)",
                         nbins, xmin, xmax,
                         nbins, xmin, xmax);

    for (int i = 1; i <= nbins; i++)
    {
        double Ei = h1->GetBinCenter(i);
        double Ci = h1->GetBinContent(i);

        if (Ci <= 0) continue;

        for (int j = 1; j <= nbins; j++)
        {
            double Ej = h1->GetBinCenter(j);
            double Cj = h1->GetBinContent(j);

            if (Cj <= 0) continue;

            // Weight = product of counts (approx correlation)
            double weight = Ci * Cj;

            GamGam->Fill(Ei, Ej, weight);
        }
    }

    TFile *fout = new TFile("gamma_gamma.root", "RECREATE");
    GamGam->Write();
    fout->Close();

    hgg->Draw("COLZ");
}
