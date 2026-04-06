TH1D* ZeroBinsInRange(TH1* h, double xmin, double xmax, double threshold)
{
    if (!h) {
        std::cerr << "Histogram pointer is null!" << std::endl;
        return nullptr;
    }

    // Clone the histogram (keeps binning, title, etc.)
    TH1D* hnew = (TH1D*)h->Clone(Form("%s_modified", h->GetName()));
    hnew->SetDirectory(0); // detach from any ROOT file

    int nbins = hnew->GetNbinsX();

    for (int i = 1; i <= nbins; i++) {
        double x = hnew->GetBinCenter(i);
        double content = hnew->GetBinContent(i);

        // Apply condition
        if (x >= xmin && x <= xmax && content >= threshold) {
            hnew->SetBinContent(i, 0.0);
            hnew->SetBinError(i, 0.0); // optional but usually correct
        }
    }
    hnew -> Draw();
    return hnew;
}
