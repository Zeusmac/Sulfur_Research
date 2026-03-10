void plot()
{
    // --- Open ROOT files ---
    TFile * f = new TFile("Mar_5_26/GammaNgated_329_326_330bckrndSub.root");
    TFile * f1 = new TFile("Mar_5_26/Gammagated_329_326_330bckrndSub.root");
    TFile * f2 = new TFile("Mar_3_26/Gammagated_peak329_326-330_321_325.root");
    TFile * f3 = new TFile("Mar_3_26/GammaNgated_peak329_326-330.root");
    TFile * f4 = new TFile("Mar_3_26/GammaNgated_bckgrnd329_326-330_321_325.root");
    TFile * f5 = new TFile("Mar_3_26/bckgrnd_no_neutron_329_321-325.root");
    TFile * f6 = new TFile("Mar_5_26/Gammagated_2789_2781_2790bckrndSub.root");
    TFile * f7 = new TFile("Mar_5_26/Gammagated_2789_2781_2790.root");
    TFile * f8 = new TFile("Mar_5_26/bckgrnd2810_2819.root");

    // --- Get histograms ---
    TH1D * bckgrndSubN = (TH1D*) f -> Get("Decay_curve_329");
    TH1D * bckgrndSub  = (TH1D*) f1 -> Get("Decay_curve_329");
    TH1D * peak        = (TH1D*) f2 -> Get("Decay_curve_329");
    TH1D * peakN       = (TH1D*) f3 -> Get("Decay_curve_329");
    TH1D * bckgrndN    = (TH1D*) f4 -> Get("Decay_curve_329_bckgrnd");
    TH1D * bckgrnd     = (TH1D*) f5 -> Get("px");
    TH1D * peak2789    = (TH1D*) f7 -> Get("Decay_curve_329");
    TH1D * bckgrnd2789 = (TH1D*) f8 -> Get("bckgrnd2789");
    TH1D * bckgrndsub2789 = (TH1D*) f6 -> Get("Decay_curve_329");

    // --- Axis formatting ---
    TH1D* hists[] = {peak, peakN, bckgrnd, bckgrndN, bckgrndSub, bckgrndSubN, peak2789};
    for(auto h : hists)
    {
	if(!h) continue; // Skip if histogram wasn't loaded
        h->GetXaxis()->SetTitle("Time (ms)");
        h->GetYaxis()->SetTitle("Counts");
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        //h->SetLineWidth(2);
        //h->SetLineColor(kBlue+1);
        //h->SetMarkerStyle(20);
       //h->SetMarkerColor(kRed+1);
    }

    bckgrnd->SetTitle("Bckgrnd_329_321-325");
    bckgrndN->SetTitle("NBckgrnd_329_321-325");
    bckgrnd2789->SetTitle("Bckgrnd_2810_2819");	
    bckgrnd->GetXaxis()->SetRangeUser(0, 1000);
    bckgrndN->GetXaxis()->SetRangeUser(0, 1000);
    bckgrnd2789->GetXaxis()->SetRangeUser(0, 1000);
		
    // --- Create canvas and legends ---
    TCanvas* ng = new TCanvas("NeutronGated","NeutronGated",1200,900);
    ng->Divide(3,3);

    TLegend* leg[9];
    for(int i = 0; i < 9; i++)
    {
	leg[i] = new TLegend(0.30, 0.72, 0.70, 0.88);
        leg[i]->SetTextSize(0.035);
        leg[i]->SetBorderSize(0);
	leg[i]->SetFillStyle(0);
        leg[i]->SetFillColorAlpha(kWhite,0.8);
        leg[i]->SetTextFont(42);
    }

    // --- Fill each pad ---
    ng->cd(1);
    peak->Draw("E");
    leg[0]->AddEntry((TObject*)0, Form("Decay: 194.761 ms #pm 12.825 ms"), "");
    leg[0]->AddEntry((TObject*)0, Form("#chi^{2}: 1.044"), "");
    leg[0]->Draw();

    ng->cd(2);
    peakN->Draw("E");
    leg[1]->AddEntry((TObject*)0, Form("Decay: 120.198 ms #pm 13.342 ms"), "");
    leg[1]->AddEntry((TObject*)0, Form("#chi^{2}: 1.116"), "");
    leg[1]->Draw();

    ng->cd(4);
    bckgrnd->Draw("E");

    ng->cd(5);
    bckgrndN->Draw("E");

    ng->cd(7);
    bckgrndSub->Draw("E");
    leg[4]->AddEntry((TObject*)0, Form("Decay: 99.253536ms #pm 4.690601ms"), "");
    leg[4]->AddEntry((TObject*)0, Form("#chi^{2}: 14.199320"), "");
    leg[4]->Draw();

    ng->cd(8);
    bckgrndSubN->Draw("E");
    leg[5]->AddEntry((TObject*)0, Form("Decay: 99.942373ms #pm 9.244286ms"), "");
    leg[5]->AddEntry((TObject*)0, Form("#chi^{2}: 3.487275"), "");
    leg[5]->Draw();

    ng->cd(3);
    peak2789->Draw("E");
    leg[6]->AddEntry((TObject*)0, Form("Decay: 126.309848ms #pm 3.248221ms "), "");
    leg[6]->AddEntry((TObject*)0, Form("#chi^{2}: 1.001477"), "");
    leg[6]->Draw();

    ng->cd(9);
	bckgrndsub2789->Draw("E");
	leg[7]->AddEntry((TObject*)0, Form("Decay: 126.309848ms #pm 3.248221ms "), "");
	leg[7]->AddEntry((TObject*)0, Form("#chi^{2}: 1.001477"), "");
	leg[7]->Draw();	
    ng->cd(6);
    bckgrnd2789->Draw("E");

    ng->Update();
}
