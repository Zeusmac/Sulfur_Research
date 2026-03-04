void plot()
{

	TFile * f = new TFile("Mar_3_26/GammaNgated_bckgrndSub329_326-330_321_325.root");
	TFile * f1 = new TFile("Mar_3_26/Gammagated_bckgrndSub329_326-330_321_325.root");
	TFile * f2 = new TFile("Mar_3_26/Gammagated_peak329_326-330_321_325.root");
	TFile * f3 = new TFile("Mar_3_26/GammaNgated_329_331_335.root");
	TFile * f4 = new TFile("Mar_3_26/GammaNgated_bckgrnd329_326-330_321_325.root");
	TFile * f5 = new TFile("Mar_3_26/bckgrnd_no_neutron_329_321-325.root");
	
	TH1D * bckgrndSubN = (TH1D*) f -> Get("Decay_curve_329");
	TH1D * bckgrndSub = (TH1D*) f1 -> Get("Decay_curve_329");
	TH1D * peak = (TH1D*) f2 -> Get("Decay_curve_329");
	TH1D * peakN = (TH1D*) f3 -> Get("Decay_curve_329");
	TH1D * bckgrndN = (TH1D*) f4 -> Get("Decay_curve_329_bckgrnd");
	TH1D * bckgrnd = (TH1D*) f5 -> Get("px");
	
	bckgrnd -> SetTitle("Bckgrnd_329_321-325");
	bckgrndN -> SetTitle("NBckgrnd_329_321-325");
	bckgrnd -> GetXaxis()->SetRangeUser(0.0, 1000);
	bckgrndN -> GetXaxis()->SetRangeUser(0.0, 1000);

	TCanvas* ng = new TCanvas("NeutronGated","NeutronGated");
	TLegend* leg[6];
	for(int i = 0; i < 6; i++)
	{
		leg[i] = new TLegend(0.1, 0.7, 0.3, 0.9);
	}	


	ng -> Divide(2,3);
	ng -> cd(1);
	peak -> Draw();
	leg[0] ->AddEntry((TObject*)0, Form("Decay: 194.760701ms +_ 12.82490ms"), "");
	leg[0] ->AddEntry((TObject*)0, Form("chi^2: 1.044156"), "");
	leg[0] ->SetTextSize(0.04); 
 	leg[0] -> Draw();
	ng -> cd(2);
	peakN -> Draw();	
	leg[1] ->AddEntry((TObject*)0, Form("Decay: 120.198383ms +_ 13.341614ms"), "");
	leg[1] ->AddEntry((TObject*)0, Form("chi^2:  1.116453"), "");
	leg[1] ->SetTextSize(0.04); 
 	leg[1] -> Draw();

	ng -> cd(3);
	bckgrnd -> Draw();

	ng -> cd(4);
	bckgrndN -> Draw();

	ng -> cd(5);
	bckgrndSub -> Draw();
	
	leg[4] ->AddEntry((TObject*)0, Form("Decay: 77.129412ms +_ 1.576034ms  "), "");
	leg[4] ->AddEntry((TObject*)0, Form("chi^2: 34.443941"), "");
	leg[4] ->SetTextSize(0.04); 
 	leg[4] -> Draw();

	ng -> cd(6);
	bckgrndSubN -> Draw();	
	leg[5] ->AddEntry((TObject*)0, Form("Decay: 102.594574ms +_ 4.2039609e-05ms"), "");
	leg[5] ->AddEntry((TObject*)0, Form("chi^2: 4.798097"), "");
	leg[5] ->SetTextSize(0.04); 
 	leg[5] -> Draw();
	ng -> Draw();

}
