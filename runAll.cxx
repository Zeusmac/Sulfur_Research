
void runAll() {
    TFile *f = new TFile("/Users/macwheeler/Sulfur_Research/Decay_Curve_2789keVFit.root");
    //gROOT->ProcessLine("TFile *f = new TFile(\"GammaSpec_44S.root\")");
   // gROOT->ProcessLine("TFile *f = new TFile(\"Decay_Curve_2789kev.root\")");
    gROOT->ProcessLine(".L MultGausTest.cxx");
    TH2D * dtime = (TH2D*) f-> Get("Decay_Curve_2789kev");


    //gROOT->ProcessLine("dtime -> ProjectionX(\"Decay_Curve_2789keV\",2782,2790)");
    //dtime -> ProjectionX("Decay_Curve",20,3500);
    //gROOT->ProcessLine("dtime -> ProjectionY(\"Gamma_spec_100ms\",1001,1101)");
    gROOT->ProcessLine("Decayp(Decay_Curve_2789keV)");
    //gROOT->ProcessLine("MultGausFit(dtime_py,1000,1000)");

    
}
