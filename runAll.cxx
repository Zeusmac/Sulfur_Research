
void runAll() {
    gROOT -> ProcessLine("TFile* f = new TFile(\"e21062_44S.root\")");
    gROOT->ProcessLine(".L MultGausTest.cxx");
    //TH2D * dt = (TH2D*) f-> Get("dtimeN");


    //dt -> ProjectionX("Decay_Curve",20,3500);
    gROOT->ProcessLine("dtime -> ProjectionX(\"Decay_curve_329\",326,330)");
    gROOT->ProcessLine("dtime -> ProjectionX(\"Decay_curve_329_bckgrnd\",321,325)");
    gROOT ->ProcessLine("Decay_curve_329 -> Add(Decay_curve_329_bckgrnd,-1)");

    //gROOT->ProcessLine("dtime -> ProjectionX(\"Decay_curve_879\",876,880)");
    //gROOT->ProcessLine("dtime -> ProjectionX(\"Decay_curve_879_bckgrnd\",860,864)");

    gROOT->ProcessLine("Decayp(Decay_curve_329,\"Gammagated_bckgrondSub329_321_325.png\",1)");
    
    //gROOT ->ProcessLine("Decay_curve_879 -> Add(Decay_curve_879_bckgrnd,-1)");

   // gROOT ->ProcessLine("Decay_curve_879 -> Add(Decay_curve_329)");
   //gROOT->ProcessLine("Decayp(Decay_curve_879)");
}

