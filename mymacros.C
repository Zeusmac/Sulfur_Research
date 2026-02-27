int gCounter  = 0;


 void Draw(TObject * object)
 {
 TCanvas * canvas = new TCanvas(Form("canvas_%i",gCounter++),"");
 object -> Draw();
 canvas -> Update();
 canvas -> Draw();
 }