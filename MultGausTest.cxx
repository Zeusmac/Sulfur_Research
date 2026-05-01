//include <iostream>
//include "TMath.h"
//include "MultGausTest.h"


double Gaus_fit(double *dim, double *par)
{
  // 5 paremeters used in function: a, x0, sigma, slope, intercept
  double a = par[0];
  double x0 = par[1];
  double sigma = par[2];
  double bkgrd = par[3] * dim[0] + par[4];
  double x = dim[0];

  return a * TMath::Exp(-((x - x0) * (x - x0)) / (2 * sigma * sigma)) + bkgrd;
}

// double LinFit(double *dim, double *par)
// {
//   return par[0] + par[1] * dim[0];
// }


// double decay_daughter(double *dim,double *par)
// { 
//     /* function for fitting decay curves 
//   a is amplitude of source,
//   b is tau of the source,
//   c is tau of the daughter,
//   bkgrd is backround,
//    batemans equation.
//   */
//   double a = par[0];
//   double b = par[1];
//   double c = par[2];
//   double d = par[3];
//   double bkgrd = par[3];
//   double x = dim[0];
//   if( x < 0 ) return bkgrd;
//   return d*((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x)))+ bkgrd;
// }

double decay_parent(double *dim, double *par)
{
  /* function for fitting decay curves 
  a is amplitude of source,
  b is tau of the source,
  bkgrd is backround,
  this first term is the exponetial decay 
  */
  double a = par[0];
  double b = par[1];
  double bkgrd = par[2];
  double x = dim[0];

  if(x < 0) return bkgrd; 
  return  a * b *  TMath::Exp(-(b) * (x)) + bkgrd;
}

Double_t decay_background(Double_t *dim, Double_t *par)
{
/* function for fitting decay curves 
  a is amplitude of source,
  b is tau of the source,
  c is tau of the daughter,
  bkgrd is backround shift,
  this first term is the exponetial decay the second term is batemans equation.
  */
  Double_t a = par[0]; //count
  Double_t b = par[1]; //lamp
  Double_t c = par[2]; //lamD
  Double_t d = par[3]; //ratio
  Double_t e = par[4]; // background counts
  Double_t f = par[5]; // lamB
  Double_t bkgrd = par[6];
  Double_t x = dim[0];

  return a * TMath::Exp(-(b) * (x)) + d*((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x))) + e * TMath::Exp(-(f) * (x)) + bkgrd;
}

Double_t decay_nbeta(Double_t *dim, Double_t *par)
{
/* function for fitting decay curves 
  a is amplitude of source,
  b is tau of the source,
  c is tau of the daughter,
  bkgrd is backround shift,
  this first term is the exponetial decay the second term is batemans equation.
  */
  Double_t a = par[0]; //count
  Double_t b = par[1]; //lamp
  Double_t c = par[2]; //lamD
  Double_t d = par[3]; //daughter ratio
  Double_t t = par[4]; //beta n ratio
  Double_t g = par[5]; // lam beta n
  Double_t e = par[6]; // background counts
  Double_t f = par[7]; // lamB
  Double_t bkgrd = par[8];
  Double_t x = dim[0];

  return a * TMath::Exp(-(b) * (x)) + d*((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x))) + t*((b*a)/(g-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(g)*(x))) + e * TMath::Exp(-(f) * (x)) + bkgrd; bkgrd;;
}


// double decayE(double *dim, double *par)
// {
//   /* function for fitting decay curves 
//   a is amplitude of source,
//   b is tau of the source,
//   c is tau of the daughter,
//   d is the ratio of daughter to source
//   e is the ratio of beta_n decay
//   f halflife of beta_n decay 
//   g is the shift
//   bkgrd is background,
//   this first term is the exponetial decay the second term is batemans equation.
//   */
//   double a = par[0];
//   double b = par[1];
//   double c = par[2];
//   double d = par[3];
//   double e = par[4];
//   double f = par[5];
//   double h = par[7];
//   double j = par[6];
//   double g = par[8];
//   double bkgrd = par[9];
//   double x = dim[0];

//   if( x < g) return bkgrd;
//   return a * TMath::Exp(-(b) * (x)) + d*((b*a)/(c-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(c)*(x))) + e*((b*a)/(f-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(f)*(x))) + h*((b*a)/(j-b))*(TMath::Exp(-(b)*(x))-TMath::Exp(-(j)*(x))) + bkgrd;
// }




//44SGated
void Decayp(TH1* hist,const char* filename, int bin)
{
  double MinX, MaxX;
  MaxX = 5000; // hist -> GetXaxis() -> GetXmax();
  MinX = 0; // hist -> GetXaxis() -> GetXmin();
  //MinX = 1;
  hist -> Rebin(bin);
  //MaxX = ;
  TCanvas * fitgraph = new TCanvas("fits");
  fitgraph -> Clear();
  Double_t Parent_halflife = TMath::Log(2)/101;
  //function to apply the fit Decayp to the root interface
  TF1 *decayp = new TF1("decayp", decay_parent, MinX, MaxX,3);
  decayp->SetParName(0, "source #     ");
  decayp->SetParName(1, "Source lambda");
  decayp->SetParName(2, "shift        ");
  decayp->SetParameters(1000000, 
                        .004,
                        2000000);
  decayp->SetParLimits(0, 0, 100e7);
  decayp->SetParLimits(1, .001, .1);
  decayp->SetParLimits(2, 0, 100e6);

  for (int i = 0; i < decayp->GetNpar(); i++)
    {
    printf("Initial p%d = %f\n", i, decayp->GetParameter(i));
    }

  hist->Fit("decayp", "RS", "", MinX, MaxX);

  Int_t n_params = decayp->GetNpar();
  for (int i = 0; i < n_params; ++i) 
  {
    Double_t par_min, par_max;
    decayp->GetParLimits(i, par_min, par_max);
    std::cout << "Parameter " << i << " (" << decayp->GetParName(i) 
                                           << ") limits: [" << par_min << ", " << par_max << "]" << std::endl;
  }

  cout<< "Min X: " << MinX << endl;
  cout<< "bin :" << bin << "ms" << endl;  

  Double_t decayhl = TMath::Log(2)/(decayp ->GetParameter(1));
  Double_t decay_error = (TMath::Log(2)/TMath::Power((decayp -> GetParameter(1)),2))*(decayp -> GetParError(1));
  printf("Decay: %7fms +_ %7fms \n", decayhl, decay_error);
  
  Double_t chi2 = decayp -> GetChisquare();
  double_t NDF = decayp -> GetNDF();
  printf("Chi^2/NDF: %7f\n", chi2/NDF);
  hist -> Print("V");

  TF1 * bck = new TF1("bck", "[0]", MinX, MaxX);
  bck -> SetParameters(decayp ->GetParameter(2));
  bck -> SetLineColor(kBlack);

  auto legend = new TLegend(0.7,0.4,.9,.7);
  legend->AddEntry(hist,"Data","E");
  legend->AddEntry(decayp,"Fit","l");
  legend->AddEntry(bck,"bck","l");

  TString title = Form("%s;Time (ms);Counts", filename);	
  hist -> SetTitle(title);
  hist -> GetXaxis() -> CenterTitle(true);
  hist -> GetYaxis() -> CenterTitle(true);
  hist -> SetLineColor(kBlue);
  hist->GetXaxis()->SetRangeUser(MinX, MaxX);
  hist -> Draw("E");
  decayp -> SetLineColor(kRed); 
  decayp -> Draw("same");
  bck -> Draw("same");
  //fitgraph -> SetLogy(); 
  fitgraph -> Draw("E");
  //legend-> Draw();
  fitgraph -> SaveAs(filename);
}

//44S
void Decay_exp(TH1 *hist)
{
  double MinX, MaxX;
  MaxX = hist -> GetXaxis() -> GetXmax();
  MinX = hist -> GetXaxis() -> GetXmin();
  TCanvas * fitgraph = new TCanvas("fits");
  fitgraph -> Clear();
  MinX = 0;
  hist = hist -> Rebin(10);
  TF1 *decay = new TF1("decay", decay_nbeta, MinX, MaxX, 9);
  decay->SetParName(0, "source #             ");
  decay->SetParName(1, "Source lambda        ");
  decay->SetParName(2, "Daughter lambda      ");
  decay->SetParName(3, "Daughter ratio       ");
  decay->SetParName(4, "beta n ratio         ");
  decay->SetParName(5, "Lam beta N           ");
  decay->SetParName(6, "Bckgrnd Source       ");
  decay->SetParName(7, "lamb                 ");
  decay->SetParName(8, "background           ");
  decay->SetParameters(287839, 
                      (TMath::Log(2)/126),
                      (TMath::Log(2)/(542)), 
                      0.7, // ratio
                      0.4, //beta n ratio
                      TMath::Log(2)/(3.13e3), //beta n lam
                      300000, //background source
                      0.0001, //lamb
                      736663); // background shift
  decay->SetParLimits(0, 100000, 1500000); // # source
  decay->SetParLimits(1, TMath::Log(2)/(126.308692 + 3.246710), TMath::Log(2)/(126.308692 - 3.246710)); // source lamp
  decay->SetParLimits(2, (TMath::Log(2)/(542)), (TMath::Log(2)/(542))); // lamd
  decay->SetParLimits(3, .4, 1); // daughter ratio 
  decay->SetParLimits(4, .00001,.45); //beta n ratio
  decay->SetParLimits(5, (TMath::Log(2)/((3.13)*1e3)), (TMath::Log(2)/((3.13)*1e3))); // lam beta n
  decay->SetParLimits(6, 100000, 500000); // background source
  decay->SetParLimits(7, .000001, .001); // lamb
  decay->SetParLimits(8, 40000, 1500000); //background

  for (int i = 0; i < decay ->GetNpar(); i++) 
    {
    printf("Initial p%d = %f\n", i, decay->GetParameter(i));
    }

  std::cout << "MinX :" << MinX << std::endl;
  hist->Fit("decay", "RS", "", MinX, MaxX);

  Int_t n_params = decay->GetNpar();
  for (int i = 0; i < n_params; ++i) 
  {
    Double_t par_min, par_max;
    decay->GetParLimits(i, par_min, par_max);
    std::cout << "Parameter " << i << " (" << decay->GetParName(i) 
                              << ") limits: [" << par_min << ", " << par_max << "]" << std::endl;
  }
  Double_t decayhl = TMath::Log(2)/(decay ->GetParameter(1));
  Double_t decay_error = (TMath::Log(2)/TMath::Power((decay -> GetParameter(1)),2))*(decay -> GetParError(1));
  printf("Decay: %7fms +_ %7fms \n", decayhl, decay_error);

  Double_t decay_daughter = TMath::Log(2)/(decay -> GetParameter(2));
  Double_t decay_daughter_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(2)),2))*(decay -> GetParError(2));
  printf("Decay_daughter: %7fms +_ %7fms \n", decay_daughter, decay_daughter_error);

  Double_t decay_betan= TMath::Log(2)/(decay -> GetParameter(5));
  Double_t decay_betan_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(5)),2))*(decay -> GetParError(5));
  printf("Decay_nbeta: %7fms +_ %7fms \n", decay_betan, decay_betan_error);

  Double_t decay_background = TMath::Log(2)/(decay -> GetParameter(7));
  Double_t decay_background_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(7)),2))*(decay -> GetParError(7));
  printf("Decay_background: %7fms +_ %7fms \n", decay_background, decay_background_error);

  Double_t chi2 = decay -> GetChisquare();
  double_t NDF = decay -> GetNDF();
  printf("Chi^2/NDF: %7f\n", chi2/NDF);
  hist -> Print("V");

  TF1 * parent = new TF1("parent", "[0] * exp(-[1]*x) + [2]", MinX, MaxX);
  parent -> SetParameters(decay ->GetParameter(0),decay -> GetParameter(1),decay -> GetParameter(8));

  TF1 * daughter = new TF1("daughter", "[3]*(([1]*[0])/([2]-[1]))*( (exp(-[1]*x)) - (exp(-[2]*x))) + [4]", MinX, MaxX);
  daughter -> SetParameters(decay ->GetParameter(0),decay -> GetParameter(1), decay -> GetParameter(2), decay -> GetParameter(3), decay -> GetParameter(8));
  daughter -> SetLineColor(kGreen);
  
  TF1 * betan = new TF1("Beta N", "[3]*(([1]*[0])/([2]-[1]))*( (exp(-[1]*x)) - (exp(-[2]*x))) + [4]", MinX, MaxX);
  betan -> SetParameters(decay ->GetParameter(0),decay -> GetParameter(1), decay -> GetParameter(5), decay -> GetParameter(4), decay -> GetParameter(8));

   TF1 * expbckgrnd = new TF1("expbckgrnd", "[0] * exp(-[1]*x) + [2]", MinX, MaxX);
  expbckgrnd -> SetParameters(decay ->GetParameter(6),decay -> GetParameter(7),decay -> GetParameter(8));

  TF1 * bckgrnd = new TF1("bckgrnd", "[0]", 0, 1000);
  bckgrnd -> SetParameters(decay -> GetParameter(8));

  // TF1 *decay_chi2_top = new TF1(
  //   "decay_chi2_top",
  //   "[0]*exp(-[1]*x) + [3]*(([1]*[0])/([2]-[1]))*(exp(-[1]*x)-exp(-[2]*x)) + [4]*exp(-[5]*x) + [6]",
  //   0, 1976);

  // decay_chi2_top-> SetParameters(decay -> GetParameter(0),
  //                             .005525,  //(.00772) //lambda /// chi2 = 2
  //                             decay -> GetParameter(2),
  //                             decay -> GetParameter(3),
  //                             decay -> GetParameter(4),
  //                             decay -> GetParameter(5),
  //                             decay -> GetParameter(6)
  //                             );
  // double chi2_sum = 0;
  // for (int i = MinX; i <= hist->GetNbinsX(); ++i) 
  //   {
  //   double x = hist->GetBinCenter(i);
  //   Double_t c = hist->GetBinContent(i);
  //   if (c > 0)  
  //     {
  //     Double_t model =  decay -> Eval(x);
  //     chi2_sum += (((model - c)*(model - c))/(c));
  //     }
  //   }

  // TF1 *decay_chi2_bottom = (TF1*)decay_chi2_top -> Clone("decay_chi2_bottom");
  // decay_chi2_bottom -> SetParameters(decay -> GetParameter(0),
  //                             (.00772), //lambda /// chi2 = 2
  //                             decay -> GetParameter(2),
  //                             decay -> GetParameter(3),
  //                             decay -> GetParameter(4),
  //                             decay -> GetParameter(5),
  //                             decay -> GetParameter(6)
  //                             );

  // TGraph *g1 = new TGraph(decay_chi2_top);
  // TGraph *g2 = new TGraph(decay_chi2_bottom);

  // int N = g1->GetN();
  // TGraph *band = new TGraph(2*N);
  // for (int i=0;i<N;i++)
  // {
  //   double x,y;
  //   g1->GetPoint(i,x,y);
  //   band->SetPoint(i,x,y);
  // }

  // for (int i=0;i<N;i++)
  // {
  //     double x,y;
  //     g2->GetPoint(N-1-i,x,y);
  //     band->SetPoint(N+i,x,y);
  // }

  
  //std::cout << "Chi2_sum: " << chi2_sum/NDF << std::endl;

  auto legend = new TLegend(0.4,0.1,.7,.3);
  legend->AddEntry(hist,"Data","E");
  legend->AddEntry(decay,"Fit","l");
  legend->AddEntry(parent,"Parent","l");
  legend->AddEntry(daughter,"Daughter","l");
  legend->AddEntry(betan,"Beta N","l");
  legend->AddEntry(expbckgrnd,"exp_bckgrnd","l");
  legend->AddEntry(bckgrnd,"bckgrnd","l");

  parent -> SetLineColor(kOrange);
  daughter -> SetLineColor(kRed-2);
  bckgrnd -> SetLineColor(kBlack);
  betan ->  SetLineColor(kGreen);
  hist -> SetLineColor(kBlack);
  decay -> SetLineColor(kRed); 
  expbckgrnd -> SetLineColor(kMagenta);  
  // band -> SetFillColorAlpha(kRed, 0.30);
  // decay_chi2_top -> SetLineColor(kRed);
  // decay_chi2_bottom -> SetLineColor(kRed);

  fitgraph -> cd();
  hist -> SetTitle("Decay of 44S Fit;Time (ms);Counts");
  hist -> GetXaxis() -> CenterTitle(true);
  hist -> GetYaxis() -> CenterTitle(true);
  //decay->GetXaxis()->SetRange(1000,2000);
  hist -> Draw( "E");
  decay -> Draw("same");
  expbckgrnd -> Draw("same");  
  //band->Draw("F same");
  //decay_chi2_bottom -> Draw("same");
  //decay_chi2_top -> Draw("same");
  parent -> Draw("same" );
  daughter -> Draw("same");
  bckgrnd -> Draw("same");
  betan -> Draw("same");
  //fitgraph -> SetLogy();
  fitgraph -> Draw( "E");
  legend -> Draw();
  fitgraph -> SaveAs("44SFit.png");
}

//46Cl
void Decayb(TH1 *hist)
{
  double MinX, MaxX;
  MaxX = hist -> GetXaxis() -> GetXmax();
  MinX = hist -> GetXaxis() -> GetXmin();
  TCanvas * fitgraph = new TCanvas("fits");
  fitgraph -> Clear();
 
  TF1 *decay = new TF1("decay", decay_background, MinX, MaxX, 7);
  decay->SetParName(0, "source #             ");
  decay->SetParName(1, "Source lambda        ");
  decay->SetParName(2, "Daughter lambda      ");
  decay->SetParName(3, "ratio                ");
  decay->SetParName(4, "Background coefficant");
  decay->SetParName(5, "LamB                 ");
  decay->SetParName(6, "background           ");
  decay->SetParameters(200, 
                      (TMath::Log(2)/232),
                      (TMath::Log(2)/8.4e3), 
                      .19, // ratio
                      10, //background coefficant
                      .0001, // lamb
                      1); // background shift
  decay->SetParLimits(0, 100, 400); // # source
  decay->SetParLimits(1, 0, 1); // source lamp
  decay->SetParLimits(2, (TMath::Log(2)/((8.4 + .5)*1e3)), (TMath::Log(2)/((8.4 - .5 )*1e3))); // lamd
  decay->SetParLimits(3, 0,.5); // ratio
  decay->SetParLimits(4, 0, 200); // background source
  decay->SetParLimits(5, 0, .0001); //lamb
  decay->SetParLimits(6, 0, 2); //backzground
  MinX = 10;
  hist->Fit("decay", "LRS", "", MinX, MaxX);
  cout<< "Min X;" << MinX << endl;
  Int_t n_params = decay->GetNpar();
  for (int i = 0; i < n_params; ++i) 
  {
    Double_t par_min, par_max;
    decay->GetParLimits(i, par_min, par_max);
    std::cout << "Parameter " << i << " (" << decay->GetParName(i) 
                              << ") limits: [" << par_min << ", " << par_max << "]" << std::endl;
  }
  Double_t decayhl = TMath::Log(2)/(decay ->GetParameter(1));
  Double_t decay_error = (TMath::Log(2)/TMath::Power((decay -> GetParameter(1)),2))*(decay -> GetParError(1));
  printf("Decay: %7fms +_ %7fms \n", decayhl, decay_error);

  Double_t decay_daughter = TMath::Log(2)/(decay -> GetParameter(2));
  Double_t decay_daughter_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(2)),2))*(decay -> GetParError(2));
  printf("Decay_daughter: %7fms +_ %7fms \n", decay_daughter, decay_daughter_error);

  Double_t decay_background = TMath::Log(2)/(decay -> GetParameter(5));
  Double_t decay_background_error = (TMath::Log(2)/TMath::Power((decay ->GetParameter(5)),2))*(decay -> GetParError(5));
  printf("Decay_background_decay: %7fms +_ %7fms \n", decay_background, decay_background_error);

  Double_t chi2 = decay -> GetChisquare();
  double_t NDF = decay -> GetNDF();
  printf("Chi^2/NDF: %7f\n", chi2/NDF);

  TF1 * parent = new TF1("parent", "[0] * exp(-[1]*x) + [2]", 0, 2000);
  parent -> SetParameters(decay ->GetParameter(0),decay -> GetParameter(1), decay -> GetParameter(3));
  parent -> SetLineColor(kOrange);

  TF1 * daughter = new TF1("daughter", "[3]*(([1]*[0])/([2]-[1]))*( (exp(-[1]*x)) - (exp(-[2]*x))) + [4]", 0, 2000);
  daughter -> SetParameters(decay ->GetParameter(0),decay -> GetParameter(1), decay -> GetParameter(2), decay -> GetParameter(3),decay -> GetParameter(6));
  daughter -> SetLineColor(kGreen);

  TF1 * bckgrnd = new TF1("bckgrnd", "[0] * exp(-[1]*x) + [2]", 0, 2000);
  bckgrnd -> SetParameters(decay ->GetParameter(4),decay -> GetParameter(5), decay -> GetParameter(6));
  bckgrnd -> SetLineColor(kBlack);

   auto legend = new TLegend(0.7,0.4,.9,.7);
  legend->AddEntry(hist,"Data","E");
  legend->AddEntry(decay,"Fit","l");
  legend->AddEntry(parent,"Parent","l");
  legend->AddEntry(daughter,"Daughter","l");
  legend->AddEntry(bckgrnd,"bckgrnd","l");

  fitgraph -> cd();
  hist -> SetTitle("Decay of 46Cl;Time (ms);Counts");
  hist -> GetXaxis() -> CenterTitle(true);
  hist -> GetYaxis() -> CenterTitle(true);
  hist -> Draw("E");
  decay -> SetLineColor(632); 
  decay -> Draw("same");
  parent -> Draw("same");
  daughter -> Draw("same");
  bckgrnd -> Draw("same");
  //fitgraph -> SetLogy();
  legend -> Draw();
  fitgraph -> Draw("E");
  fitgraph -> SaveAs("Cl46_Fit_MLM2.png");
}


void MultGausFit(TH1 *hist, double lower, double upper){
     //for (int i = 0; i < PeakNo; i++) {

      TCanvas * fitgraph = new TCanvas("fits");
      TF1 *F = new TF1("F", Gaus_fit, lower, upper, 5);
      F->SetParName(0, "a");
      F->SetParName(1, "x0");
      F->SetParName(2, "sigma");
      F->SetParName(3, "slope");
      F->SetParName(4, "intercept");
      F->SetParLimits(0, 0,2e6);
      F->SetParLimits(1, lower, upper);
      F->SetParLimits(2, 0, 10);
      F->SetParLimits(4, 0, 4e6);
      F->SetParameters(103558,2786.81,1.63737,2.9846e+06,-270.236);
      //hist->GetListOfFunctions()->Add(F);
      hist->Fit(F, "SR","", lower, upper);
      hist -> Print("V");
      hist -> Draw();
      F -> Draw("same");
      fitgraph -> Draw("E");
      fitgraph -> SaveAs("Peak1155KeV.png");
      
 // }
}
