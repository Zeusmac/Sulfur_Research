#define DecayFinder_cxx

#include "DecayFinder.h"
#include <TH2.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <iostream>
#include "Tools.h"
#include <vector>


//===================== Settings
Bool_t debug = false;

TString pidCutFileName = "PIDCuts.root";

double distThreshold = 0.1; //TODO correct for postion
double forwardTime = 800; //ms 
double backwardTime = 200; //ms
double lowcutgamEn = 20;
double highcutgamEn = 6000;
double addback_time = 20;
double beta_time = 0.6;

int These_isotopes[] = {10,11,12,13,14}; // elements that you want to look at
int size = 2; // size of These_isotopes
//===================== histograms

TH2F * hPID;
TH2F ** hPIDdecay;
TH2F ** hPIDImp;
TH2F * py;
TH2F ** hGamGam;

TH1F ** hDecay;
TH2F * hDist; 
TH1F ** hvetoF;

TH2F ** hRawGamEn;
TH2F ** hGamEnT;
TH2F ** hAddback;
TH2F ** hAddbackdecay;
TH2F ** hGamEnDecay;
TH2F ** hImpGamEn;
TH2F ** hBetaGamEn;
TH2F ** hGamBkgnd;
TH1F ** Imp_gam_dT;
TH1F ** beta_gam_dT;
TH1F * Add_dT;

TH2F * hIonsPos;
TH2F * hBetaPos;
TH2F * hBetabgnd;
TH2F * hIonbgnd;

//===================== Parameters

int tick2ns = 8;
double tick2ms = tick2ns / 1e6;
double tick2us = tick2ns / 1e3;
double tick2s = tick2ns / 1e9;


TFile * cutFile;
TCutG * cut = NULL;
TObjArray * cutList;
int numCut = 0;
int gamT_size = TMath::QuietNaN();

int numMatch = 0;
std::vector<double> null;
std::vector<std::vector<double>> decay_data;

std::vector<std::vector<int>> isotopes = {
                            {19,  7}, {20,  7}, {21,  7},
                            {22,  8}, {23,  8}, {24,  8},
                  {24,  9}, {25,  9}, {26,  9}, {27,  9},           {29,  9},
                  {27, 10}, {28, 10}, {29, 10}, {30, 10}, {31, 10}, {32, 10},
        {29, 11}, {30, 11}, {31, 11}, {32, 11}, {33, 11}, {34, 11}, {35, 11},
        {32, 12}, {33, 12}, {34, 12}, {35, 12}, {36, 12}, {37, 12}, {38, 12},
        {35, 13}, {36, 13}, {37, 13}, {38, 13}, {39, 13}, {40, 13}, {41, 13},
                  {39, 14}, {40, 14}, {41, 14}, {42, 14}, 
                            {43, 15}, {44, 15}
                          }; 

std::string isoSym(int Z){
  switch (Z) {
    case  3 : return "Li"; break;
    case  4 : return "Be"; break;
    case  5 : return "B"; break;
    case  6 : return "C"; break;
    case  7 : return "N"; break;
    case  8 : return "O"; break;
    case  9 : return "F"; break;
    case 10 : return "Ne"; break;
    case 11 : return "Na"; break;
    case 12 : return "Mg"; break;
    case 13 : return "Al"; break;
    case 14 : return "Si"; break;
    case 15 : return "P"; break;
    case 16 : return "S"; break;
  }
  return "XX";
}

double dist(double x1, double x2, double y1, double y2)
{
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
bool look_at(int Z, int these_isotopes[], int size)
{
  for( int i =0; i < size; i++)
  {
    if ( these_isotopes[i] == Z ) return true;
  }
  return false;
}
//===================== Begins
void DecayFinder::Begin(TTree * /*tree*/){


  TString option = GetOption();

  hPID = new TH2F("hPID", "PID (beam & Implant & back veto); A/Q; Z",160, -175, -150, 1500, 1000, 2500);

  hDist = new TH2F("hDist", "dist; dx; dy", 100, -distThreshold, distThreshold, 100, -distThreshold, distThreshold);
  
  hIonsPos = new TH2F("hIonsPos", "Ions Pos; x; y", 200, 0, 1, 200, 0, 1);
  hBetaPos = new TH2F("hBetaPos", "Beta Pos; x; y", 200, 0, 1, 200, 0, 1);
  hBetabgnd = new TH2F("hBetaBgnd", "Beta Pos; x; y", 200, 0, 1, 200, 0, 1);
  hIonbgnd = new TH2F("hIonBgnd", "Beta Pos; x; y", 200, 0, 1, 200, 0, 1);
   
  Add_dT = new TH1F("Add_dT", "time difference between Corelated gammas of different crystals in a clover", 800, -200, 200);
  
  //-------------- GetCut;
  cutFile = new TFile(pidCutFileName, "READ");
  bool listExist = cutFile->GetListOfKeys()->Contains("cutList");
  if( !listExist ) {
    cutList = new TObjArray();
    numCut = (int) isotopes.size();
    
    printf("============= no Cut found, creating new default Cuts\n");
    
    for( int i = 0 ; i < numCut ; i++){
      
      cut = new TCutG();
      
      printf("%2d | A: %2d, Z: %d \n", i,  isotopes[i][0], isotopes[i][1]); 
        
      TString name; name.Form("cut%02d%s", isotopes[i][0], isoSym(isotopes[i][1]).c_str());
      cut->SetName(name);
      cut->SetVarX("AoQ");
      cut->SetVarY("Z");
      cut->SetTitle(Form("cut%2d", i));
      cut->SetLineColor(i+1);
      
      cut->SetPoint(0, (isotopes[i][0] - 0.45)/isotopes[i][1], isotopes[i][1] - 0.45);
      cut->SetPoint(1, (isotopes[i][0] - 0.45)/isotopes[i][1], isotopes[i][1] + 0.45);
      cut->SetPoint(2, (isotopes[i][0] + 0.45)/isotopes[i][1], isotopes[i][1] + 0.45);
      cut->SetPoint(3, (isotopes[i][0] + 0.45)/isotopes[i][1], isotopes[i][1] - 0.45);
      cut->SetPoint(4, (isotopes[i][0] - 0.45)/isotopes[i][1], isotopes[i][1] - 0.45);

      cutList->AddAtAndExpand(cut, i);
    }
  }
  else
  {
    cutList = (TObjArray*) cutFile->FindObjectAny("cutList");
    numCut = cutList->GetLast()+1;
    printf("------------ found %d cuts  in %s \n", numCut, pidCutFileName.Data());
    for( int k = 0; k < numCut; k++){
      if( cutList->At(k) != NULL )
      {
        printf("found a cut at %2d \n", k);
      }
      else
      {
        printf("     No cut at %2d \n", k);
      }
    }
  }

  hDecay = new TH1F *[numCut];
  hvetoF = new TH1F *[numCut];
  Imp_gam_dT = new TH1F*[numCut];
  beta_gam_dT = new TH1F*[numCut];
  hPIDdecay = new TH2F *[numCut];
  hPIDImp = new TH2F *[numCut];
  hRawGamEn = new TH2F *[numCut];
  hImpGamEn = new TH2F *[numCut];
  hBetaGamEn = new TH2F *[numCut];
  hGamBkgnd = new TH2F *[numCut];
  hGamEnT = new TH2F *[numCut];
  hAddback = new TH2F *[numCut];
  hAddbackdecay = new TH2F *[numCut];
  hGamGam = new TH2F *[numCut];
  hGamEnDecay = new TH2F *[numCut];
  
  for( int i = 0; i < numCut; i++)
  {
    hPIDdecay[i] = new TH2F(Form("hPIDDecay%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "PID (Valid beta decay); A/Q; Z", 160, -175, -150, 1500, 1000, 2500);
    hPIDImp[i] = new TH2F(Form("hPIDImp%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "PID (Valid beta decay); A/Q; Z",160, -175, -150, 1500, 1000, 2500);
    hRawGamEn[i] = new TH2F(Form("hRawGamEn%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "RawGamma Energy vs crystal", 52, 0, 52, 4000, 0, 4000);
    hImpGamEn[i] = new TH2F(Form("hImpGamEn%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Gamma Energy vs time difference between gamma and implant", 10000, -5000, 5000, 1000, 0, 1000);
    hBetaGamEn[i] = new TH2F(Form("hBetaGamEn%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Gamma Energy vs crystal", 52, 0, 52, 4000, 0, 4000);
    hGamBkgnd[i] = new TH2F(Form("hGamBkgnd%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "BkgndGamma Energy vs crystal", 52, 0, 52, 4000, 0, 4000);
    hGamEnT[i] = new TH2F(Form("hGamEnT%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Gamma Energy vs time",3000, 0, 3000, 1200, 0, 1200000);
    hAddback[i] = new TH2F(Form("hAddback%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Addback vs Clovers",13, 1, 13, 4000, 0, 4000);
    hAddbackdecay[i] = new TH2F(Form("hAddbackdecay%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Addback vs Clovers",1000,0,forwardTime,4000, 0, 4000);
    hGamGam[i] = new TH2F(Form("hGamGam%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Gamma Energy vs Gamma Energy",4000, 0, 4000, 4000, 0,4000);
    hGamEnDecay[i] = new TH2F(Form("hGamEnDecay%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "Gamma Energy vs Decay",400, 0, 400, 4000, 0,4000);
    Imp_gam_dT[i] = new TH1F(Form("imp_gam_dT%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "time difference between Corelated Ions and gammas", 10100, -100, 10000);
    beta_gam_dT[i] = new TH1F(Form("beta_gam_dT%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), "time difference between Corelated betas and gammas", 800, -200, 200);
    hDecay[i] = new TH1F(Form("hDecay%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), Form("Decay cut-%02d ; [ms]; count", i), 1000, -backwardTime, forwardTime); 
    hvetoF[i] = new TH1F(Form("hvetoF%02d%s",isotopes[i][0], isoSym(isotopes[i][1]).c_str()), Form("veto-F cut-%02d ; [ms]; count", i), 100, 0, 6000); 
    decay_data.push_back(null);
  }
  printf("=====================================\n");
}


//===================== Process
Bool_t DecayFinder::Process(Long64_t entry)
{
  Int_t cloverID[52] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13};//for crystalId -> cloverID (ryans cloverID mapping will always be 1-52s
  //if( entry > (Long64_t)1000 ) return kTRUE;

  b_Z->GetEntry(entry);
  b_AoQ->GetEntry(entry);
  b_crossTime->GetEntry(entry);
  b_crossEnergy->GetEntry(entry);
  b_flag->GetEntry(entry);
  b_vetoFlag->GetEntry(entry);
  b_xIons->GetEntry(entry);
  b_yIons->GetEntry(entry);
  b_dyIonsTime->GetEntry(entry);
  b_veto_r->GetEntry(entry);
  b_veto_f->GetEntry(entry);
  b_gamEn -> GetEntry(entry);
  b_gamID -> GetEntry(entry);
  b_gamT -> GetEntry(entry);
  b_gamNum -> GetEntry(entry);

  ULong64_t ImplantEntry = 0;
  ULong64_t ImplantTime = 0;
  int foundbeta = 0;
  double ImplantX = TMath::QuietNaN();
  double ImplantY = TMath::QuietNaN();
  bool found = false; // for finding the number fo correlated implants
  
  if( debug && entry % 10 == 0 ) printf("------------- %llu\n", entry);
  
  if( (flag == 3 || flag == 7) && (vetoFlag & 2) == 0 ) hPID->Fill(AoQ, Z); // with Beam + Ions, and has no veto_rear

  int cutID = -1;
  for(int i = 0; i < numCut; i++){
    cut = (TCutG*) cutList->At(i);
    if( cut->IsInside(AoQ, Z) ) cutID = i;
  }
  if( cutID == -1 ) return kTRUE;

  //if(look_at(isotopes[cutID][1],These_isotopes, size) == false) return kTRUE;

  // std::cout << "gamNum" << gamNum << std::endl;
  if( gamNum > 0 ) {
    for ( int l =0; l < gamNum;l++ )
      {
        if( (gamEn[l] < lowcutgamEn) || (gamEn[l] > highcutgamEn) ) continue;
        hRawGamEn[cutID] -> Fill(gamID[l],gamEn[l]);
        hGamEnT[cutID] -> Fill(gamT[l],gamEn[l]);
      }
  }

  
  if( (flag == 3 || flag == 7) && (vetoFlag & 2) == 0 ) // has beam +Ions , has no veto_rear
  {
    hPIDImp[cutID] -> Fill(AoQ,Z);
  }



  if( flag >= 4)  return kTRUE; /// has beta event
  ///if( flag != 3 ) return kTRUE; ///

  //if( (flag & 1) == 0 ) return kTRUE; /// no beam
  //if( (flag & 2) == 0 ) return kTRUE; /// no Ions
  //if( (flag & 4) == 4 ) return kTRUE; /// has only beta

  //if( veto_f > 0 ) return kTRUE;
  //if( veto_r > 0 ) return kTRUE;

  
  if( !((vetoFlag & 2) == 0)) return kTRUE; // events that have at rear veto hit
  //if( dyIonsTime[0] == 0 ) return kTRUE;
  //if (xIons < .2 || xIons > .8 ) return kTRUE;
  //if (yIons < .2 || yIons > .8 ) return kTRUE;

  ImplantEntry = entry;
  ImplantTime = dyIonsTime[0];
  ImplantX = xIons;
  ImplantY = yIons;
  
  if( TMath::IsNaN(xIons) ) return kTRUE;
  if( TMath::IsNaN(yIons) ) return kTRUE;
  
    
  if( debug ) printf("===========%8lld, %13llu, %f, %f\n", ImplantEntry, ImplantTime, ImplantX, ImplantY); 
    
  //searching from past 10k events to all.
  for( ULong64_t k = ImplantEntry - 10000; k < totNumEntry ; k++)
  { 
    b_flag->GetEntry(k);
    b_vetoFlag->GetEntry(k);
    b_dyBetaTime->GetEntry(k);
    b_crossTime->GetEntry(k);
    b_crossEnergy->GetEntry(k);
    b_veto_f->GetEntry(k);
    b_veto_r->GetEntry(k);
    b_xBeta->GetEntry(k);
    b_yBeta->GetEntry(k);
    b_gamT -> GetEntry(k);
    b_gamEn -> GetEntry(k);
    b_gamID -> GetEntry(k);
    b_gamNum -> GetEntry(k);
    
    if ( k == ImplantEntry ) continue; // skip itself

    // has only front veto hits
    // if( (vetoFlag & 2) == 2 ) continue; /// only front veto can hit
    if( (vetoFlag & 4) == 4 ) continue; 

    //if (flag == 4) continue; 
    if( (flag & 4) == 0 ) continue; /// Has beta  
    if( (flag & 2) == 2 ) continue; /// Has no implant   
    if( !TMath::IsNaN(veto_r) ) continue; /// no veto_rear

    //if( dyBetaTime[0] == 0) continue;

    if( dist(ImplantX, xBeta, ImplantY, yBeta) > distThreshold ) continue; // when dist between ions and beta > distThread
    
    if( k < ImplantEntry && ( ImplantTime - dyBetaTime[0] ) * tick2ms > backwardTime) continue; /// when past event, more than backwardTime
    if( k > ImplantEntry && ( dyBetaTime[0] - ImplantTime ) * tick2ms > forwardTime ) break; /// when future event, more then forwardTime
    
    

    //===== print status
    clock.Stop("timer");
    Double_t time = clock.GetRealTime("timer");
    clock.Start("timer");
    printf(" %llu[%4.2f%%], searched next %lld events | match case %d | expect : %5.2f min\r", 
         ImplantEntry, ImplantEntry*100./totNumEntry, k - ImplantEntry, 
         numMatch, totNumEntry*time/(ImplantEntry+1.)/60.);
      
    if( debug ) printf("           %8lld, %13llu, %f, %f\n", k, dyBetaTime[0], xBeta, yBeta);
    if( debug ) printf("+++++++++++%8s, %13llu | %.f ms\n", "+", dyBetaTime[0], (dyBetaTime[0] - ImplantTime)*tick2ms);
    
    hDist->Fill(ImplantX - xBeta, ImplantY - yBeta);
      
    numMatch ++;           
    double dT = 0;
    if( dyBetaTime[0] > ImplantTime )
    {
      dT = (dyBetaTime[0] - ImplantTime)*tick2ms;
      found = true; // found a valid implant
      foundbeta++;
    }
    else
    {
      dT = (ImplantTime - dyBetaTime[0])*tick2ms;
      dT = -dT;
    }
    hDecay[cutID]->Fill(dT);

    decay_data[cutID].push_back(dT);

    hPIDdecay[cutID] -> Fill(AoQ,Z);

    //=========== Gamma rays coicidence with implants
    for(int j = 0; j < gamNum; j++)
    {
      if( (gamEn[j] < lowcutgamEn) || (gamEn[j] > highcutgamEn) ) continue; 

      if(dT < 0 )
      { 
      hGamBkgnd[cutID] -> Fill(gamID[j],gamEn[j]);   //=========== Backround Gamma rays
      hBetabgnd -> Fill(xBeta,yBeta);
      hIonbgnd -> Fill(xBeta,yBeta);
      }

      double Imp_dT = 0;
      if (gamT[j] > ImplantTime)
      {
        Imp_dT = (gamT[j] - ImplantTime)*tick2ns; 
      }
      else
      {
        Imp_dT = (ImplantTime - gamT[j])*tick2ns; 
        Imp_dT = -Imp_dT;
      }
      hImpGamEn[cutID] -> Fill(Imp_dT,gamEn[j]);
      if(foundbeta == 1 )
      {
      Imp_gam_dT[cutID] -> Fill(Imp_dT);
      }
    }
      
    if((dT < 0) || (foundbeta != 1)) continue;
    hBetaPos->Fill(xBeta, yBeta);
    hIonsPos->Fill( xIons, yIons);
    // ============== Gamma rays in coicdence with beta
    for(ULong64_t p = k - 1000; p < (k + 1500); p++) // looking at gammas 1k(15ms) events ahead of this beta
    {  
      b_flag->GetEntry(p);
      b_vetoFlag->GetEntry(p);
      b_gamT -> GetEntry(p);
      b_gamEn -> GetEntry(p);
      b_gamID -> GetEntry(p);
      b_gamNum -> GetEntry(p);
      b_veto_r->GetEntry(p);


      int temp_gam[100] = {0};
      for(int j = 0; j < gamNum; j++)
      {
        Int_t tempint1 = cloverID[gamID[j]];
        if( (gamEn[j] < lowcutgamEn ) || (gamEn[j] > highcutgamEn) ) continue;

        double beta_dT = 0;
        if( gamT[j] > dyBetaTime[0] )
        {
          beta_dT = (gamT[j] - dyBetaTime[0])*tick2ns; 
        }
        else
        {
          beta_dT = (dyBetaTime[0] - gamT[j])*tick2ns; 
          beta_dT = -beta_dT;
        }
        beta_gam_dT[cutID] -> Fill(beta_dT);

        if( (abs(beta_dT) < beta_time) )
        {
          hBetaGamEn[cutID] -> Fill(gamID[j],gamEn[j]);
          hGamEnDecay[cutID] -> Fill(dT,gamEn[j]);    

          for(int k = 0; k < gamNum; k++)
          {
            Int_t tempint2 = cloverID[gamID[k]];
            if(k == j) continue;
            if((gamEn[k] < lowcutgamEn) || (gamEn[k] > highcutgamEn)) continue;
            double beta_dT = 0;
            if (gamT[k] > dyBetaTime[0])
            {
              beta_dT = (gamT[k] - dyBetaTime[0])*tick2us; 
            }
            else
            {
              beta_dT = (dyBetaTime[0] - gamT[k])*tick2us; 
              beta_dT = -beta_dT;
            }
            if( abs(beta_dT) < beta_time )
            {
              hGamGam[cutID] -> Fill(gamEn[j],gamEn[k]);
              if ( (tempint1 != tempint2) ) continue;
              if ( (temp_gam[j] == 1) || (temp_gam[k] == 1) ) continue;
              double gam_dT;
              if(gamT[k] >= gamT[j])
              {
                gam_dT = gamT[k] - gamT[j];
              }
              else
              {
                gam_dT = gamT[j] - gamT[k];
                gam_dT = -dT;
              }
              Add_dT -> Fill(gam_dT);
              if(abs(gam_dT) < addback_time)
              {
                temp_gam[j] += 1;
                temp_gam[k] += 1;
                double sum = gamEn[j] + gamEn[k];
                hAddback[cutID] -> Fill(tempint1,sum);
                hAddbackdecay[cutID] -> Fill(dT,sum);
              }
            }
          }
          if( temp_gam[j] == 0 )
          {
          hAddback[cutID] -> Fill(tempint1,gamEn[j]);
          hAddbackdecay[cutID] -> Fill(dT,gamEn[j]);
          }
        }
      }
    }
  }
    //hvetoF[cutID]->Fill(crossEnergy);
  return kTRUE;
}


void DecayFinder::Terminate(){
  
  for(int t = 0; t < numCut; t++)
  {
    //if(!(isotopes[t][1] == 11)) continue;
    store_decayXY(Form("/mnt/data_8TB/Mac_Stuffs/Decaytxt/hDecay%02d%s",isotopes[t][0], isoSym(isotopes[t][1]).c_str()),hDecay[t]);
  }
  //============== Save the results
  TFile * haha = new TFile("results.root", "recreate");
  haha->cd();
  for( int i = 0; i < numCut; i++){
    if(!(isotopes[i][1] == 11)) continue;
    hDecay[i]->Write();
  }
  haha->Close();
  
  //============== Plot result
  gStyle->SetOptStat("neiou");
  TCanvas * cDecay = new TCanvas("cDecay", "Decay", 1000, 1000);
  cDecay->Divide(2,2);
  
  int padID = 1; cDecay->cd(padID); 
  hPID->Draw("colz");
  for( int i = 0; i < numCut; i++){
    if(!(isotopes[i][1] == 11)) continue;
    cut = (TCutG*) cutList->At(i);
    cut->SetLineColor(i+2);
    cut->Draw("same");
  }
  
  padID++; cDecay->cd(padID); cDecay->cd(padID)->SetLogy();
  double yMax = 0;
  for( int i = 0; i < numCut; i++){
    if( hDecay[i]->GetMaximum() > yMax ) yMax = hDecay[i]->GetMaximum();
  }
  for( int i = 0; i < numCut; i++){
    if(!(isotopes[i][1] == 11)) continue;
    hDecay[i]->SetLineColor(i+2);
    hDecay[i]->SetMaximum(yMax *1.2);
    hDecay[i]->Draw(i == 0 ? "" : "same");
  }
  
  //padID++; cDecay->cd(padID);
  //hDist->Draw("colz");
  //
  //padID++; cDecay->cd(padID);
  //for( int i = 0; i < numCut; i++){
  //  hvetoF[i]->SetLineColor(i+2);
  //	hvetoF[i]->Draw(i == 0 ? "" : "same");
  //}

  padID++; cDecay->cd(padID);
  hIonsPos->Draw("colz");

  padID++; cDecay->cd(padID);
  hBetaPos->Draw("colz");

  //================ Save historgram 
  TString fHistRootName = "/mnt/data_8TB/Mac_Stuffs/DecayResults/con5_decay.root";
  if( true ){
    TFile * fHist = new TFile(fHistRootName, "recreate");
    fHist->cd();
    
    for( int i = 0 ; i < numCut; i++){

      //if( look_at(isotopes[i][1],These_isotopes,size) == false ) continue;

      hDecay[i]->Write();
      hPIDdecay[i]->Write();
      hPIDImp[i] -> Write();
      Imp_gam_dT[i] -> Write();
      beta_gam_dT[i] -> Write();
      hImpGamEn[i] -> Write();
      hBetaGamEn[i] -> Write();
      hRawGamEn[i] -> Write();
      hGamBkgnd[i] -> Write();
      hGamEnT[i] -> Write();
      hAddback[i] -> Write();
      hGamEnDecay[i] ->Write();
      hGamGam[i] -> Write();
      hAddbackdecay[i] -> Write();
    }
    Add_dT -> Write();
    hPID-> Write();
    hIonsPos -> Write();
    hBetaPos -> Write();
    hBetabgnd -> Write();
    hIonbgnd -> Write();
    fHist->Close();
    
    printf("---- Save PID histogram as %s \n", fHistRootName.Data());
  }

}
