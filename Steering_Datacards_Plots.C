/*
Maco to run production of data cards or drawing plots from data cards serially.
Written by Merijn van de klundert, mvandekl@cern.ch
*/


#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;

void Steering_Datacards_Plots(){

  //compile and load all macros needed. Configurable_mt.h MUST be compiled and loaded first of course
  gROOT->ProcessLine(".L Configurable_mt.h++");
  gROOT->ProcessLine(".L ProduceDataCards.C++");
  gROOT->ProcessLine(".L PlotDataCards.C++");
  
  char cmdphi[1000];

  //here set the values for the oarguments of the macro.
  bool ProduceDataCard=true;
  bool ProducePlot=true;
  
  int CattoReview=0; //don't run plot macro with cat=-1, it only runs for 1 cat..
  int era_=2017;
  int decaymode=0;//if initialise to -1, no selection on decay mode made (mu-tau)!
  int CPMethod=0;
  int Observable=21;
  TString inputdir="/Users/klundert/DESY/HiggsCPProjectSoftware/Outputs/2019_11_26_DeepTau/predictions_2017/";

  for(int i=0;i<2;i++){ //here iterate over the observable of interest
    CPMethod=i;
    decaymode=i;
    //  CattoReview=i;
    
    if(ProduceDataCard){
    sprintf(cmdphi,"ProduceDataCards(%i,%i,%i,%i,%i,\"%s\")",CattoReview,era_,decaymode,CPMethod,Observable,inputdir.Data());
    gROOT->ProcessLine(cmdphi);}

    if(ProducePlot){
      sprintf(cmdphi,"PlotDataCards(%i,%i,%i,%i,%i)",CattoReview,era_,decaymode,CPMethod,Observable);    
      gROOT->ProcessLine(cmdphi);
    }  
  }
   
  return;
}

