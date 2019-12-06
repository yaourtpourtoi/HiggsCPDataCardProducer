/*
Maco to add data cards for different configurations. It works in tandem with Configurable_mt
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
#include "Configurable_mt.h"

using namespace std;

//int era, decaymode, CPMethod,
void InitStatics();

//specify in TargetFileName where the combined data card is stored..
void AddDataCards(TString TargetFileName="CombinedDataCard.root"){

  //  gROOT->ProcessLine(".L Configurable_mt.h++"); //uncomment f you want to compile a macro dependent on Configurable class..

  //check if the data card exist
  if(!gSystem->AccessPathName(TargetFileName)) throw std::invalid_argument("Target file given to hadd exists already, exiting");
  
  TString CombineDataCard="hadd ";
  CombineDataCard+=TargetFileName;
  CombineDataCard+=" ";

  //initialise the static members, such that we can obtain smoothly the paths; no objects needed
  Configurable_mt::era=2017;
  Configurable_mt::decaymode=0;
  Configurable_mt::CPMethod=0;
  Configurable_mt::Observable=21;
  vector <TString> DatacardPaths;

  for(int index=0;index<2;index++){
  Configurable_mt::decaymode=index;
  Configurable_mt::CPMethod=index;    
  InitStatics();
  
  CombineDataCard+=Configurable_mt::OutputFileDataCard;
  CombineDataCard+=" ";
  }
  
  cout<<"Before applying following merging command: "<<CombineDataCard <<endl;
  gSystem->Exec(CombineDataCard);
									  
  cout<<"Combined data card created in "<<TargetFileName<<endl;

  return;
}
void InitStatics(){

  Configurable_mt::InitObservableNames(); //It initialised what is drawn in the histogram, the obs. name for output path, and sets axes titles
  Configurable_mt::InitCatNames(); ////This will initialise a vector mt_* era_* <category> *_channel* *_CP calc method (decay plane, IP method)*. Needed to get the histogram in a dir recognisable by CombineHarvester later. We store a histo finally in catname/samplename. mt_ is hardcoded in configurable_mt. We clear the vector catnames first, since we use push_back  
  Configurable_mt::InitOutputBaseDir(); //All statics must have been set then correctly. This sets the base directory OutputDir
  Configurable_mt::InitOutputFileDataCard(); //initialises  OutputFileDataCard (it uses OutputDir)
  return;
}
