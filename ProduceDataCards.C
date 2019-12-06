/*
Macro to produce datacards for Higgs CP measurement. It works in tandem with Configurable_mt
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

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"

#include <iostream>
#include "Configurable_mt.h"

using namespace std;

void ProduceDataCards(
				int CattoReview=-1,//-1: loop over all samples available. N: only loop over sample N => allows to run things parallel
				int era_=2017,
				int decaymode=0,//if initialise to -1, no selection made!
				int CPMethod=0, //initialised the charged particle cutoff and the appendix for the cat.
				int Observable=21,
				//TString directory = "/Users/klundert/DESY/HiggsCPProjectSoftware/Outputs/2019_11_17_Unfiltered_NewTopNorm/predictions_2017/", //latest versionbefore switchng to deeptau
				TString directory = "/Users/klundert/DESY/HiggsCPProjectSoftware/Outputs/2019_11_26_DeepTau/predictions_2017/"
				){

  //  gROOT->Reset();//should reset all globals, but it doesn't.
  //Initialise all static configurables, before working with objects (or static members) per cat-sample pair. This must be done before any other action, otherwise objects won't initialise properly. This steps neeed to be done in any macro working with configurable_mt 

  bool CleanAfterWritingDatacard=true; //important: this deletes everything declared in the plot macro. Absolutely necessary for running serially in a root session. If like to review some control canvas, put to fale, but only put to false when run in single iteration
  bool HighVerbose=false;
  
  Configurable_mt::InputFileDir=directory;
  Configurable_mt::era=era_;
  Configurable_mt::decaymode=decaymode;
  Configurable_mt::CPMethod=CPMethod;
  Configurable_mt::Observable=Observable;   
  Configurable_mt::InitObservableNames(); //It initialised what is drawn in the histogram, the obs. name for output path, and sets axes titles
  Configurable_mt::InitCatNames(); ////This will initialise a vector mt_* era_* <category> *_channel* *_CP calc method (decay plane, IP method)*. Needed to get the histogram in a dir recognisable by CombineHarvester later. We store a histo finally in catname/samplename mt_ is hardcoded in configurable_mt. We clear the vector catnames first, since we use push_back

  Configurable_mt::InitOutputBaseDir(); //All statics must have been set then correctly. This sets the base directory OutputDir
  Configurable_mt::InitOutputFileDataCard(); //initialises  OutputFileDataCard (it uses OutputDir)

  //here we initialise the bound for number cats we want to loop over
  int CATSINCLASS=Configurable_mt::CatNames.size();
  int NCats, CatCounterAddition;//NCats should contain the number of cats to loop over; either 1 or Configurable_mt::CatNames.size();
  if(CattoReview>=CATSINCLASS){throw std::invalid_argument("For the index given, no category exist for this era, exiting");}
  else if(CattoReview==-1){
    NCats=Configurable_mt::CatNames.size();
    CatCounterAddition=0;}
  else{
    NCats=1;
    CatCounterAddition=CattoReview;}
  
  int ctr=0;
  std::vector< std::vector< Configurable_mt* > > Configs(NCats, std::vector< Configurable_mt* >( Configurable_mt::SampleNames.size() ) );//also create a QCD storage object
  int EvalCtr;
  for(int catctr=0; catctr<NCats;catctr++){
    EvalCtr=catctr;
    EvalCtr+=CatCounterAddition; //either add 0, or the right index, and then loop will terminate anyway for NCats=1

    for(uint samplectr=0; samplectr<Configurable_mt::SampleNames.size();samplectr++){
      Configs[catctr][samplectr]=new Configurable_mt(EvalCtr, samplectr);
      Configs[catctr][samplectr]->Initialise();//Sets everything cat-sample dependent. assigns Cuts (also pt cuts), selection (cat, decay, calc method), weights, sample and cat name, intialises axes binning and histograms, and specific root input file name.
    }
  }
  //end Initialisation loop

  //make the output dir. Comment: for reasons not yet sorted, we must create the root file here, not earlier in the program
  gSystem->mkdir(Configurable_mt::OutputFileDataCard,kTRUE); //it will ask now to create a dir like a/b/file.root, but this seems ok
  TFile * f =new TFile(Configurable_mt::OutputFileDataCard,"recreate");
  cout<<"filling "<<Configurable_mt::OutputFileDataCard<<endl;


  //for control purposes
  TCanvas* dummy =new TCanvas("dummy","dummy",500,500);   //a dummy and control canvas. Draw method draws per default in a canvas, with dummy we avoid to overwrite the canvas for control purposes 
  TCanvas* ControlCanvas =new TCanvas("ControlCanvas","ControlCanvas",500,500);
  ControlCanvas->Divide(NCats,Configurable_mt::SampleNames.size()*2); //factor 2 is to draw os and ss histogram
  ctr=0;

  // --- start 2d loop, here we fill 2-D histograms 
  for(int catctr=0; catctr<NCats;catctr++){
    TDirectory *CatDir = f->mkdir(Configs[catctr][0]->CatName);
    //later here may initialise a loop over systematics; per systematic want to derive QCD, if loop over sys, enclose the sample loop..

    for(uint samplectr=0; samplectr<Configurable_mt::SampleNames.size()-1;samplectr++){//don't include QCD, we'll do after sample loop completed!

      cout<<"Opening Input File Name "<<Configs[catctr][samplectr]->InputFileName<<endl;//works 
      TFile* file=new TFile(Configs[catctr][samplectr]->InputFileName,"open");

      TTree * tree = (TTree*)file->Get("TauCheck"); 

      //obtain the OS histograms
      gROOT->cd(); //this will put us to the directory where the histograms are located. it is needed to draw for control purposes..
      dummy->cd();//to fetch the drawing of the histogram

      tree->Draw(Configs[catctr][samplectr]->ObservableDraw+">>"+Configs[catctr][samplectr]->HistoName2D,Configs[catctr][samplectr]->FinalCut);

      if(samplectr==0) cout<<"Configs[catctr][samplectr]->HistoName2D "<<Configs[catctr][samplectr]->HistoName2D <<"\n Configs[catctr][samplectr]->ObservableDraw "<< Configs[catctr][samplectr]->ObservableDraw<<"\n Configs[catctr][samplectr]->FinalCut "<<Configs[catctr][samplectr]->FinalCut <<"\n Configs[catctr][samplectr]->TwoDimHisto->GetSumOfWeights()  "<<Configs[catctr][samplectr]->TwoDimHisto->GetSumOfWeights()<<"\n" <<endl;
      else if(HighVerbose) cout<<"Configs[catctr][samplectr]->HistoName2D "<<Configs[catctr][samplectr]->HistoName2D <<"\n Configs[catctr][samplectr]->ObservableDraw "<< Configs[catctr][samplectr]->ObservableDraw<<"\n Configs[catctr][samplectr]->FinalCut) "<<Configs[catctr][samplectr]->FinalCut <<"\n Configs[catctr][samplectr]->TwoDimHisto->GetSumOfWeights()  "<<Configs[catctr][samplectr]->TwoDimHisto->GetSumOfWeights() <<"\n"<<endl;
      
      ctr++; 
      ControlCanvas->cd(ctr);
     if(!((catctr==0&&samplectr==0)||(catctr==1&&samplectr==0))) Configs[catctr][samplectr]->TwoDimHisto->Draw();
      ControlCanvas->Update();

      CatDir->cd();
      Configs[catctr][samplectr]->Write2DHistoToDir(CatDir);

      //identical syntax, now for SS region
      gROOT->cd();
      ctr++;
      dummy->cd();
      tree->Draw(Configs[catctr][samplectr]->ObservableDraw+">>"+Configs[catctr][samplectr]->HistoNameSS2D,Configs[catctr][samplectr]->FinalCutSS);
      ControlCanvas->cd(ctr);      
      Configs[catctr][samplectr]->TwoDimHistoSS->Draw();
      ControlCanvas->Update();

       //here add the bg samples to qcd object. Note we identify bg as not data and not 125. This MUST always be in the sample names in configurable!
      if(!(Configs[catctr][samplectr]->HistoName2D.Contains("data_obs")||Configs[catctr][samplectr]->HistoName2D.Contains("125")||Configs[catctr][samplectr]->HistoName2D.Contains("OldSignals"))){
	Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHistoSS->Add(Configs[catctr][samplectr]->TwoDimHistoSS);}      
    }

    //now the ss backgrounds have been filled. copy ss data, subtract the backgrounds
    Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHisto->Add(Configs[catctr][0]->TwoDimHistoSS);    
    Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHisto->Add(Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHistoSS,-1);
    cout<<"Constructed the QCD bg"<<endl;
    
    CatDir->cd();    
    Configs[catctr][Configurable_mt::SampleNames.size()-1]->Write2DHistoToDir(CatDir);
    //prompt QCD for quick check if nothing suspicious
    if(HighVerbose){
      cout<<"Configs[catctr][Configurable_mt::SampleNames.size()-1]->HistoName2D "<<Configs[catctr][Configurable_mt::SampleNames.size()-1]->HistoName2D <<" Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHisto->GetSumOfWeights()  "<<Configs[catctr][Configurable_mt::SampleNames.size()-1]->TwoDimHisto->GetSumOfWeights() <<endl;}
    
  }
  // --- end 2d loop  
  
  TCanvas* ControlCanvas2 =new TCanvas("ControlCanvas2","ControlCanvas2",500,500);
  ControlCanvas2->Divide(NCats,Configurable_mt::SampleNames.size());
  
  ctr=0;

  //3th loop: unroll or project to 1d histo, dependent on observables
  for(uint catctr=0; catctr<NCats;catctr++){

    TDirectory *CatDir=f->GetDirectory(Configs[catctr][0]->CatName);    
    for(int samplectr=0; samplectr<Configurable_mt::SampleNames.size();samplectr++){
      gROOT->cd(); //this will put us to the directory where the histograms are located. it is needed to draw for control purposes..

      ctr++;
      if(Observable<21) Configurable_mt::ProjectonXaxis(Configs[catctr][samplectr]->TwoDimHisto)->Copy(*Configs[catctr][samplectr]->OneDimHisto);
      else{Configurable_mt::Unfold(Configs[catctr][samplectr]->TwoDimHisto)->Copy(*Configs[catctr][samplectr]->OneDimHisto);}

      Configs[catctr][samplectr]->Write1DHistoToDir(CatDir);      
      ControlCanvas2->cd(ctr);
      //draw, but don't unblind data in signal cats
      if(!((catctr==0&&samplectr==0)||(catctr==1&&samplectr==0))) Configs[catctr][samplectr]->OneDimHisto->Draw();      
      ControlCanvas2->Update();           
    }
  }
  //end 3th loop 
  f->Write();

  //control loop if we didn't pick too fine binning for the CP angl ehistograms -> then first need to add all backgrounds.. exclude data and signals
  //We currently do this only if we run over a single category for the CP angle distributions..


  TCanvas * SumBackgroundsBinRatioCanvas =new TCanvas("SumBackgroundsBinRatioCanvas","SumBackgroundsBinRatioCanvas");
    
  if(NCats==1&&Observable==21){
    TH1F* SumBackgrounds=new TH1F();
    Configs[0][0]->OneDimHisto->Copy(*SumBackgrounds);
    SumBackgrounds->Reset();
    
    for(uint catctr=0; catctr<NCats;catctr++){
      for(uint samplectr=0; samplectr<Configurable_mt::SampleNames.size();samplectr++){
	
	if(!(Configs[catctr][samplectr]->HistoName2D.Contains("data_obs")||Configs[catctr][samplectr]->HistoName2D.Contains("125")||Configs[catctr][samplectr]->HistoName2D.Contains("OldSignals"))){	
	  SumBackgrounds->Add(Configs[catctr][samplectr]->OneDimHisto);
	}  
      }  
    }
    
    cout<<"********* Checking for bins with bad sigma/content ratio *********"<<endl;
    TH1F* SumBackgroundsBinRatio=new TH1F("SumBackgroundsBinRatio","SumBackgroundsBinRatio",200,-1,1);
    
    int NBinsX=SumBackgrounds->GetXaxis()->GetNbins();
  for(int index=1;index<NBinsX;index++){
    double Ratio=SumBackgrounds->GetBinError(index)/SumBackgrounds->GetBinContent(index);
    SumBackgroundsBinRatio->Fill(Ratio);
    if(Ratio>0.5){
      cout<<"BIN PROBLEM:Ratio "<<Ratio<<endl;
      cout<<"index "<<index<<endl;
    }
    //look if empty bin in data. Only do for non-signal categories 
    if(Configs[0][0]->OneDimHisto->GetBinContent(index)==0 && CattoReview>1) cout<<"data seems to have empty bin, index "<<index <<endl;
  }
  
  cout<<endl;
  cout<<"********* Finished for bins with bad sigma/content ratio *********"<<endl;
  
  SumBackgrounds->SetName("SumBackgrounds");
  SumBackgrounds->SetTitle("SumBackgrounds");
  f->cd();
  SumBackgrounds->Write();
  SumBackgroundsBinRatio->Write();

  SumBackgroundsBinRatio->SetDirectory(0);
  SumBackgroundsBinRatio->Draw();
  SumBackgroundsBinRatioCanvas->Update();
  
  }

  f->Close();
  delete dummy;


  if(CleanAfterWritingDatacard){    
    //Now clean up the objects, if run mutliple times they're not cleared automatically
    for(int catctr=0; catctr<NCats;catctr++){
      for(int samplectr=0; samplectr<Configurable_mt::SampleNames.size();samplectr++){
	delete Configs[catctr][samplectr];
      }}    
    delete ControlCanvas;
    delete ControlCanvas2;
    delete SumBackgroundsBinRatioCanvas;
  }        

  cout<<"create file in Configurable_mt::OutputFileDataCard  "<<Configurable_mt::OutputFileDataCard<<endl;
  return;
}

