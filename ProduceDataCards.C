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

#include "Configurable_mt.h"

using namespace std;

void ProduceDataCards(
				int PredictedClass=-1,//-1: loop over all NN classes available. N: compute datacard only of chosen class
				int era_=2017,
				int decaymode=1,//if initialise to -1, no selection made!
				int CPMethod=1, //initialised the charged particle cutoff and the appendix for the cat.
				int Observable=21,
				//TString directory = "/Users/klundert/DESY/HiggsCPProjectSoftware/Outputs/2019_11_17_Unfiltered_NewTopNorm/predictions_2017/", //latest versionbefore switchng to deeptau
				TString directory = "/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HiggsCP/Outputs/test_1S_3B/NTuples_mt_2017/predictions_2017/"
				){

  vector<TString> SampleList={"data","Zll","TT","ST","VV","W","ggHSM", "ggHPS", "ggHMM","qqHSM", "qqHPS", "qqHMM"};
  vector<int> catIndices={0,1,2,3,-1};
  const int nCategories=4;
  const int nSignCategories=1;
  const bool useEmbedded=false;
  const bool FFmethod=false;
  
  if(useEmbedded)SampleList.push_back("embedded");
  else SampleList.push_back("Ztt");
  if(FFmethod)SampleList.push_back("fakes");
  else SampleList.push_back("QCD");
  const int nSamples=SampleList.size();
  

  //  gROOT->Reset();//should reset all globals, but it doesn't.
  //Initialise all static configurables, before working with objects (or static members) per cat-sample pair. This must be done before any other action, otherwise objects won't initialise properly. This steps neeed to be done in any macro working with configurable_mt 

  bool CleanAfterWritingDatacard=true; //important: this deletes everything declared in the plot macro. Absolutely necessary for running serially in a root session. If like to review some control canvas, put to fale, but only put to false when run in single iteration
  bool HighVerbose=false;


  Configurable_mt::useEmbedded=useEmbedded;
  Configurable_mt::FFmethod=FFmethod;


  Configurable_mt::nSignalCategories=nSignCategories;
  Configurable_mt::nCategories=nCategories;
  Configurable_mt::CommonCuts="(byVVLooseDeepTau2017v2p1VSe_2 >0.5 && byTightDeepTau2017v2p1VSmu_2 >0.5 && puppimt_1<50)";
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
  //int nCategories, CatCounterAddition;//nCategories should contain the number of cats to loop over; either 1 or nCategories;

  if(PredictedClass>=nCategories){throw std::invalid_argument("For the index given, no category exist for this era, exiting");}

  int ctr=0;
  //  std::vector< std::vector< Configurable_mt* > > Configs(nCategories, std::vector< Configurable_mt* >(nSamples));//also create a QCD storage object
  map<pair<int,TString>,Configurable_mt*> Configs;
  
  int EvalCtr;

  for(int catIndex : catIndices){
    unsigned int ISample=0;
    for(TString sample : SampleList){
      Configs.insert({std::make_pair(catIndex,sample),new Configurable_mt(catIndex, sample)});
    Configs.at(std::make_pair(catIndex,sample))->Initialise();//Sets everything cat-sample dependent. assigns Cuts (also pt cuts), selection (cat, decay, calc method), weights, sample and cat name, intialises axes binning and histograms, and specific root input file name.
      ISample++;
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
  ControlCanvas->Divide(nCategories,nSamples*2); //factor 2 is to draw os and ss histogram
  ctr=0;

  // --- start 2d loop, here we fill 2-D histograms 
  for(int catIndex : catIndices){
    TDirectory *CatDir = f->mkdir(Configs.at(std::make_pair(catIndex,"data"))->CatName);
    //later here may initialise a loop over systematics; per systematic want to derive QCD, if loop over sys, enclose the sample loop..
    for(TString sample : SampleList){
      cout << sample << endl;

      if(sample=="QCD"||sample=="fakes")continue;
      //for(uint samplectr=0; samplectr<nSamples-1;samplectr++){//don't include QCD, we'll do after sample loop completed!

      cout<<"Opening Input File Name "<<Configs.at(std::make_pair(catIndex,sample))->InputFileName<<endl;//works 
      TFile* file=new TFile(Configs.at(std::make_pair(catIndex,sample))->InputFileName,"open");

      TTree * tree = (TTree*)file->Get("TauCheck"); 

      //obtain the OS histograms
      gROOT->cd(); //this will put us to the directory where the histograms are located. it is needed to draw for control purposes..
      dummy->cd();//to fetch the drawing of the histogram

      tree->Draw(Configs.at(std::make_pair(catIndex,sample))->ObservableDraw+">>"+Configs.at(std::make_pair(catIndex,sample))->HistoName2D,Configs.at(std::make_pair(catIndex,sample))->SRCut);

      if(sample=="data") cout<<"Configs.at(std::make_pair(catIndex,sample))->HistoName2D "<<Configs.at(std::make_pair(catIndex,sample))->HistoName2D <<"\n Configs.at(std::make_pair(catIndex,sample))->ObservableDraw "<< Configs.at(std::make_pair(catIndex,sample))->ObservableDraw<<"\n Configs.at(std::make_pair(catIndex,sample))->SRCut "<<Configs.at(std::make_pair(catIndex,sample))->SRCut <<"\n Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto->GetSumOfWeights()  "<<Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto->GetSumOfWeights()<<"\n" <<endl;
      else if(HighVerbose) cout<<"Configs.at(std::make_pair(catIndex,sample))->HistoName2D "<<Configs.at(std::make_pair(catIndex,sample))->HistoName2D <<"\n Configs.at(std::make_pair(catIndex,sample))->ObservableDraw "<< Configs.at(std::make_pair(catIndex,sample))->ObservableDraw<<"\n Configs.at(std::make_pair(catIndex,sample))->SRCut) "<<Configs.at(std::make_pair(catIndex,sample))->SRCut <<"\n Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto->GetSumOfWeights()  "<<Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto->GetSumOfWeights() <<"\n"<<endl;
      
      ctr++; 
      ControlCanvas->cd(ctr);
      if(!(sample=="data"&&(catIndex<nSignCategories||catIndex==-1))) Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto->Draw();
      ControlCanvas->Update();

      CatDir->cd();
      Configs.at(std::make_pair(catIndex,sample))->Write2DHistoToDir(CatDir);

      //identical syntax, now for SS region
      gROOT->cd();
      ctr++;
      dummy->cd();
      tree->Draw(Configs.at(std::make_pair(catIndex,sample))->ObservableDraw+">>"+Configs.at(std::make_pair(catIndex,sample))->HistoNameSS2D,Configs.at(std::make_pair(catIndex,sample))->SSCut);
      if(FFmethod)tree->Draw(Configs.at(std::make_pair(catIndex,sample))->ObservableDraw+">>"+Configs.at(std::make_pair(catIndex,sample))->HistoNameAR2D,Configs.at(std::make_pair(catIndex,sample))->ARCut);
      ControlCanvas->cd(ctr);      
      Configs.at(std::make_pair(catIndex,sample))->TwoDimHistoSS->Draw();
      if(FFmethod)Configs.at(std::make_pair(catIndex,sample))->TwoDimHistoAR->Draw();
      ControlCanvas->Update();

       //here add the bg samples to qcd object. Note we identify bg as not data and not 125. This MUST always be in the sample names in configurable!
      if(!(Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("data_obs")||Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("125")||Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("OldSignals"))){
	Configs.at(std::make_pair(catIndex,"QCD"))->TwoDimHistoSS->Add(Configs.at(std::make_pair(catIndex,sample))->TwoDimHistoSS);
        if(FFmethod)Configs.at(std::make_pair(catIndex,"fakes"))->TwoDimHistoAR->Add(Configs.at(std::make_pair(catIndex,sample))->TwoDimHistoAR);
      }
    }

    //now the ss backgrounds have been filled. copy ss data, subtract the backgrounds
    Configs.at(std::make_pair(catIndex,"QCD"))->TwoDimHisto->Add(Configs.at(std::make_pair(catIndex,"data"))->TwoDimHistoSS);    
    Configs.at(std::make_pair(catIndex,"QCD"))->TwoDimHisto->Add(Configs.at(std::make_pair(catIndex,"QCD"))->TwoDimHistoSS,-1);
    cout<<"Constructed the QCD bg"<<endl;
    
    if(FFmethod){
      Configs.at(std::make_pair(catIndex,"fakes"))->TwoDimHisto->Add(Configs.at(std::make_pair(catIndex,"data"))->TwoDimHistoAR);    
      Configs.at(std::make_pair(catIndex,"fakes"))->TwoDimHisto->Add(Configs.at(std::make_pair(catIndex,"fakes"))->TwoDimHistoAR,-1);
    cout<<"Constructed the fakes bg"<<endl;
   

    }


    CatDir->cd();    
    Configs.at(std::make_pair(catIndex,"QCD"))->Write2DHistoToDir(CatDir);
    if(FFmethod)Configs.at(std::make_pair(catIndex,"fakes"))->Write2DHistoToDir(CatDir);
    //prompt QCD for quick check if nothing suspicious
    if(HighVerbose){
      cout<<"Configs[catIndex][nSamples-1]->HistoName2D "<<Configs.at(std::make_pair(catIndex,"QCD"))->HistoName2D <<" Configs[catIndex][nSamples-1]->TwoDimHisto->GetSumOfWeights()  "<<Configs.at(std::make_pair(catIndex,"QCD"))->TwoDimHisto->GetSumOfWeights() <<endl;}
    
  }
  // --- end 2d loop  
  
  TCanvas* ControlCanvas2 =new TCanvas("ControlCanvas2","ControlCanvas2",500,500);
  ControlCanvas2->Divide(nCategories,nSamples);
  
  ctr=0;

  //3th loop: unroll or project to 1d histo, dependent on observables
  for(uint catIndex : catIndices){

    TDirectory *CatDir=f->GetDirectory(Configs.at(std::make_pair(catIndex,"data"))->CatName);    
    for(TString sample : SampleList){
      gROOT->cd(); //this will put us to the directory where the histograms are located. it is needed to draw for control purposes..

      ctr++;
      if(Observable<21) Configurable_mt::ProjectonXaxis(Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto)->Copy(*Configs.at(std::make_pair(catIndex,sample))->OneDimHisto);
      else{Configurable_mt::Unfold(Configs.at(std::make_pair(catIndex,sample))->TwoDimHisto)->Copy(*Configs.at(std::make_pair(catIndex,sample))->OneDimHisto);}

      Configs.at(std::make_pair(catIndex,sample))->Write1DHistoToDir(CatDir);      
      ControlCanvas2->cd(ctr);
      //draw, but don't unblind data in signal cats
      if(!(sample=="data"&&(catIndex<nSignCategories||catIndex==-1))) Configs.at(std::make_pair(catIndex,sample))->OneDimHisto->Draw();      
      ControlCanvas2->Update();           
    }
  }
  //end 3th loop 
  f->Write();

  //control loop if we didn't pick too fine binning for the CP angl ehistograms -> then first need to add all backgrounds.. exclude data and signals
  //We currently do this only if we run over a single category for the CP angle distributions..


  TCanvas * SumBackgroundsBinRatioCanvas =new TCanvas("SumBackgroundsBinRatioCanvas","SumBackgroundsBinRatioCanvas");
    
  if(PredictedClass==-1&&Observable==21){
    TH1F* SumBackgrounds=new TH1F();
    Configs.at(std::make_pair(0,"data"))->OneDimHisto->Copy(*SumBackgrounds);
    SumBackgrounds->Reset();
    
    for(int catIndex : catIndices){
      for(TString sample : SampleList){
	
	if(!(Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("data_obs")||Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("125")||Configs.at(std::make_pair(catIndex,sample))->HistoName2D.Contains("OldSignals"))){	
	  SumBackgrounds->Add(Configs.at(std::make_pair(catIndex,sample))->OneDimHisto);
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
    if(Configs.at(std::make_pair(0,"data"))->OneDimHisto->GetBinContent(index)==0 && PredictedClass>1) cout<<"data seems to have empty bin, index "<<index <<endl;
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
    for(int catIndex : catIndices){
      for(TString sample : SampleList){
	delete Configs.at(std::make_pair(catIndex,sample));
      }}    
    delete ControlCanvas;
    delete ControlCanvas2;
    delete SumBackgroundsBinRatioCanvas;
  }        

  cout<<"create file in Configurable_mt::OutputFileDataCard  "<<Configurable_mt::OutputFileDataCard<<endl;
  return;
}

