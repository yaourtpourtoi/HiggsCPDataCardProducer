
/*
Plotting tool to plot histograms from datacards. It uses Configurable_mt. 
Written by Merijn van de klundert, mvandekl@cern.ch
*/

#include <iostream>
#include <vector>
#include <map>
#include <iomanip>
#include "Plotting.h"
#include "Plotting_Style.h"
#include "HttStylesNew.cc"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TEfficiency.h"
#include "TMath.h"

#include "TH2F.h"
#include "TSystem.h"
#include "Configurable_mt.h"


void   PlotDataCards(
		       int CattoReview=8,
		       int era_=2017,
		       int decaymode=0,//if initialise to -1, no selection made!
		       int CPMethod=0,
		       int Observable=21){

  if(CattoReview==-1||CattoReview>8){throw std::invalid_argument("CattoReview cannot be -1 or larger than 8 when making a plot");}

  // Initialise a few settings
  bool CleanafterWritingPDF=true; //important: this deletes everything declared with new in the plot macro. Good practise for running serially, but can't inspect the plot in root afterwards anymore. only put to false when run the plotmacro single iteration
  
  bool blindData = true;
  bool logY = true;
  if(Observable<20) logY=false;
  
  if(CattoReview>1) blindData = false;//Currently, fully blind for ggH and qqH category
  int nbMin =0;    //bins used for blinding
  int nbMax = 500;

  bool plotLeg = true;
  bool showSignal = true;

  //here set static configurables. This must be done before any other action, otherwise objects won't initialise properly
  cout<<"init from macro"<< endl;
  Configurable_mt::era=era_;
  Configurable_mt::decaymode=decaymode;
  Configurable_mt::CPMethod=CPMethod;
  Configurable_mt::Observable=Observable; //should find a trick to init obsevable string..
  
  Configurable_mt::InitObservableNames(); //call this to initialise the observable name for the path and for draw method. Per default this observable is plotted against DNN score. will init static ObservableDraw and ObservablePath;
  Configurable_mt::InitCatNames(); //call this first to add decay mode and calc method to the category string
  Configurable_mt::InitOutputBaseDir(); //All statics must have been set then correctly. This sets the base directory OutputDir
  Configurable_mt::InitOutputFileDataCard(); //initialises  OutputFileDataCard (needed to fetch data card).
  Configurable_mt::InitOutputPrefitPlotDir(); //inits OutputDirPreFit needed to write the plot somewhere
  cout<<"Configurable_mt::OutputFileDataCard "<<Configurable_mt::OutputFileDataCard<<endl;
  
  //open the root file for the static parameters set above
  TFile * f =new TFile(Configurable_mt::OutputFileDataCard,"open");
  //we now retrieve the histograms from the data cards and assign the the 
  TH1D * histData=new TH1D(); 
  TH1D * QCD=new TH1D();     
  TH1D * Ztt=new TH1D();      
  TH1D * Zll=new TH1D();      
  TH1D * W=new TH1D();        
  TH1D * TT=new TH1D();       
  TH1D * VV=new TH1D();      
  TH1D * EWK=new TH1D();//ewk is propagated throughout the macro, may want to merge vv and W later into EWK       
  TH1D * SMH=new TH1D();     
  TH1D * SMHVBF=new TH1D();     
  TH1D * SMHOdd=new TH1D();     
  TH1D * SMHVBFOdd=new TH1D();    
  TH1D * SMHMixed=new TH1D();
  TH1D * SMHVBFMixed=new TH1D();
  
  int nBins=0;
  int NSamples=Configurable_mt::SampleNames.size();
  
  TH2F* TwoDHistforBinInfo;
  for(int samplectr=0; samplectr<NSamples;samplectr++){
    TString CatName=Configurable_mt::CatNames[CattoReview];
    TString SampleName=Configurable_mt::SampleNames[samplectr];
    TString HistoName=CatName+"/"+SampleName;

    f->cd(CatName); //go to correct directory in the file

    cout<<"fetching: "<< HistoName <<endl;    
    TH1F* temp=(TH1F*)f->Get(HistoName);
    nBins=temp->GetXaxis()->GetNbins();
    
    //here the histograms still have identical name, so norm to binwidth for NN score histograms easiest here.
    if(Observable==20){
      for(int binctr=1; binctr<nBins+1;binctr++){
	temp->SetBinContent(binctr,temp->GetBinContent(binctr)/temp->GetBinWidth(binctr));
	temp->SetBinError(binctr,temp->GetBinError(binctr)/temp->GetBinWidth(binctr));}}
    
    //here match now, safest is by name:    
    if(SampleName=="data_obs") temp->Copy(*histData);
    if(SampleName=="Ztt") temp->Copy(*Ztt);
    if(SampleName=="Zll") temp->Copy(*Zll);
    if(SampleName=="W") temp->Copy(*W);
    if(SampleName=="TT") temp->Copy(*TT);
    if(SampleName=="VV") temp->Copy(*VV);
    if(SampleName=="ggH_sm_htt125") temp->Copy(*SMH);
    if(SampleName=="ggH_ps_htt125") temp->Copy(*SMHOdd);
    if(SampleName=="ggH_mm_htt125") temp->Copy(*SMHMixed);    
    if(SampleName=="qqH_sm_htt125") temp->Copy(*SMHVBF);
    if(SampleName=="qqH_ps_htt125") temp->Copy(*SMHVBFOdd);
    if(SampleName=="qqH_mm_htt125") temp->Copy(*SMHVBFMixed);
    if(SampleName=="QCD") temp->Copy(*QCD);

    //this is useful for obtaining the NN bin edges later, when plotting tlines
    TString ExampleHisto2D=HistoName+"_2D";
    if(samplectr==0) TwoDHistforBinInfo=(TH2F*)f->Get(ExampleHisto2D);
  }    
 

  // ************************************
  // ***** Summarize backgrounds  *******
  // ************************************
  
  cout << endl;
  cout << "DATA : Sum of weights = " << histData->GetSumOfWeights() << " : Integral = " << histData->Integral(1,nBins+1) << endl;
  cout << "QCD : Sum of weights = " << QCD->GetSumOfWeights() << " : Integral = " << QCD->Integral(1,nBins+1) << endl;
  cout << "VV  : Sum of weights = " << VV->GetSumOfWeights()  << " : Integral = " << VV->Integral(1,nBins+1)  << endl;
  cout << "W   : Sum of weights = " << W->GetSumOfWeights()   << " : Integral = " << W->Integral(1,nBins+1)   << endl;
  //  cout << "EWK  : Sum of weights = " << EWK->GetSumOfWeights()  << " : Integral = " << EWK->Integral(1,nBins+1)  << endl;
  cout << "TT  : Sum of weights = " << TT->GetSumOfWeights()  << " : Integral = " << TT->Integral(1,nBins+1)  << endl;
  cout << "Ztt : Sum of weights = " << Ztt->GetSumOfWeights() << " : Integral = " << Ztt->Integral(1,nBins+1) << endl;
  cout << "Zll : Sum of weights = " << Zll->GetSumOfWeights() << " : Integral = " << Zll->Integral(1,nBins+1) << endl;


  // ***********************************
  // **** Initialise ad-hoc systematic uncertainties *****
  // ***********************************

  TH1D * dummy = (TH1D*)Ztt->Clone("dummy");
  float errQCD = 0.15; // ad-hoc sys uncertainty of QCD background
  float errVV  = 0.15; // ad-hoc sys uncertainty of VV background
  float errW   = 0.10; // ad-hoc sys uncertainty of W+Jets background
  float errTT  = 0.07; // ad-hoc sys uncertainty of TT background
  float errEWK   = 0.12; // ad-hoc sys uncertainty of EWK background
  
  for (int iB=1; iB<=nBins; ++iB) {   // Add general systematic uncertainties to each bin as error
    float eQCD   = errQCD*QCD->GetBinContent(iB);
    float eVV    = errVV*VV->GetBinContent(iB);
    float eW     = errW*W->GetBinContent(iB);
    //float eEWK     = errEWK*EWK->GetBinContent(iB);
    float eTT    = errTT*TT->GetBinContent(iB);
    float err2   = eQCD*eQCD+eVV*eVV + eW*eW + eTT*eTT;     //TO DO: WAS IST MIT DEM FEHLER AUF DY?
    //float err2   = eQCD*eQCD+eEWK*eEWK + eTT*eTT;     //TO DO: WAS IST MIT DEM FEHLER AUF DY?
    float errTot = TMath::Sqrt(err2);
    //    cout << "eQCD: " << eQCD << "  eEWK: " << eEWK << "  eTT: " << eTT << "  eTotal: " << errTot << endl;
    dummy -> SetBinError(iB,errTot);
    SMH   -> SetBinError(iB,0);
    SMHVBF-> SetBinError(iB,0);
    SMHOdd-> SetBinError(iB,0);
    SMHVBFOdd-> SetBinError(iB,0);
    SMHMixed-> SetBinError(iB,0);
    SMHVBFMixed-> SetBinError(iB,0);
  }
  cout << endl;
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // FINAL PLOTTING
  ///////////////////////////////////////////////////////////////////////////////////////////////

  ModTDRStyle();
  bool LargeScale = true;  
  TCanvas* canv1;
  if(LargeScale)canv1 = new TCanvas("c1", "c1", 2000,800);
  else canv1 = new TCanvas("c1", "c1", 1200,1000);
  canv1->cd();
  vector<TPad*> pads = TwoPadSplit(0.29, 0.00, 0.00);
  pads[0]->SetLogy(logY);

  vector<TH1*> h = CreateAxisHists(2, histData, histData->GetXaxis()->GetXmin(), histData->GetXaxis()->GetXmax()-0.01);
  h[0]->GetYaxis()->SetMaxDigits(3);
  h[0] -> GetYaxis()->SetTitle(Configurable_mt::ObservableYaxis);
  pads[0]->cd();
  h[0] -> Draw();
  pads[0]->Update();
  pads[1] -> cd();
  h[1]    -> Draw();
  SetupTwoPadSplitAsRatio(pads, "Obs/Exp", true, 0.65, 1.35);

  h[1] -> GetXaxis()->SetTitle(Configurable_mt::ObservableXaxis);
  h[1] -> GetYaxis()->SetRangeUser(0.5,2);
  h[1] -> GetYaxis()->SetNdivisions(4);
  h[1] -> GetXaxis()->SetTitleOffset(0.95);
  h[1] -> GetXaxis()->SetNdivisions(505);
  h[1] -> GetYaxis()->SetTitleOffset(1.1);
  if(LargeScale)  h[1] -> GetYaxis()->SetTitleOffset(0.9);
  pads[0] -> cd();
  
  h[0] -> GetYaxis()->SetTitleOffset(1.6);
  if(LargeScale)h[0] -> GetYaxis()->SetTitleOffset(1);
  pads[1] -> SetGrid(0,1);
  if(logY) h[0] -> SetMinimum(1);
  pads[0] -> cd();

  TLegend *legend = PositionedLegend(0.5, 0.15, 3, 0.03); //widht, height, position, offset, see plotting.h
  legend -> SetTextFont(42);
  legend-> SetNColumns(3);

  pads[1] -> cd();
  TLegend *legendRatio = PositionedLegend(0.5, 0.04, 2, 0.03); //widht, height, position, offset, see plotting.h
  legendRatio -> SetTextFont(42);
  legendRatio-> SetNColumns(3);

  pads[0] -> cd();
  histData -> SetMarkerColor(1);
  histData -> SetLineColor(1);
  histData -> SetFillColor(1);
  histData -> SetFillStyle(0);
  histData -> SetLineWidth(2);
  histData -> SetMarkerStyle(20);
  histData -> SetMarkerSize(1.1);

  //assign colours, use colours coherent with IC: https://github.com/danielwinterbottom/ICHiggsTauTau/blob/master/python/plotting.py#L2158-L2183
  InitHist(QCD,TColor::GetColor(250,202,255));
  InitHist(Ztt,TColor::GetColor(248,206,104));
  InitHist(Zll,TColor::GetColor(100,192,232));
  InitHist(TT,TColor::GetColor(155,152,204));
  // InitHist(EWK,TColor::GetColor(222,90,106));
  InitHist(W,TColor::GetColor(222,90,106));
  InitHist(VV,TColor::GetColor("#6F2D35")); //no colur in IC list, for them its merged into EWK
  
  // Add all bkg contributions to one stack plot. Note that this order is coherent with SM analysis material
  THStack *stack = new THStack("Background",""); //setting the HIST option avoids that the individual errors are drawn..
  stack -> Add(QCD,"HIST");
  stack -> Add(VV,"HIST");
  //  stack -> Add(EWK,"HIST");
  stack -> Add(W,"HIST");
  stack -> Add(TT,"HIST");
  stack -> Add(Zll,"HIST");
  stack -> Add(Ztt,"HIST");
  stack -> Draw("hsame"); //it was hsame
  canv1->Update();
  
  InitSignal(SMH ,2);
  InitSignal(SMHVBF ,3);
  InitSignal(SMHOdd,4);
  InitSignal(SMHVBFOdd,5);
  InitSignal(SMHMixed,6);
  InitSignal(SMHVBFMixed,7); 
  
  if (showSignal){
    pads[0] -> cd();
    
    SMH->Draw("hsame");
    SMHVBF->Draw("hsame");
    SMHOdd->Draw("hsame");
    SMHVBFOdd->Draw("hsame");
    
    /* Mixed sampels. Per default won't display for visibility
       SMHMixed->Draw("hsame");
       SMHVBFMixed->Draw("hsame"); */  
  }
  
  canv1->Update();
  canv1->cd();

  //blind data here if needed. Must be done here, since data is used below in manipulations
  if (blindData){
    for (int iB=nbMin; iB<=nbMax; ++iB){
      histData->SetBinContent(iB,-1);
      histData->SetBinError(iB,0);}}
  
  // Initialize a histogram which adds all error up
  TH1D * bkgdErr = (TH1D*)stack->GetStack()->Last()->Clone("bkgdErr");
  float errLumi = 0.03;
  float errMuon = 0.03;
  float errElectron = 0.04;
  for (int iB=1; iB<=nBins; ++iB) {
    QCD->SetBinError(iB,0);
    VV->SetBinError(iB,0);
    TT->SetBinError(iB,0);
    W->SetBinError(iB,0);
    //  EWK->SetBinError(iB,0);
    Ztt->SetBinError(iB,0);
    Zll->SetBinError(iB,0);

    float eStat =  bkgdErr->GetBinError(iB);
    float X = bkgdErr->GetBinContent(iB);
    float eLumi = errLumi * X;
    float eMuon = errMuon * X;
    float eElectron = errElectron * X;
    float eBkg = dummy->GetBinError(iB);
    float Err = TMath::Sqrt(eStat*eStat+eLumi*eLumi+eBkg*eBkg+eMuon*eMuon+eElectron*eElectron);
    bkgdErr->SetBinError(iB,Err);
    cout<<"BIN: "<<iB << " eStat = " << eStat << " : eLumi = "<< eLumi <<" : eBkg = " << eBkg << " : eTOT = "<<Err<<endl;
    }
  
  canv1->cd();
  pads[0] -> cd();

  bkgdErr -> SetMarkerSize(0);
  int new_idx = 923;
  bkgdErr -> SetFillColor(new_idx);
  bkgdErr -> SetFillStyle(3002);//3004 goes also, but it isn't very nice in the ratio plot
  bkgdErr -> SetLineWidth(1);
  bkgdErr -> Draw("e2same");
  canv1   -> Update();

  TH1D * ratioH    = (TH1D*)histData -> Clone("ratioH");
  TH1D * ratioErrH = (TH1D*)bkgdErr  -> Clone("ratioErrH");
  ratioH -> Divide((TH1D*)stack->GetStack()->Last()); //Divide by the sum of the THStack
  
  // Set error of MC bkg correctly in ratio
  for (int iB=1; iB<=nBins; ++iB) {
    ratioErrH -> SetBinContent(iB,1.0);
    ratioErrH -> SetBinError(iB,0.0);
    float xBkg   = bkgdErr -> GetBinContent(iB);
    float errBkg = bkgdErr -> GetBinError(iB);
    if (xBkg>0) {
      float relErr = errBkg/xBkg;
      ratioErrH->SetBinError(iB,relErr);
    }
  }

  pads[1]->cd();

  ratioErrH -> SetFillColor(new_idx);
  ratioErrH -> SetFillStyle(3002);
  ratioErrH -> SetLineWidth(1);
  ratioErrH->Draw("e2same");
  ratioH->Draw("pe0same");

  //Merijn: include the signal ratio w.r.t bg, i.e. s+b/b. Won't propagate or display errors in rescaled signal.
  TH1D * ratioS_ggH_SM    = (TH1D*)SMH -> Clone("ratioS_ggH_SM");
  TH1D * ratioS_VBF_SM    = (TH1D*)SMHVBF -> Clone("ratioS_VBF_SM");
  TH1D * ratioS_ggH_BSM   = (TH1D*)SMHOdd -> Clone("ratioS_ggH_BSM");
  TH1D * ratioS_VBF_BSM   = (TH1D*)SMHVBFOdd -> Clone("ratioS_VBF_BSM");
  
  ratioS_ggH_SM->Add((TH1D*)stack->GetStack()->Last());
  ratioS_ggH_BSM->Add((TH1D*)stack->GetStack()->Last());
  ratioS_VBF_SM->Add((TH1D*)stack->GetStack()->Last());
  ratioS_VBF_BSM->Add((TH1D*)stack->GetStack()->Last());
  
  ratioS_ggH_SM->Divide((TH1D*)stack->GetStack()->Last());
  ratioS_ggH_BSM->Divide((TH1D*)stack->GetStack()->Last());
  ratioS_VBF_SM->Divide((TH1D*)stack->GetStack()->Last());
  ratioS_VBF_BSM->Divide((TH1D*)stack->GetStack()->Last());
  
  if(Observable!=21&&Observable!=22){
    ratioS_ggH_SM->Draw("hist same");
    ratioS_VBF_SM->Draw("hist same");
    legendRatio->AddEntry(ratioS_ggH_SM,"SM ggH (125)+bkg","l");
    legendRatio->AddEntry(ratioS_VBF_SM,"SM VBF (125)+bkg","l");}
  
  if((Observable==21||Observable==22)&&CattoReview!=1){
    ratioS_ggH_SM->Draw("hist same");
    ratioS_ggH_BSM->Draw("hist same");
    legendRatio->AddEntry(ratioS_ggH_SM,"SM ggH (125)+bkg","l");
    legendRatio->AddEntry(ratioS_ggH_BSM,"BSM ggH (125)+bkg","l");}
  
  if((Observable==21||Observable==22)&&CattoReview==1){
    ratioS_VBF_SM->Draw("hist same");
    ratioS_VBF_BSM->Draw("hist same");
    legendRatio->AddEntry(ratioS_VBF_SM,"SM VBF (125)+bkg","l");
    legendRatio->AddEntry(ratioS_VBF_BSM,"BSM VBF (125)+bkg","l");}
  
  legendRatio->Draw();
  pads[0]->cd();
  histData->Draw("pesame");

  //Routine to plot Tlines and NN score bins for observable 21
  if(Observable==21){
    double xcoor;
    int NBinsphiCP=TwoDHistforBinInfo->GetXaxis()->GetNbins();
    for(int bindex=0;bindex<TwoDHistforBinInfo->GetYaxis()->GetNbins();bindex++){
      
      
    }        
  }

  //some cosmetics
  FixTopRange(pads[0], GetPadYMax(pads[0])*1, 0.3); //Merijn: factor 2 fives better results in some cases..
  if (era_== 2016) DrawTitle(pads[0], "35.9 fb^{-1} (13 TeV, 2017)", 3);
  if (era_== 2017) DrawTitle(pads[0], "41.9 fb^{-1} (13 TeV, 2017)", 3);
  if (era_== 2018) DrawTitle(pads[0], "59.7 fb^{-1} (13 TeV, 2017)", 3);
  DrawTitle(pads[0], "#scale[1.2]{         #bf{CMS} Work in progress}", 1);
  FixBoxPadding(pads[0], legend, 0.05);

  pads[0]->cd(); //Draw here legend. do it here, since all configs need to be set first. THis order is coherent with SM results
  legend -> AddEntry(Ztt,"Ztt","f");
  legend -> AddEntry(Zll,"Zll","f");  
  legend -> AddEntry(TT,"t#bar{t}","f");
  //  legend -> AddEntry(EWK,"EWK","f");
  legend -> AddEntry(VV,"VV","f");
  legend -> AddEntry(W,"W","f");
  legend -> AddEntry(QCD,"QCD","f");
  legend  -> AddEntry(bkgdErr, "Bkg. uncertainty" , "F" );
  
  legend->AddEntry(SMH,"SM ggH (125)","l");
  legend->AddEntry(SMHVBF,"SM VBF (125)","l");
  legend->AddEntry(SMHOdd,"Pseudo ggH (125)","l");
  legend->AddEntry(SMHVBFOdd,"Pseudo VBF (125)","l");
  
  /*    if desired add legends for mixed samples
	legend->AddEntry(SMHMixed,"MaxMix(125) #times 1","f");
	legend->AddEntry(SMHVBFMixed,"MaxMix(125) VBF #times 1","f");*/  
  legend -> AddEntry(histData, "Observed", "ple");

  
  /* if desired, draw sum of backgrounds. For NN, first normalise to binwidth of course
  TH1F* SumBackgrounds=(TH1F*)f->Get("SumBackgrounds");
  SumBackgrounds->Draw("same"); //for controll
  legend->AddEntry(SumBackgrounds,"Sum backgrounds","f");
  */

  legend->Draw();
  FixOverlay();
  canv1->Update();
  pads[0]->GetFrame()->Draw();

  //Set final range y-axis
  double MaximalYExpected=bkgdErr->GetBinContent(bkgdErr->GetMaximumBin());
  h[0]->GetYaxis()->SetRangeUser(0.1,1000*MaximalYExpected);
  if(logY==false) h[0]->GetYaxis()->SetRangeUser(0.1,1.5*MaximalYExpected);
  pads[0]->Update(); 

  //draw a string such as mt_2017_ggH_mupi_IP, such that we can identify what is drawn..  
  TString category=Configurable_mt::CatNames[CattoReview];
  TLegend* justcatindexlegend=new TLegend(0.2,0.75,0.4,0.9);
  justcatindexlegend->SetHeader(category);
  justcatindexlegend->Draw();

  //make a directory if not yet exisiting, draw save the canvas
  gSystem->mkdir(Configurable_mt::OutputDirPreFit,kTRUE); //it will ask now recursively to create a dir like a/b/file.root. Nonetheless of the .root extension it works
  canv1 -> Print( Configurable_mt::OutputDirPreFit + "MuTau_" + Configurable_mt::ObservablePath +"_Cat_" + Configurable_mt::BaseCatNames[CattoReview] + ".pdf" );

  //clean up; THIS IS OBLIGATORY IF WE RUN THE MACRO SERIALLY! Deleting everything below will remove the screen output: https://www.youtube.com/watch?v=BvYJB8NOyBo&t=10s
  if(CleanafterWritingPDF){
    
    for(uint i=0;i<pads.size();i++) delete pads[i];
    for(uint i=0;i<h.size();i++) delete h[i];
    
    delete legend;   
    delete legendRatio;
    delete justcatindexlegend;    
    delete stack;
    delete bkgdErr;
    delete ratioH;
    delete ratioErrH;
    delete ratioS_ggH_SM;
    delete ratioS_VBF_SM;
    delete ratioS_ggH_BSM;
    delete ratioS_VBF_BSM;
    
    delete canv1;
    
    delete histData;
    delete QCD;     
    delete Ztt;      
    delete Zll;      
    delete W;        
    delete TT;       
    delete VV;      
    delete EWK;       
    delete SMH;     
    delete SMHVBF;     
    delete SMHOdd;     
    delete SMHVBFOdd;    
    delete SMHMixed;     
    delete SMHVBFMixed;
    f->Close();
  }
  
  return;
}
