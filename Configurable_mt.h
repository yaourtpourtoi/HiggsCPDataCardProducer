/*
Helper class to produce data cards and make plots of the data cards.
Written by Merijn van de klundert, mvandekl@cern.ch
*/


#ifndef Configurable_mt_H
#define Configurable_mt_H


#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include "TDirectory.h"
#include "TMath.h"
#include "TROOT.h"


using namespace std;

class Configurable_mt : public TObject { // this is needed in case we'd like to add the class to root

   public :
  // Default CTOR, needed by ROOT
  Configurable_mt(){};
  
  //constructor. We only initialise objects for a specific sample-cat combination, anything else is specified via the static members
  Configurable_mt(int category_, int sample_){
    cat=category_;
    sample=sample_;
  }
  
  ~Configurable_mt(){    
    delete xbins;
    delete ybins;

    
    cout<<"OneDimHisto->GetName() "<<OneDimHisto->GetName()<<endl;
    cout<<"OneDimHistoSS->GetName() "<<OneDimHistoSS->GetName()<<endl;
    cout<<"TwoDimHisto->GetName() "<<TwoDimHisto->GetName()<<endl;
    cout<<"TwoDimHistoSS->GetName() "<<TwoDimHistoSS->GetName()<<endl;
    
    delete OneDimHisto;
    delete OneDimHistoSS;    
    delete TwoDimHisto;
    delete TwoDimHistoSS;
    //alternatively, use something like delete gROOT->FindObject("mt_2017_NoCat_mutau_IPdata_obsSS");    
  }
  
  //STATIC DATA MEMBERS
  //SampleNames: always Data first (call data_obs), qcd last. Signals need to have 125 in the sample name in order not to be subtracted in qcd calc. These signal names are format as expected by CH
  inline static const vector<TString> SampleNames{"data_obs","Ztt","Zll","TT","VV","W","ggH_sm_htt125", "ggH_ps_htt125", "ggH_mm_htt125","qqH_sm_htt125", "qqH_ps_htt125", "qqH_mm_htt125","QCD"};
  
  //Important: BaseCatNames ALWAYS expects NoCat at the end! for the last entry no category cut will be applied.
  inline static const vector<TString> BaseCatNames{"ggH","qqH","Ztt","Zll","W","tt","QCD","misc","NoCat"};
  inline static vector<TString> CatNames; //CatNames will just be BaseCatNames, prepended with mt_ , et_. This vector is cleared and push_back is used to fill
  inline static TString InputFileDir; 
  inline static TString OutputFileDir;

  inline static TString OutputDir; //tricky to make static if we run in parallel
  inline static TString OutputFileDataCard;
  inline static TString OutputDirPreFit; //tricky to make static if we run in parallel


 //these things we WON'T be able probably to run multi-threaded over, these configs must be done via separate processes
 inline static int era;
 inline static int decaymode; //0 is pion, 1 is rho. -1: no cut on decay mode! note: decay mode+method will specify the pt cuts in the selection cuts
 inline static int CPMethod; //0: IP method. 1: DP method. 
 inline static int Observable;
 inline static TString ObservableDraw; //declared in  InitObservableNames. This is used in main program from now onwards in the tree:: Draw!
 inline static TString ObservablePath; //declared in InitObservableNames. This is used to put a label to the plot and generate a path not containing semicolons
 inline static TString ObservableXaxis; //declared in InitObservableNames. used in plottng macro to assign x-axis label
 inline static TString ObservableYaxis; //declared in InitObservableNames. used in plottng macro to assign y-axis label
 
 //Default initialised. These may need to be made era-dependent
 inline static  double TTnorm = 1.0;   //scale factors for normalization of ttbar and wjets
 inline static  double Wnorm  = 1.0;  
 inline static  TString topweight="1*";   
 inline static  TString qcdweight="1*";
 inline static  TString zptmassweight="1*";                
 
 //configurable data members of the object
 uint cat; //cat is also the handle to run code in parallel
 int sample;
 
 //data members assign via initialise() function, after the static configurables have been set
 TString CatName, SampleName;
 TString InputFileName;
 int muonptcut; //currently always 20, but NOTE for decay plane method it may be lower!
 int TauChConstCut; //Note: it 
 TString selcut;
 TString FinalCut;// will contain the final cut applied to ->Draw() Method.
 TString FinalCutSS;// will contain the final cut applied to ->Draw() Method.


 //for the histogram
 TH1F* OneDimHisto=NULL;
 TH1F* OneDimHistoSS=NULL;
 TH2F* TwoDimHisto=NULL;
 TH2F* TwoDimHistoSS=NULL;

 // these are (must!) initialised in InitHistogram
 TString HistoName, HistoTitle;
 TString HistoName2D, HistoTitle2D;
 TString HistoNameSS, HistoTitleSS;
 TString HistoNameSS2D, HistoTitleSS2D;

 //fixed in InitialiseXandYAxis
 float xmin, xmax;
 int nbinsx;
 double * xbins=NULL; //initialise to NULL, avoids crash when delete and not assigned. 
 float ymin, ymax;
 int nbinsy;
 double *  ybins=NULL;

 TString EventsWeight;
 //later may extend with asking for tree name or vector with systematics..

 inline static void InitOutputBaseDir(){
   if(era>2018||era<2016){throw std::invalid_argument("Era given not valid, please choose 2016, 2017 or 2018");}
   if(decaymode>11||era<-1){throw std::invalid_argument("decaymode given not valid, please enter -1, 0,1, or 11");}
   if(CPMethod>1||era<0){throw std::invalid_argument("CPMethod given not valid, please choose 0 or 1");}

   OutputDir="";
   OutputDir+=era;
   OutputDir+="/";
   OutputDir+="DM_";
   OutputDir+=decaymode;
   //OutputDir+="/"; OutputDir+="Method_"; OutputDir+=CPMethod; we don't set the CP method. It it STICTLY only relevant for observable 21, and then it enters the path via the obsevable name.
   OutputDir+="/";
   //we don't add ObservablePath, since for plots, typically you want all plots for a certain configuration in one directory
 }
 
 inline static void InitOutputFileDataCard(){
   //add observable for path where data card is looked up
   OutputFileDataCard=OutputDir;
   OutputFileDataCard+=ObservablePath;
   OutputFileDataCard+="/";
   OutputFileDataCard+="Datacard.root";
 }

 inline static void InitOutputPrefitPlotDir(){
   //for plots, likely it is easiest to have them all in one directory, since many (~20). thus don't add ObservablePath, we add this to the pdf name instead
   OutputDirPreFit="Plots_Prefit/";
   OutputDirPreFit+=OutputDir;
 }
 
 void Initialise(){
   //here define define cuts and weights. If different for 2016/2018, redeclar with simple if statement
   SampleName=SampleNames[sample];
   InputFileName+=InputFileDir;
   InputFileName+=SampleNametoInputFile(SampleName);
   CatName=CatNames[cat];
   
   //define the event weight
   EventsWeight = "xsec_lumi_weight*puweight*effweight*mcweight*zptweight*topptweight*";
   TString TSLeaf=""; //base: ggH_sm_htt125
   TString SampleNameDummy=SampleName;//the step may not be needed..
   if(SampleName.Contains("ggH")){TSLeaf=SampleNameDummy.ReplaceAll("ggH_","");TSLeaf.Prepend("gen_");TSLeaf.Append("*");} 
   if(SampleName.Contains("qqH")){TSLeaf=SampleNameDummy.ReplaceAll("qqH_","");TSLeaf.Prepend("gen_");TSLeaf.Append("*");}
   EventsWeight+=TSLeaf;
     
   //define pt and category cuts
   muonptcut=20; 
  if(decaymode==-1||decaymode==11) TauChConstCut=0;  //
  if(decaymode==0&&CPMethod==0) TauChConstCut=0;  //
  if((decaymode==1&&CPMethod==0)||(decaymode==2&&CPMethod==0)) TauChConstCut=40; //ip method for rho
  if((decaymode==1&&CPMethod==1)||(decaymode==2&&CPMethod==1)) TauChConstCut=0;  //dp method for rho
  
  selcut="";
  if(cat<(Configurable_mt::CatNames.size()-1)){//only add cat in case it is not "no cat"
    selcut="&&predicted_class==";
    selcut+=cat;}
  
  if(decaymode!=-1){
    selcut+="&&tau_decay_mode_2==";
    selcut+=decaymode;}
  
  selcut+="&&pt_1>";
  selcut+=muonptcut;
  selcut+="&&pt_2>";
  selcut+=TauChConstCut;
   
  selcut+="&&byMediumDeepTau2017v2p1VSjet_2 >0.5 &&mt_1<50"; 

  //here may want to load sample-specific cuts import from a cut directory

  if(sample==0){
    FinalCut="(os>0.5"+selcut+")";
    FinalCutSS=qcdweight+"(os<0.5"+selcut+")";}
  
  else if(sample==1){//Ztt
    TString isZTT="(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))*";
    FinalCut=isZTT+EventsWeight+zptmassweight+"(os>0.5"+selcut+")";
    FinalCutSS=isZTT+EventsWeight+zptmassweight+qcdweight+"(os<0.5"+selcut+")";}
  
  else if(sample==2){//Zll
    TString isZLL="!(((gen_match_1==3||gen_match_1==4)&&gen_match_2==5)||((gen_match_2==3||gen_match_2==4)&&gen_match_1==5))*";
    FinalCut=isZLL+EventsWeight+zptmassweight+"(os>0.5"+selcut+")";
    FinalCutSS=isZLL+EventsWeight+zptmassweight+qcdweight+"(os<0.5"+selcut+")";
  }    
  else if(sample==3){
    FinalCut=EventsWeight+topweight+"(os>0.5"+selcut+")";
    FinalCutSS=EventsWeight+topweight+qcdweight+"(os<0.5"+selcut+")";}
  else{
    FinalCut=EventsWeight+"(os>0.5"+selcut+")";
    FinalCutSS=EventsWeight+qcdweight+"(os<0.5"+selcut+")";}
    
  InitialiseXandYAxis(cat);
  InitHistogram();
}
  
//HELPER FUNCTIONS
 TString SampleNametoInputFile(TString sample){
   TString out="mt-NOMINAL_ntuple_";
   TString SampleAppendix="";
   if(sample.Contains("data_obs")){
     SampleAppendix="Data";}
   else if(sample.Contains("ggH")){
     SampleAppendix="ggH125";}
   else if(sample.Contains("qqH")){
   SampleAppendix="qqH125";}
   else if(sample.Contains("Ztt")){
     SampleAppendix="DY";}
   else if(sample.Contains("Zll")){
     SampleAppendix="DY";}
   
   else{
     SampleAppendix+=sample;}
   out+=SampleAppendix;
   out+=".root";
   return out;
 }

 void InitHistogram(){    
   HistoName=CatName+SampleName;//we put this way to remain able to call with unique name. When unfoldig we rename to sample name
   HistoNameSS=CatName+SampleName+"SS";

   HistoName2D=CatName+SampleName+"_2D";
   HistoNameSS2D=CatName+SampleName+"SS"+"_2D";

   OneDimHisto=new TH1F(HistoName,HistoName,nbinsx,xbins);
   OneDimHisto->Sumw2();
   OneDimHisto->SetXTitle(ObservableXaxis);
   OneDimHisto->SetYTitle(ObservableYaxis);
   
   OneDimHistoSS=new TH1F(HistoNameSS,HistoNameSS,nbinsx,xbins);
   OneDimHistoSS->Sumw2();
   OneDimHistoSS->SetXTitle(ObservableXaxis);
   OneDimHistoSS->SetYTitle(ObservableYaxis);
   
   TwoDimHisto=new TH2F(HistoName2D,HistoName2D,nbinsx,xbins,nbinsy,ybins);
   TwoDimHisto->Sumw2();
   TwoDimHisto->SetXTitle(ObservableXaxis);
   TwoDimHisto->SetYTitle(ObservableYaxis);
   
   TwoDimHistoSS=new TH2F(HistoNameSS2D,HistoNameSS2D,nbinsx,xbins,nbinsy,ybins);
   TwoDimHistoSS->Sumw2();
   TwoDimHistoSS->SetXTitle(ObservableXaxis);
   TwoDimHistoSS->SetYTitle(ObservableYaxis); }

 void Write1DHistoToDir(TDirectory *CatDir){
   CatDir->cd();    
    TH1F* tmphisto=new TH1F();
    OneDimHisto->Copy(*tmphisto);
    tmphisto->SetName(SampleName);
    tmphisto->Write();
    delete tmphisto;}
    
  void Write2DHistoToDir(TDirectory *CatDir){
    CatDir->cd();    
    TH2F* tmphisto=new TH2F();
    TwoDimHisto->Copy(*tmphisto);
    tmphisto->SetName(SampleName+"_2D");
    tmphisto->Write();
    delete tmphisto;}

inline static TH1F * Unfold(TH2F * histInput) {
  int nBinsX = histInput->GetNbinsX();
  int nBinsY = histInput->GetNbinsY();
  int nBins = nBinsX * nBinsY;

  TString NameNew = TString(histInput->GetName())+TString("_unfolded");
  TH1F * histOutput = new TH1F(NameNew,"",nBins,0,float(nBins));
  int iBin = 1;
  for (int j=1; j<=nBinsY; ++j) {
    for (int i=1; i<=nBinsX; ++i) {
      histOutput->SetBinContent(iBin,histInput->GetBinContent(i,j));
      histOutput->SetBinError(iBin,histInput->GetBinError(i,j));
      iBin++;}}
  return histOutput;}
 
 inline static TH1D * ProjectonXaxis(TH2F * histInput){
   TH1D * histOutput=histInput->ProjectionX();
   return histOutput;}
 
 inline static void InitCatNames(){
   CatNames.clear();
   TString BaseCat="";  
   for(uint i=0; i<BaseCatNames.size();i++){
     BaseCat="mt_";
     BaseCat+=era;
     BaseCat+="_";
     BaseCat+=BaseCatNames[i];
     if(decaymode==-1) BaseCat+="_mutau";
     if(decaymode==0) BaseCat+="_mupi";
     if(decaymode==1) BaseCat+="_murho";
     if(decaymode==11) BaseCat+="_muother";
     if(CPMethod==0) BaseCat+="_IP";
     if(CPMethod==1) BaseCat+="_DP";
     CatNames.push_back(BaseCat);
   }}

 inline static void InitObservableNames(){
   if(Observable>21||Observable<0){throw std::invalid_argument("For index given, no observable is in the lookup table, exiting");}
   ObservableYaxis="Occ.";//correct for nearly all
   ObservableDraw="predicted_prob:"; //we always draw two dimensionally w.r.t. the nn score, below we add the specific observable for x-axis. Adjust ObservableDraw below if you want to unroll in a different observable   
   if(Observable==0){ ObservablePath="pt_1";   ObservableXaxis="pt #mu [GeV]";}
   if(Observable==1){ ObservablePath="pt_2";   ObservableXaxis="pt #tau [GeV]";}
   if(Observable==2){ ObservablePath="jpt_1";   ObservableXaxis="Leading jet pT [GeV]"; }
   if(Observable==3){ ObservablePath="jpt_2";   ObservableXaxis="Trailing jet pT [GeV]"; }
   if(Observable==4){ ObservablePath="bpt_1";   ObservableXaxis="Leading b-jet pT [GeV]"; }
   if(Observable==5){ ObservablePath="bpt_2";   ObservableXaxis="Trailing b-jet pT [GeV]"; }   
   if(Observable==6){ ObservablePath="njets";   ObservableXaxis="Jet multiplicity"; }
   if(Observable==7){ ObservablePath="nbtag";   ObservableXaxis="B-jet multiplicity"; }
   if(Observable==8){ ObservablePath="mt_1";   ObservableXaxis="#mu transverse mass [GeV]"; }
   if(Observable==9){ ObservablePath="pt_tt";   ObservableXaxis="Di-#tau pT [GeV]"; }
   if(Observable==10){ ObservablePath="mjj";   ObservableXaxis="Jet pair invar. mass [GeV]"; }
   if(Observable==11){ ObservablePath="jdeta";   ObservableXaxis="Jet pair #Delta#eta"; }
   if(Observable==12){ ObservablePath="dijetpt";   ObservableXaxis="Di-jet pT [GeV]"; }
   if(Observable==13){ ObservablePath="m_vis";   ObservableXaxis="Visible mass [GeV]"; }
   if(Observable==14){ ObservablePath="eta_1";   ObservableXaxis="#mu #eta"; }
   if(Observable==15){ ObservablePath="eta_2";   ObservableXaxis="#tau #eta"; }
   if(Observable==16){ ObservablePath="met_rcmr";   ObservableXaxis="met [GeV]"; }
   if(Observable==17){ ObservablePath="pt_sv";   ObservableXaxis="pT SVfit [GeV]"; }
   if(Observable==18){ ObservablePath="";   ObservableXaxis=""; }//empty slots!
   if(Observable==19){ ObservablePath="";   ObservableXaxis=""; }
   //observables > 20 reserved for NN score and CP angles
   if(Observable==20){ ObservablePath="predicted_prob";   ObservableXaxis="NN output"; ObservableYaxis="dN/d(NN output)";}
   if(Observable==21){ ObservablePath="acotautau_0"; ObservablePath+=CPMethod;  ObservableXaxis="#phi CP"; }
   
   ObservableDraw+=ObservablePath;
 }

 //axes initialiser. Can only be called with an object, since need a category (i.e. not static!)
 void InitialiseXandYAxis(int cat){

   //we only declare here the number of bins and ranges for obs. 0-19. Below, axes are initialised quidistantly in case we don't do something special for obs. 20 and 21.
   if(Observable==0){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==1){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==2){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==3){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==4){ nbinsx  = 16; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==5){ nbinsx  = 16; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==6){ nbinsx  = 5; xmin=0;xmax=5; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==7){ nbinsx  = 5; xmin=0;xmax=5; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==8){ nbinsx  = 22; xmin=0;xmax=60; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==9){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==10){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==11){ nbinsx  = 8; xmin=0;xmax=5.6; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==12){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==13){ nbinsx  = 30; xmin=0;xmax=160;ymin=0; ymax=1; nbinsy=1;}
   if(Observable==14){ nbinsx  = 50; xmin=-2.5;xmax=2.5; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==15){ nbinsx  = 50; xmin=-2.5;xmax=2.5; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==16){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==17){ nbinsx  = 32; xmin=0;xmax=160; ymin=0; ymax=1; nbinsy=1;}
   if(Observable==18){ ymin=0; ymax=1; nbinsy=1;}
   if(Observable==19){ ymin=0; ymax=1; nbinsy=1;}
    
   //below declare binnings for the NN axes for different categories.       
   if(Observable==20){
     ymin=0; ymax=1; nbinsy=1; //for plotting nn score, we just need one bin on y-axis
     AxesLUT(cat, nbinsx, xbins, decaymode, CPMethod);
     
     /* //uncomment these lines if instead want NN score in flat binning in 10 bins from 0 to 1
	nbinsx=10; xmax=1; xmin=0;
	xbins=NULL;*/
     
   }//sets the axes xbins for the category, puts nbinsx to size-1        
   
   if(Observable==21){//note the difference: we flip here the axis binning w.r.t. observable 20
     nbinsx  = 5; xmin=0;xmax=2*TMath::Pi();
     AxesLUT(cat,nbinsy, ybins, decaymode, CPMethod);}
   
   //here we make equidistant binnings in case an axis has not been defined yet
   if(xbins==NULL){
     xbins=new double[nbinsx+1];
     double step=(xmax-xmin)/nbinsx;
     for(int bc=0;bc<nbinsx+1;bc++){
       xbins[bc]=xmin+bc*step;}}
   
   if(ybins==NULL){
     ybins=new double[nbinsy+1];
     double step=(ymax-ymin)/nbinsy;     
     for(int bc=0;bc<nbinsy+1;bc++){
       ybins[bc]=ymin+bc*step;}}
 }

 //this is NOT static; it depends on CAT, thus only call it for an object
 void AxesLUT(int CAT, int & NBINS, double*& axis, int DM, int CM){
   //Here declare everything as function of CAT and whatever dependence you like to add; pt Higgs, ...!!
   vector<double> tmp;

   //Provisional binning for mu+pi case, which we adopt for dm=-1 also
   if(DM==0||DM==-1){
     if(CAT==0) tmp={0, 0.25, 0.3, 1};
     if(CAT==1) tmp={0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1}; //last bin must be 0.8 to 1, otherwise ratio goes bad.. cutting bins further in two gives issues
     if(CAT==2) tmp={0, 0.35, 0.45, 0.55, 0.65, 1}; 
     if(CAT==3) tmp={0, 0.25, 0.3, 0.4, 0.5, 0.6, 1}; 
     if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45, 1}; 
     if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     if(CAT==8) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something
   }

   //mu+rho in dp method
   if(DM==1&&CM==1){
     if(CAT==0) tmp={0, 0.25, 0.3,0.35, 1}; 
     if(CAT==1) tmp={0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; 
     if(CAT==2) tmp={0, 0.35, 0.45, 0.55, 0.65, 1}; 
     if(CAT==3) tmp={0, 0.25, 0.3, 0.4, 0.45, 0.55, 1}; 
     if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45, 1}; 
     if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     if(CAT==8) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
   }

   if(DM==11){//other category. We only need signal cats, rest is to avoid crash..
     if(CAT==0) tmp={0,0.3,1}; 
     if(CAT==1) tmp={0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1}; 
     if(CAT==2) tmp={0, 0.35, 0.45, 0.55, 0.65, 1}; 
     if(CAT==3) tmp={0, 0.25, 0.3, 0.4, 0.45, 0.55, 1}; 
     if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45,0.5,0.6,1}; 
     if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     if(CAT==8) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
   }
   
   NBINS=tmp.size()-1;
   axis=new double[tmp.size()];
   for(uint index=0; index<tmp.size();index++){
     axis[index]=tmp[index];}
 }
 
 // this is needed in case we'd like to add the class to root
  ClassDef (Configurable_mt,1)        
    };

#endif
