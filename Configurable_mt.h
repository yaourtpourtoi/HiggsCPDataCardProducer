/*
Helper class to produce data cards and make plots of the data cards.
Written by Merijn van de klundert, mvandekl@cern.ch
*/


#ifndef Configurable_mt_H
#define Configurable_mt_H


#include "TH2F.h"
#include "TSystem.h"
#include <iostream>
#include <map>
#include <iomanip>

#include "TDirectory.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"

template<typename KeyTemp, typename ValTemp>
std::vector<KeyTemp> extract_keys(std::map<KeyTemp, ValTemp> const& input_map) {
  std::vector<KeyTemp> keys;
  for (auto const& element : input_map) {
    keys.push_back(element.first);
  }
  return keys;
}


using namespace std;

class Configurable_mt : public TObject { // this is needed in case we'd liek to add the class to root

   public :
  // Default CTOR, needed by ROOT
  Configurable_mt(){};
  
  //constructor. We only initialise objects for a specific sample-cat combination, anything else is specified via the static members
  Configurable_mt(int category_, TString sample_){
    cat=category_;
    sample=sample_;
  }
  
  ~Configurable_mt(){    
    delete xbins;
    delete ybins;

    
    cout<<"OneDimHisto->GetName() "<<OneDimHisto->GetName()<<endl;
    cout<<"OneDimHistoSS->GetName() "<<OneDimHistoSS->GetName()<<endl;
    cout<<"OneDimHistoAR->GetName() "<<OneDimHistoAR->GetName()<<endl;
    cout<<"TwoDimHisto->GetName() "<<TwoDimHisto->GetName()<<endl;
    cout<<"TwoDimHistoSS->GetName() "<<TwoDimHistoSS->GetName()<<endl;
    cout<<"TwoDimHistoAR->GetName() "<<TwoDimHistoAR->GetName()<<endl;
    
    delete OneDimHisto;
    delete OneDimHistoSS;    
    delete OneDimHistoAR;    
    delete TwoDimHisto;
    delete TwoDimHistoSS;
    delete TwoDimHistoAR;
    //alternatively, use something like delete gROOT->FindObject("mt_2017_NoCat_mutau_IPdata_obsSS");    
  }
  
  //STATIC DATA MEMBERS
  //SampleNames: always Data first (call data_obs), qcd last. Signals need to have 125 in the sample name in order not to be subtracted in qcd calc. These signal names are format as expected by CH
  //inline static const vector<TString> SampleNames{"data_obs","Ztt","Zll","TT","VV","W","ggH_sm_htt125", "ggH_ps_htt125", "ggH_mm_htt125","qqH_sm_htt125", "qqH_ps_htt125", "qqH_mm_htt125","QCD"};
  inline static const map<TString,TString> SampleNames{
    {"data", "data_obs"},
    {"Ztt", "Ztt"},
    {"Zll", "Zll"},
    {"TT", "TT"},
    {"ST", "ST"},
    {"VV", "VV"},
    {"W", "W"},
    {"ggHSM", "ggH_sm_htt125"},
    {"ggHPS", "ggH_ps_htt125"},
    {"ggHMM", "ggH_mm_htt125"},
    {"qqHSM", "qqH_sm_htt125"},
    {"qqHPS", "qqH_ps_htt125"},
    {"qqHMM", "qqH_mm_htt125"},
    {"QCD", "QCD"},
    {"fakes", "fakes"},
    {"embedded", "EMB"}
  };
  inline static const map<TString,TString> FileNames{
    {"data", "mt-NOMINAL_ntuple_Data"},
    {"Ztt", "mt-NOMINAL_ntuple_DY"},
    {"Zll", "mt-NOMINAL_ntuple_DY"},
    {"TT", "mt-NOMINAL_ntuple_TT"},
    {"ST", "mt-NOMINAL_ntuple_ST"},
    {"VV", "mt-NOMINAL_ntuple_VV"},
    {"W", "mt-NOMINAL_ntuple_W"},
    {"ggHSM", "mt-NOMINAL_ntuple_ggH125"},
    {"ggHPS", "mt-NOMINAL_ntuple_ggH125"},
    {"ggHMM", "mt-NOMINAL_ntuple_ggH125"},
    {"qqHSM", "mt-NOMINAL_ntuple_qqH125"},
    {"qqHPS", "mt-NOMINAL_ntuple_qqH125"},
    {"qqHMM", "mt-NOMINAL_ntuple_qqH125"},
    {"embedded", "mt-NOMINAL_ntuple_EMB"}
  };
  
  //Important: BaseCatNames ALWAYS expects NoCat at the end! for the last entry no category cut will be applied.
  //  inline static const map<int,TString> CategoryNames{{0,"ggH"},{1,"qqH"},{2,"Ztt"},{3,"Zll"},{4,"W"},{5,"tt"},{6,"QCD"},{7,"misc"},{-1,"NoCat"}};


  inline static const map<int,TString> CategoryNames{{0,"higgs"},{1,"ztt"},{2,"top"},{3,"fakes"},{-1,"NoCat"}};


  inline static map<int,TString> CatNames; //CatNames will just be CategoryNames, prepended with mt_ , et_. This vector is cleared and push_back is used to fill
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
 inline static TString ObservableName;
 inline static TString ObservableDraw; //declared in  InitObservableNames. This is used in main program from now onwards in the tree:: Draw!
 inline static TString ObservablePath; //declared in InitObservableNames. This is used to put a label to the plot and generate a path not containing semicolons
 inline static TString ObservableXaxis; //declared in InitObservableNames. used in plottng macro to assign x-axis label
 inline static TString ObservableYaxis; //declared in InitObservableNames. used in plottng macro to assign y-axis label

 inline static TString CommonCuts;
 inline static int nSignalCategories;
 inline static int nCategories;
 inline static bool useEmbedded;
 inline static bool FFmethod;


 //Default initialised. These may need to be made era-dependent
 inline static  double TTnorm = 1.0;   //scale factors for normalization of ttbar and wjets
 inline static  double Wnorm  = 0.97;  
 inline static  TString qcdweight="*1.06";               
 
 //configurable data members of the object
 int cat; //cat is also the handle to run code in parallel
 TString sample;
 
 //data members assign via initialise() function, after the static configurables have been set
 TString CatName, SampleName;
 TString InputFileName;
 int muonptcut; //currently always 20, but NOTE for decay plane method it may be lower!
 int TauChConstCut; //Note: it 
 TString selcut;
 TString SRCut;// will contain the final cut applied to ->Draw() Method.
 TString SSCut;// will contain the final cut applied to ->Draw() Method.
 TString ARCut;

 //for the histogram
 TH1F* OneDimHisto=NULL;
 TH1F* OneDimHistoSS=NULL;
 TH1F* OneDimHistoAR=NULL;
 TH2F* TwoDimHisto=NULL;
 TH2F* TwoDimHistoSS=NULL;
 TH2F* TwoDimHistoAR=NULL;

 // these are (must!) initialised in InitHistogram
 TString HistoName, HistoTitle;
 TString HistoName2D, HistoTitle2D;
 TString HistoNameSS, HistoTitleSS;
 TString HistoNameSS2D, HistoTitleSS2D;
 TString HistoNameAR, HistoTitleAR;
 TString HistoNameAR2D, HistoTitleAR2D;

 //fixed in InitialiseXandYAxis
 float xmin, xmax;
 int nbinsx;
 double * xbins=NULL; //initialise to NULL, avoids crash when delete and not assigned. 
 float ymin, ymax;
 int nbinsy;
 double *  ybins=NULL;

 TString EventWeight;
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

   SampleName=SampleNames.at(sample);
   InputFileName+=InputFileDir;
   if(sample!="QCD"&&sample!="fakes")InputFileName+=SampleNametoInputFile(sample);
   CatName=CatNames.at(cat);

   //define the event weight
   EventWeight = "*xsec_lumi_weight*weight";
   TString TSLeaf=""; //base: ggH_sm_htt125
   //TString SampleNameDummy=SampleName;//the step may not be needed..
   if(sample.Contains("SM"))EventWeight+="*gen_sm_htt125";
   else if(sample.Contains("PS"))EventWeight+="*gen_ps_htt125";
   else if(sample.Contains("MM"))EventWeight+="*gen_mm_htt125";
     
   //define pt and category cuts
  
   //selcut="(byMediumDeepTau2017v2p1VSjet_2 >0.5 && mt_1<50)"; 
  selcut=CommonCuts; 
  if(cat!=-1){//only add cat in case it is not "no cat"
    selcut+="*(predicted_class=="+TString::Itoa(cat,10)+")";
  }
  
  if(decaymode!=-1){
    selcut+="*(tau_decay_mode_2=="+TString::Itoa(decaymode,10)+")";
  }
  
  selcut+="*(pt_1>20)";
  selcut+="*(pt_2>30)";
  if(CPMethod==0&&(decaymode==1||decaymode==2)) selcut+="*(chpt_2>40)";//ip method for rho

  //+=TauChConstCut;
   
  SRCut=selcut+"*(os>0.5 && byMediumDeepTau2017v2p1VSjet_2>0.5 )";
  SSCut=selcut+qcdweight+"*(os<0.5 && byMediumDeepTau2017v2p1VSjet_2>0.5 )";
  ARCut=selcut+"*ff_nom*(os>0.5 && byMediumDeepTau2017v2p1VSjet_2<0.5 )";

  //here may want to load sample-specific cuts import from a cut directory

  if(sample!="data"&&sample!="embedded"){//weights applied only on MC
    SRCut+=EventWeight;
    ARCut+=EventWeight;
    SSCut+=EventWeight;
    if(useEmbedded&&!(sample.Contains("ggH")||sample.Contains("qqH"))){
      SRCut+="*(!(gen_match_2==5 &&gen_match_1==4))";
      SSCut+="*(!(gen_match_2==5 &&gen_match_1==4))";
      ARCut+="*(!(gen_match_2==5 &&gen_match_1==4))";
    }
    if(FFmethod&&!(sample.Contains("ggH")||sample.Contains("qqH"))){
      SRCut+="*(gen_match_2!=6)";
      SSCut+="*(gen_match_2!=6)";
      ARCut+="*(gen_match_2!=6)";
    }

  }  
  else if(sample=="embedded"){
    SRCut+="*embweight*effweight*mcweight";
    ARCut+="*embweight*effweight*mcweight";
    SSCut+="*embweight*effweight*mcweight";
  }

  if(sample=="Ztt"){//Ztt
    TString isZTT="*(gen_match_2==5)";
    SRCut+=isZTT+"*zptweight";
    ARCut+=isZTT+"*zptweight";
    SSCut+=isZTT+"*zptweight";}
  
  else if(sample=="Zll"){//Zll
    TString isZLL="*!(gen_match_2==5)";
    SRCut+=isZLL+"*zptweight";
    SSCut+=isZLL+"*zptweight";
    ARCut+=isZLL+"*zptweight";
  }    
  else if(sample=="TT"){
    SRCut+="*topptweight";
    ARCut+="*topptweight";
    SSCut+="*topptweight";}
    
  InitialiseXandYAxis(cat);
  InitHistogram();
}
  
//HELPER FUNCTIONS
 TString SampleNametoInputFile(TString sample){
   TString filename=FileNames.at(sample);
   /*
   TString filename="mt-NOMINAL_ntuple_";
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
     filename+=SampleAppendix;*/
   filename+=".root";
   return filename;
 }

 void InitHistogram(){    
   HistoName=CatName+SampleName;//we put this way to remain able to call with unique name. When unfoldig we rename to sample name
   HistoNameSS=CatName+SampleName+"SS";
   HistoNameAR=CatName+SampleName+"AR";

   HistoName2D=CatName+SampleName+"_2D";
   HistoNameSS2D=CatName+SampleName+"SS"+"_2D";
   HistoNameAR2D=CatName+SampleName+"AR"+"_2D";

   OneDimHisto=new TH1F(HistoName,HistoName,nbinsx,xbins);
   OneDimHisto->Sumw2();
   OneDimHisto->SetXTitle(ObservableXaxis);
   OneDimHisto->SetYTitle(ObservableYaxis);
   
   OneDimHistoSS=new TH1F(HistoNameSS,HistoNameSS,nbinsx,xbins);
   OneDimHistoSS->Sumw2();
   OneDimHistoSS->SetXTitle(ObservableXaxis);
   OneDimHistoSS->SetYTitle(ObservableYaxis);
   
   OneDimHistoAR=new TH1F(HistoNameAR,HistoNameAR,nbinsx,xbins);
   OneDimHistoAR->Sumw2();
   OneDimHistoAR->SetXTitle(ObservableXaxis);
   OneDimHistoAR->SetYTitle(ObservableYaxis);
   
   TwoDimHisto=new TH2F(HistoName2D,HistoName2D,nbinsx,xbins,nbinsy,ybins);
   TwoDimHisto->Sumw2();
   TwoDimHisto->SetXTitle(ObservableXaxis);
   TwoDimHisto->SetYTitle(ObservableYaxis);
   
   TwoDimHistoSS=new TH2F(HistoNameSS2D,HistoNameSS2D,nbinsx,xbins,nbinsy,ybins);
   TwoDimHistoSS->Sumw2();
   TwoDimHistoSS->SetXTitle(ObservableXaxis);
   TwoDimHistoSS->SetYTitle(ObservableYaxis);

   TwoDimHistoAR=new TH2F(HistoNameAR2D,HistoNameAR2D,nbinsx,xbins,nbinsy,ybins);
   TwoDimHistoAR->Sumw2();
   TwoDimHistoAR->SetXTitle(ObservableXaxis);
   TwoDimHistoAR->SetYTitle(ObservableYaxis); }

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
   vector<int> categoryIndices=extract_keys(CategoryNames);
   for(int categoryIndex : categoryIndices){
     BaseCat="mt_";
     BaseCat+=era;
     BaseCat+="_";
     BaseCat+=CategoryNames.at(categoryIndex);
     if(decaymode==-1) BaseCat+="_mutau";
     else if(decaymode==0) BaseCat+="_mupi";
     else if(decaymode==1) BaseCat+="_murho";
     else BaseCat+="_muother";
     if(CPMethod==0) BaseCat+="_IP";
     else if(CPMethod==1) BaseCat+="_DP";
     CatNames.insert({categoryIndex,BaseCat});
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
   if(Observable==21){ ObservablePath="acotautau_helix_0"; ObservablePath+=CPMethod;  ObservableXaxis="#phi CP"; }
   
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
     if(CAT==0) tmp={0, 0.5, 0.8, 1};
     else if(CAT==1) tmp={0, 0.55, 0.75, 1}; //last bin must be 0.8 to 1, otherwise ratio goes bad.. cutting bins further in two gives issues
     else if(CAT==2) tmp={0, 0.55, 0.75, 1}; 
     else if(CAT==3) tmp={0, 0.5, 0.7, 1}; 
     else if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45, 1}; 
     else if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     else if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     else if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     else if(CAT==-1) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something
   }

   //mu+rho in dp method
   else if(DM==1&&CM==1){
     if(CAT==0) tmp={0, 0.5, 0.8, 1}; 
     else if(CAT==1) tmp={0, 0.55, 0.75, 1}; 
     else if(CAT==2) tmp={0, 0.55, 0.75, 1}; 
     else if(CAT==3) tmp={0, 0.5, 0.7, 1}; 
     else if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45, 1}; 
     else if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     else if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     else if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     else if(CAT==-1) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
   }

   else if(DM==11||DM==10){//other category. We only need signal cats, rest is to avoid crash..
     if(CAT==0) tmp={0, 0.5, 0.8, 1}; 
     else if(CAT==1) tmp={0, 0.55, 0.75, 1}; 
     else if(CAT==2) tmp={0, 0.55, 0.75, 1}; 
     else if(CAT==3) tmp={0, 0.5, 0.7, 1}; 
     else if(CAT==4) tmp={0, 0.25, 0.35,0.4, 0.45,0.5,0.6,1}; 
     else if(CAT==5) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 1}; 
     else if(CAT==6) tmp={0, 0.3, 0.4, 0.5, 0.6, 1}; 
     else if(CAT==7) tmp={0, 0.35, 0.45, 0.55, 0.65, 0.8, 1}; 
     else if(CAT==-1) tmp={0, 0.3, 0.6, 1}; //DNN score range in case we pick no category, i.e. it is void observable in for the DNN score or CP angles, but need to provide something.
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
