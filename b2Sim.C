// time is in microseconds
// model with absorption
#include "modelFit.hh"
using namespace TMath;

// typical expect 5 pe/kev and 2000 kev source
double nprim = 5.E3; //primary photons created 
double PPM=0;
double absorption = 0.62;
double kprime =1.0;
double kp=1.0;
double area = 0.6*0.6;
enum {NSIPM=12}; 
int theSet; // the dopant set
int ievent=0;
TRandom3* ran3;
TDirectory* dEvent; 
int levelColor[4];


modelFit* model[NTYPES];
TF1* fit[NTYPES];
TF1* fAbs;
TH1D *hWaveEv[NSETS];
TH1D *hWaveSum[NSETS];
TH1D *hWaveInt[NSETS];
TH1D *hAccept;
vector<TVector3*> position;
vector<double> gain;
vector<double> acceptEff;


TFile *fout;

double getNorm() {
  return 1.;
}


static double absFunc(double *xx, double *par) {
  double c = xx[0];
  double a0 = 0.62;
  double c0 = par[0];
  if(c>0.1)  return a0;
  double a = 1.-1.*Exp(-c/c0);
  return a;
};

void makeWaves() {
  TH1D* h = (TH1D*) fit[4]->GetHistogram();
  for (int isipm = 0; isipm  < NSIPM; ++isipm) {
    hWaveSum[isipm] = (TH1D*) h->Clone();
    int level= isipm/3;
    hWaveSum[isipm]->Reset("ICESM");
    hWaveSum[isipm]->SetName(Form("WaveSumSipm-%i-level-%i",isipm,level));
    hWaveSum[isipm]->SetTitle(Form("Wave sum Sipm %i level %i",isipm,level));
    hWaveSum[isipm]->GetXaxis()->SetTitle("time nanoseconds");
    hWaveSum[isipm]->GetYaxis()->SetTitle("PE");
    hWaveSum[isipm]->SetLineColor(levelColor[level]);
    hWaveSum[isipm]->SetMarkerColor(levelColor[level]);
    printf(" made %s \n", hWaveSum[isipm]->GetName());
  }
  hAccept = new TH1D("Accept"," acceptEffance by sipm ",NSIPM,0,NSIPM);
  hAccept->GetXaxis()->SetTitle("sipm number");
  hAccept->GetYaxis()->SetTitle("PE");


}

void makeIntegrals() {
  // integral plots
  int startBin = hWaveSum[0]->FindBin(900.);
  for (int i = 0; i <  NSIPM; ++i)
  {
    double runPPM = fit[4]->GetParameter(2);
    double runKp = fit[4]->GetParameter(4);
    double runKprime = fit[4]->GetParameter(7);
    double yieldSum=0;
    double yieldSumErr=0;
    TString htitle; htitle.Form("IntegralToFitSet%i-PPM%.2f-kP%.2f-kPrime-%.2f", i,runPPM,runKp,runKprime);
    hWaveInt[i] = (TH1D *)hWaveSum[i]->Clone(htitle);
    hWaveInt[i]->Reset("ICESM");
    hWaveInt[i]->SetTitle(Form("integral set %i PPM %i ", i, int(runPPM)));
    printf(" set %i %s bins %i \n",i,hWaveInt[i]->GetName(),hWaveInt[i]->GetNbinsX());
    for (int ibin = 0; ibin < hWaveInt[i]->GetNbinsX(); ++ibin)
    {
      hWaveInt[i]->SetBinContent(ibin,0);
      hWaveInt[i]->SetBinError(ibin, 0);
      if(ibin<startBin) continue;
      double val = hWaveSum[i]->GetBinContent(ibin);
      yieldSum += val;
      yieldSumErr = hWaveSum[i]->GetBinError(ibin);
      hWaveInt[i]->SetBinContent(ibin, yieldSum);
      hWaveInt[i]->SetBinError(ibin, yieldSumErr);
    }
  }
}


void plotSet() {

  TString canTitle;
  TCanvas *can; 
  canTitle.Form("ModelSet-%.3f",PPM);
  can = new TCanvas(canTitle, canTitle);
  can->SetLogy(1);

  for (int i = 0; i < NTYPES; ++i){
    printf(" fit name %s ", fit[i]->GetName());
    fit[i]->GetYaxis()->SetRangeUser(1.E-3,1E3);
    fit[i]->SetLineColor(model[i]->setColor[i]);
    if(i==0) fit[i]->Draw("");
    else fit[i]->Draw("same");
  }
  gPad->SetLogy();
  can->BuildLegend(.3,.2,.3,.2,canTitle);
  can->Modified(); can->Update();
  can->Print(".pdf");
}



void makeModel()
{
  /*
    param 0 binw 4.000E-02 +/- 0.000E+00 
	  param 1 norm 1.000E+04 +/- 0.000E+00 
	  param 2 PPM 1.000E+03 +/- 0.000E+00 
	  param 3 tau3 1.600E+03 +/- 0.000E+00 
	  param 4 kp 1.000E+00 +/- 0.000E+00 
	  param 5 sfrac 1.400E-01 +/- 0.000E+00 
	  param 6 ab 0.000E+00 +/- 0.000E+00 
	  param 7 kprime 0.000E+00 +/- 0.000E+00 
	  param 8 type 4.000E+00 +/- 0.000E+00 
  ** 
  */
  if(theSet>=NSETS) return;
  int fitColor[NTYPES] = {kYellow, kBlue, kRed, kGreen, kBlack};
  TString typeName[NTYPES];
  typeName[0]=TString("singlet");
  typeName[1]=TString("triplet");
  typeName[2]=TString("mixed");
  typeName[3]=TString("xenon");
  typeName[4]=TString("total");
  for (int itype = 0; itype < NTYPES; ++itype){
    model[itype]= new modelFit(itype,theSet);
    fit[itype] = model[itype]->fp;
    fit[itype]->SetName(Form("fit-%s-%.2f",typeName[itype].Data(),PPM));
    PPM =  model[itype]->XPPM[theSet];
    double afactor = absorption;
    if(PPM<1&&absorption>0.001)  afactor = fAbs->Eval( PPM );
    fit[itype]->SetParameter(2, PPM);
    fit[itype]->SetParameter(4,kp);
    //absorption = fAbs->Eval(PPM);
    fit[itype]->SetParameter(6,afactor);
    fit[itype]->SetParameter(7,kprime);
    //model[itype][theSet]->show();
    printf(" fit %i set %i ppm %.2f abs %.4f fit %s \n",itype,theSet,PPM, afactor, fit[itype]->GetName() );
    fout->Append(fit[itype]);
  }
 
  return;
}

void makeAbsFunction() 
{
  fAbs = new TF1("fAbs",absFunc,0,1,1);   // create TF1 class.
  fAbs->SetParName(0,"C0");
  fAbs->SetParameter(0,.103);
  fAbs->SetNpx(1000000); 
  /*
  TCanvas* canabs = new TCanvas("absFitFunc","absFitFunc");
  fabs->Draw("");
  canabs->Print(".pdf");
  */
}


void genEvent() 
{
  ++ievent;
  dEvent->cd();

  for(int isipm=0; isipm<NSIPM; ++ isipm) {
    int level = isipm/3;
    hWaveEv[isipm] = (TH1D*) hWaveSum[0]->Clone();
    hWaveEv[isipm]->Reset("ICESM");
    hWaveEv[isipm]->SetName(Form("WaveEv-%i-Ev%i",isipm,ievent));
    hWaveEv[isipm]->SetTitle(Form("Wave Sipm %i event %i ",isipm,ievent));
    hWaveEv[isipm]->GetXaxis()->SetTitle("time nanoseconds");
    hWaveEv[isipm]->GetYaxis()->SetTitle("PE");
    hWaveEv[isipm]->SetMarkerStyle(20);
    hWaveEv[isipm]->SetLineColor(levelColor[level]);
    hWaveEv[isipm]->SetMarkerColor(levelColor[level]);
    hWaveEv[isipm]->SetMarkerSize(0.2);

    for(int ibin=0; ibin< hWaveSum[isipm]->GetNbinsX() ; ++ ibin) {
      double val = model[4]->fp->Eval(hWaveSum[isipm]->GetBinCenter(ibin));
      double rval = ran3->PoissonD(val*acceptEff[isipm]);
      // accumulate
      hWaveEv[isipm]->SetBinContent(ibin,rval);
      hWaveSum[isipm]->SetBinContent(ibin,rval +  hWaveSum[isipm]->GetBinContent(ibin));
      for(unsigned i=0; i<position.size(); ++i) hAccept->SetBinContent(isipm+1, rval + hAccept->GetBinContent(isipm+1));
    }
  }
  fout->cd();
}

void makeBacon()
{
  // geometry
  double polar = 30.*Pi()/180.; // rad
  double cosPolar = Cos(polar);
  double height[4]={1,10,20,30}; // cm
  int iphi = -1;
  for(int isipm = 0; isipm<NSIPM; ++ isipm) {
    int level = isipm/3;
    ++iphi;
    if(iphi>2) iphi=0;
    double phi = double(iphi)*2*Pi()/3.;
    TVector3 * pos = new TVector3(1.,1.,1.);
    double r = height[level]/cosPolar;
    pos->SetMag(r);
    pos->SetTheta(polar);
    pos->SetPhi(phi);
    //printf(" %i -- r %f phi %i  \n", isipm, r, iphi);
    pos->Print();
    position.push_back(pos);
    gain.push_back(1.);
    acceptEff.push_back(0.15*area/r*4*Pi());
  }
}

void b2Sim(int ngen = 2, int iset=0) 
{
  ran3 = new TRandom3(65539);
  levelColor[0]=kBlack;
  levelColor[1]=kBlue;
  levelColor[2]=kRed;
  levelColor[3]=kGreen;

  theSet=iset;
  TString fileTitle;
  fileTitle.Form("bacon-model-kp-%.2f-kPrime-%.2f-Abs-%.2f.root",kp,kprime,absorption);
  cout << fileTitle << endl;
  fout = new TFile(fileTitle, "RECREATE");
  dEvent = fout->mkdir("events");
  makeAbsFunction(); 
  makeModel();
  makeWaves();
  fout->ls();
  plotSet();
  makeBacon();
  printf(" made  %lu \n",position.size());
  for(unsigned i=0; i<position.size(); ++i) printf(" %i) r %f theta %f phi %f acceptEff %f \n",i,position[i]->Mag(), position[i]->Theta(),position[i]->Phi(), acceptEff[i]);
  
  // gen events 
  for(int ig=0; ig<ngen; ++ig) genEvent();

  TString canTitle;
  canTitle.Form("bacon2-%i-events",ngen);
  TCanvas *canB = new TCanvas(canTitle,canTitle);
  canB->SetLogy();
  canB->SetGridy();
  canB->SetGridx();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int iphi=0;
  for (int isipm = 0; isipm < NSIPM ; ++isipm)
  {
    int level = isipm/3;
    ++iphi;
    if(iphi>2) iphi=0;
    printf(" isipm %i icount %i \n",isipm,iphi);

    if(iphi>0) continue;
    hWaveSum[isipm]->SetStats(1);
    hWaveSum[isipm]->SetMarkerStyle(20);
    //hWaveSum[isipm]->SetLineStyle(isipm);
    hWaveSum[isipm]->SetMarkerSize(0.2);
    if(isipm==0) hWaveSum[isipm]->Draw("");
    else hWaveSum[isipm]->Draw("same");
  }
  canB->BuildLegend();
  canB->Print(".pdf");

  TCanvas* canA = new TCanvas("acceptEff","acceptEff");
  canA->SetTickx(); canA->SetTicky();
  canA->SetGridx(); canA->SetGridy();
  hAccept->Draw();



  fout->Write();
}
