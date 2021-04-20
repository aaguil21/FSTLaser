#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include "Fit/Fitter.h"
#include "Fit/Chi2FCN.h"
#include <Math/Functor.h>
#include "Math/WrappedParamFunction.h"
#include "Math/AdaptiveIntegratorMultiDim.h"

#include "../../src/FstUtil/FstCons.h"
#include "../../src/FstUtil/FstRawHit.h"
#include "../../src/FstUtil/FstCluster.h"
#include "../../src/FstUtil/FstTrack.h"
#include "../../src/FstUtil/FstEvent.h"
#include "./draw.h"

R__LOAD_LIBRARY(../../lib/libFstEvent.dylib)

typedef std::vector<double> dVec;
typedef std::tuple<double, double, double, double, double, double, int> tPars; // phi_rot_ist2, x2_shift, y2_shift, z0_shift, xCut, yCut, nOffset

tPars minuitAlignment(dVec x0_orig, dVec y0_orig, dVec x1_orig, dVec y1_orig, dVec x3_orig, dVec y3_orig, tPars cParameters);
tPars updateFitParameters(tPars fitPars, double xCut, double yCut, int nOffset);

double gaussian(double *var, double *par)
{
  double residual = var[0];

  double norm  = par[0];
  double mean  = par[1];
  double sigma = par[2];
  double binwidth = par[3];

  double sigmaSquare = sigma*sigma;
  double Norm = norm/(sqrt(2.0*TMath::Pi())*sigma);
  double Power = -0.5*(residual-mean)*(residual-mean)/sigmaSquare;

  double y = Norm*exp(Power)*binwidth;

  return y;
}

int getAlignment_FSTClusters_3Layer(float nFstHitsCut = 4.0, int numOfUsedTimeBins = 2, int sensorId = 0, string mod = "Mod04")
{
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1001);
  gStyle->SetStatX(0.95); gStyle->SetStatY(0.90);
  gStyle->SetStatW(0.20); gStyle->SetStatH(0.20);

  dVec r0_fst;
  dVec phi0_fst;
  dVec x0_fst;
  dVec y0_fst;
  dVec x1_ist;
  dVec y1_ist;
  dVec x3_ist;
  dVec y3_ist;

  // clear vector
  r0_fst.clear();
  phi0_fst.clear();
  x0_fst.clear();
  y0_fst.clear();
  x1_ist.clear();
  y1_ist.clear();
  x3_ist.clear();
  y3_ist.clear();

  const double phi_rot_ist1 = -0.00113551; // aligned IST1 & IST3 to IST2
  const double x1_shift = 0.530444;
  const double y1_shift = 1.17133;
  const double phi_rot_ist3 = -0.00379431;
  const double x3_shift = 0.45811;
  const double y3_shift = 1.0116;

  std::string inputlist = Form("../../list/FST/alignment/%s/FstAlignment_%s.list",mod.c_str(),mod.c_str());
  cout << "input list set to: " << inputlist.c_str() << endl;
  TChain *mChainInPut = new TChain("mTree_FstEvent");
  if (!inputlist.empty())   // if input file is ok
  {
    cout << "Open input probability file list" << endl;
    ifstream in(inputlist.c_str());  // input stream
    if(in)
    {
      cout << "input file probability list is ok" << endl;
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  string addfile;
	  addfile = str;
	  mChainInPut->AddFile(addfile.c_str(),-1,"mTree_FstEvent");
	  Long64_t file_entries = mChainInPut->GetEntries();
	  cout << "File added to data chain: " << addfile.c_str() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: input probability file input is problemtic" << endl;
      return false;
    }
  }

  FstEvent *mFstEvent = new FstEvent();
  mChainInPut->SetBranchAddress("FstEvent",&mFstEvent);
  long NumOfEvents = (long)mChainInPut->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;

  // if(NumOfEvents > 1000) NumOfEvents = 1000;
  // NumOfEvents = 500;
  mChainInPut->GetEntry(0);

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mChainInPut->GetEntry(i_event);

    std::vector<FstCluster *> clusterVec_fst;
    std::vector<FstCluster *> clusterVec_ist1;
    std::vector<FstCluster *> clusterVec_ist3;
    clusterVec_fst.clear();
    clusterVec_ist1.clear();
    clusterVec_ist3.clear();

    for(int i_cluster = 0; i_cluster < mFstEvent->getNumClusters(); ++i_cluster)
    { // get Clusters info for IST1 IST2 and IST3
      FstCluster *fstCluster = mFstEvent->getCluster(i_cluster);
      if(fstCluster->getLayer() == 0 && fstCluster->getClusterType() == 2 && fstCluster->getIsSeed() && fstCluster->getSensor() == sensorId)
      { // use Scan cluster
	clusterVec_fst.push_back(fstCluster);
      }
      if(fstCluster->getLayer() == 1 && fstCluster->getClusterType() == 1)
      {
	clusterVec_ist1.push_back(fstCluster);
      }
      if(fstCluster->getLayer() == 3 && fstCluster->getClusterType() == 1)
      {
	clusterVec_ist3.push_back(fstCluster);
      }
    }

    for(int i_ist1 = 0; i_ist1 < clusterVec_ist1.size(); ++i_ist1)
    {
      double x1 = clusterVec_ist1[i_ist1]->getMeanX();
      double y1 = clusterVec_ist1[i_ist1]->getMeanY();
      const double z1_ist = FST::pitchLayer12 + FST::pitchLayer23;
      double row1 = clusterVec_ist1[i_ist1]->getMeanRow();

      double x1_corr = x1*TMath::Cos(phi_rot_ist1) + y1*TMath::Sin(phi_rot_ist1) + x1_shift;
      double y1_corr = y1*TMath::Cos(phi_rot_ist1) - x1*TMath::Sin(phi_rot_ist1) + y1_shift;

      for(int i_ist3 = 0; i_ist3 < clusterVec_ist3.size(); ++i_ist3)
      {
	double x3 = clusterVec_ist3[i_ist3]->getMeanX();
	double y3 = clusterVec_ist3[i_ist3]->getMeanY();
	const double z3_ist = 0.0;
	double row3 = clusterVec_ist3[i_ist3]->getMeanRow();

	double x3_corr = x3*TMath::Cos(phi_rot_ist3) + y3*TMath::Sin(phi_rot_ist3) + x3_shift;
	double y3_corr = y3*TMath::Cos(phi_rot_ist3) - x3*TMath::Sin(phi_rot_ist3) + y3_shift;

	if( abs(row3-row1) < 17)
	{
	  if(clusterVec_fst.size() > 0) // FST has at least one hit
	  {
	    for(int i_fst = 0; i_fst < clusterVec_fst.size(); ++i_fst)
	    {
	      double r0 = clusterVec_fst[i_fst]->getMeanX();
	      double phi0 = clusterVec_fst[i_fst]->getMeanY();
	      double rMinTest = FST::mFstRMin[sensorId];
	      double rMaxTest = FST::mFstRMax[sensorId];
	      if(r0 > rMinTest && r0 < rMaxTest)
	      {
		double x0 = r0*TMath::Cos(phi0);
		double y0 = r0*TMath::Sin(phi0);
		r0_fst.push_back(r0);
		phi0_fst.push_back(phi0);
		x0_fst.push_back(x0);
		y0_fst.push_back(y0);
		x1_ist.push_back(x1_corr); // aligned w.r.t. IST2
		y1_ist.push_back(y1_corr);
		x3_ist.push_back(x3_corr);
		y3_ist.push_back(y3_corr);
	      }
	    }
	  }
	}
      }
    }
  }

  // mFile_InPut->Close();

  cout << "Start Minuit Fit for alignment => " << endl;
  // -------------------------------------
  // for plotting
  int numBinX   = 50;
  double maxX   = 160.0;
  double minX   = -160.0;
  int numBinY   = 200;
  double maxY   = 100.0;
  double minY   = -100.0;
  int numBinR   = 50;
  double maxR   = 160.0;
  double minR   = -160.0;
  int numBinPhi = 1570;
  double maxPhi = 0.5*TMath::Pi();
  double minPhi = -0.5*TMath::Pi();
  if(sensorId > 0)
  {
    numBinX   = 75;
    maxX      = 320.0;
    minX      = -160.0;
    numBinY   = 625;
    maxY      = 100.0;
    minY      = -100.0;
    numBinR   = 75;
    maxR      = 320.0;
    minR      = -160.0;
    numBinPhi = 1570;
    maxPhi    = 0.5*TMath::Pi();
    minPhi    = -0.5*TMath::Pi();
  }
  TH1F *h_mXResidual      = new TH1F("h_mXResidual","h_mXResidual",numBinX,minX,maxX);
  TH1F *h_mYResidual      = new TH1F("h_mYResidual","h_mYResidual",numBinY,minY,maxY);
  TH1F *h_mRResidual      = new TH1F("h_mRResidual","h_mRResidual",numBinR,minR,maxR);
  TH1F *h_mPhiResidual    = new TH1F("h_mPhiResidual","h_mPhiResidual",numBinPhi,minPhi,maxPhi);

  TH1F *h_mXCutResidual      = new TH1F("h_mXCutResidual","h_mXCutResidual",numBinX,minX,maxX);
  TH1F *h_mYCutResidual      = new TH1F("h_mYCutResidual","h_mYCutResidual",numBinY,minY,maxY);
  TH1F *h_mRCutResidual      = new TH1F("h_mRCutResidual","h_mRCutResidual",numBinR,minR,maxR);
  TH1F *h_mPhiCutResidual    = new TH1F("h_mPhiCutResidual","h_mPhiCutResidual",numBinPhi,minPhi,maxPhi);

  h_mXResidual->SetTitle("Corrected X-residual FST Scan Clusters");
  h_mXResidual->GetXaxis()->SetTitle("x-residual (mm)");
  h_mXResidual->GetXaxis()->SetTitleSize(0.06);
  h_mXResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mXResidual->GetYaxis()->SetTitleSize(0.06);
  h_mXResidual->SetLineColor(1);
  h_mXResidual->SetLineWidth(1);
  h_mXCutResidual->SetMarkerStyle(24);
  h_mXCutResidual->SetMarkerColor(2);
  h_mXCutResidual->SetMarkerSize(0.9);

  h_mYResidual->SetTitle("Corrected Y-residual FST Scan Clusters");
  h_mYResidual->GetXaxis()->SetTitle("Y-residual (mm)");
  h_mYResidual->GetXaxis()->SetTitleSize(0.06);
  h_mYResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mYResidual->GetYaxis()->SetTitleSize(0.06);
  h_mYResidual->SetLineColor(1);
  h_mYResidual->SetLineWidth(1);
  h_mYCutResidual->SetMarkerStyle(24);
  h_mYCutResidual->SetMarkerColor(2);
  h_mYCutResidual->SetMarkerSize(0.9);

  h_mRResidual->SetTitle("Corrected R-residual FST Scan Clusters");
  h_mRResidual->GetXaxis()->SetTitle("r-residual (mm)");
  h_mRResidual->GetXaxis()->SetTitleSize(0.06);
  h_mRResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mRResidual->GetYaxis()->SetTitleSize(0.06);
  h_mRResidual->SetLineColor(1);
  h_mRResidual->SetLineWidth(1);
  h_mRCutResidual->SetMarkerStyle(24);
  h_mRCutResidual->SetMarkerColor(2);
  h_mRCutResidual->SetMarkerSize(0.9);

  h_mPhiResidual->SetTitle("Corrected #phi-residual FST Scan Clusters");
  h_mPhiResidual->GetXaxis()->SetTitle("#phi-residual (rad)");
  h_mPhiResidual->GetXaxis()->SetTitleSize(0.06);
  h_mPhiResidual->GetXaxis()->SetRangeUser(-0.5*TMath::Pi(),0.2*TMath::Pi());
  h_mPhiResidual->GetYaxis()->SetTitle("No. Tracks");
  h_mPhiResidual->GetYaxis()->SetTitleSize(0.06);
  h_mPhiResidual->SetLineColor(1);
  h_mPhiResidual->SetLineWidth(1);
  h_mPhiCutResidual->SetMarkerStyle(24);
  h_mPhiCutResidual->SetMarkerColor(2);
  h_mPhiCutResidual->SetMarkerSize(0.9);

  TLegend *leg = new TLegend(0.18,0.75,0.5,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(h_mXResidual,"Current Pars","l");
  leg->AddEntry(h_mXCutResidual,"Next Fit","p");
  // -------------------------------------

  TCanvas *c_Residual = new TCanvas("c_Residual","c_Residual",10,10,800,800);
  c_Residual->Divide(2,2);
  for(int i_pad = 0; i_pad < 4; ++i_pad)
  {
    c_Residual->cd(i_pad+1)->SetLeftMargin(0.15);
    c_Residual->cd(i_pad+1)->SetBottomMargin(0.15);
    c_Residual->cd(i_pad+1)->SetTicks(1,1);
    c_Residual->cd(i_pad+1)->SetGrid(0,0);
  }
  string outputname = Form("./figures/Sensor%dResidual_FSTClusters_3Layer.pdf",sensorId);
  string output_start = Form("./figures/Sensor%dResidual_FSTClusters_3Layer.pdf[",sensorId);
  c_Residual->Print(output_start.c_str()); // open pdf file

  //------------------------------
  // initial fit
  int nFit = 0;
  cout << "Iteration " << nFit << " ==> " << endl;
  tPars initPars = std::make_tuple(0.0,0.0,0.0,0.0,1000.0,1000.0,0); // set initial fit parameters
  tPars fitPars = minuitAlignment(x0_fst, y0_fst, x1_ist, y1_ist, x3_ist, y3_ist, initPars);
  cout << endl;
  //------------------------------
  // plotting initial distribution
  {
    double phi_ist2_temp = std::get<0>(initPars); // aligment parameters from previous iteration
    double x2_shift_temp = std::get<1>(initPars);
    double y2_shift_temp = std::get<2>(initPars);
    double z0_shift_temp = std::get<3>(initPars);
    double xCut_temp     = std::get<4>(initPars);
    double yCut_temp     = std::get<5>(initPars);

    // cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;
    h_mXResidual->Reset();
    h_mYResidual->Reset();
    h_mRResidual->Reset();
    h_mPhiResidual->Reset();
    h_mXCutResidual->Reset();
    h_mYCutResidual->Reset();
    h_mRCutResidual->Reset();
    h_mPhiCutResidual->Reset();
    for(int i_cluster = 0; i_cluster < x1_ist.size(); i_cluster++)
    {
      double x1_corr_temp = x1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
      double y1_corr_temp = y1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

      double x3_corr_temp = x3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
      double y3_corr_temp = y3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

      double x0_proj_temp   = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
      double y0_proj_temp   = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
      double r0_proj_temp   = TMath::Sqrt(x0_proj_temp*x0_proj_temp+y0_proj_temp*y0_proj_temp);
      double phi0_proj_temp = TMath::ATan2(y0_proj_temp,x0_proj_temp);

      double dx   = x0_fst[i_cluster] - x0_proj_temp;
      double dy   = y0_fst[i_cluster] - y0_proj_temp;
      double dr   = r0_fst[i_cluster] - r0_proj_temp;
      double dphi = phi0_fst[i_cluster] - phi0_proj_temp;
      // cout << "dx = " << dx << ", dy = " << dy << endl;

      h_mXResidual->Fill(dx);
      h_mYResidual->Fill(dy);
      h_mRResidual->Fill(dr);
      h_mPhiResidual->Fill(dphi);
      if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
      {
	h_mXCutResidual->Fill(dx);
	h_mYCutResidual->Fill(dy);
	h_mRCutResidual->Fill(dr);
	h_mPhiCutResidual->Fill(dphi);
      }
    }

    c_Residual->cd();
    string TitleName = "Without Alignment";
    TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.1,1.0,TitleName.c_str());
    canvas_title->SetBorderSize(1);
    canvas_title->Draw();
    c_Residual->cd(1);
    h_mXResidual->DrawCopy("hE");

    c_Residual->cd(2);
    h_mYResidual->DrawCopy("hE");

    c_Residual->cd(3);
    h_mRResidual->DrawCopy("hE");

    c_Residual->cd(4);
    h_mPhiResidual->DrawCopy("hE");
    c_Residual->Update();
    c_Residual->Print(outputname.c_str());
  }

  // reject trakcs with xCut = 50 and yCut = 5
  tPars itPars = updateFitParameters(fitPars, 50, 10, 1);
  int nOffset_temp = std::get<6>(itPars);
  while(nOffset_temp > 0)
  { 
    // plotting
    {
      double phi_ist2_temp = std::get<0>(itPars); // aligment parameters from previous iteration
      double x2_shift_temp = std::get<1>(itPars);
      double y2_shift_temp = std::get<2>(itPars);
      double z0_shift_temp = std::get<3>(itPars);
      double xCut_temp     = std::get<4>(itPars);
      double yCut_temp     = std::get<5>(itPars);

      // cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;
      h_mXResidual->Reset();
      h_mYResidual->Reset();
      h_mRResidual->Reset();
      h_mPhiResidual->Reset();
      h_mXCutResidual->Reset();
      h_mYCutResidual->Reset();
      h_mRCutResidual->Reset();
      h_mPhiCutResidual->Reset();
      for(int i_cluster = 0; i_cluster < x1_ist.size(); i_cluster++)
      {
	double x1_corr_temp = x1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y1_corr_temp = y1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x3_corr_temp = x3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y3_corr_temp = y3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x0_proj_temp   = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double y0_proj_temp   = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double r0_proj_temp   = TMath::Sqrt(x0_proj_temp*x0_proj_temp+y0_proj_temp*y0_proj_temp);
	double phi0_proj_temp = TMath::ATan2(y0_proj_temp,x0_proj_temp);

	double dx   = x0_fst[i_cluster] - x0_proj_temp;
	double dy   = y0_fst[i_cluster] - y0_proj_temp;
	double dr   = r0_fst[i_cluster] - r0_proj_temp;
	double dphi = phi0_fst[i_cluster] - phi0_proj_temp;
	// cout << "dx = " << dx << ", dy = " << dy << endl;

	h_mXResidual->Fill(dx);
	h_mYResidual->Fill(dy);
	h_mRResidual->Fill(dr);
	h_mPhiResidual->Fill(dphi);
	if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
	{
	  h_mXCutResidual->Fill(dx);
	  h_mYCutResidual->Fill(dy);
	  h_mRCutResidual->Fill(dr);
	  h_mPhiCutResidual->Fill(dphi);
	}
      }

      c_Residual->cd();
      string TitleName = Form("Iteration %d",nFit);
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.1,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Residual->cd(1);
      h_mXResidual->DrawCopy("h");
      h_mXCutResidual->DrawCopy("p same");
      PlotLine(xCut_temp, xCut_temp, 0.0, 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-xCut_temp, -xCut_temp, 0., 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      leg->Draw("same");

      c_Residual->cd(2);
      h_mYResidual->GetXaxis()->SetRangeUser(-15.0,15.0);
      h_mYResidual->DrawCopy("h");
      h_mYCutResidual->DrawCopy("p same");
      PlotLine(yCut_temp, yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-yCut_temp, -yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);

      c_Residual->cd(3);
      h_mPhiResidual->GetXaxis()->SetRangeUser(-0.1,0.1);
      h_mRResidual->DrawCopy("h");
      h_mRCutResidual->DrawCopy("p same");

      c_Residual->cd(4);
      h_mPhiResidual->DrawCopy("h");
      h_mPhiCutResidual->DrawCopy("p same");
      c_Residual->Update();
      c_Residual->Print(outputname.c_str());
    }

    // fit
    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x0_fst, y0_fst, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<6>(fitPars);
    itPars = updateFitParameters(fitPars,50,10,nOffset_temp);
    cout << endl;
  }

  // reject trakcs with xCut = 25 and yCut = 5 
  itPars = updateFitParameters(fitPars, 25, 5, 1);
  nOffset_temp = std::get<6>(itPars);
  while(nOffset_temp > 0)
  {
    // plotting
    {
      double phi_ist2_temp = std::get<0>(itPars); // aligment parameters from previous iteration
      double x2_shift_temp = std::get<1>(itPars);
      double y2_shift_temp = std::get<2>(itPars);
      double z0_shift_temp = std::get<3>(itPars);
      double xCut_temp     = std::get<4>(itPars);
      double yCut_temp     = std::get<5>(itPars);

      // cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;
      h_mXResidual->Reset();
      h_mYResidual->Reset();
      h_mRResidual->Reset();
      h_mPhiResidual->Reset();
      h_mXCutResidual->Reset();
      h_mYCutResidual->Reset();
      h_mRCutResidual->Reset();
      h_mPhiCutResidual->Reset();
      for(int i_cluster = 0; i_cluster < x1_ist.size(); i_cluster++)
      {
	double x1_corr_temp = x1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y1_corr_temp = y1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x3_corr_temp = x3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y3_corr_temp = y3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x0_proj_temp   = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double y0_proj_temp   = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double r0_proj_temp   = TMath::Sqrt(x0_proj_temp*x0_proj_temp+y0_proj_temp*y0_proj_temp);
	double phi0_proj_temp = TMath::ATan2(y0_proj_temp,x0_proj_temp);

	double dx   = x0_fst[i_cluster] - x0_proj_temp;
	double dy   = y0_fst[i_cluster] - y0_proj_temp;
	double dr   = r0_fst[i_cluster] - r0_proj_temp;
	double dphi = phi0_fst[i_cluster] - phi0_proj_temp;
	// cout << "dx = " << dx << ", dy = " << dy << endl;

	h_mXResidual->Fill(dx);
	h_mYResidual->Fill(dy);
	h_mRResidual->Fill(dr);
	h_mPhiResidual->Fill(dphi);
	if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
	{
	  h_mXCutResidual->Fill(dx);
	  h_mYCutResidual->Fill(dy);
	  h_mRCutResidual->Fill(dr);
	  h_mPhiCutResidual->Fill(dphi);
	}
      }

      c_Residual->cd();
      string TitleName = Form("Iteration %d",nFit);
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.1,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Residual->cd(1);
      h_mXResidual->DrawCopy("h");
      h_mXCutResidual->DrawCopy("p same");
      PlotLine(xCut_temp, xCut_temp, 0.0, 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-xCut_temp, -xCut_temp, 0., 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      leg->Draw("same");

      c_Residual->cd(2);
      h_mYResidual->DrawCopy("h");
      h_mYCutResidual->DrawCopy("p same");
      PlotLine(yCut_temp, yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-yCut_temp, -yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);

      c_Residual->cd(3);
      h_mRResidual->DrawCopy("h");
      h_mRCutResidual->DrawCopy("p same");

      c_Residual->cd(4);
      h_mPhiResidual->DrawCopy("h");
      h_mPhiCutResidual->DrawCopy("p same");
      c_Residual->Update();
      c_Residual->Print(outputname.c_str());
    }

    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x0_fst, y0_fst, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<6>(fitPars);
    itPars = updateFitParameters(fitPars,25,5,nOffset_temp);
    cout << endl;
  }

  itPars = updateFitParameters(fitPars, 15, 2.5, 1);
  nOffset_temp = std::get<6>(itPars);
  while(nOffset_temp > 0)
  {
    // plotting
    {
      double phi_ist2_temp = std::get<0>(itPars); // aligment parameters from previous iteration
      double x2_shift_temp = std::get<1>(itPars);
      double y2_shift_temp = std::get<2>(itPars);
      double z0_shift_temp = std::get<3>(itPars);
      double xCut_temp     = std::get<4>(itPars);
      double yCut_temp     = std::get<5>(itPars);

      // cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;
      h_mXResidual->Reset();
      h_mYResidual->Reset();
      h_mRResidual->Reset();
      h_mPhiResidual->Reset();
      h_mXCutResidual->Reset();
      h_mYCutResidual->Reset();
      h_mRCutResidual->Reset();
      h_mPhiCutResidual->Reset();
      for(int i_cluster = 0; i_cluster < x1_ist.size(); i_cluster++)
      {
	double x1_corr_temp = x1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y1_corr_temp = y1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x3_corr_temp = x3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
	double y3_corr_temp = y3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

	double x0_proj_temp   = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double y0_proj_temp   = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
	double r0_proj_temp   = TMath::Sqrt(x0_proj_temp*x0_proj_temp+y0_proj_temp*y0_proj_temp);
	double phi0_proj_temp = TMath::ATan2(y0_proj_temp,x0_proj_temp);

	double dx   = x0_fst[i_cluster] - x0_proj_temp;
	double dy   = y0_fst[i_cluster] - y0_proj_temp;
	double dr   = r0_fst[i_cluster] - r0_proj_temp;
	double dphi = phi0_fst[i_cluster] - phi0_proj_temp;
	// cout << "dx = " << dx << ", dy = " << dy << endl;

	h_mXResidual->Fill(dx);
	h_mYResidual->Fill(dy);
	h_mRResidual->Fill(dr);
	h_mPhiResidual->Fill(dphi);
	if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
	{
	  h_mXCutResidual->Fill(dx);
	  h_mYCutResidual->Fill(dy);
	  h_mRCutResidual->Fill(dr);
	  h_mPhiCutResidual->Fill(dphi);
	}
      }

      c_Residual->cd();
      string TitleName = Form("Iteration %d",nFit);
      TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.1,1.0,TitleName.c_str());
      canvas_title->SetBorderSize(1);
      canvas_title->Draw();
      c_Residual->cd(1);
      h_mXResidual->DrawCopy("h");
      h_mXCutResidual->DrawCopy("p same");
      PlotLine(xCut_temp, xCut_temp, 0.0, 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-xCut_temp, -xCut_temp, 0., 0.8*h_mXResidual->GetMaximum(), 4, 2, 2);
      leg->Draw("same");

      c_Residual->cd(2);
      h_mYResidual->DrawCopy("h");
      h_mYCutResidual->DrawCopy("p same");
      PlotLine(yCut_temp, yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);
      PlotLine(-yCut_temp, -yCut_temp, 0.0, 0.8*h_mYResidual->GetMaximum(), 4, 2, 2);

      c_Residual->cd(3);
      h_mRResidual->DrawCopy("h");
      h_mRCutResidual->DrawCopy("p same");

      c_Residual->cd(4);
      h_mPhiResidual->DrawCopy("h");
      h_mPhiCutResidual->DrawCopy("p same");
      c_Residual->Update();
      c_Residual->Print(outputname.c_str());
    }

    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x0_fst, y0_fst, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<6>(fitPars);
    itPars = updateFitParameters(fitPars,15,2.5,nOffset_temp);
    cout << endl;
  }

  cout << "Finish Minuit Fit for alignment => " << endl;
  // plotting final alignment
  {
    double phi_ist2_temp = std::get<0>(fitPars); // aligment parameters from previous iteration
    double x2_shift_temp = std::get<1>(fitPars);
    double y2_shift_temp = std::get<2>(fitPars);
    double z0_shift_temp = std::get<3>(fitPars);
    double xCut_temp     = std::get<4>(fitPars);
    double yCut_temp     = std::get<5>(fitPars);

    // cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;
    h_mXResidual->Reset();
    h_mYResidual->Reset();
    h_mRResidual->Reset();
    h_mPhiResidual->Reset();
    h_mXCutResidual->Reset();
    h_mYCutResidual->Reset();
    h_mRCutResidual->Reset();
    h_mPhiCutResidual->Reset();
    for(int i_cluster = 0; i_cluster < x1_ist.size(); i_cluster++)
    {
      double x1_corr_temp = x1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
      double y1_corr_temp = y1_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

      double x3_corr_temp = x3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
      double y3_corr_temp = y3_ist[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_ist[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

      double x0_proj_temp   = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
      double y0_proj_temp   = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
      double r0_proj_temp   = TMath::Sqrt(x0_proj_temp*x0_proj_temp+y0_proj_temp*y0_proj_temp);
      double phi0_proj_temp = TMath::ATan2(y0_proj_temp,x0_proj_temp);

      double dx   = x0_fst[i_cluster] - x0_proj_temp;
      double dy   = y0_fst[i_cluster] - y0_proj_temp;
      double dr   = r0_fst[i_cluster] - r0_proj_temp;
      double dphi = phi0_fst[i_cluster] - phi0_proj_temp;
      // cout << "dx = " << dx << ", dy = " << dy << endl;

      h_mXResidual->Fill(dx);
      h_mYResidual->Fill(dy);
      h_mRResidual->Fill(dr);
      h_mPhiResidual->Fill(dphi);
      if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
      {
	h_mXCutResidual->Fill(dx);
	h_mYCutResidual->Fill(dy);
	h_mRCutResidual->Fill(dr);
	h_mPhiCutResidual->Fill(dphi);
      }
    }

    c_Residual->cd();
    string TitleName = "Final Alignment";
    TPaveLabel* canvas_title = new TPaveLabel(0.0,0.98,0.1,1.0,TitleName.c_str());
    canvas_title->SetBorderSize(1);
    canvas_title->Draw();
    c_Residual->cd(1);
    if(sensorId > 0) h_mXResidual->GetXaxis()->SetRangeUser(-160.0,160.0);
    h_mXResidual->DrawCopy("hE");
    TF1 *f_gausX = new TF1("f_gausX",gaussian,-150.0,150.0,4);
    f_gausX->SetParameter(0,100.0);
    f_gausX->SetParameter(1,0.0);
    f_gausX->SetParameter(2,10.0);
    f_gausX->FixParameter(3,h_mXResidual->GetBinWidth(1));
    f_gausX->SetRange(-30,30);
    h_mXResidual->Fit(f_gausX,"R");
    f_gausX->Draw("l same");

    c_Residual->cd(2);
    h_mYResidual->DrawCopy("hE");
    TF1 *f_gausY = new TF1("f_gausY",gaussian,-150.0,150.0,4);
    f_gausY->SetParameter(0,100.0);
    f_gausY->SetParameter(1,0.0);
    f_gausY->SetParameter(2,10.0);
    f_gausY->FixParameter(3,h_mYResidual->GetBinWidth(1));
    f_gausY->SetRange(-5.0,5.0);
    h_mYResidual->Fit(f_gausY,"R");
    f_gausY->Draw("l same");

    c_Residual->cd(3);
    if(sensorId > 0) h_mRResidual->GetXaxis()->SetRangeUser(-160.0,160.0);
    h_mRResidual->DrawCopy("hE");
    TF1 *f_gausR = new TF1("f_gausR",gaussian,-150.0,150.0,4);
    f_gausR->SetParameter(0,100.0);
    f_gausR->SetParameter(1,0.0);
    f_gausR->SetParameter(2,10.0);
    f_gausR->FixParameter(3,h_mRResidual->GetBinWidth(1));
    f_gausR->SetRange(-30,30);
    h_mRResidual->Fit(f_gausR,"R");
    f_gausR->Draw("l same");

    c_Residual->cd(4);
    h_mPhiResidual->DrawCopy("hE");
    TF1 *f_gausPhi = new TF1("f_gausPhi",gaussian,-1.0,1.0,4);
    f_gausPhi->SetParameter(0,100.0);
    f_gausPhi->SetParameter(1,0.0);
    f_gausPhi->SetParameter(2,10.0);
    f_gausPhi->FixParameter(3,h_mPhiResidual->GetBinWidth(1));
    f_gausPhi->SetRange(-0.02,0.02);
    h_mPhiResidual->Fit(f_gausPhi,"R");
    f_gausPhi->Draw("l same");

    c_Residual->Update();
    c_Residual->Print(outputname.c_str());
  }

  cout << "Minuit minimization: phi_rot_ist2 = " << std::get<0>(fitPars) << ", x2_shift = " << std::get<1>(fitPars) << ", y2_shift = " << std::get<2>(fitPars) << ", z0_shift = " << std::get<3>(fitPars) << endl;

  string output_stop = Form("./figures/Sensor%dResidual_FSTClusters_3Layer.pdf]",sensorId); 
  c_Residual->Print(output_stop.c_str()); // close pdf file

  return 1;
}

tPars minuitAlignment(dVec x0_orig, dVec y0_orig, dVec x1_orig, dVec y1_orig, dVec x3_orig, dVec y3_orig, tPars cParameters)
{
  dVec x0_temp, x0_fit;
  dVec y0_temp, y0_fit;
  dVec x1_temp, x1_fit;
  dVec y1_temp, y1_fit;
  dVec x3_temp, x3_fit;
  dVec y3_temp, y3_fit;

  // clear vector
  x0_temp.clear();
  y0_temp.clear();
  x1_temp.clear();
  y1_temp.clear();
  x3_temp.clear();
  y3_temp.clear();

  x0_fit.clear();
  y0_fit.clear();
  x1_fit.clear();
  y1_fit.clear();
  x3_fit.clear();
  y3_fit.clear();

  for(int i_cluster = 0; i_cluster < x1_orig.size(); ++i_cluster)
  { // set hit containers
    x0_temp.push_back(x0_orig[i_cluster]);
    y0_temp.push_back(y0_orig[i_cluster]);
    x1_temp.push_back(x1_orig[i_cluster]);
    y1_temp.push_back(y1_orig[i_cluster]);
    x3_temp.push_back(x3_orig[i_cluster]);
    y3_temp.push_back(y3_orig[i_cluster]);
  }

  // reject offset trakcs
  const double phi_ist2_temp = std::get<0>(cParameters); // aligment parameters from previous iteration
  const double x2_shift_temp = std::get<1>(cParameters);
  const double y2_shift_temp = std::get<2>(cParameters);
  const double z0_shift_temp = std::get<3>(cParameters);
  const double xCut_temp     = std::get<4>(cParameters);
  const double yCut_temp     = std::get<5>(cParameters);
  cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;

  int numOfUsedHits = 0;
  for(int i_cluster = 0; i_cluster < x1_temp.size(); i_cluster++)
  {
    double x1_corr_temp = x1_temp[i_cluster]*TMath::Cos(phi_ist2_temp) + y1_temp[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
    double y1_corr_temp = y1_temp[i_cluster]*TMath::Cos(phi_ist2_temp) - x1_temp[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

    double x3_corr_temp = x3_temp[i_cluster]*TMath::Cos(phi_ist2_temp) + y3_temp[i_cluster]*TMath::Sin(phi_ist2_temp) + x2_shift_temp;
    double y3_corr_temp = y3_temp[i_cluster]*TMath::Cos(phi_ist2_temp) - x3_temp[i_cluster]*TMath::Sin(phi_ist2_temp) + y2_shift_temp;

    double x0_proj_temp = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);
    double y0_proj_temp = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*(FST::pitchLayer03+z0_shift_temp)/(FST::pitchLayer12+FST::pitchLayer23);

    double dx = x0_temp[i_cluster] - x0_proj_temp;
    double dy = y0_temp[i_cluster] - y0_proj_temp;
    // cout << "dx = " << dx << ", dy = " << dy << endl;

    if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
    {
      x0_fit.push_back(x0_temp[i_cluster]);
      y0_fit.push_back(y0_temp[i_cluster]);
      x1_fit.push_back(x1_temp[i_cluster]);
      y1_fit.push_back(y1_temp[i_cluster]);
      x3_fit.push_back(x3_temp[i_cluster]);
      y3_fit.push_back(y3_temp[i_cluster]);
      numOfUsedHits++;
    }
  }
  cout << "numOfUsedTracks = " << numOfUsedHits << endl;

  cout << "Start the Minuit minimization!" << endl;
  auto chi2Function = [&](const Double_t *par)
  { //minimisation function computing the sum of squares of residuals
    double phi_rot_ist2 = par[0];
    double x2_shift     = par[1];
    double y2_shift     = par[2];
    double z0_shift     = par[3];

    const double x_weight = 1.0;
    const double y_weight = 400.0;
    double f = 0;
    for (int i_cluster = 0; i_cluster < numOfUsedHits; i_cluster++)
    {
      double x1_corr = x1_fit[i_cluster]*TMath::Cos(phi_rot_ist2) + y1_fit[i_cluster]*TMath::Sin(phi_rot_ist2) + x2_shift;
      double y1_corr = y1_fit[i_cluster]*TMath::Cos(phi_rot_ist2) - x1_fit[i_cluster]*TMath::Sin(phi_rot_ist2) + y2_shift;

      double x3_corr = x3_fit[i_cluster]*TMath::Cos(phi_rot_ist2) + y3_fit[i_cluster]*TMath::Sin(phi_rot_ist2) + x2_shift;
      double y3_corr = y3_fit[i_cluster]*TMath::Cos(phi_rot_ist2) - x3_fit[i_cluster]*TMath::Sin(phi_rot_ist2) + y2_shift;

      double x0_proj = x3_corr + (x1_corr-x3_corr)*(FST::pitchLayer03+z0_shift)/(FST::pitchLayer12+FST::pitchLayer23);
      double y0_proj = y3_corr + (y1_corr-y3_corr)*(FST::pitchLayer03+z0_shift)/(FST::pitchLayer12+FST::pitchLayer23);

      double dx = x0_fit[i_cluster] - x0_proj;
      double dy = y0_fit[i_cluster] - y0_proj;

      f += x_weight*dx*dx + y_weight*dy*dy;
    }
    return f;
  };

  // wrap chi2 funciton in a function object for the fit
  // 3 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(chi2Function,4);
  ROOT::Fit::Fitter fitter;

  double pStart[4] = {phi_ist2_temp,x2_shift_temp,y2_shift_temp,z0_shift_temp};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("phi_rot_ist2");
  fitter.Config().ParSettings(0).SetLimits(0.5*TMath::Pi(),1.5*TMath::Pi());
  fitter.Config().ParSettings(1).SetName("x2_shift");
  fitter.Config().ParSettings(2).SetName("y2_shift");
  fitter.Config().ParSettings(3).SetName("z0_shift");

  // do the fit
  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("find alignment","Minuit Fit failed");
    return std::make_tuple(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1);
  }

  const ROOT::Fit::FitResult & result = fitter.Result();
  // result.Print(std::cout);
  const double *fitpar = result.GetParams();
  double phi_rot_ist2_MF = fitpar[0];
  double x2_shift_MF     = fitpar[1];
  double y2_shift_MF     = fitpar[2];
  double z0_shift_MF     = fitpar[3];

  // calculate nOffset with current fit parameters
  int nOffset = 0;
  for (int i_cluster = 0; i_cluster < numOfUsedHits; i_cluster++)
  {
    double x1_corr_MF = x1_fit[i_cluster]*TMath::Cos(phi_rot_ist2_MF) + y1_fit[i_cluster]*TMath::Sin(phi_rot_ist2_MF) + x2_shift_MF;
    double y1_corr_MF = y1_fit[i_cluster]*TMath::Cos(phi_rot_ist2_MF) - x1_fit[i_cluster]*TMath::Sin(phi_rot_ist2_MF) + y2_shift_MF;

    double x3_corr_MF = x3_fit[i_cluster]*TMath::Cos(phi_rot_ist2_MF) + y3_fit[i_cluster]*TMath::Sin(phi_rot_ist2_MF) + x2_shift_MF;
    double y3_corr_MF = y3_fit[i_cluster]*TMath::Cos(phi_rot_ist2_MF) - x3_fit[i_cluster]*TMath::Sin(phi_rot_ist2_MF) + y2_shift_MF;

    double x0_proj_MF = x3_corr_MF + (x1_corr_MF-x3_corr_MF)*(FST::pitchLayer03+z0_shift_MF)/(FST::pitchLayer12+FST::pitchLayer23);
    double y0_proj_MF = y3_corr_MF + (y1_corr_MF-y3_corr_MF)*(FST::pitchLayer03+z0_shift_MF)/(FST::pitchLayer12+FST::pitchLayer23);

    double dx = x0_fit[i_cluster] - x0_proj_MF;
    double dy = y0_fit[i_cluster] - y0_proj_MF;

    if( !(abs(dx) < xCut_temp && abs(dy) < yCut_temp) )
    {
      nOffset++;
    }
  }

  tPars tParameters = std::make_tuple(phi_rot_ist2_MF, x2_shift_MF, y2_shift_MF, z0_shift_MF, xCut_temp, yCut_temp, nOffset);

  // cout << "Minuit minimization: phi_rot_ist1 = " << phi_rot_ist1_MF << ", phi_rot_ist3 = " << phi_rot_ist3_MF << ", x_shift = " << x_shift_MF << ", y_shift = " << y_shift_MF << endl;
  cout << "Minuit minimization: phi_rot_ist2 = " << std::get<0>(tParameters) << ", x2_shift = " << (double)std::get<1>(tParameters) << ", y2_shift = " << std::get<2>(tParameters) << ", z0_shift = " << std::get<3>(tParameters) << ", nOffset = " << std::get<6>(tParameters) << endl;

  return tParameters;
}

tPars updateFitParameters(tPars fitPars, double xCut, double yCut, int nOffset)
{
  const double phi_ist2_old = std::get<0>(fitPars);
  const double x2_shift_old = std::get<1>(fitPars);
  const double y2_shift_old = std::get<2>(fitPars);
  const double z0_shift_old = std::get<3>(fitPars);
  const double xCut_old     = std::get<4>(fitPars);
  const double yCut_old     = std::get<5>(fitPars);
  const double nOffset_old  = std::get<6>(fitPars);

  tPars newPars = std::make_tuple(phi_ist2_old, x2_shift_old, y2_shift_old, z0_shift_old, xCut, yCut, nOffset);

  return newPars;
}
