#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>

#include <TMath.h>
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

typedef std::vector<double> dVec;
typedef std::tuple<double, double, double, double, double, double, double, double, int> tPars; // phi_rot_ist1, x1_shift, y1_shift, phi_rot_ist3, x3_shift, y3_shift, xCut, yCut, nOffset

tPars minuitAlignment(dVec x2_orig, dVec y2_orig, dVec x1_orig, dVec y1_orig, dVec x3_orig, dVec y3_orig, tPars cParameters);
tPars updateFitParameters(tPars fitPars, double xCut, double yCut, int nOffset);

R__LOAD_LIBRARY(../../lib/libFstEvent.dylib)

int getAlignment_ISTClusters()
{
  dVec x1_ist;
  dVec y1_ist;
  dVec x2_ist;
  dVec y2_ist;
  dVec x3_ist;
  dVec y3_ist;

  // clear vector
  x1_ist.clear();
  y1_ist.clear();
  x2_ist.clear();
  y2_ist.clear();
  x3_ist.clear();
  y3_ist.clear();

  FstEvent *mFstEvent = new FstEvent();
  std::string inputfile = "../../output/FstClusters_HV140V_woPed.root";
  std::cout << "inputfile = " << inputfile.c_str() << std::endl;
  TFile *mFile_InPut = TFile::Open(inputfile.c_str());
  TTree *mTree_FstEvent = (TTree*)mFile_InPut->Get("mTree_FstEvent");
  mTree_FstEvent->SetBranchAddress("FstEvent",&mFstEvent);

  long NumOfEvents = (long)mTree_FstEvent->GetEntries();
  cout << "total number of events: " << NumOfEvents << endl;
  // if(NumOfEvents > 1000) NumOfEvents = 1000;
  // NumOfEvents = 500;
  mTree_FstEvent->GetEntry(0);

  for(int i_event = 0; i_event < NumOfEvents; ++i_event)
  {
    if(i_event%1000==0) cout << "processing events:  " << i_event << "/" << NumOfEvents << endl;
    mTree_FstEvent->GetEntry(i_event);

    std::vector<FstCluster *> clusterVec_ist1;
    std::vector<FstCluster *> clusterVec_ist2;
    std::vector<FstCluster *> clusterVec_ist3;
    clusterVec_ist1.clear();
    clusterVec_ist2.clear();
    clusterVec_ist3.clear();

    for(int i_cluster = 0; i_cluster < mFstEvent->getNumClusters(); ++i_cluster)
    { // get Hits info for IST1 IST2 and IST3
      FstCluster *fstCluster = mFstEvent->getCluster(i_cluster);
      if(fstCluster->getLayer() == 1)
      {
	clusterVec_ist1.push_back(fstCluster);
      }
      if(fstCluster->getLayer() == 2)
      {
	clusterVec_ist2.push_back(fstCluster);
      }
      if(fstCluster->getLayer() == 3)
      {
	clusterVec_ist3.push_back(fstCluster);
      }
    }

    const double xMax = 24.0*FST::pitchColumn; // center of 23 column + 0.5*pitchColumn
    const double xMin = 20.0*FST::pitchColumn; // center of 20 column - 0.5*pitchColumn
    if(clusterVec_ist2.size() > 0) // IST2 has at least one hit
    {
      for(int i_ist1 = 0; i_ist1 < clusterVec_ist1.size(); ++i_ist1)
      {
	double x1 = clusterVec_ist1[i_ist1]->getMeanX();
	double y1 = clusterVec_ist1[i_ist1]->getMeanY();
	const double z1_ist = FST::pitchLayer12 + FST::pitchLayer23;
	double col1 = clusterVec_ist1[i_ist1]->getMeanColumn();
	double row1 = clusterVec_ist1[i_ist1]->getMeanRow();
	for(int i_ist3 = 0; i_ist3 < clusterVec_ist3.size(); ++i_ist3)
	{
	  double x3 = clusterVec_ist3[i_ist3]->getMeanX();
	  double y3 = clusterVec_ist3[i_ist3]->getMeanY();
	  const double z3_ist = 0.0;
	  double col3 = clusterVec_ist3[i_ist3]->getMeanColumn();
	  double row3 = clusterVec_ist3[i_ist3]->getMeanRow();

	  double x2_proj = x3 + (x1-x3)*FST::pitchLayer23/(FST::pitchLayer12+FST::pitchLayer23);
	  double y2_proj = y3 + (y1-y3)*FST::pitchLayer23/(FST::pitchLayer12+FST::pitchLayer23);

	  for(int i_ist2 = 0; i_ist2 < clusterVec_ist2.size(); ++i_ist2)
	  {
	    double x2 = clusterVec_ist2[i_ist2]->getMeanX();
	    double y2 = clusterVec_ist2[i_ist2]->getMeanY();
	    const double z2_ist = FST::pitchLayer23;

	    if( abs(row3-row1) < 17 && x2_proj >= xMin && x2_proj <= xMax)
	    {
	      x1_ist.push_back(x1);
	      y1_ist.push_back(y1);
	      x2_ist.push_back(x2);
	      y2_ist.push_back(y2);
	      x3_ist.push_back(x3);
	      y3_ist.push_back(y3);
	    }
	  }
	}
      }
    }
  }

  mFile_InPut->Close();

  cout << "Start Minuit Fit for alignment => " << endl;

  int nFit = 0;
  cout << "Iteration " << nFit << " ==> " << endl;
  tPars initPars = std::make_tuple(0.0,0.0,0.0,0.0,0.0,0.0,1000.0,1000.0,0); // set initial fit parameters
  tPars fitPars = minuitAlignment(x2_ist, y2_ist, x1_ist, y1_ist, x3_ist, y3_ist, initPars);
  cout << endl;

  // reject trakcs with xCut = 50 and yCut = 5
  tPars itPars = updateFitParameters(fitPars, 50, 5, 1);
  int nOffset_temp = std::get<8>(itPars);
  while(nOffset_temp > 0)
  { 
    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x2_ist, y2_ist, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<8>(fitPars);
    itPars = updateFitParameters(fitPars,50,5,nOffset_temp);
    cout << endl;
  }

  // reject trakcs with xCut = 25 and yCut = 2.5 
  itPars = updateFitParameters(fitPars, 25, 2.5, 1);
  nOffset_temp = std::get<8>(itPars);
  while(nOffset_temp > 0)
  {
    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x2_ist, y2_ist, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<8>(fitPars);
    itPars = updateFitParameters(fitPars,25,2.5,nOffset_temp);
    cout << endl;
  }

  // reject trakcs with xCut = 12 and yCut = 1.5 
  itPars = updateFitParameters(fitPars, 12, 1.5, 1);
  nOffset_temp = std::get<8>(itPars);
  while(nOffset_temp > 0)
  {
    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x2_ist, y2_ist, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<8>(fitPars);
    itPars = updateFitParameters(fitPars,12,1.5,nOffset_temp);
    cout << endl;
  }

  // reject trakcs with xCut = 6 and yCut = 0.6 
  itPars = updateFitParameters(fitPars, 6, 0.6, 1);
  nOffset_temp = std::get<8>(itPars);
  while(nOffset_temp > 0)
  {
    nFit++;
    cout << "Iteration " << nFit << " ==> " << endl;
    fitPars = minuitAlignment(x2_ist, y2_ist, x1_ist, y1_ist, x3_ist, y3_ist, itPars);
    nOffset_temp = std::get<8>(fitPars);
    itPars = updateFitParameters(fitPars,6,0.6,nOffset_temp);
    cout << endl;
  }

  cout << "Finish Minuit Fit for alignment => " << endl;

  cout << "Minuit minimization: phi_rot_ist1 = " << std::get<0>(fitPars) << ", x1_shift = " << std::get<1>(fitPars) << ", y1_shift = " << std::get<2>(fitPars) << endl;
  cout << "Minuit minimization: phi_rot_ist3 = " << std::get<3>(fitPars) << ", x3_shift = " << std::get<4>(fitPars) << ", y3_shift = " << std::get<5>(fitPars) << endl;

  return 1;
}

tPars minuitAlignment(dVec x2_orig, dVec y2_orig, dVec x1_orig, dVec y1_orig, dVec x3_orig, dVec y3_orig, tPars cParameters)
{
  dVec x2_temp, x2_fit;
  dVec y2_temp, y2_fit;
  dVec x1_temp, x1_fit;
  dVec y1_temp, y1_fit;
  dVec x3_temp, x3_fit;
  dVec y3_temp, y3_fit;

  // clear vector
  x2_temp.clear();
  y2_temp.clear();
  x1_temp.clear();
  y1_temp.clear();
  x3_temp.clear();
  y3_temp.clear();

  x2_fit.clear();
  y2_fit.clear();
  x1_fit.clear();
  y1_fit.clear();
  x3_fit.clear();
  y3_fit.clear();

  for(int i_cluster = 0; i_cluster < x1_orig.size(); ++i_cluster)
  { // set hit containers
    x2_temp.push_back(x2_orig[i_cluster]);
    y2_temp.push_back(y2_orig[i_cluster]);
    x1_temp.push_back(x1_orig[i_cluster]);
    y1_temp.push_back(y1_orig[i_cluster]);
    x3_temp.push_back(x3_orig[i_cluster]);
    y3_temp.push_back(y3_orig[i_cluster]);
  }

  // const double pitchLayer03 = 134.9375; // mm
  const double pitchLayer12 = 34.925; //mm, distances between the 1st & 2nd Layer => Added by Xu Sun @ 02/13/2020
  const double pitchLayer23 = 30.1625; //mm, distances between the 2nd & 3rd Layer

  // reject offset trakcs
  const double phi_ist1_temp = std::get<0>(cParameters); // aligment parameters from previous iteration
  const double x1_shift_temp = std::get<1>(cParameters);
  const double y1_shift_temp = std::get<2>(cParameters);
  const double phi_ist3_temp = std::get<3>(cParameters); // aligment parameters from previous iteration
  const double x3_shift_temp = std::get<4>(cParameters);
  const double y3_shift_temp = std::get<5>(cParameters);
  const double xCut_temp     = std::get<6>(cParameters);
  const double yCut_temp     = std::get<7>(cParameters);
  cout << "reject trakcs out of xCut = " << xCut_temp << " and yCut: " << yCut_temp << endl;

  int numOfUsedHits = 0;
  for(int i_cluster = 0; i_cluster < x1_temp.size(); i_cluster++)
  {
    double x1_corr_temp = x1_temp[i_cluster]*TMath::Cos(phi_ist1_temp) + y1_temp[i_cluster]*TMath::Sin(phi_ist1_temp) + x1_shift_temp;
    double y1_corr_temp = y1_temp[i_cluster]*TMath::Cos(phi_ist1_temp) - x1_temp[i_cluster]*TMath::Sin(phi_ist1_temp) + y1_shift_temp;

    double x3_corr_temp = x3_temp[i_cluster]*TMath::Cos(phi_ist3_temp) + y3_temp[i_cluster]*TMath::Sin(phi_ist3_temp) + x3_shift_temp;
    double y3_corr_temp = y3_temp[i_cluster]*TMath::Cos(phi_ist3_temp) - x3_temp[i_cluster]*TMath::Sin(phi_ist3_temp) + y3_shift_temp;

    double x2_proj_temp = x3_corr_temp + (x1_corr_temp-x3_corr_temp)*pitchLayer23/(pitchLayer12+pitchLayer23);
    double y2_proj_temp = y3_corr_temp + (y1_corr_temp-y3_corr_temp)*pitchLayer23/(pitchLayer12+pitchLayer23);

    double dx = x2_temp[i_cluster] - x2_proj_temp;
    double dy = y2_temp[i_cluster] - y2_proj_temp;
    // cout << "dx = " << dx << ", dy = " << dy << endl;

    if(abs(dx) < xCut_temp && abs(dy) < yCut_temp)
    {
      x2_fit.push_back(x2_temp[i_cluster]);
      y2_fit.push_back(y2_temp[i_cluster]);
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
    double phi_rot_ist1 = par[0];
    double x1_shift     = par[1];
    double y1_shift     = par[2];
    double phi_rot_ist3 = par[3];
    double x3_shift     = par[4];
    double y3_shift     = par[5];

    const double x_weight = 1.0;
    const double y_weight = 5.0;
    double f = 0;
    for (int i_cluster = 0; i_cluster < numOfUsedHits; i_cluster++)
    {
      double x1_corr = x1_fit[i_cluster]*TMath::Cos(phi_rot_ist1) + y1_fit[i_cluster]*TMath::Sin(phi_rot_ist1) + x1_shift;
      double y1_corr = y1_fit[i_cluster]*TMath::Cos(phi_rot_ist1) - x1_fit[i_cluster]*TMath::Sin(phi_rot_ist1) + y1_shift;

      double x3_corr = x3_fit[i_cluster]*TMath::Cos(phi_rot_ist3) + y3_fit[i_cluster]*TMath::Sin(phi_rot_ist3) + x3_shift;
      double y3_corr = y3_fit[i_cluster]*TMath::Cos(phi_rot_ist3) - x3_fit[i_cluster]*TMath::Sin(phi_rot_ist3) + y3_shift;

      double x2_proj = x3_corr + (x1_corr-x3_corr)*pitchLayer23/(pitchLayer12+pitchLayer23);
      double y2_proj = y3_corr + (y1_corr-y3_corr)*pitchLayer23/(pitchLayer12+pitchLayer23);

      double dx = x2_fit[i_cluster] - x2_proj;
      double dy = y2_fit[i_cluster] - y2_proj;

      f += x_weight*dx*dx + y_weight*dy*dy;
    }
    return f;
  };

  // wrap chi2 funciton in a function object for the fit
  // 3 is the number of fit parameters (size of array par)
  ROOT::Math::Functor fcn(chi2Function,6);
  ROOT::Fit::Fitter fitter;

  // double pStart[3] = {0.01,140,-20};
  double pStart[6] = {phi_ist1_temp,x1_shift_temp,y1_shift_temp,phi_ist3_temp,x3_shift_temp,y3_shift_temp};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("phi_rot_ist1");
  fitter.Config().ParSettings(1).SetName("x1_shift");
  fitter.Config().ParSettings(2).SetName("y1_shift");
  fitter.Config().ParSettings(3).SetName("phi_rot_ist3");
  fitter.Config().ParSettings(4).SetName("x3_shift");
  fitter.Config().ParSettings(5).SetName("y3_shift");

  // do the fit
  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("find alignment","Minuit Fit failed");
    return std::make_tuple(-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1);
  }

  const ROOT::Fit::FitResult & result = fitter.Result();
  // result.Print(std::cout);
  const double *fitpar = result.GetParams();
  double phi_rot_ist1_MF = fitpar[0];
  double x1_shift_MF     = fitpar[1];
  double y1_shift_MF     = fitpar[2];
  double phi_rot_ist3_MF = fitpar[3];
  double x3_shift_MF     = fitpar[4];
  double y3_shift_MF     = fitpar[5];

  // calculate nOffset with current fit parameters
  int nOffset = 0;
  for (int i_cluster = 0; i_cluster < numOfUsedHits; i_cluster++)
  {
    double x1_corr_MF = x1_fit[i_cluster]*TMath::Cos(phi_rot_ist1_MF) + y1_fit[i_cluster]*TMath::Sin(phi_rot_ist1_MF) + x1_shift_MF;
    double y1_corr_MF = y1_fit[i_cluster]*TMath::Cos(phi_rot_ist1_MF) - x1_fit[i_cluster]*TMath::Sin(phi_rot_ist1_MF) + y1_shift_MF;

    double x3_corr_MF = x3_fit[i_cluster]*TMath::Cos(phi_rot_ist3_MF) + y3_fit[i_cluster]*TMath::Sin(phi_rot_ist3_MF) + x3_shift_MF;
    double y3_corr_MF = y3_fit[i_cluster]*TMath::Cos(phi_rot_ist3_MF) - x3_fit[i_cluster]*TMath::Sin(phi_rot_ist3_MF) + y3_shift_MF;

    double x2_proj_MF = x3_corr_MF + (x1_corr_MF-x3_corr_MF)*pitchLayer23/(pitchLayer12+pitchLayer23);
    double y2_proj_MF = y3_corr_MF + (y1_corr_MF-y3_corr_MF)*pitchLayer23/(pitchLayer12+pitchLayer23);

    double dx = x2_fit[i_cluster] - x2_proj_MF;
    double dy = y2_fit[i_cluster] - y2_proj_MF;

    if( !(abs(dx) < xCut_temp && abs(dy) < yCut_temp) )
    {
      nOffset++;
    }
  }

  tPars tParameters = std::make_tuple(phi_rot_ist1_MF, x1_shift_MF, y1_shift_MF, phi_rot_ist3_MF, x3_shift_MF, y3_shift_MF, xCut_temp, yCut_temp, nOffset);

  // cout << "Minuit minimization: phi_rot_ist1 = " << phi_rot_ist1_MF << ", phi_rot_ist3 = " << phi_rot_ist3_MF << ", x_shift = " << x_shift_MF << ", y_shift = " << y_shift_MF << endl;
  cout << "Minuit minimization: phi_rot_ist1 = " << std::get<0>(tParameters) << ", x1_shift = " << (double)std::get<1>(tParameters) << ", y1_shift = " << std::get<2>(tParameters) <<  ", phi_rot_ist3 = " << std::get<3>(tParameters) << ", x3_shift = " << (double)std::get<4>(tParameters) << ", y1_shift = " << std::get<5>(tParameters) << ", nOffset = " << std::get<8>(tParameters) << endl;

  return tParameters;
}

tPars updateFitParameters(tPars fitPars, double xCut, double yCut, int nOffset)
{
  const double phi_ist1_old = std::get<0>(fitPars);
  const double x1_shift_old = std::get<1>(fitPars);
  const double y1_shift_old = std::get<2>(fitPars);
  const double phi_ist3_old = std::get<3>(fitPars);
  const double x3_shift_old = std::get<4>(fitPars);
  const double y3_shift_old = std::get<5>(fitPars);
  const double xCut_old     = std::get<6>(fitPars);
  const double yCut_old     = std::get<7>(fitPars);
  const double nOffset_old  = std::get<8>(fitPars);

  tPars newPars = std::make_tuple(phi_ist1_old, x1_shift_old, y1_shift_old, phi_ist3_old, x3_shift_old, y3_shift_old, xCut, yCut, nOffset);

  return newPars;
}
