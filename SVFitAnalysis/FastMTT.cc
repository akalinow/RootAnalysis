#include "FastMTT.h"

#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TMatrixD.h"

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::Likelihood(){

  covMET.ResizeTo(2,2);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
Likelihood::~Likelihood(){}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setLeptonInputs(const TLorentzVector & aLeg1P4,
                                 const TLorentzVector & aLeg2P4,
                                 int aLeg1DecayType,
                                 int aLeg2DecayType){
  leg1P4 = aLeg1P4;
  leg2P4 = aLeg2P4;

  mVis = (leg1P4 + leg2P4).M();
  mVisLeg1 = leg1P4.M();
  mVisLeg2 = leg2P4.M();
  mVisLeg2 = 0.3;///TEST

  leg1DecayType = aLeg1DecayType;
  leg2DecayType = aLeg2DecayType;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setMETInputs(const TLorentzVector & aMET,
                              const TMatrixD& aCovMET){
  recoMET = aMET;
  covMET = aCovMET;

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Likelihood::setParameters(const std::vector<double> & aPars){

  parameters = aPars;

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double Likelihood::massLikelihood(const double & m) const{

  double coeff1 = parameters[0];
  double coeff2 = parameters[1];
  double mShift = m*coeff2;

  double mTau = 1.77685;

  double x1Min = std::min(1.0, std::pow(mVisLeg1/mTau,2));
  double x2Min = std::max(std::pow(mVisLeg2/mTau,2), std::pow(mVis/mShift,2));
  double x2Max = std::min(1.0, std::pow(mVis/mShift,2)/x1Min);

  double value = 2.0*std::pow(mVis,2)*std::pow(mShift,-coeff1)*(log(x2Max)-log(x2Min) + std::pow(mVis/mShift,2)*(1 - std::pow(x2Min,-1)));
  if(mShift<mVis) return 0.0;

  return value*1E12;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////x
double Likelihood::metTF(const TLorentzVector & metP4,
                         const TLorentzVector & nuP4,
                         const TMatrixD& covMET) const{

  double  aMETx = metP4.X();
  double  aMETy = metP4.Y();

  double invCovMETxx = covMET(1,1);
  double invCovMETxy = -covMET(0,1);
  double invCovMETyx = -covMET(1,0);
  double invCovMETyy = covMET(0,0);
  double covDet = invCovMETxx*invCovMETyy - invCovMETxy*invCovMETyx;

  if( std::abs(covDet)<1E-10){
    std::cerr << "Error: Cannot invert MET covariance Matrix (det=0) !!"
	      <<"METx: "<<aMETy<<" METy: "<<aMETy
	      << std::endl;
    return 0;
  }
  double const_MET = 1./(2.*TMath::Pi()*TMath::Sqrt(covDet));

  double residualX = aMETx - (nuP4.X());
  double residualY = aMETy - (nuP4.Y());

  double pull2 = residualX*(invCovMETxx*residualX + invCovMETxy*residualY) +
    residualY*(invCovMETyx*residualX + invCovMETyy*residualY);
  pull2/=covDet;

  return const_MET*TMath::Exp(-0.5*pull2);
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
double Likelihood::value(const double *x) const{

  TLorentzVector testP4 = leg1P4*(1/x[0]) + leg2P4*(1/x[1]);
  TLorentzVector testMET = leg1P4*(1.0/x[0] - 1) + leg2P4*(1.0/x[1] - 1);

  double value = metTF(recoMET, testMET, covMET);
  value *= massLikelihood(testP4.M());

  return -value;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
FastMTT::FastMTT(){

  minimizerName = "Minuit2";
  minimizerAlgorithm = "Migrad";
  initialize();
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
FastMTT::~FastMTT(){

  delete minimizer;
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::initialize(){

  minimizer = ROOT::Math::Factory::CreateMinimizer(minimizerName, minimizerAlgorithm);
  minimizer->SetMaxFunctionCalls(10000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);

  std::vector<std::string> varNames = {"x1", "x2"};
  nVariables = varNames.size();
  std::vector<double> initialValues(nVariables,0.5);
  std::vector<double> stepSizes(nVariables, 0.01);

  for(unsigned int iVar=0; iVar<nVariables; ++iVar){
    minimizer->SetVariable(iVar, varNames[iVar].c_str(), initialValues[iVar], stepSizes[iVar]);
  }

  std::vector<double> shapeParams = {6, 1.0/1.15};
  setLikelihoodParams(shapeParams);
  likelihoodFunctor = new ROOT::Math::Functor(&myLikelihood, &Likelihood::value, nVariables);
  minimizer->SetFunction(*likelihoodFunctor);
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::setLikelihoodParams(const std::vector<double> & aPars){

   myLikelihood.setParameters(aPars);

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::minimize(const TLorentzVector & aLeg1P4,
                       const TLorentzVector & aLeg2P4,
                       const TLorentzVector & aMET,
                       const TMatrixD & aCovMET,
                       int aLeg1DecayType,
                       int aLeg2DecayType){

  clock.Reset();
  clock.Start("minimize");

  for(int iTest = 0; iTest<1; ++iTest){

  myLikelihood.setLeptonInputs(aLeg1P4, aLeg2P4, aLeg1DecayType, aLeg2DecayType);
  myLikelihood.setMETInputs(aMET, aCovMET);

  minimizer->SetVariableLimits(0, 0.0 ,1.0);
  minimizer->SetVariableLimits(1, 0.0 ,1.0);
  minimizer->Minimize();

  const double *theMinimum = minimizer->X();

  bestP4 = aLeg1P4*(1.0/theMinimum[0]) + aLeg2P4*(1.0/theMinimum[1]);

   if(false){
  std::cout<<" minimizer "
	   <<" nCalls: "<<minimizer->NCalls()
    	   <<" nIterations: "<<minimizer->NIterations()
           <<" x1Max: "<<theMinimum[0]
           <<" x2Max: "<<theMinimum[1]
	   <<" x3Max: "<<theMinimum[2]
           <<" max LLH: "<<minimizer->MinValue()
	   <<" m: "<<bestP4.M()
           <<std::endl;
  }
}
  clock.Stop("minimize");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void FastMTT::scan(const TLorentzVector & aLeg1P4,
                       const TLorentzVector & aLeg2P4,
                       const TLorentzVector & aMET,
                       const TMatrixD & aCovMET,
                       int aLeg1DecayType,
                       int aLeg2DecayType){

  clock.Reset();
  clock.Start("scan");

  for(int iTest = 0; iTest<1; ++iTest){

  myLikelihood.setLeptonInputs(aLeg1P4, aLeg2P4, aLeg1DecayType, aLeg2DecayType);
  myLikelihood.setMETInputs(aMET, aCovMET);

  double lh = 0.0;
  double maxLH = -99.0;

  double x[2] = {0.5, 0.5};
  double theMinimum[2] = {0.5, 0.5};
  int nGridPoints = 100;
  int nCalls = 0;
    for(int iX2 = 1; iX2<nGridPoints;++iX2){
      x[1] = (double)iX2/nGridPoints;
      for(int iX1 = 1; iX1<nGridPoints;++iX1){
        x[0] = (double)iX1/nGridPoints;

         lh = - myLikelihood.value(x);
        ++nCalls;

    if(lh>maxLH){
      maxLH = lh;
      theMinimum[0] = x[0];
      theMinimum[1] = x[1];
    }
  }
}

  bestP4 = aLeg1P4*(1.0/theMinimum[0]) + aLeg2P4*(1.0/theMinimum[1]);
}
  clock.Stop("scan");
}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getCpuTime(const std::string & method){

  return clock.GetCpuTime(method.c_str());

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
double FastMTT::getRealTime(const std::string & method){

  return clock.GetRealTime(method.c_str());

}
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
