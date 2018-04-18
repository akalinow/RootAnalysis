
#ifndef FastMTT_FastMTT_H
#define FastMTT_FastMTT_H

#include <string>
#include <vector>

namespace ROOT{
  namespace Math{
    class Minimizer;
    class Functor;
  }
}
class TF1;
class TVector2;

#include "TMatrixD.h"
#include "TLorentzVector.h"
#include "TBenchmark.h"

namespace fastMTT {
  double likelihoodFunc(double *x, double *par);
}

class Likelihood{

public:

  Likelihood();

  ~Likelihood();

  double value(const double *x) const;

  void setLeptonInputs(const TLorentzVector & aLeg1P4,
                       const TLorentzVector & aLeg2P4,
                       int aLeg1DecayType,
                       int aLeg2DecayType);

  void setMETInputs(const TLorentzVector & aMET,
                    const TMatrixD& aCovMET);

  void setParameters(const std::vector<double> & parameters);

  double massLikelihood(const double & m) const;

  double metTF(const TLorentzVector & metP4,
               const TLorentzVector & nuP4,
               const TMatrixD& covMET) const;

private:

  TLorentzVector leg1P4, leg2P4;
  TLorentzVector recoMET;
  TMatrixD covMET;

  double mVis, mVisLeg1, mVisLeg2;

  int leg1DecayType, leg2DecayType;

  std::vector<double> parameters;

};
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
class FastMTT {

 public:

  FastMTT();

  ~FastMTT();

  void initialize();

  ///Run the minimalization procedure for given inputs.
  ///Results are stored in internal variables accesed by
  ///relevant get methods.
  void minimize(const TLorentzVector & aLeg1P4,
                const TLorentzVector & aLeg2P4,
                const TLorentzVector & aMET,
                const TMatrixD & aCovMET,
                int leg1DecayType,
                int leg2DecayType);

  ///Run a scan over x1 and x2 [0,1] rectangle for given inputs.
  ///Results are stored in internal variables accesed by
  ///relevant get methods.
  void scan(const TLorentzVector & aLeg1P4,
                const TLorentzVector & aLeg2P4,
                const TLorentzVector & aMET,
                const TMatrixD & aCovMET,
                int leg1DecayType,
                int leg2DecayType);

  ///Set likelihood shape parameters.
  void setLikelihoodParams(const std::vector<double> & aPars);

  ///Retrieve the four momentum corresponding to the likelihood maximum
  const TLorentzVector & getBestP4() const { return bestP4; }

  ///Retrieve the CPU timing for given methods
  ///Possible values:
  /// scan
  /// minimize
  double getCpuTime(const std::string & method);

  ///Retrieve the CPU timing for given methods
  ///Possible values:
  /// scan
  /// minimize
  double getRealTime(const std::string & method);

 private:

   // Minimizer types and algorithms.
   // minimizerName               minimizerAlgorithm
   // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
   //  Minuit2                     Fumili2
   //  Fumili
   //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
   //                              BFGS2, SteepestDescent
   //  GSLMultiFit
   //   GSLSimAn
   //   Genetic
   std::string minimizerName;
   std::string  minimizerAlgorithm;

   ROOT::Math::Minimizer* minimizer;

   ///Minimum location
   std::vector<double> minimalizationResult;

  ///Dimenstion of minimalization space
  unsigned int nVariables;

  ///Names of variables to be minimized
  std::vector<std::string> varNames;

  ///Values of variables to be minimized
  std::vector<double> variables;

  ///Step sizes for each minimized variable
  std::vector<double> stepSizes;

  ROOT::Math::Functor *likelihoodFunctor;

  Likelihood myLikelihood;

  TLorentzVector bestP4;

  TBenchmark clock;

};

#endif
