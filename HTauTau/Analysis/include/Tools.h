#ifndef RootAnalysis_Tools_H
#define RootAnalysis_Tools_H

#include <string>
#include "AnalysisEnums.h"
#include "EventProxyHTT.h"

namespace HTTAnalysis {

  ///Return current luminosity estimate.
  float getLumi();

  ///Return sample cross section in [pb]
  float getCrossSection(const std::string & sampleName);

  ///Return generator weight. Most samples have large values of weights
  ///which are constant up to + or - sign. We normalise those weights to +-1.
  float getGenWeight(const EventProxyHTT & myEventProxy);
  
  ///Get jets separated by deltaR from tau an muon.
  std::vector<HTTParticle> getSeparatedJets(const EventProxyHTT & myEventProxy,
					    const HTTParticle & aLeg1,
					    const HTTParticle & aLeg2, 
					    float deltaR);

  ///Return human readable systematic effect name for given systematic effect number number.
  std::string systEffectName(unsigned int iSystEffect);

  ///Return human readable systematic effect name for given systematic effect number number.
  ///Replace CAT pattern by correct category name
  std::string systEffectName(unsigned int iCategory, unsigned int iSystEffect,
			     const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester);

  ///Return human readable sample name (Data, WJets, etc).
  std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name (Data, WJets, etc).
  ///Make the methos static, so other modules can use it.
  ///Method used when sample coding in TTree is not present.
  ///In this case a ROOT file name is used to decode the sample type.
  std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy);

  ///Return sample name for DY. Name encoded jet bin, and decay mode.
  std::string getDYSampleName(const EventProxyHTT & myEventProxy);

  //Return name sample name suffix for different particles matched to reconstructed tau
  std::string getMatchingName(const EventProxyHTT & myEventProxy); 

  ///Return string encoding di-tau decay mode.
  ///The event can belong to more than one category
  std::vector<std::string> getTauDecayName(int decModeMinus, int decModePlus);

  ///Check if the decMode points to single prong hadronic tau decay
  bool isOneProng(int decMode);

  ///Check if the decMode points to leptonic tau decay
  bool isLepton(int decMode);
}
#endif
