#ifndef RootAnalysis_Tools_H
#define RootAnalysis_Tools_H

#include <string>
#include "AnalysisEnums.h"
#include "EventProxyHTT.h"

namespace TauLFVAnalysis {

  ///Return current luminosity estimate.
  float getLumi();

  ///Return sample cross section in [pb]
  float getCrossSection(const std::string & sampleName);

  ///Return human readable sample name.
  std::string getSampleName(const EventProxyHTT & myEventProxy);

  ///Return human readable sample name. 
  ///Make the methos static, so other modules can use it.
  ///Method used when sample coding in TTree is not present.
  ///In this case a ROOT file name is used to decode the sample type.
  std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy);

}
#endif
