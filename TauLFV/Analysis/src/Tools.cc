#include <string>
#include <vector>
#include "Tools.h"
#include "AnalysisEnums.h"

namespace TauLFVAnalysis {

float getLumi(){
 
        float run2016 = 36; //fb-1 Placeholder. Put correct value here
	return run2016;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
float getCrossSection(const std::string & sampleName){

  double crossSection = 1.0;               
                
  if(sampleName=="DYLowM") {
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
    crossSection = 1.0;
  }
  
  return crossSection;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getSampleName(const EventProxyHTT & myEventProxy){

        return getSampleNameFromFileName(myEventProxy);

}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
std::string getSampleNameFromFileName(const EventProxyHTT & myEventProxy){

        std::string fileName = myEventProxy.getTTree()->GetCurrentFile()->GetName();
        std::string sampleName = "Unknown";

        if(fileName.find("DsToTauTo3Mu")!=std::string::npos) sampleName = "DsToTauTo3Mu";
        else if(fileName.find("BuToTauTo3Mu")!=std::string::npos) sampleName = "BuToTauTo3Mu";
	else if(fileName.find("BdToTauTo3Mu")!=std::string::npos) sampleName = "BdToTauTo3Mu";
	else if(fileName.find("DoubleMuonLowMass")!=std::string::npos) sampleName = "Data";
 
        return sampleName;
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
}
