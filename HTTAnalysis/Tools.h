#ifndef RootAnalysis_Tools_H
#define RootAnalysis_Tools_H

#include <string>

namespace HTTAnalysis {

///Return current luminosity estimate.
float getLumi();

///Return sample cross section in [pb]
float getCrossSection(const std::string & sampleName);

///Return human readable category name for given category number.
std::string categoryName(unsigned int iCategory);

///Return human readable systematic effect name for given systematic effect number number.
std::string systEffectName(unsigned int iSystEffect);

///Return string encoding di-tau decay mode.
///The event can belong to more than one category
std::vector<std::string> getTauDecayName(int decModeMinus, int decModePlus);

///Check if the decMode points to single prong hadronic tau decay
bool isOneProng(int decMode);

///Check if the decMode points to leptonic tau decay
bool isLepton(int decMode);

}
#endif
