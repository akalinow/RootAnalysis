#ifndef RootAnalysis_Tools_H
#define RootAnalysis_Tools_H

#include <string>

namespace HTTAnalysis {

std::string categoryName(unsigned int iCategory);

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
