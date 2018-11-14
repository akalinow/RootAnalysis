#ifndef WarsawAnalysis_HSCPDataFormats_HSCPEvent_h
#define WarsawAnalysis_HSCPDataFormats_HSCPEvent_h

#include "TLorentzVector.h"
#include <map>
#include <vector>


class HSCPEvent;

class HSCPParticle{

  friend class HSCPEvent;
  friend class EventProxyHSCP;

  public:

  HSCPParticle(){ clear();}

  ~HSCPParticle(){}

  void clear();

  const float & getPt() const {return pt;}

  const float & getEta() const {return eta;}

  const float & getPhi() const {return phi;}

 protected:

  float pt, eta, phi;
  
};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HSCPEvent{

  friend class EventProxyHSCP;

 public:

  HSCPEvent();

  ~HSCPEvent(){}

  ///Data member setters.
  void setRun(unsigned int x){runId = x;}

  void setEvent(unsigned long int x){eventId = x;}

  ///Reset class data members
  void clear();

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}

  const HSCPParticle & getCandidate(unsigned int i) const;

 protected:

  ///Event run and number
  int runId;
  int eventId;

  std::vector<HSCPParticle> candidates;
  HSCPParticle dummyCandidate;

};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
#endif
