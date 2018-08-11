//! ObjectMessenger derived class for Machine Learning analyses.
/*!
  \author Rafal Maselek
  \date May 2018
  
  ObjectMessenger is a class that implements the basic interface for storing objects of different types.
  It allows to store single objects as well as collections (std::vector) of objects.
  The base class has implemented storage, putting and getting for primitive types: bool, unsigned, int, float, double.
  MLObjectMessenger adds implementation for HTTParticle and HTTAnalysis::sysEffects types.
*/
#ifndef RootAnalysis_MLObjectMessenger_H
#define RootAnalysis_MLObjectMessenger_H

#include "ObjectMessenger.h"
// #include "EventProxyBase.h"
#include "EventProxyHTT.h"

class MLObjectMessenger : public ObjectMessenger
{
 public:
  
  //! Constructor with object's name. The name should start with "ML";
  inline MLObjectMessenger(const std::string& aName) : ObjectMessenger(aName)
  {
    // std::string name = aName;
    // myName_ = name;
  }

  //! Destructor.
  virtual ~MLObjectMessenger(){;}
  
  //! Getter for the name.
  // virtual inline const std::string & name(){return myName_;};

  //! Method for putting an object into messenger
  template<typename T> inline void putObject(const T* obj, const std::string& flavour)
  {
    try {
      ObjectMessenger::putObject(obj, flavour);
    }
    catch(const std::exception& e) {
        std::throw_with_nested(std::runtime_error("[ERROR] PUTTING OBJECT INTO MLObjectMessenger FAILED!"));
    }
  };

  template<typename T> inline void putObject(const T obj, const std::string& flavour)
    {putObject(ptr(obj),flavour);};


  //! Method for putting an object into vector inside a messenger
  template<typename T> inline void putObjectVector(const T* obj, const std::string& flavour)
  {
    try {
      ObjectMessenger::putObjectVector(obj, flavour);
    }
    catch(const std::exception& e) {
        std::throw_with_nested(std::runtime_error("[ERROR] PUTTING VECTOR INTO MLObjectMessenger FAILED!"));
    }
  };
  template<typename T> inline void putObjectVector(const T obj, const std::string& flavour)
    {putObjectVector(ptr(obj),flavour);};


  //! Method for overwriting existing objects
  template<typename T> inline bool overwriteObject(const T* obj, const std::string& flavour)
  {
    try {
      return ObjectMessenger::overwriteObject(obj, flavour);
    }
    catch(const std::exception& e) {
        std::throw_with_nested(std::runtime_error("[ERROR] OVERWRITING OBJECT IN MLObjectMessenger FAILED!"));
    }
  };
  template<typename T> inline bool overwriteObject(const T obj, const std::string& flavour)
    {return overwriteObject(ptr(obj),flavour);};


  //! Method for accessing an object from the messenger
  template<typename T> inline const T* getObject(T* obj, const std::string& flavour) const
  {
    try {
      return ObjectMessenger::getObject(obj, flavour);
    }
    catch(const std::exception& e) {
        std::throw_with_nested(std::runtime_error("[ERROR] GETTING OBJECT FROM MLObjectMessenger FAILED!"));
    }
  };
  template<typename T> inline const T* getObject(T obj, const std::string& flavour) const
    {getObject(ptr(obj),flavour);};


  //! Method for accessing an object from a vector in the messenger
  template<typename T> inline const std::vector<const T*>* getObjectVector(T* obj, const std::string& flavour) const 
  {
    try {
      return ObjectMessenger::getObject(obj, flavour);
    }
    catch(const std::exception& e) {
        std::throw_with_nested(std::runtime_error("[ERROR] GETTING VECTOR FROM MLObjectMessenger FAILED!"));
    }
  };
  template<typename T> inline const std::vector<const T*>* getObjectVector(T obj, const std::string& flavour) const
    {return getObjectVector(ptr(obj), flavour);};


  //! Method cleaning the messenger contents
  inline void clear()
  {
    ObjectMessenger::clear();
    mapHTTParticle_.clear();
    mapSystEffect_.clear();
    mapHTTParticleVector_.clear();
  }

 private:

  std::string myName_; /*!< Name of MLObjectMesssenger. It should start with "ML" in order to work properly. */

  ///Containers for carried objects to be defined
  std::map<std::string, const HTTAnalysis::sysEffects*> mapSystEffect_;
  std::map<std::string, const HTTParticle *> mapHTTParticle_;
  std::map<std::string, std::vector<const HTTParticle *>> mapHTTParticleVector_;
 
};

/*******************************************************************************************************************/
/*** HTTPARTICLE *** HTTPARTICLE *** HTTPARTICLE *** HTTPARTICLE *** HTTPARTICLE *** HTTPARTICLE *** HTTPARTICLE ***/
/*******************************************************************************************************************/
//! Putter for particle objects. Used to put legs and jets into shared containers.
template<> inline void MLObjectMessenger::putObjectVector(const HTTParticle* obj, const std::string & flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapHTTParticleVector_.find(name);
  if (iter == mapHTTParticleVector_.end() )
  {
    // key does not exists, create the vector
     std::vector<const HTTParticle*> v;
     mapHTTParticleVector_[name] = v;
  }
  mapHTTParticleVector_[name].push_back(obj);

}

//! Getter for vector of particle objects. Used to access legs and jets.
template<> inline const std::vector<const HTTParticle*>* MLObjectMessenger::getObjectVector(HTTParticle* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapHTTParticleVector_.at(name);
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapHTTParticleVector_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLObjectMessenger::getObject!"));
  } 
}

//! Putter for particle objects. Doesn't copy the object.
template<> inline void MLObjectMessenger::putObject(const HTTParticle* obj, const std::string & flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapHTTParticle_.find(name);
  if (iter == mapHTTParticle_.end())
  {
    mapHTTParticle_[name]=obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" !"));
}

//! Getter for particle objects.
template<> inline const HTTParticle* MLObjectMessenger::getObject(HTTParticle* obj, const std::string & flavour) const
{
  try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return mapHTTParticle_.at(name);
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapHTTParticle_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLObjectMessenger::getObject!"));
  }
}

/***************************************************************************************************************************/
/*** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS *** SYSEFFECTS ***/
/***************************************************************************************************************************/
//! Putter for systematic errors objects. Doesn't copy the object.
template<> inline void MLObjectMessenger::putObject(const HTTAnalysis::sysEffects* obj, const std::string & flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapSystEffect_.find(name);
  if (iter == mapSystEffect_.end())
  {
    mapSystEffect_[name]=obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE MLObjectMessenger::overwriteObject INSTEAD."));
}

//! Getter for systematic errors objects.
template<> inline const HTTAnalysis::sysEffects* MLObjectMessenger::getObject(HTTAnalysis::sysEffects* obj, const std::string & flavour) const
{
  try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return mapSystEffect_.at(name);
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapSystEffect_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN MLObjectMessenger::getObject!"));
  }
}

#endif

