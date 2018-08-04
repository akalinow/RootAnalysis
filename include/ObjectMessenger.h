//! Base class for storing objects passed along Analyzers.
/*!
  \author Rafal Maselek
  \date May 2018

  ObjectMessenger is a class that implements the basic interface for storing objects of different types.
  It allows to store single objects as well as collections (std::vector) of objects.
  The base class has implemented storage, putting and getting for primitive types: bool, unsigned, int, float, double.
*/

#ifndef RootAnalysis_ObjectMessenger_H
#define RootAnalysis_ObjectMessenger_H

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <exception>
#include <typeinfo>



class ObjectMessenger{

 public:

  //! Default no-parameter constructor
  ObjectMessenger(){}

  //! Constructor with a name
  ObjectMessenger(const std::string & aName);

  //! Destructor
  virtual ~ObjectMessenger(){;}
  
  //! Getter for name
  inline const std::string& name(){return myName_;};

  //! Clever way to extend functionality from pointers to references as well 
  template<typename T> inline const T* ptr(T& obj) {return const_cast<const T*>(&obj); } //turn reference into pointer!
  template<typename T> inline const T* ptr(T * obj) {return const_cast<const T*>(obj); } //obj is already pointer, return it!

  //! Method for putting an object into messenger
  template<typename T> inline void putObject(const T* obj, const std::string& flavour)
    {throw std::logic_error(std::string("[ERROR] UNDEFINED ACTION FOR: ObjectMessenger::putObject, for type: ")+typeid(T).name());};
  template<typename T> inline void putObject(const T obj, const std::string& flavour)
    {putObject(ptr(obj),flavour);};

  //! Method for putting an object into vector inside a messenger
  template<typename T> inline void putObjectVector(const T* obj, const std::string& flavour)
    {throw std::logic_error(std::string("[ERROR] UNDEFINED ACTION FOR: ObjectMessenger::putObjectVector, for type: ")+typeid(T).name());};
  template<typename T> inline void putObjectVector(const T obj, const std::string& flavour)
    {putObjectVector(ptr(obj),flavour);};

  //! Method for overwriting existing objects
  template<typename T> inline bool overwriteObject(const T* obj, const std::string& flavour)
    {throw std::logic_error(std::string("[ERROR] UNDEFINED ACTION FOR: ObjectMessenger::overwriteObject, for type: ")+typeid(T).name());};
  template<typename T> inline bool overwriteObject(const T obj, const std::string& flavour)
    {return overwriteObject(ptr(obj),flavour);};

  //! Method for accessing an object from the messenger
  template<typename T> inline const T* getObject(T* obj, const std::string& flavour) const
    {throw std::logic_error(std::string("[ERROR] UNDEFINED ACTION FOR: ObjectMessenger::getObject, for type: ")+typeid(T).name());};
  template<typename T> inline const T* getObject(T obj, const std::string& flavour) const
    {getObject(ptr(obj),flavour);};

  //! Method for accessing an object from a vector in the messenger
  template<typename T> inline const std::vector<const T*>* getObjectVector(T* obj, const std::string& flavour) const
    {throw std::logic_error(std::string("[ERROR] UNDEFINED ACTION FOR: ObjectMessenger::getObjectVector, for type: ")+typeid(T).name());};
  template<typename T> inline const std::vector<const T*>* getObjectVector(T obj, const std::string& flavour) const
    {return getObjectVector(ptr(obj), flavour);};

  //! Method cleaning the messenger contents
  void clear();

 private:

  std::string myName_; /*!< Name of the ObjectMessenger */

  //! Containers for carried objects to be defined in class
  /*! [!IMPORTANT] In the specialization below, varaibles of primitive types will be copied to maps instead of copying pointers. */
  std::map<std::string, bool> mapBool_;
  std::map<std::string, unsigned> mapUnsigned_;
  std::map<std::string, int> mapInt_;
  std::map<std::string, float> mapFloat_;
  std::map<std::string, double> mapDouble_;
 
};

/************************************************************************************************************************/
/*** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL *** BOOL ***/
/************************************************************************************************************************/
//! Putter for bools. Copies variables.
template<> inline void ObjectMessenger::putObject(const bool* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapBool_.find(name);
  if (iter == mapBool_.end())
  {
    mapBool_[name]=*obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE ObjectMessenger::overwriteObject INSTEAD."));
}

//! Function to overwrite a bool parameter previously put. 
template<> inline bool ObjectMessenger::overwriteObject(const bool* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapBool_.find(name);
  if (iter == mapBool_.end())
  {
    return false; /*!< key not found */
  }
  else
  { 
    mapBool_[name]=*obj;
    return true;
  }
}

//! Getter for bool objects.
template<> inline const bool* ObjectMessenger::getObject(bool* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapBool_.at(name);  
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapBool_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN ObjectMessenger::getObject!"));
  }
}

/************************************************************************************************************************/
/*** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED *** UNSIGNED ***/
/************************************************************************************************************************/
//! Putter for unsigned. Copies variables.
template<> inline void ObjectMessenger::putObject(const unsigned* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapUnsigned_.find(name);
  if (iter == mapUnsigned_.end())
  {
    mapUnsigned_[name]=*obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE ObjectMessenger::overwriteObject INSTEAD."));
}

//! Function to overwrite an unsigned parameter previously put. 
template<> inline bool ObjectMessenger::overwriteObject(const unsigned* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapUnsigned_.find(name);
  if (iter == mapUnsigned_.end())
  {
    return false; /*!< key not found */
  }
  else
  { 
    mapUnsigned_[name]=*obj;
    return true;
  }
}

//! Getter for unsigned objects.
template<> inline const unsigned* ObjectMessenger::getObject(unsigned* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapUnsigned_.at(name);  
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapUnsigned_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN ObjectMessenger::getObject!"));
  }  
}

/************************************************************************************************************************/
/*** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT *** INT ***/
/************************************************************************************************************************/
//! Putter for ints. Copies variables.
template<> inline void ObjectMessenger::putObject(const int* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapInt_.find(name);
  if (iter == mapInt_.end())
  {
    mapInt_[name]=*obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE ObjectMessenger::overwriteObject INSTEAD."));
}

//! Function to overwrite an int parameter previously put. 
template<> inline bool ObjectMessenger::overwriteObject(const int* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapInt_.find(name);
  if (iter == mapInt_.end())
  {
    return false; /*!< key not found */
  }
  else
  { 
    mapInt_[name]=*obj;
    return true;
  }
}

//! Getter for int objects.
template<> inline const int* ObjectMessenger::getObject(int* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapInt_.at(name);  
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapInt_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN ObjectMessenger::getObject!"));
  }  
}

/************************************************************************************************************************/
/*** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT *** FLOAT ***/
/************************************************************************************************************************/
//! Putter for floats. Copies variables.
template<> inline void ObjectMessenger::putObject(const float* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapFloat_.find(name);
  if (iter == mapFloat_.end())
  {
    mapFloat_[name]=*obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE ObjectMessenger::overwriteObject INSTEAD."));
}

//! Function to overwrite a float parameter previously put. 
template<> inline bool ObjectMessenger::overwriteObject(const float* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapFloat_.find(name);
  if (iter == mapFloat_.end())
  {
    return false; /*!< key not found */
  }
  else
  { 
    mapFloat_[name]=*obj;
    return true;
  }
}

//! Getter for float objects.
template<> inline const float* ObjectMessenger::getObject(float* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapFloat_.at(name);  
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapFloat_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN  ObjectMessenger::getObject!"));
  }  
}

/************************************************************************************************************************/
/*** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE *** DOUBLE ***/
/************************************************************************************************************************/
//! Putter for doubles. Copies variables.
template<> inline void ObjectMessenger::putObject(const double* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapDouble_.find(name);
  if (iter == mapDouble_.end())
  {
    mapDouble_[name]=*obj;
  }
  else throw std::invalid_argument(std::string("[ERROR] TRYING TO OVERWRITE: \"")+name+std::string("\" USE ObjectMessenger::overwriteObject INSTEAD."));
}

//! Function to overwrite a double parameter previously put. 
template<> inline bool ObjectMessenger::overwriteObject(const double* obj, const std::string& flavour)
{
  std::string name = flavour;
  std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
  auto iter = mapDouble_.find(name);
  if (iter == mapDouble_.end())
  {
    return false; /*!< key not found */
  }
  else
  { 
    mapDouble_[name]=*obj;
    return true;
  }
}

//! Getter for double objects.
template<> inline const double* ObjectMessenger::getObject(double* obj, const std::string & flavour) const
{
 try 
  {
    std::string name = flavour;
    std::for_each(name.begin(), name.end(), [](char & c){c = ::tolower(c);});
    return &mapDouble_.at(name);  
  } 
  catch(const std::out_of_range& e)
  {
    std::throw_with_nested(std::runtime_error(std::string("[ERROR] NO ELEMENT \"")+flavour+std::string("\" IN mapDouble_!")));
  }
  catch(const std::exception& e)
  {
     std::throw_with_nested(std::runtime_error("[ERROR] UNKNOWN ERROR IN ObjectMessenger::getObject!"));
  } 
}

#endif

