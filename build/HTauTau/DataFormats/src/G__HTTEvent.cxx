// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__HTTEvent
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "HTTEvent.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *HTTEvent_Dictionary();
   static void HTTEvent_TClassManip(TClass*);
   static void *new_HTTEvent(void *p = nullptr);
   static void *newArray_HTTEvent(Long_t size, void *p);
   static void delete_HTTEvent(void *p);
   static void deleteArray_HTTEvent(void *p);
   static void destruct_HTTEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HTTEvent*)
   {
      ::HTTEvent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HTTEvent));
      static ::ROOT::TGenericClassInfo 
         instance("HTTEvent", "HTTEvent.h", 30,
                  typeid(::HTTEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HTTEvent_Dictionary, isa_proxy, 4,
                  sizeof(::HTTEvent) );
      instance.SetNew(&new_HTTEvent);
      instance.SetNewArray(&newArray_HTTEvent);
      instance.SetDelete(&delete_HTTEvent);
      instance.SetDeleteArray(&deleteArray_HTTEvent);
      instance.SetDestructor(&destruct_HTTEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HTTEvent*)
   {
      return GenerateInitInstanceLocal((::HTTEvent*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HTTEvent*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HTTEvent_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HTTEvent*)nullptr)->GetClass();
      HTTEvent_TClassManip(theClass);
   return theClass;
   }

   static void HTTEvent_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HTTParticle_Dictionary();
   static void HTTParticle_TClassManip(TClass*);
   static void *new_HTTParticle(void *p = nullptr);
   static void *newArray_HTTParticle(Long_t size, void *p);
   static void delete_HTTParticle(void *p);
   static void deleteArray_HTTParticle(void *p);
   static void destruct_HTTParticle(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HTTParticle*)
   {
      ::HTTParticle *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HTTParticle));
      static ::ROOT::TGenericClassInfo 
         instance("HTTParticle", "HTTEvent.h", 204,
                  typeid(::HTTParticle), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HTTParticle_Dictionary, isa_proxy, 4,
                  sizeof(::HTTParticle) );
      instance.SetNew(&new_HTTParticle);
      instance.SetNewArray(&newArray_HTTParticle);
      instance.SetDelete(&delete_HTTParticle);
      instance.SetDeleteArray(&deleteArray_HTTParticle);
      instance.SetDestructor(&destruct_HTTParticle);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HTTParticle*)
   {
      return GenerateInitInstanceLocal((::HTTParticle*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HTTParticle*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HTTParticle_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HTTParticle*)nullptr)->GetClass();
      HTTParticle_TClassManip(theClass);
   return theClass;
   }

   static void HTTParticle_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HTTPair_Dictionary();
   static void HTTPair_TClassManip(TClass*);
   static void *new_HTTPair(void *p = nullptr);
   static void *newArray_HTTPair(Long_t size, void *p);
   static void delete_HTTPair(void *p);
   static void deleteArray_HTTPair(void *p);
   static void destruct_HTTPair(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HTTPair*)
   {
      ::HTTPair *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HTTPair));
      static ::ROOT::TGenericClassInfo 
         instance("HTTPair", "HTTEvent.h", 300,
                  typeid(::HTTPair), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HTTPair_Dictionary, isa_proxy, 4,
                  sizeof(::HTTPair) );
      instance.SetNew(&new_HTTPair);
      instance.SetNewArray(&newArray_HTTPair);
      instance.SetDelete(&delete_HTTPair);
      instance.SetDeleteArray(&deleteArray_HTTPair);
      instance.SetDestructor(&destruct_HTTPair);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HTTPair*)
   {
      return GenerateInitInstanceLocal((::HTTPair*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HTTPair*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HTTPair_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HTTPair*)nullptr)->GetClass();
      HTTPair_TClassManip(theClass);
   return theClass;
   }

   static void HTTPair_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_HTTEvent(void *p) {
      return  p ? new(p) ::HTTEvent : new ::HTTEvent;
   }
   static void *newArray_HTTEvent(Long_t nElements, void *p) {
      return p ? new(p) ::HTTEvent[nElements] : new ::HTTEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_HTTEvent(void *p) {
      delete ((::HTTEvent*)p);
   }
   static void deleteArray_HTTEvent(void *p) {
      delete [] ((::HTTEvent*)p);
   }
   static void destruct_HTTEvent(void *p) {
      typedef ::HTTEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HTTEvent

namespace ROOT {
   // Wrappers around operator new
   static void *new_HTTParticle(void *p) {
      return  p ? new(p) ::HTTParticle : new ::HTTParticle;
   }
   static void *newArray_HTTParticle(Long_t nElements, void *p) {
      return p ? new(p) ::HTTParticle[nElements] : new ::HTTParticle[nElements];
   }
   // Wrapper around operator delete
   static void delete_HTTParticle(void *p) {
      delete ((::HTTParticle*)p);
   }
   static void deleteArray_HTTParticle(void *p) {
      delete [] ((::HTTParticle*)p);
   }
   static void destruct_HTTParticle(void *p) {
      typedef ::HTTParticle current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HTTParticle

namespace ROOT {
   // Wrappers around operator new
   static void *new_HTTPair(void *p) {
      return  p ? new(p) ::HTTPair : new ::HTTPair;
   }
   static void *newArray_HTTPair(Long_t nElements, void *p) {
      return p ? new(p) ::HTTPair[nElements] : new ::HTTPair[nElements];
   }
   // Wrapper around operator delete
   static void delete_HTTPair(void *p) {
      delete ((::HTTPair*)p);
   }
   static void deleteArray_HTTPair(void *p) {
      delete [] ((::HTTPair*)p);
   }
   static void destruct_HTTPair(void *p) {
      typedef ::HTTPair current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HTTPair

namespace ROOT {
   static TClass *vectorlElongsPdoublegR_Dictionary();
   static void vectorlElongsPdoublegR_TClassManip(TClass*);
   static void *new_vectorlElongsPdoublegR(void *p = nullptr);
   static void *newArray_vectorlElongsPdoublegR(Long_t size, void *p);
   static void delete_vectorlElongsPdoublegR(void *p);
   static void deleteArray_vectorlElongsPdoublegR(void *p);
   static void destruct_vectorlElongsPdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<long double>*)
   {
      vector<long double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<long double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<long double>", -2, "vector", 389,
                  typeid(vector<long double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlElongsPdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<long double>) );
      instance.SetNew(&new_vectorlElongsPdoublegR);
      instance.SetNewArray(&newArray_vectorlElongsPdoublegR);
      instance.SetDelete(&delete_vectorlElongsPdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlElongsPdoublegR);
      instance.SetDestructor(&destruct_vectorlElongsPdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<long double> >()));

      ::ROOT::AddClassAlternate("vector<long double>","std::vector<long double, std::allocator<long double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<long double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlElongsPdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<long double>*)nullptr)->GetClass();
      vectorlElongsPdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlElongsPdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlElongsPdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long double> : new vector<long double>;
   }
   static void *newArray_vectorlElongsPdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<long double>[nElements] : new vector<long double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlElongsPdoublegR(void *p) {
      delete ((vector<long double>*)p);
   }
   static void deleteArray_vectorlElongsPdoublegR(void *p) {
      delete [] ((vector<long double>*)p);
   }
   static void destruct_vectorlElongsPdoublegR(void *p) {
      typedef vector<long double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<long double>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = nullptr);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 389,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)nullptr)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETVector2gR_Dictionary();
   static void vectorlETVector2gR_TClassManip(TClass*);
   static void *new_vectorlETVector2gR(void *p = nullptr);
   static void *newArray_vectorlETVector2gR(Long_t size, void *p);
   static void delete_vectorlETVector2gR(void *p);
   static void deleteArray_vectorlETVector2gR(void *p);
   static void destruct_vectorlETVector2gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector2>*)
   {
      vector<TVector2> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector2>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector2>", -2, "vector", 389,
                  typeid(vector<TVector2>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector2gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector2>) );
      instance.SetNew(&new_vectorlETVector2gR);
      instance.SetNewArray(&newArray_vectorlETVector2gR);
      instance.SetDelete(&delete_vectorlETVector2gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector2gR);
      instance.SetDestructor(&destruct_vectorlETVector2gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector2> >()));

      ::ROOT::AddClassAlternate("vector<TVector2>","std::vector<TVector2, std::allocator<TVector2> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TVector2>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector2gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TVector2>*)nullptr)->GetClass();
      vectorlETVector2gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector2gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector2gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector2> : new vector<TVector2>;
   }
   static void *newArray_vectorlETVector2gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector2>[nElements] : new vector<TVector2>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector2gR(void *p) {
      delete ((vector<TVector2>*)p);
   }
   static void deleteArray_vectorlETVector2gR(void *p) {
      delete [] ((vector<TVector2>*)p);
   }
   static void destruct_vectorlETVector2gR(void *p) {
      typedef vector<TVector2> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TVector2>

namespace ROOT {
   static TClass *vectorlETLorentzVectorgR_Dictionary();
   static void vectorlETLorentzVectorgR_TClassManip(TClass*);
   static void *new_vectorlETLorentzVectorgR(void *p = nullptr);
   static void *newArray_vectorlETLorentzVectorgR(Long_t size, void *p);
   static void delete_vectorlETLorentzVectorgR(void *p);
   static void deleteArray_vectorlETLorentzVectorgR(void *p);
   static void destruct_vectorlETLorentzVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TLorentzVector>*)
   {
      vector<TLorentzVector> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TLorentzVector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TLorentzVector>", -2, "vector", 389,
                  typeid(vector<TLorentzVector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETLorentzVectorgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TLorentzVector>) );
      instance.SetNew(&new_vectorlETLorentzVectorgR);
      instance.SetNewArray(&newArray_vectorlETLorentzVectorgR);
      instance.SetDelete(&delete_vectorlETLorentzVectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlETLorentzVectorgR);
      instance.SetDestructor(&destruct_vectorlETLorentzVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TLorentzVector> >()));

      ::ROOT::AddClassAlternate("vector<TLorentzVector>","std::vector<TLorentzVector, std::allocator<TLorentzVector> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TLorentzVector>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETLorentzVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TLorentzVector>*)nullptr)->GetClass();
      vectorlETLorentzVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETLorentzVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETLorentzVectorgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector> : new vector<TLorentzVector>;
   }
   static void *newArray_vectorlETLorentzVectorgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TLorentzVector>[nElements] : new vector<TLorentzVector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETLorentzVectorgR(void *p) {
      delete ((vector<TLorentzVector>*)p);
   }
   static void deleteArray_vectorlETLorentzVectorgR(void *p) {
      delete [] ((vector<TLorentzVector>*)p);
   }
   static void destruct_vectorlETLorentzVectorgR(void *p) {
      typedef vector<TLorentzVector> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TLorentzVector>

namespace ROOT {
   static TClass *vectorlEHTTParticlegR_Dictionary();
   static void vectorlEHTTParticlegR_TClassManip(TClass*);
   static void *new_vectorlEHTTParticlegR(void *p = nullptr);
   static void *newArray_vectorlEHTTParticlegR(Long_t size, void *p);
   static void delete_vectorlEHTTParticlegR(void *p);
   static void deleteArray_vectorlEHTTParticlegR(void *p);
   static void destruct_vectorlEHTTParticlegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HTTParticle>*)
   {
      vector<HTTParticle> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HTTParticle>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HTTParticle>", -2, "vector", 389,
                  typeid(vector<HTTParticle>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHTTParticlegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HTTParticle>) );
      instance.SetNew(&new_vectorlEHTTParticlegR);
      instance.SetNewArray(&newArray_vectorlEHTTParticlegR);
      instance.SetDelete(&delete_vectorlEHTTParticlegR);
      instance.SetDeleteArray(&deleteArray_vectorlEHTTParticlegR);
      instance.SetDestructor(&destruct_vectorlEHTTParticlegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HTTParticle> >()));

      ::ROOT::AddClassAlternate("vector<HTTParticle>","std::vector<HTTParticle, std::allocator<HTTParticle> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HTTParticle>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHTTParticlegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HTTParticle>*)nullptr)->GetClass();
      vectorlEHTTParticlegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHTTParticlegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHTTParticlegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HTTParticle> : new vector<HTTParticle>;
   }
   static void *newArray_vectorlEHTTParticlegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HTTParticle>[nElements] : new vector<HTTParticle>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHTTParticlegR(void *p) {
      delete ((vector<HTTParticle>*)p);
   }
   static void deleteArray_vectorlEHTTParticlegR(void *p) {
      delete [] ((vector<HTTParticle>*)p);
   }
   static void destruct_vectorlEHTTParticlegR(void *p) {
      typedef vector<HTTParticle> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HTTParticle>

namespace ROOT {
   static TClass *vectorlEHTTPairgR_Dictionary();
   static void vectorlEHTTPairgR_TClassManip(TClass*);
   static void *new_vectorlEHTTPairgR(void *p = nullptr);
   static void *newArray_vectorlEHTTPairgR(Long_t size, void *p);
   static void delete_vectorlEHTTPairgR(void *p);
   static void deleteArray_vectorlEHTTPairgR(void *p);
   static void destruct_vectorlEHTTPairgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<HTTPair>*)
   {
      vector<HTTPair> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<HTTPair>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<HTTPair>", -2, "vector", 389,
                  typeid(vector<HTTPair>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHTTPairgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<HTTPair>) );
      instance.SetNew(&new_vectorlEHTTPairgR);
      instance.SetNewArray(&newArray_vectorlEHTTPairgR);
      instance.SetDelete(&delete_vectorlEHTTPairgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHTTPairgR);
      instance.SetDestructor(&destruct_vectorlEHTTPairgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<HTTPair> >()));

      ::ROOT::AddClassAlternate("vector<HTTPair>","std::vector<HTTPair, std::allocator<HTTPair> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<HTTPair>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHTTPairgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<HTTPair>*)nullptr)->GetClass();
      vectorlEHTTPairgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHTTPairgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHTTPairgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HTTPair> : new vector<HTTPair>;
   }
   static void *newArray_vectorlEHTTPairgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<HTTPair>[nElements] : new vector<HTTPair>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHTTPairgR(void *p) {
      delete ((vector<HTTPair>*)p);
   }
   static void deleteArray_vectorlEHTTPairgR(void *p) {
      delete [] ((vector<HTTPair>*)p);
   }
   static void destruct_vectorlEHTTPairgR(void *p) {
      typedef vector<HTTPair> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<HTTPair>

namespace {
  void TriggerDictionaryInitialization_libHTTEvent_Impl() {
    static const char* headers[] = {
"HTTEvent.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.26.04-edd28/x86_64-centos7-gcc11-opt/include",
"/scratch/Magisterka/RootAnalysis/HTauTau/DataFormats/src/../include",
"/scratch/Magisterka/RootAnalysis/HTauTau/DataFormats/src",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/6.26.04-edd28/x86_64-centos7-gcc11-opt/include/",
"/scratch/Magisterka/RootAnalysis/build/HTauTau/DataFormats/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libHTTEvent dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TLorentzVector.h")))  __attribute__((annotate("$clingAutoload$HTTEvent.h")))  TLorentzVector;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$TVector2.h")))  __attribute__((annotate("$clingAutoload$HTTEvent.h")))  TVector2;
class __attribute__((annotate("$clingAutoload$HTTEvent.h")))  HTTPair;
class __attribute__((annotate("$clingAutoload$HTTEvent.h")))  HTTParticle;
class __attribute__((annotate("$clingAutoload$HTTEvent.h")))  HTTEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libHTTEvent dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "HTTEvent.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HTTEvent", payloadCode, "@",
"HTTPair", payloadCode, "@",
"HTTParticle", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libHTTEvent",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libHTTEvent_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libHTTEvent_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libHTTEvent() {
  TriggerDictionaryInitialization_libHTTEvent_Impl();
}
