
// Arnaud Desitter
// $Id: cxSmartPtr.h,v 1.1 2002-06-18 14:49:58 arnaud Exp $

// Basic "smart pointer" class to handle an Iris Explorer reference counted 
// structure (the design is inspired by the CORBA smart pointer).
// Make it exception safe (no memory leak).

#ifndef _CXSMARTPTR_H
#define _CXSMARTPTR_H

#include <cx/DataOps.h>   // to get cxDataRefDec()

template<class T>
class cxSmartPtrVar
{
  T* ptr_;

  void destroy(){ if (ptr_) cxDataRefDec(ptr_); };

  cxSmartPtrVar(const cxSmartPtrVar<T>& r);
  cxSmartPtrVar<T>& operator=(const cxSmartPtrVar<T>&);
public:
  cxSmartPtrVar() : ptr_(0) { }
  cxSmartPtrVar(T* p) : ptr_(p) { }
  ~cxSmartPtrVar() { destroy(); }

  cxSmartPtrVar<T>& operator=(T*);
  T* operator->() { return ptr_; }
  const T& in() const { return *ptr_; }
//  T& inout() { return *ptr_; }
//  T*& out() { delete ptr_; ptr_ = 0; return ptr_; }
  T* _retn() { T* ret = ptr_; ptr_ = 0; return ret; }
#ifdef HAVE_NO_CONST_TYPE_CONVERSION_OVERLOAD
  operator T&() const
    { return *ptr_; }
#else
# ifndef OB_ONLY_IN_TYPE_CONVERSION
  operator T&()
    { return *ptr_; }
# endif
  operator const T&() const
    { return *ptr_; }
#endif
  operator T*&() { return ptr_; }
};

template<class T>
cxSmartPtrVar<T>&
cxSmartPtrVar<T>::operator=(T* p)
{
  destroy();
  ptr_ = p;
  return *this;
}

//
// Instantiate some Explorer data type
//

#include<cx/DataAccess.h>
#define CXSMARTPOINTER(TYPE)                 \
    typedef cxSmartPtrVar< TYPE > TYPE##_var;

CXSMARTPOINTER(cxLattice)
CXSMARTPOINTER(cxData)
CXSMARTPOINTER(cxCoord)
CXSMARTPOINTER(cxPyramid)
CXSMARTPOINTER(cxPyramidDictionary)
CXSMARTPOINTER(cxConnection)
CXSMARTPOINTER(cxGeometry)

#endif
