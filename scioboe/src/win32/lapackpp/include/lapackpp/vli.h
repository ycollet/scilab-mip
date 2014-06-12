//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector LongInt Class
//
//      A lightweight vector class with minimal overhead.
//
//      shallow assignment
//      unit stride
//      inlined access A(i)
//      optional (compile-time) array bounds checking through 
//              VECTOR_LONG_INT_BOUNDS_CHECK
//      A(i) is the same as A[i]
//      auto conversion to long int*
//      a null vector has size of 0, but has the ref_count structure
//              has been initalized
//

#ifndef _VECTOR_LONG_INT_H_
#define _VECTOR_LONG_INT_H_    

#include <iostream>       // for formatted printing of matrices

#ifndef __ASSERT_H
#include <cassert>     // cheap "error" protection used in checking
#endif                  // checking array bounds.

#include "arch.h"

typedef struct vrefLongInt {
      typedef long int value_type;
      int        sz;                                        
      value_type*    data;                                       
      int        ref_count;
      int        vref_ref_count;
      vrefLongInt(value_type *_data, int _sz)
	 : sz(_sz)
	 , data(_data)
	 , ref_count(2)
	 , vref_ref_count(1)
      {};
      vrefLongInt(int _sz)
	 : sz(_sz)
	 , data(new value_type[sz])
	 , ref_count(1)
	 , vref_ref_count(1)
      {};
} vrefLongInt;
                        


class DLLIMPORT VectorLongInt
{
   public:
      /// The type of the values in this vector
      typedef long int value_type;
   private:
      /// The type of the internal management structure
      typedef vrefLongInt vref_type;
      vref_type *p;
      value_type *data;            // performance hack, avoid COMPLEX
      // indirection to data.

      /** Dereferences the vref management structure,
       * deleting it if that was the last reference. */
      inline void unref_vref();
      /** Make this class a reference to the given vref
       * management structure. Make sure to call unref_vref()
       * beforehand, if suitable. */
      inline void ref_vref(vref_type* other);

    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
    //inline VectorLongInt();     // this should behave as VectorLongInt(0)
    VectorLongInt(unsigned);                             
    VectorLongInt(unsigned, long int);   // can't be inlined because of 'for'
                                       // statement.
    VectorLongInt(long int*, unsigned);
    VectorLongInt(long int*, unsigned, unsigned, bool);
    VectorLongInt(const VectorLongInt&); 
    ~VectorLongInt() ;                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline long int&        operator[](int); 
    inline long int&        operator[](int) const;  // read only
    inline long int&        operator()(int); 
    inline long int&        operator()(int) const; // read only
    inline              operator    long int*(void); 
    inline int          size() const;
    inline int          null() const;
           int          resize(unsigned d);
    inline int          ref_count() const;  // return the number of ref counts
    inline long int*        addr() const;
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    inline  VectorLongInt& operator=(const VectorLongInt&);
            VectorLongInt& operator=(long int);
    inline  VectorLongInt& ref(const VectorLongInt &);
            VectorLongInt& inject(const VectorLongInt&);
            VectorLongInt& copy(const VectorLongInt&);

    /* I/O */                                                      
    friend std::ostream&   operator<<(std::ostream&, const VectorLongInt&);       

};                                                                     


    // operators and member functions

inline int VectorLongInt::null()    const
{
    return (size() == 0) ;
}

inline int VectorLongInt::size() const
{
    return   p-> sz;
}


inline int VectorLongInt::ref_count() const
{
    return p->ref_count;
}

inline long int* VectorLongInt::addr() const
{
    return data;
}

inline VectorLongInt::operator long int*(void) 
{
    return data;
}


inline long int& VectorLongInt::operator()(int i)
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif 
    return data[i];
}

inline long int& VectorLongInt::operator()(int i) const
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline long int& VectorLongInt::operator[](int i)
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

//  *CHANGE*  [] is the same as ()
inline long int& VectorLongInt::operator[](int i) const
{
#ifdef VECTOR_LONG_INT_BOUNDS_CHECK
    assert(0<=i && i<size());
#endif  
    return data[i];
}

inline void VectorLongInt::ref_vref(vref_type* other)
{
   p = other;
   data = p->data;
   p->ref_count++;
   p->vref_ref_count++;
}

inline void VectorLongInt::unref_vref()
{
   if (--(p->ref_count) == 0)              // perform garbage col.
   {
      delete [] p->data;
      delete p;
   } 
   else
   {
      // Check whether the internal management structure needs to
      // be deleted
      if (--(p->vref_ref_count) == 0)
	 delete p;
   }
}

inline VectorLongInt& VectorLongInt::ref(const VectorLongInt& m)
{
   if (&m != this)
   {
      unref_vref();
      ref_vref(m.p);
   }
   return *this;
}

inline VectorLongInt& VectorLongInt::operator=(const VectorLongInt& m)
{

    return  ref(m);
}




#endif 
// _VECTOR_LONG_INT_H_

