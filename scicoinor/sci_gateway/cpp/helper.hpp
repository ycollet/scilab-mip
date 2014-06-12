#ifndef HELPER_H
#define HELPER_H

#include <CoinMessageHandler.hpp>
#include <CbcEventHandler.hpp>
#include <ClpEventHandler.hpp>

#define SCICOINOR_ERROR if(_SciErr.iErr)	\
    {						\
      printError(&_SciErr, 0);			\
      return _SciErr.iErr;			\
    }

class DerivedHandler : public CoinMessageHandler 
{
public:
  virtual int print();
};

#undef DEBUG

// scilab master declares the basbrk symbol which allows to handle the ctrl+c event in scilab
#define HANDLE_CTRLC 1

#ifdef HANDLE_CTRLC
class MyClpEventHandler : public ClpEventHandler 
{
public:
  MyClpEventHandler();
  MyClpEventHandler(ClpSimplex * model);
  virtual ~MyClpEventHandler();
  MyClpEventHandler(const MyClpEventHandler & rhs);
  MyClpEventHandler& operator=(const MyClpEventHandler & rhs);
  virtual int event(Event whichEvent);
  void setInCbc(bool _inCbc) {inCbc = _inCbc;}
  bool getInCbc() const {return inCbc;}
  virtual ClpEventHandler * clone() const ;
protected:
  bool inCbc;
};

class MyCbcEventHandler : public CbcEventHandler 
{
public:
  MyCbcEventHandler();
  MyCbcEventHandler(CbcModel * model);
  virtual ~MyCbcEventHandler();
  MyCbcEventHandler(const MyCbcEventHandler & rhs);
  MyCbcEventHandler & operator=(const MyCbcEventHandler & rhs);
  virtual CbcAction event(CbcEvent whichEvent);
  void setStopOnFirstSol(bool _stopOnFirstSol) {stopOnFirstSol = _stopOnFirstSol;}
  bool getStopOnFirstSol() const {return stopOnFirstSol;}
  virtual CbcEventHandler * clone() const ;
protected:
  bool stopOnFirstSol;
};
#endif
#endif
