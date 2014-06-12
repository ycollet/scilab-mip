#include <string.h>

#include <CoinMessageHandler.hpp>

extern "C" {
#include <stack-c.h>
#include <sciprint.h>
}

#undef min
#undef max

#include <helper.hpp>

/////////////////////
// Message handler //
/////////////////////

int DerivedHandler::print()
{
  sciprint("%s\n", messageBuffer());
  return 0;
}

///////////////////////
// Clp event handler //
///////////////////////

MyClpEventHandler::MyClpEventHandler() : ClpEventHandler(), inCbc(false)
{
}

MyClpEventHandler::MyClpEventHandler(const MyClpEventHandler & rhs) : ClpEventHandler(rhs), inCbc(rhs.inCbc)
{  
}

MyClpEventHandler::MyClpEventHandler(ClpSimplex * model) : ClpEventHandler(model), inCbc(false)
{
}

MyClpEventHandler::~MyClpEventHandler()
{
}

MyClpEventHandler & MyClpEventHandler::operator=(const MyClpEventHandler& rhs)
{
  if (this != &rhs) 
    {
      ClpEventHandler::operator=(rhs);
      inCbc = rhs.inCbc;
    }
  return *this;
}

ClpEventHandler * MyClpEventHandler::clone() const
{
  return new MyClpEventHandler(*this);
}

int MyClpEventHandler::event(Event whichEvent)
{
  switch(whichEvent)
    {
    case endOfIteration:
    case endOfFactorization:
    case endOfValuesPass:
    case node:
    case treeStatus:
    case solution:
    case theta:
      /* .true. fortran is -1 :'( */
      if (C2F(basbrk).iflag == -1)
	{
#ifdef DEBUG
	  sciprint("inCbc = %d\n", (int)inCbc);
	  sciprint("MyClpEventHandler: iflag = %d interruptible = %d\n", C2F(basbrk).iflag, C2F(basbrk).interruptible);
#endif
	  if (!inCbc) C2F(basbrk).iflag = 0;
	  //C2F(basbrk).interruptible = 0;
	  // If return code -1 then carries on 
	  // if 0 sets ClpModel::status() to 5 (stopped by event) and will return to user. 
	  // At present if <-1 carries on and if >0 acts as if 0 - this may change 
	  return 0;
	}
      return -1;
      break;
    default:
      return -1;
      break;
    }
}


///////////////////////
// Cbc event handler //
///////////////////////

MyCbcEventHandler::MyCbcEventHandler() : CbcEventHandler(), stopOnFirstSol(false)
{
}

MyCbcEventHandler::MyCbcEventHandler(const MyCbcEventHandler & rhs) : CbcEventHandler(rhs), stopOnFirstSol(rhs.stopOnFirstSol)
{  
}

MyCbcEventHandler::MyCbcEventHandler(CbcModel * model) : CbcEventHandler(model), stopOnFirstSol(false)
{
}

MyCbcEventHandler::~MyCbcEventHandler()
{
}

MyCbcEventHandler & MyCbcEventHandler::operator=(const MyCbcEventHandler& rhs)
{
  if (this != &rhs)
    {
      CbcEventHandler::operator=(rhs);
    }
  stopOnFirstSol = rhs.stopOnFirstSol;

  return *this;
}

CbcEventHandler * MyCbcEventHandler::clone() const
{
  return new MyCbcEventHandler(*this);
}

CbcEventHandler::CbcAction MyCbcEventHandler::event(CbcEvent whichEvent)
{
  switch(whichEvent)
    {
    case solution:
      if (stopOnFirstSol) return stop;
    case heuristicSolution:
      if (stopOnFirstSol) return stop;
    case node:
    case treeStatus:
    default:
      /* .true. fortran is -1 :'( */
      if (C2F(basbrk).iflag == -1)
	{
#ifdef DEBUG
	  sciprint("MyCbcEventHandler: iflag = %d interruptible = %d\n", C2F(basbrk).iflag, C2F(basbrk).interruptible);
#endif
	  C2F(basbrk).iflag = 0;
	  //C2F(basbrk).interruptible = 0;
	  // If return code -1 then carries on 
	  // if 0 sets ClpModel::status() to 5 (stopped by event) and will return to user. 
	  // At present if <-1 carries on and if >0 acts as if 0 - this may change 
	  return stop;
	}

      return noAction;
    }
}
