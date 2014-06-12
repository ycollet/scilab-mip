// Copyright (C) 2009-2010 Yann COLLETTE. All Rights Reserved.
// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef SCILABIPOPTJOURNAL_HPP
#define SCILABIPOPTJOURNAL_HPP

#include <IpJournalist.hpp>

using namespace Ipopt;

// Class ScilabJournal.
// ---------------------------------------------------------------
// This class encapsulates journal output to the Scilab console.
class ScilabJournal : public Journal
{
public:
  
  // The constructor.
  ScilabJournal (EJournalLevel default_level);
  
  // The destructor.
  virtual ~ScilabJournal() { };

protected:
  virtual void PrintImpl      (EJournalCategory category, EJournalLevel level, const char* str);
  virtual void PrintfImpl     (EJournalCategory category, EJournalLevel level, const char* pformat, va_list ap);
  virtual void FlushBufferImpl();
};
#endif
