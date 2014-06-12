// Copyright (C) 2009-2010 Yann COLLETTE. All Rights Reserved.
// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include <stdio.h>

#include <scilabjournal.hpp>
#include <scilabexception.hpp>

extern "C"
{
#include <stack-c.h>
#include <sciprint.h>
}

// Function definitions for class ScilabJournal.

ScilabJournal::ScilabJournal (EJournalLevel default_level) : Journal("scilab", default_level) { }

void ScilabJournal::PrintImpl (EJournalCategory category, EJournalLevel level, const char* str) 
{
  sciprint((char *)str);
}

void ScilabJournal::PrintfImpl (EJournalCategory category, EJournalLevel level, const char* pformat, va_list ap) 
{
  const int maxStrLen = 1024;
  char      s[maxStrLen];

  if (vsnprintf(s,maxStrLen,pformat,ap) >= maxStrLen)
    throw ScilabException("String buffer it too short for all the characters to be printed to Scilab console");
  sciprint(s);
}

void ScilabJournal::FlushBufferImpl() { }
