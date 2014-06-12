// Copyright (C) 2009-2010 Yann COLLETTE. All Rights Reserved.
// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#include <scilabexception.hpp>

// Function definitions for class ScilabException

ScilabException::ScilabException (const char* message) throw() : exception() 
{ 
  this->message = message;
}

ScilabException::ScilabException (const ScilabException& source) throw() : exception() 
{
  message = source.message;
}

ScilabException& ScilabException::operator=(const ScilabException& source) 
{ 
  message = source.message; 
  return *this;
}
