// Copyright (C) 2009-2010 Yann COLLETTE. All Rights Reserved.
// Copyright (C) 2007 Peter Carbonetto. All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Peter Carbonetto
//         Dept. of Computer Science
//         University of British Columbia
//         May 19, 2007

#ifndef SCILABEXCEPTION_HPP
#define SCILABEXCEPTION_HPP

#include <exception>

// Class ScilabException
// -----------------------------------------------------------------
// It is assumed that the argument passed to the constructor persists
// as long as the ScilabException object is in scope. Usually, this
// means that it should persist for the duration of the entire
// program. This is always the case if the input "message" is a literal.

class ScilabException : public std::exception 
{
public:
  ScilabException (const char* message) throw();
  ~ScilabException()                    throw() { };
  
  // The copy constructor makes a shallow copy.
  ScilabException (const ScilabException& source) throw();
  
  // The copy assignment operator makes a shallow copy as well.
  ScilabException& operator= (const ScilabException& source);
  
  // Return the message string.
  virtual const char* what () const throw() { return message; };
  
private:
    const char* message;  // The error message.
};
#endif
