#ifndef ERROR_H
#define ERROR_H

#include <assert.h>
#include "gestchar.h"

enum CODE_ERROR { NONE, COMPILE_ERROR, EXEC_ERROR, MEM_ERROR, UNKNOWN };

class Error
{
  GestChar * ErrorMessage;
  CODE_ERROR code;
 public:
  Error() {}
  Error(char * Text) 
    { 
      ErrorMessage = new GestChar(Text);
      assert(ErrorMessage);
    }
  virtual void WriteError() const { cout << (*ErrorMessage).Data(); }
  virtual char * getErrorMessage() { return (*ErrorMessage).Data(); }
  void InitCodeError(CODE_ERROR c) { code = c; }
  CODE_ERROR CodeError() { return code; }
};

class ErrorCompile : public Error
{
  GestChar * s;
 public:
 ErrorCompile(char * Text) : Error(Text) 
  { InitCodeError(COMPILE_ERROR); }

  void WriteError() const 
  { 
    cout << "\n*** Compile error ***\n";
    Error::WriteError();
    cout << endl; 
  }

  virtual char * getErrorMessage()
  { 
    s = new GestChar("*** FreeFem compile error ***\r\n");
    assert(s);
    *s = *s + Error::getErrorMessage();
    *s = *s + "\r\n";
    return (*s).Data();
  }
};


class ErrorExec : public Error
{
  GestChar * s;
 public:
 ErrorExec(char * Text) : Error(Text) { InitCodeError(EXEC_ERROR); }
  void WriteError() const 
  { 
    cout << "\n*** Execute error ***\n";
    Error::WriteError();
    cout << endl; 
  }

  virtual char * getErrorMessage()
  { 
    s = new GestChar("*** FreeFem execute error ***\r\n");
    assert(s);
    *s = *s + Error::getErrorMessage();
    *s = *s + "\r\n";
    return (*s).Data();
  }
};

class ErrorMemory : public Error
{
  GestChar * s;
 public:
 ErrorMemory(char * Text) : Error(Text) { InitCodeError(MEM_ERROR); }
  void WriteError() const 
  { 
    cout << "\n-> Memory allocation problem\n";
    Error::WriteError();
    cout << endl; 
  }

  virtual char * getErrorMessage()
  { 
    s = new GestChar("*** FreeFem memory allocation problem ***\r\n");
    assert(s);
    *s = *s + Error::getErrorMessage();
    *s = *s + "\r\n";
    return (*s).Data();
  }
};

#endif
