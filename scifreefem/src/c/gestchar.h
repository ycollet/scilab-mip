#ifndef GESTCHAR_H
#define GESTCHAR_H

#include <stdlib.h>

class GestChar
{
  char * Text;
 public:
  GestChar() {}
  GestChar(const char * );
  GestChar(const GestChar & );
  ~GestChar() { free(Text); }
  char * Data() const { return Text; };

  GestChar & operator = (const GestChar & s);
};

GestChar operator + (const GestChar & s1, const GestChar & s2);

void InitStr(char * & Dest, char * Source);

#endif
