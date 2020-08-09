#include <string>
#include <iostream>

using namespace std;

#include <assert.h>
#include <string.h>

#include "gestchar.h"

GestChar::GestChar(const char * s)
{ 
  Text = (char *) malloc((unsigned) (strlen(s) + 1) * sizeof(char));
  assert(Text);
  strcpy(Text,s); 
}

GestChar::GestChar(const GestChar & s)
{
  strcpy(Text,s.Text);
}

GestChar & GestChar::operator = (const GestChar & s)
{
  if (Text != s.Text)
    {
      if (Text) free(Text);
      Text = (char *) malloc((unsigned) (strlen(s.Text) + 1) * 
			     sizeof(char));
      assert(Text);
      strcpy(Text,s.Text);
    }

  return *this;
}

GestChar operator + (const GestChar & s1, const GestChar & s2)
{
  char * s;
  s = (char *) malloc((strlen(s1.Data()) + strlen(s2.Data()) + 2) * 
		      sizeof(char));
  strcpy(s,s1.Data());
  strcat(s,s2.Data());

  return GestChar(s);
}

void InitStr(char * & Dest, char * Source)
{
  int Len = strlen(Source);
  Dest = (char *) malloc ((unsigned) (Len + 1) * sizeof(char));
  strcpy(Dest,Source);
}
