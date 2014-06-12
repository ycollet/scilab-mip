#include <stack-c.h>
#include <sciprint.h>
#include <string.h>

#include <helper.h>

//////////////////////////////
// Various helper functions //
//////////////////////////////

int find_label_partial(char **LabelList, int nblabels, const char * LabelToFind)
{
  int Pos = -1, i;

  if (LabelList!=NULL)
    {
      for(i=0; i<nblabels; i++)
	{
	  // A bug in scilab: if the mlist contains only the type, the C API returns nblabels==2 !!
	  if (LabelList[i]!=NULL)
	    {
	      if (strncmp(LabelList[i],LabelToFind,strlen(LabelToFind))==0)
		{
		  Pos = i;
		  return Pos;
		} 
	    }
	} 
    }
  return Pos;
} 

int find_label(char **LabelList, int nblabels, const char * LabelToFind)
{
  int Pos = -1, i;

  if (LabelList!=NULL)
    {
      for(i=0; i<nblabels; i++)
	{
	  // A bug in scilab: if the mlist contains only the type, the C API returns nblabels==2 !!
	  if (LabelList[i]!=NULL)
	    {
	      if (strncmp(LabelList[i],LabelToFind,strlen(LabelToFind))==0)
		{
		  Pos = i;
		  return Pos;
		} 
	    }
	} 
    }

  return Pos;
} 
