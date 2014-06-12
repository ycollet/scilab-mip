#ifndef HELPER_H
#define HELPER_H

#define SCICOINOR_ERROR if(_SciErr.iErr)	\
    {						\
      printError(&_SciErr, 0);			\
      return _SciErr.iErr;			\
    }

#undef DEBUG

// scilab master declares the basbrk symbol which allows to handle the ctrl+c event in scilab
#define HANDLE_CTRLC 1

#endif
