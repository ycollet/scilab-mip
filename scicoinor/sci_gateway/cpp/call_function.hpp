#ifndef CALL_FUNCTION_HPP
#define CALL_FUNCTION_HPP

#include <vector>

using namespace std;

typedef void (*voidf)();

typedef struct {
  char *name;
  voidf f;
} FTAB;

struct param_fobj {
  int * param_addr;
  int   param_type;
};

voidf GetFunctionPtr(char *, int, FTAB *, voidf, int *, int*, int*, int*, vector<param_fobj> &);
#endif
