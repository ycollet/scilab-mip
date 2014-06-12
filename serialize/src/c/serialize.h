#ifndef SERIALIZE_H
#define SERIALIZE_H

#include <api_scilab.h>

union scilab_data_t;

struct serialize_list
{
  int nb_elems;
  int type;
  int precision;
  union scilab_data_t * data;
};

// TODO:
// - poly
// - [u]int64
// - handle
// - u_function
// - c_function
// - lib
// - implicit_poly
// - intrinsic_function

struct serialize_sci_matrix
{
  int piRows, piCols;
  int isComplex;
  double * real_data;
  double * imag_data;
};

struct serialize_sci_boolean
{
  int piRows, piCols;
  int * bool_data;
};

struct serialize_sci_sparse
{
  int piRows, piCols;
  int piNbItems, isComplex;
  int * piNbItemRow, * piColsPos;
  double * real_data, * imag_data;
};

struct serialize_sci_boolean_sparse
{
  int piRows, piCols;
  int piNbItems;
  int * piNbItemRow, * piColPos;
};

struct serialize_sci_int8
{
  int piRows, piCols;
  char * int8_data;
};

struct serialize_sci_int16
{
  int piRows, piCols;
  short * int16_data;
};

struct serialize_sci_int32
{
  int piRows, piCols;
  int * int32_data;
};

struct serialize_sci_uint8
{
  int piRows, piCols;
  unsigned char * uint8_data;
};

struct serialize_sci_uint16
{
  int piRows, piCols;
  unsigned short * uint16_data;
};

struct serialize_sci_uint32
{
  int piRows, piCols;
  unsigned int * int32_data;
};

struct serialize_sci_strings
{
  int piRows, piCols;
  char ** pstStrings;
};

struct serialize_sci_list
{
  int nbItem;
  struct serialize_list * list;
};

struct serialize_sci_mlist
{
  int nbItem;
  struct serialize_list * mlist;
};

struct serialize_sci_tlist
{
  int nbItem;
  struct serialize_list * tlist;
};

struct serialize_sci_pointer
{
  void * ptr;
};

union scilab_data_t
{
  struct serialize_sci_matrix         matrix;
  struct serialize_sci_boolean        boolean;
  struct serialize_sci_sparse         sparse;
  struct serialize_sci_boolean_sparse boolean_sparse;
  struct serialize_sci_int8           int8;
  struct serialize_sci_int16          int16;
  struct serialize_sci_int32          int32;
  struct serialize_sci_uint8          uint8;
  struct serialize_sci_uint16         uint16;
  struct serialize_sci_uint32         uint32;
  struct serialize_sci_strings        strings;
  struct serialize_sci_list           list;
  struct serialize_sci_mlist          mlist;
  struct serialize_sci_tlist          tlist;
  struct serialize_sci_pointer        pointer;
  struct serialize_list               element;
};

void serialize(int * piAddress, struct serialize_list * _list);
void unserialize(int iVar, struct serialize_list * _list, int do_free, int in_list, int index_item, int * parent);
#endif
