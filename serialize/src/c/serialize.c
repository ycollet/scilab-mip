#include <string.h>

#include <api_scilab.h>
#include <stack-c.h>
#include <MALLOC.h>
#include <Scierror.h>

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
  int * piNbItemRow, * piColPos;
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
  unsigned int * uint32_data;
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

void serialize(int * piAddress, struct serialize_list * _list)
{
  int type, precision, nb_elems, i;
  SciErr _SciErr;
  double * tmp_real = NULL;
  double * tmp_imag = NULL;
  int * tmp_bool = NULL;
  int * tmp_nbitemrow = NULL;
  int * tmp_picolpos  = NULL;
  char * tmp_int8 = NULL;
  short * tmp_int16 = NULL;
  int * tmp_int32 = NULL;
  unsigned char * tmp_uint8 = NULL;
  unsigned short * tmp_uint16 = NULL;
  unsigned int * tmp_uint32 = NULL;
  int * piItemAddress = NULL;

  _SciErr = getVarType(pvApiCtx, piAddress, &type);

  switch(type)
    {
    case sci_matrix:
      _list->nb_elems = 1;
      _list->type = sci_matrix;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      _list->data[0].matrix.isComplex = isVarComplex(pvApiCtx, piAddress);
      if (_list->data[0].matrix.isComplex)
        {
          _SciErr = getComplexMatrixOfDouble(pvApiCtx, piAddress,&_list->data[0].matrix.piRows, &_list->data[0].matrix.piCols, &tmp_real, &tmp_imag);
          _list->data[0].matrix.real_data = (double *)MALLOC(_list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
          memcpy(_list->data[0].matrix.real_data, tmp_real, _list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
          _list->data[0].matrix.imag_data = (double *)MALLOC(_list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
          memcpy(_list->data[0].matrix.imag_data, tmp_imag, _list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
        }
      else
        {
          _SciErr = getMatrixOfDouble(pvApiCtx, piAddress,&_list->data[0].matrix.piRows, &_list->data[0].matrix.piCols, &tmp_real);
          _list->data[0].matrix.real_data = (double *)MALLOC(_list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
          memcpy(_list->data[0].matrix.real_data, tmp_real, _list->data[0].matrix.piRows*_list->data[0].matrix.piCols*sizeof(double));
        }
      break;
    case sci_boolean:
      _list->nb_elems = 1;
      _list->type = sci_boolean;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      _SciErr = getMatrixOfBoolean(pvApiCtx, piAddress,&_list->data[0].boolean.piRows, &_list->data[0].boolean.piCols, &tmp_bool);
      _list->data[0].boolean.bool_data = (int *)MALLOC(_list->data[0].boolean.piRows*_list->data[0].boolean.piCols*sizeof(int));
      memcpy(_list->data[0].boolean.bool_data, tmp_bool, _list->data[0].boolean.piRows*_list->data[0].boolean.piCols*sizeof(int));
      break;
    case sci_sparse:
      _list->nb_elems = 1;
      _list->type = sci_sparse;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      _list->data[0].sparse.isComplex = isVarComplex(pvApiCtx, piAddress);
      if (_list->data[0].sparse.isComplex)
        {
          _SciErr = getComplexSparseMatrix(pvApiCtx, piAddress, &_list->data[0].sparse.piRows, &_list->data[0].sparse.piCols, 
                                           &_list->data[0].sparse.piNbItems, 
                                           &tmp_nbitemrow, &tmp_picolpos, 
                                           &tmp_real, &tmp_imag);
          _list->data[0].sparse.piNbItemRow = (int *)MALLOC(_list->data[0].sparse.piRows*sizeof(int));
          _list->data[0].sparse.piColPos    = (int *)MALLOC(_list->data[0].sparse.piNbItems*sizeof(int));
          _list->data[0].sparse.real_data   = (double *)MALLOC(_list->data[0].sparse.piNbItems*sizeof(int));
          _list->data[0].sparse.imag_data   = (double *)MALLOC(_list->data[0].sparse.piNbItems*sizeof(int));
          memcpy(_list->data[0].sparse.piNbItemRow, tmp_nbitemrow, _list->data[0].sparse.piRows*sizeof(int));
          memcpy(_list->data[0].sparse.piColPos,    tmp_picolpos,  _list->data[0].sparse.piNbItems*sizeof(int));
          memcpy(_list->data[0].sparse.real_data,   tmp_real,      _list->data[0].sparse.piNbItems*sizeof(int));
          memcpy(_list->data[0].sparse.imag_data,   tmp_imag,      _list->data[0].sparse.piNbItems*sizeof(int));
        }
      else
        {
          _SciErr = getSparseMatrix(pvApiCtx, piAddress, &_list->data[0].sparse.piRows, &_list->data[0].sparse.piCols, 
                                    &_list->data[0].sparse.piNbItems, 
                                    &tmp_nbitemrow, &tmp_picolpos, &tmp_real);
          _list->data[0].sparse.piNbItemRow = (int *)MALLOC(_list->data[0].sparse.piRows*sizeof(int));
          _list->data[0].sparse.piColPos    = (int *)MALLOC(_list->data[0].sparse.piNbItems*sizeof(int));
          _list->data[0].sparse.real_data   = (double *)MALLOC(_list->data[0].sparse.piNbItems*sizeof(int));
          memcpy(_list->data[0].sparse.piNbItemRow, tmp_nbitemrow, _list->data[0].sparse.piRows*sizeof(int));
          memcpy(_list->data[0].sparse.piColPos,    tmp_picolpos,  _list->data[0].sparse.piNbItems*sizeof(int));
          memcpy(_list->data[0].sparse.real_data,   tmp_real,      _list->data[0].sparse.piNbItems*sizeof(int));
        }
      break;
    case sci_boolean_sparse:
      _list->nb_elems = 1;
      _list->type = sci_boolean_sparse;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      _SciErr = getBooleanSparseMatrix(pvApiCtx, piAddress, &_list->data[0].boolean_sparse.piRows, &_list->data[0].boolean_sparse.piCols, 
                                       &_list->data[0].boolean_sparse.piNbItems,
                                       &tmp_nbitemrow, &tmp_picolpos);
      _list->data[0].boolean_sparse.piNbItemRow = (int *)MALLOC(_list->data[0].boolean_sparse.piRows*sizeof(int));
      _list->data[0].boolean_sparse.piColPos    = (int *)MALLOC(_list->data[0].boolean_sparse.piNbItems*sizeof(int));
      memcpy(_list->data[0].boolean_sparse.piNbItemRow, tmp_nbitemrow, _list->data[0].boolean_sparse.piRows*sizeof(int));
      memcpy(_list->data[0].boolean_sparse.piColPos,    tmp_picolpos,  _list->data[0].boolean_sparse.piNbItems*sizeof(int));
      break;
    case sci_ints:
      _list->nb_elems = 1;
      _SciErr = getMatrixOfIntegerPrecision(pvApiCtx, piAddress, &precision);
      _list->type      = sci_ints;
      _list->precision = precision;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      switch(precision)
        {
        case SCI_INT8:
          _SciErr = getMatrixOfInteger8(pvApiCtx, piAddress, &_list->data[0].int8.piRows, &_list->data[0].int8.piCols, &tmp_int8);
          _list->data[0].int8.int8_data = (char *)MALLOC(_list->data[0].int8.piRows*_list->data[0].int8.piCols*sizeof(char));
          memcpy(_list->data[0].int8.int8_data, tmp_int8, _list->data[0].int8.piRows*_list->data[0].int8.piCols*sizeof(char));
          break;
        case SCI_INT16:
          _SciErr = getMatrixOfInteger16(pvApiCtx, piAddress, &_list->data[0].int16.piRows, &_list->data[0].int16.piCols, &tmp_int16);
          _list->data[0].int16.int16_data = (short *)MALLOC(_list->data[0].int16.piRows*_list->data[0].int16.piCols*sizeof(short));
          memcpy(_list->data[0].int16.int16_data, tmp_int16, _list->data[0].int16.piRows*_list->data[0].int16.piCols*sizeof(short));
          break;
        case SCI_INT32:
          _SciErr = getMatrixOfInteger32(pvApiCtx, piAddress, &_list->data[0].int32.piRows, &_list->data[0].int32.piCols, &tmp_int32);
          _list->data[0].int32.int32_data = (int *)MALLOC(_list->data[0].int32.piRows*_list->data[0].int32.piCols*sizeof(int));
          memcpy(_list->data[0].int32.int32_data, tmp_int32, _list->data[0].int32.piRows*_list->data[0].int32.piCols*sizeof(int));
          break;
        case SCI_UINT8:
          _SciErr = getMatrixOfUnsignedInteger8(pvApiCtx, piAddress, &_list->data[0].uint8.piRows, &_list->data[0].uint8.piCols, &tmp_uint8);
          _list->data[0].uint8.uint8_data = (unsigned char *)MALLOC(_list->data[0].uint8.piRows*_list->data[0].uint8.piCols*sizeof(unsigned char));
          memcpy(_list->data[0].uint8.uint8_data, tmp_uint8, _list->data[0].uint8.piRows*_list->data[0].uint8.piCols*sizeof(unsigned char));
          break;
        case SCI_UINT16:
          _SciErr = getMatrixOfUnsignedInteger16(pvApiCtx, piAddress, &_list->data[0].uint16.piRows, &_list->data[0].uint16.piCols, &tmp_uint16);
          _list->data[0].uint16.uint16_data = (unsigned short *)MALLOC(_list->data[0].uint16.piRows*_list->data[0].uint16.piCols*sizeof(unsigned short));
          memcpy(_list->data[0].uint16.uint16_data, tmp_uint16, _list->data[0].uint16.piRows*_list->data[0].uint16.piCols*sizeof(unsigned short));
          break;
        case SCI_UINT32:
          _SciErr = getMatrixOfUnsignedInteger32(pvApiCtx, piAddress, &_list->data[0].uint32.piRows, &_list->data[0].uint32.piCols, &tmp_uint32);
          _list->data[0].uint32.uint32_data = (unsigned int *)MALLOC(_list->data[0].uint32.piRows*_list->data[0].uint32.piCols*sizeof(unsigned int));
          memcpy(_list->data[0].uint32.uint32_data, tmp_int32, _list->data[0].uint32.piRows*_list->data[0].uint32.piCols*sizeof(unsigned int));
          break;
        }
      break;
    case sci_strings:
      _list->nb_elems = 1;
      _list->type = sci_strings;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));
      getAllocatedMatrixOfString(pvApiCtx, piAddress, &_list->data[0].strings.piRows, &_list->data[0].strings.piCols, &_list->data[0].strings.pstStrings);
      break;
    case sci_list:
      _SciErr = getListItemNumber(pvApiCtx, piAddress, &nb_elems);
      _list->nb_elems = nb_elems;
      _list->type = sci_list;
      _list->data = (union scilab_data_t *)MALLOC(nb_elems*sizeof(union scilab_data_t));
      for(i=0;i<nb_elems;i++)
        {
          _SciErr = getListItemAddress(pvApiCtx, piAddress, i+1, &piItemAddress);
          serialize(piItemAddress, &_list->data[i].element);
        }
      break;
    case sci_tlist:
      _SciErr = getListItemNumber(pvApiCtx, piAddress, &nb_elems);
      _list->nb_elems = nb_elems;
      _list->type = sci_tlist;
      _list->data = (union scilab_data_t *)MALLOC(nb_elems*sizeof(union scilab_data_t));
      for(i=0;i<nb_elems;i++)
        {
          _SciErr = getListItemAddress(pvApiCtx, piAddress, i+1, &piItemAddress);
          serialize(piItemAddress, &_list->data[i].element);
        }
      break;
    case sci_mlist:
      _SciErr = getListItemNumber(pvApiCtx, piAddress, &nb_elems);
      _list->nb_elems = nb_elems;
      _list->type = sci_mlist;
      _list->data = (union scilab_data_t *)MALLOC(nb_elems*sizeof(union scilab_data_t));
      for(i=0;i<nb_elems;i++)
        {
          _SciErr = getListItemAddress(pvApiCtx, piAddress, i+1, &piItemAddress);
          serialize(piItemAddress, &_list->data[i].element);
        }
      break;
    case sci_pointer:
      _list->nb_elems = 1;
      _list->type = sci_pointer;
      _list->data = (union scilab_data_t *)MALLOC(1*sizeof(union scilab_data_t));

      _SciErr = getPointer(pvApiCtx, piAddress, &_list->data[0].pointer.ptr);
      break;
    default:
      Scierror(999,"unsupported data type - %d", type);
      break;
    }
}

void unserialize(int iVar, struct serialize_list * _list, int do_free, int in_list, int index_item, int * parent)
{
  int i;
  int * p_list_addr = NULL;
  SciErr _SciErr;

  switch(_list->type)
    {
    case sci_matrix:
      if (_list->data[0].matrix.isComplex)
        {
          if (in_list)
            {
              _SciErr = createComplexMatrixOfDoubleInList(pvApiCtx, iVar, parent, index_item, 
                                                          _list->data[0].matrix.piRows, _list->data[0].matrix.piCols, 
                                                          _list->data[0].matrix.real_data, 
                                                          _list->data[0].matrix.imag_data);
            }
          else
            {
              _SciErr = createComplexMatrixOfDouble(pvApiCtx, iVar, _list->data[0].matrix.piRows, _list->data[0].matrix.piCols, 
                                                    _list->data[0].matrix.real_data, 
                                                    _list->data[0].matrix.imag_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].matrix.real_data);
              FREE(_list->data[0].matrix.imag_data);
            }
        }
      else
        {
          if (in_list)
            {
              _SciErr = createMatrixOfDoubleInList(pvApiCtx, iVar, parent, index_item, 
                                                   _list->data[0].matrix.piRows, _list->data[0].matrix.piCols, 
                                                   _list->data[0].matrix.real_data);
            }
          else
            {
              _SciErr = createMatrixOfDouble(pvApiCtx, iVar, _list->data[0].matrix.piRows, _list->data[0].matrix.piCols, 
                                             _list->data[0].matrix.real_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].matrix.real_data);
            }
        }
      
      if (do_free) FREE(_list->data);
      break;
    case sci_boolean:
      if (in_list)
        {
          _SciErr = createMatrixOfBooleanInList(pvApiCtx, iVar, parent, index_item,
                                                _list->data[0].boolean.piRows, _list->data[0].boolean.piCols, 
                                                _list->data[0].boolean.bool_data);
        }
      else
        {
          _SciErr = createMatrixOfBoolean(pvApiCtx, iVar, _list->data[0].boolean.piRows, _list->data[0].boolean.piCols, 
                                          _list->data[0].boolean.bool_data);
        }

      if (do_free)
        {
          FREE(_list->data[0].boolean.bool_data);
          FREE(_list->data);
        }
      break;
    case sci_sparse:
      if (_list->data[0].sparse.isComplex)
        {
          if (in_list)
            {
              _SciErr = createComplexSparseMatrixInList(pvApiCtx, iVar, parent, index_item,
                                                        _list->data[0].sparse.piRows, _list->data[0].sparse.piCols, 
                                                        _list->data[0].sparse.piNbItems, 
                                                        _list->data[0].sparse.piNbItemRow, _list->data[0].sparse.piColPos, 
                                                        _list->data[0].sparse.real_data, _list->data[0].sparse.imag_data);
            }
          else
            {
              _SciErr = createComplexSparseMatrix(pvApiCtx, iVar, _list->data[0].sparse.piRows, _list->data[0].sparse.piCols, 
                                                  _list->data[0].sparse.piNbItems, 
                                                  _list->data[0].sparse.piNbItemRow, _list->data[0].sparse.piColPos, 
                                                  _list->data[0].sparse.real_data, _list->data[0].sparse.imag_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].sparse.piNbItemRow);
              FREE(_list->data[0].sparse.piColPos);
              FREE(_list->data[0].sparse.real_data);
              FREE(_list->data[0].sparse.imag_data);
              FREE(_list->data);
            }
        }
      else
        {
          if (in_list)
            {
              _SciErr = createSparseMatrixInList(pvApiCtx, iVar, parent, index_item, 
                                                 _list->data[0].sparse.piRows, _list->data[0].sparse.piCols, 
                                                 _list->data[0].sparse.piNbItems, 
                                                 _list->data[0].sparse.piNbItemRow, _list->data[0].sparse.piColPos, _list->data[0].sparse.real_data);
            }
          else
            {
              _SciErr = createSparseMatrix(pvApiCtx, iVar, _list->data[0].sparse.piRows, _list->data[0].sparse.piCols, 
                                           _list->data[0].sparse.piNbItems, 
                                           _list->data[0].sparse.piNbItemRow, _list->data[0].sparse.piColPos, _list->data[0].sparse.real_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].sparse.piNbItemRow);
              FREE(_list->data[0].sparse.piColPos);
              FREE(_list->data[0].sparse.real_data);
              FREE(_list->data);
            }
        }
      break;
    case sci_boolean_sparse:
      if (in_list)
        {
          _SciErr = createBooleanSparseMatrixInList(pvApiCtx, iVar, parent, index_item,
                                                    _list->data[0].boolean_sparse.piRows, _list->data[0].boolean_sparse.piCols, 
                                                    _list->data[0].boolean_sparse.piNbItems,
                                                    _list->data[0].boolean_sparse.piNbItemRow, _list->data[0].boolean_sparse.piColPos);
        }
      else
        {
          _SciErr = createBooleanSparseMatrix(pvApiCtx, iVar, _list->data[0].boolean_sparse.piRows, _list->data[0].boolean_sparse.piCols, 
                                              _list->data[0].boolean_sparse.piNbItems,
                                              _list->data[0].boolean_sparse.piNbItemRow, _list->data[0].boolean_sparse.piColPos);
        }

      if (do_free)
        {
          FREE(_list->data[0].boolean_sparse.piNbItemRow);
          FREE(_list->data[0].boolean_sparse.piColPos);
          FREE(_list->data);
        }
      break;
    case sci_ints:
      switch(_list->precision)
        {
        case SCI_INT8:
          if (in_list)
            {
              _SciErr = createMatrixOfInteger8InList(pvApiCtx, iVar, parent, index_item,
                                                     _list->data[0].int8.piRows, _list->data[0].int8.piCols, 
                                                     _list->data[0].int8.int8_data);
            }
          else
            {
              _SciErr = createMatrixOfInteger8(pvApiCtx, iVar, _list->data[0].int8.piRows, _list->data[0].int8.piCols, 
                                               _list->data[0].int8.int8_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].int8.int8_data);
              FREE(_list->data);
            }
          break;
        case SCI_INT16:
          if (in_list)
            {
              _SciErr = createMatrixOfInteger16InList(pvApiCtx, iVar, parent, index_item,
                                                      _list->data[0].int16.piRows, _list->data[0].int16.piCols, 
                                                      _list->data[0].int16.int16_data);
            }
          else
            {
              _SciErr = createMatrixOfInteger16(pvApiCtx, iVar, _list->data[0].int16.piRows, _list->data[0].int16.piCols, 
                                                _list->data[0].int16.int16_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].int16.int16_data);
              FREE(_list->data);
            }
          break;
        case SCI_INT32:
          if (in_list)
            {
              _SciErr = createMatrixOfInteger32InList(pvApiCtx, iVar, parent, index_item,
                                                      _list->data[0].int32.piRows, _list->data[0].int32.piCols, 
                                                      _list->data[0].int32.int32_data);
            }
          else
            {
              _SciErr = createMatrixOfInteger32(pvApiCtx, iVar, _list->data[0].int32.piRows, _list->data[0].int32.piCols, 
                                                _list->data[0].int32.int32_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].int32.int32_data);
              FREE(_list->data);
            }
          break;
        case SCI_UINT8:
          if (in_list)
            {
              _SciErr = createMatrixOfUnsignedInteger8InList(pvApiCtx, iVar, parent, index_item, 
                                                             _list->data[0].uint8.piRows, _list->data[0].uint8.piCols, 
                                                             _list->data[0].uint8.uint8_data);
            }
          else
            {
              _SciErr = createMatrixOfUnsignedInteger8(pvApiCtx, iVar, _list->data[0].uint8.piRows, _list->data[0].uint8.piCols, 
                                                       _list->data[0].uint8.uint8_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].uint8.uint8_data);
              FREE(_list->data);
            }
          break;
        case SCI_UINT16:
          if (in_list)
            {
              _SciErr = createMatrixOfUnsignedInteger16InList(pvApiCtx, iVar, parent, index_item,
                                                              _list->data[0].uint16.piRows, _list->data[0].uint16.piCols, 
                                                              _list->data[0].uint16.uint16_data);
            }
          else
            {
              _SciErr = createMatrixOfUnsignedInteger16(pvApiCtx, iVar, _list->data[0].uint16.piRows, _list->data[0].uint16.piCols, 
                                                        _list->data[0].uint16.uint16_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].uint16.uint16_data);
              FREE(_list->data);
            }
          break;
        case SCI_UINT32:
          if (in_list)
            {
              _SciErr = createMatrixOfUnsignedInteger32InList(pvApiCtx, iVar, parent, index_item,
                                                              _list->data[0].uint32.piRows, _list->data[0].uint32.piCols, 
                                                              _list->data[0].uint32.uint32_data);
            }
          else
            {
              _SciErr = createMatrixOfUnsignedInteger32(pvApiCtx, iVar, _list->data[0].uint32.piRows, _list->data[0].uint32.piCols, 
                                                        _list->data[0].uint32.uint32_data);
            }

          if (do_free)
            {
              FREE(_list->data[0].uint32.uint32_data);
              FREE(_list->data);
            }
          break;
        }
      break;
    case sci_strings:
      if (in_list)
        {
          _SciErr = createMatrixOfStringInList(pvApiCtx, iVar, parent, index_item, 
                                               _list->data[0].strings.piRows, _list->data[0].strings.piCols, 
                                               _list->data[0].strings.pstStrings);
        }
      else
        {
          _SciErr = createMatrixOfString(pvApiCtx, iVar, _list->data[0].strings.piRows, _list->data[0].strings.piCols, 
                                         _list->data[0].strings.pstStrings);
        }
      
      if (do_free)
        {
          FREE(_list->data[0].strings.pstStrings);
          FREE(_list->data);
        }
      break;
    case sci_list:
      _SciErr = createList(pvApiCtx, iVar, _list->nb_elems, &p_list_addr);

      for(i=0;i<_list->nb_elems;i++)
        {
          unserialize(iVar, &_list->data[i].element, do_free, 1, i+1, p_list_addr);
        }

      if (do_free)
        {
          FREE(_list->data);
        }
      break;
    case sci_tlist:
      _SciErr = createTList(pvApiCtx, iVar, _list->nb_elems, &p_list_addr);

      for(i=0;i<_list->nb_elems;i++)
        {
          unserialize(iVar, &_list->data[i].element, do_free, 1, i+1, p_list_addr);
        }

      if (do_free)
        {
          FREE(_list->data);
        }
      break;
    case sci_mlist:
      _SciErr = createMList(pvApiCtx, iVar, _list->nb_elems, &p_list_addr);

      for(i=0;i<_list->nb_elems;i++)
        {
          unserialize(iVar, &_list->data[i].element, do_free, 1, i+1, p_list_addr);
        }

      if (do_free)
        {
          FREE(_list->data);
        }
      break;
    case sci_pointer:
      if (in_list)
        {
          _SciErr = createPointerInList(pvApiCtx, iVar, parent, index_item, _list->data[0].pointer.ptr);
        }
      else
        {
          _SciErr = createPointer(pvApiCtx, iVar, _list->data[0].pointer.ptr);
        }

      if (do_free)
        {
          FREE(_list->data);
        }
      break;
    default:
      Scierror(999,"unserialize: unsupported data type - %d", _list->type);
      break;
    }
}
