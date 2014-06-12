#ifndef HELPER_H
#define HELPER_H

#include <stdio.h>

#define DBGPRINTF printf

int find_label_partial(char **LabelList, int nblabels, const char * LabelToFind);
int find_label(char **LabelList, int nblabels, const char * LabelToFind);

#ifdef DEBUG
#define INIT_PARAM(POS)							\
  int     m_tmp,   n_tmp,   l_tmp;					\
  int     m_param, n_param, l_param;					\
  int     m_label, n_label, pos_label;					\
  char ** LabelList = NULL;						\
  GetRhsVar(POS, MATRIX_ORIENTED_TYPED_LIST_DATATYPE, &m_param, &n_param, &l_param); \
  DBGPRINTF("PARAM_IN: m_param = %d n_param = %d\n", m_param, n_param);	\
  if ((m_param!=0)&&(n_param!=0))					\
    {									\
      GetListRhsVar(POS, 1, MATRIX_OF_STRING_DATATYPE, &m_label, &n_label, &LabelList); \
      DBGPRINTF("PARAM_IN: m_label = %d n_label = %d\n", m_label, n_label); \
    }									\
  else									\
    {									\
      m_label = 0;							\
      n_label = 0;							\
      LabelList = NULL;							\
    }
#else
#define INIT_PARAM(POS)							\
  int     m_tmp,   n_tmp,   l_tmp;					\
  int     m_param, n_param, l_param;					\
  int     m_label, n_label, pos_label;					\
  char ** LabelList = NULL;						\
  GetRhsVar(POS, MATRIX_ORIENTED_TYPED_LIST_DATATYPE, &m_param, &n_param, &l_param); \
  if ((m_param!=0)&&(n_param!=0))					\
    {									\
      GetListRhsVar(POS, 1, MATRIX_OF_STRING_DATATYPE, &m_label, &n_label, &LabelList); \
    }									\
  else									\
    {									\
      m_label = 0;							\
      n_label = 0;							\
      LabelList = NULL;							\
    }
#endif

#ifdef DEBUG
#define GET_PARAM_INT(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_INTEGER_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = *istk(l_tmp);						\
      DBGPRINTF("DEBUG: Int - Parameter %s = %d\n", NAME, RESULT);	\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
      DBGPRINTF("DEBUG: Int - Parameter %s = %d\n", NAME, RESULT);	\
    }

#define GET_PARAM_DOUBLE(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_DOUBLE_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = *stk(l_tmp);						\
      DBGPRINTF("DEBUG: Real - Parameter %s = %f\n", NAME, RESULT);	\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
      DBGPRINTF("DEBUG: Real - Parameter %s = %f\n", NAME, RESULT);	\
    }

#define GET_PARAM_VEC_DOUBLE(NAME, RESULT, SIZE, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_DOUBLE_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = stk(l_tmp);						\
      SIZE = m_tmp*n_tmp;						\
      DBGPRINTF("DEBUG: Real - Parameter %s != NULL\n", NAME);		\
    }									\
  else									\
    {									\
      RESULT = NULL;							\
      DBGPRINTF("DEBUG: Real - Parameter %s == NULL \n", NAME);		\
    }

#define GET_PARAM_STRING(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, "c", &m_tmp, &n_tmp, &l_tmp); \
      RESULT = cstk(l_tmp);						\
      DBGPRINTF("DEBUG: String - Parameter %s = %s\n", NAME, RESULT);	\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
      DBGPRINTF("DEBUG: String - Parameter %s = %s\n", NAME, RESULT);	\
    }
#else
#define GET_PARAM_INT(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_INTEGER_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = *istk(l_tmp);						\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
    }

#define GET_PARAM_DOUBLE(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_DOUBLE_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = *stk(l_tmp);						\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
    }

#define GET_PARAM_VEC_DOUBLE(NAME, RESULT, SIZE, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, MATRIX_OF_DOUBLE_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = stk(l_tmp);						\
      SIZE = m_tmp*n_tmp;						\
    }									\
  else									\
    {									\
      RESULT = NULL;							\
    }

#define GET_PARAM_STRING(NAME, RESULT, DEFAULT, FOUND)			\
  pos_label = find_label(LabelList, m_param, NAME);			\
  FOUND = pos_label;							\
  if (pos_label!=-1)							\
    {									\
      GetListRhsVar(PARAM_IN, pos_label+1, STRING_DATATYPE, &m_tmp, &n_tmp, &l_tmp); \
      RESULT = cstk(l_tmp);						\
    }									\
  else									\
    {									\
      RESULT = DEFAULT;							\
    }
#endif
#endif
