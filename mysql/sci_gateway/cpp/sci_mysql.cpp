extern "C" {

//#define USE_OLD_FUNCTIONS 1

#ifdef _MSC_VER
#include <Windows.h>
#include <mysql.h>
#include <my_global.h>
#include <my_sys.h>
#else
#include <mysql/my_global.h>
#include <mysql/my_sys.h>
#include <mysql/mysql.h>
#endif
}

#include <string.h>

extern "C" {
#include <stack-c.h>
#include <sciprint.h>
#include <Scierror.h>
#include <MALLOC.h>
#include <freeArrayOfString.h>
}

extern "C" int sci_my_init(char * fname)
{
  my_init();
  return 0;
}

extern "C" int sci_mysql_affected_rows(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_affected_rows = 1, n_affected_rows = 1, l_affected_rows;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs + 1, "i", &m_affected_rows, &n_affected_rows, &l_affected_rows);

  *istk(l_affected_rows) =  mysql_affected_rows(mysql_ptr);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_change_user(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_bool = 1, n_bool = 1, l_bool;
  int m_user, n_user, l_user;
  int m_password, n_password, l_password;
  int m_db, n_db, l_db;
  char * db = NULL;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(4,4);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_user,     &n_user,     &l_user);
  GetRhsVar(3, "c", &m_password, &n_password, &l_password);
  if (Rhs>=4)
    {
      GetRhsVar(4, "c", &m_db,       &n_db,       &l_db);
      if ((m_db!=0)&&(n_db!=0)) db = cstk(l_db);
      else                      db = NULL;
    }
  else
    {
      db = NULL;
    }

  CreateVar(Rhs+1,"i", &m_bool, &n_bool, &l_bool);

  *istk(l_bool) = mysql_change_user(mysql_ptr, cstk(l_user), cstk(l_password), db);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_character_set_name(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_charset = 1, n_charset = 1, l_charset;
  char * res_string = NULL;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  res_string = (char *)mysql_character_set_name(mysql_ptr);

  m_charset = 1;
  n_charset = strlen(res_string);

  CreateVar(Rhs+1,"c", &m_charset, &n_charset, &l_charset);
  strncpy(cstk(l_charset),res_string,strlen(res_string));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_close(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  if (mysql_ptr) mysql_close(mysql_ptr);

  return 0;
}

#ifdef USE_OLD_FUNCTIONS
extern "C" int sci_mysql_connect(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_sql_pointer_out = 1, n_sql_pointer_out = 1,  l_sql_pointer_out;
  int m_host, n_host, l_host;
  int m_user, n_user, l_user;
  int m_passwd, n_passwd, l_passwd;
  MYSQL * mysql_ptr_in = NULL;
  MYSQL * mysql_ptr_out = NULL;

  CheckRhs(4,4);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  if ((m_sql_pointer_in!=0)&&(n_sql_pointer_in!=0)) mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  else                                              mysql_ptr_in = NULL;

  GetRhsVar(2, "c", &m_host,   &n_host,   &l_host);
  GetRhsVar(3, "c", &m_user,   &n_user,   &l_user);
  GetRhsVar(4, "c", &m_passwd, &n_passwd, &l_passwd);

  mysql_ptr_out = mysql_connect(mysql_ptr_in, cstk(l_host), cstk(l_user), cstk(l_passwd));

  if (mysql_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_pointer_out, &n_sql_pointer_out, mysql_ptr_out);
    }
  else
    {
      m_sql_pointer_out = 0; n_sql_pointer_out = 0;
      CreateVar(Rhs + 1, "d", &m_sql_pointer_out, &n_sql_pointer_out, &l_sql_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}
#endif

#ifdef USE_OLD_FUNCTIONS
extern "C" int sci_mysql_create_db(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_db, n_db, l_db;
  int m_out = 1, n_out = 1, l_out;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_db, &n_db, &l_db);

  CreateVar(Rhs+1, "i", &m_out, &n_out, &l_out);

  *istk(l_out) =  mysql_create_db(mysql_ptr, cstk(l_db));

  LhsVar(1) = Rhs + 1;

  return 0;
}
#endif

extern "C" int sci_mysql_data_seek(char * fname)
{
  int m_sql_res_pointer_in,  n_sql_res_pointer_in,  l_sql_res_pointer_in;
  int m_offset, n_offset, l_offset;
  int m_out = 1, n_out = 1, l_out;
  MYSQL_RES * mysql_res_ptr = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_offset, &n_offset, &l_offset);

  mysql_data_seek(mysql_res_ptr, *istk(l_offset));
  
  return 0;
}

extern "C" int sci_mysql_debug(char * fname)
{
  int m_debug, n_debug, l_debug;

  GetRhsVar(1, "c", &m_debug, &n_debug, &l_debug);

  mysql_debug(cstk(l_debug));

  return 0;
}

#ifdef USE_OLD_FUNCTIONS
extern "C" int sci_mysql_drop_db(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_db, n_db, l_db;
  int m_out = 1, n_out = 1, l_out;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_db, &n_db, &l_db);
  CreateVar(Rhs+1, "i", &m_out, &n_out, &l_out);

  *istk(l_out) =  mysql_drop_db(mysql_ptr, cstk(l_db));

  LhsVar(1) = Rhs+1;

  return 0;
}
#endif

extern "C" int sci_mysql_dump_debug_info(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_out = 1, n_out = 1, l_out;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1,"i", &m_out, &n_out, &l_out);

  *istk(l_out) =  mysql_dump_debug_info(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_eof(char * fname)
{
  int m_sql_res_pointer_in,  n_sql_res_pointer_in,  l_sql_res_pointer_in;
  int m_bool = 1, n_bool = 1, l_bool;
  MYSQL_RES * mysql_res_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1,"i", &m_bool, &n_bool, &l_bool);

  *istk(l_bool) =  mysql_eof(mysql_res_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_errno(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_errno = 1, n_errno = 1, l_errno;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1,"i", &m_errno, &n_errno, &l_errno);

  *istk(l_errno) =  mysql_errno(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_error(char * fname)
{
  int m_sql_pointer_in,  n_sql_pointer_in,  l_sql_pointer_in;
  int m_error = 1, n_error = 1, l_error;
  MYSQL * mysql_ptr = NULL;
  char * result_string = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  result_string = (char *)mysql_error(mysql_ptr);

  if (result_string)
    {
      m_error = 1; 
      n_error = strlen(result_string);
      CreateVar(Rhs+1, "c", &m_error, &n_error, &l_error);
      
      strncpy(cstk(l_error),result_string,strlen(result_string));
    }
  else
    {
      m_error = 0; 
      n_error = 0;
      CreateVar(Rhs+1, "d", &m_error, &n_error, &l_error);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_fetch_field(char * fname)
{
  int m_sql_res_pointer_in,     n_sql_res_pointer_in,     l_sql_res_pointer_in;
  int m_extra       = 14, n_extra       = 1, l_extra;
  int m_list_labels = 14, n_list_labels = 1, l_list_list_labels, l_list_labels;
  int m_name        = 1,  n_name        = 1, l_list_name,        l_name;
  int m_org_name    = 1,  n_org_name    = 1, l_list_org_name,    l_org_name;
  int m_table       = 1,  n_table       = 1, l_list_table,       l_table;
  int m_org_table   = 1,  n_org_table   = 1, l_list_org_table,   l_org_table;
  int m_db          = 1,  n_db          = 1, l_list_db,          l_db;
  int m_catalog     = 1,  n_catalog     = 1, l_list_catalog,     l_catalog;
  int m_def         = 1,  n_def         = 1, l_list_def,         l_def;
  int m_length      = 1,  n_length      = 1, l_list_length,      l_length;
  int m_max_length  = 1,  n_max_length  = 1, l_list_max_length,  l_max_length;
  int m_flags       = 1,  n_flags       = 1, l_list_flags,       l_flags;
  int m_decimals    = 1,  n_decimals    = 1, l_list_decimals,    l_decimals;
  int m_charsetnr   = 1,  n_charsetnr   = 1, l_list_charsetnr,   l_charsetnr;
  int m_type        = 1,  n_type        = 1, l_list_type,        l_type;

  const char * FieldNames[] = {"mysql_field", 
			       "name",                 /* Name of column */
			       "org_name",             /* Original column name, if an alias */
			       "table",                /* Table of column if column was a field */
			       "org_table",            /* Org table name, if table was an alias */
			       "db",                   /* Database for table */
			       "catalog",              /* Catalog for table */
			       "def",                  /* Default value (set by mysql_list_fields) */
			       "length",               /* Width of column (create length) */
			       "max_length",           /* Max width for selected set */
			       "flags",                /* Div flags */
			       "decimals",             /* Number of decimals in field */
			       "charsetnr",            /* Character set */
			       "type"};                /* Type of field. See mysql_com.h for types */

  MYSQL_RES   * mysql_res_ptr   = NULL;
  MYSQL_FIELD * mysql_field_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  mysql_field_ptr = mysql_fetch_field(mysql_res_ptr);
  
  if (mysql_field_ptr)
    {
      if (mysql_field_ptr->name)
	{
	  m_name = 1; n_name = strlen(mysql_field_ptr->name);
	  CreateVar(Rhs+1, "c", &m_name, &n_name, &l_name);
	  strncpy(cstk(l_name), mysql_field_ptr->name, n_name);
	}
      else
	{
	  m_name = 0; n_name = 0;
	  CreateVar(Rhs+1, "d", &m_name, &n_name, &l_name);
	}

      if (mysql_field_ptr->org_name)
	{
	  m_org_name = 1; n_org_name = strlen(mysql_field_ptr->org_name);
	  CreateVar(Rhs+2, "c", &m_org_name, &n_org_name, &l_org_name);
	  strcpy(cstk(l_org_name),  mysql_field_ptr->org_name);
	}
      else
	{
	  m_org_name = 0; n_org_name = 0;
	  CreateVar(Rhs+2, "d", &m_org_name, &n_org_name, &l_org_name);
	}

      if (mysql_field_ptr->table)
	{
	  m_table = 1; n_table = strlen(mysql_field_ptr->table);
	  CreateVar(Rhs+3, "c", &m_table, &n_table, &l_table);
	  strcpy(cstk(l_table), mysql_field_ptr->table);
	}
      else
	{
	  m_table = 0; n_table = 0;
	  CreateVar(Rhs+3, "d", &m_table, &n_table, &l_table);
	}
	 
      if (mysql_field_ptr->org_table)
	{
	  m_org_table = 1; n_org_table = strlen(mysql_field_ptr->org_table);
	  CreateVar(Rhs+4, "c", &m_org_table, &n_org_table, &l_org_table);
	  strcpy(cstk(l_org_table), mysql_field_ptr->org_table);
	}
      else
	{
	  m_org_table = 0; n_org_table = 0;
	  CreateVar(Rhs+4, "d", &m_org_table, &n_org_table, &l_org_table);
	}

      if (mysql_field_ptr->db)
	{
	  m_db = 1; n_db = strlen(mysql_field_ptr->db);
	  CreateVar(Rhs+5, "c", &m_db, &n_db, &l_db);
	  strcpy(cstk(l_db), mysql_field_ptr->db);
	}
      else
	{
	  m_db = 0; n_db = 0;
	  CreateVar(Rhs+5, "d", &m_db, &n_db, &l_db);
	}

      if (mysql_field_ptr->catalog)
	{
	  m_catalog = 1; n_catalog = strlen(mysql_field_ptr->catalog);
	  CreateVar(Rhs+6, "c", &m_catalog, &n_catalog, &l_catalog);
	  strcpy(cstk(l_catalog), mysql_field_ptr->catalog);
	}
      else
	{
	  m_catalog = 0; n_catalog = 0;
	  CreateVar(Rhs+6, "d", &m_catalog, &n_catalog, &l_catalog);
	}

      if (mysql_field_ptr->def)
	{
	  m_def = 1; n_def = strlen(mysql_field_ptr->def);      
	  CreateVar(Rhs+7, "c", &m_def, &n_def, &l_def);
	  strcpy(cstk(l_def), mysql_field_ptr->def);
	}
      else
	{
	  m_def = 0; n_def = 0;
	  CreateVar(Rhs+7, "d", &m_def, &n_def, &l_def);
	}

      CreateVar(Rhs+8,  "i", &m_length,     &n_length,     &l_length);
      CreateVar(Rhs+9,  "i", &m_max_length, &n_max_length, &l_max_length);
      CreateVar(Rhs+10, "i", &m_flags,      &n_flags,      &l_flags);
      CreateVar(Rhs+11, "i", &m_decimals,   &n_decimals,   &l_decimals);
      CreateVar(Rhs+12, "i", &m_charsetnr,  &n_charsetnr,  &l_charsetnr);
      CreateVar(Rhs+13, "i", &m_type,       &n_type,       &l_type);
            
      *istk(l_length)     = mysql_field_ptr->length;
      *istk(l_max_length) = mysql_field_ptr->max_length;
      *istk(l_flags)      = mysql_field_ptr->flags;
      *istk(l_decimals)   = mysql_field_ptr->decimals;
      *istk(l_charsetnr)  = mysql_field_ptr->charsetnr;
      *istk(l_type)       = (int)mysql_field_ptr->type;
      
      CreateVar(Rhs+14, "m", &m_extra, &n_extra, &l_extra);
      CreateListVarFromPtr(Rhs+14, 1, "S",  &m_list_labels, &n_list_labels, FieldNames);

      if (mysql_field_ptr->name)
	{
	  CreateListVarFrom(Rhs+14, 2, "c", &n_name, &m_name, &l_list_name, &l_name);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 2, "d", &n_name, &m_name, &l_list_name, &l_name);
	}

      if (mysql_field_ptr->org_name)
	{
	  CreateListVarFrom(Rhs+14, 3, "c", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 3, "d", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
	}

      if (mysql_field_ptr->table)
	{
	  CreateListVarFrom(Rhs+14, 4, "c", &n_table, &m_table, &l_list_table, &l_table);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 4, "d", &n_table, &m_table, &l_list_table, &l_table);
	}

      if (mysql_field_ptr->org_table)
	{
	  CreateListVarFrom(Rhs+14, 5, "c", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 5, "d", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
	}

      if (mysql_field_ptr->db)
	{
	  CreateListVarFrom(Rhs+14, 6, "c", &n_db, &m_db, &l_list_db, &l_db);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 6, "d", &n_db, &m_db, &l_list_db, &l_db);
	}

      if (mysql_field_ptr->catalog)
	{
	  CreateListVarFrom(Rhs+14, 7, "c", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 7, "d", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
	}
	
      if (mysql_field_ptr->def)
	{
	  CreateListVarFrom(Rhs+14, 8, "c", &n_def, &m_def, &l_list_def, &l_def);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 8, "d", &n_def, &m_def, &l_list_def, &l_def);
	}

      CreateListVarFrom(Rhs+14, 9,  "i", &m_length,     &n_length,     &l_list_length,     &l_length);
      CreateListVarFrom(Rhs+14, 10, "i", &m_max_length, &n_max_length, &l_list_max_length, &l_max_length);
      CreateListVarFrom(Rhs+14, 11, "i", &m_flags,      &n_flags,      &l_list_flags,      &l_flags);
      CreateListVarFrom(Rhs+14, 12, "i", &m_decimals,   &n_decimals,   &l_list_decimals,   &l_decimals);
      CreateListVarFrom(Rhs+14, 13, "i", &m_charsetnr,  &n_charsetnr,  &l_list_charsetnr,  &l_charsetnr);
      CreateListVarFrom(Rhs+14, 14, "i", &m_type,       &n_type,       &l_list_type,       &l_type);
      
      LhsVar(1) = Rhs + 14;
    }
  else
    {
      m_extra = 0; n_extra = 0;

      CreateVar(Rhs + 1, "d", &m_extra, &n_extra, &l_extra);

      LhsVar(1) = Rhs + 1;
    }

  return 0;
}

// YC: DOC de cette fonction a faire
// extern "C" int sci_mysql_fetch_fields(char * fname)
// {
//   int m_sql_res_pointer_in,     n_sql_res_pointer_in,     l_sql_res_pointer_in;
//   int m_extra       = 14, n_extra       = 1, l_extra;
//   int m_list_labels = 14, n_list_labels = 1, l_list_list_labels, l_list_labels;
//   int m_name        = 1,  n_name        = 1, l_list_name,        l_name;
//   int m_org_name    = 1,  n_org_name    = 1, l_list_org_name,    l_org_name;
//   int m_table       = 1,  n_table       = 1, l_list_table,       l_table;
//   int m_org_table   = 1,  n_org_table   = 1, l_list_org_table,   l_org_table;
//   int m_db          = 1,  n_db          = 1, l_list_db,          l_db;
//   int m_catalog     = 1,  n_catalog     = 1, l_list_catalog,     l_catalog;
//   int m_def         = 1,  n_def         = 1, l_list_def,         l_def;
//   int m_length      = 1,  n_length      = 1, l_list_length,      l_length;
//   int m_max_length  = 1,  n_max_length  = 1, l_list_max_length,  l_max_length;
//   int m_flags       = 1,  n_flags       = 1, l_list_flags,       l_flags;
//   int m_decimals    = 1,  n_decimals    = 1, l_list_decimals,    l_decimals;
//   int m_charsetnr   = 1,  n_charsetnr   = 1, l_list_charsetnr,   l_charsetnr;
//   int m_type        = 1,  n_type        = 1, l_list_type,        l_type;

//   const char * FieldNames[] = {"mysql_field", 
// 			       "name",                 /* Name of column */
// 			       "org_name",             /* Original column name, if an alias */
// 			       "table",                /* Table of column if column was a field */
// 			       "org_table",            /* Org table name, if table was an alias */
// 			       "db",                   /* Database for table */
// 			       "catalog",              /* Catalog for table */
// 			       "def",                  /* Default value (set by mysql_list_fields) */
// 			       "length",               /* Width of column (create length) */
// 			       "max_length",           /* Max width for selected set */
// 			       "flags",                /* Div flags */
// 			       "decimals",             /* Number of decimals in field */
// 			       "charsetnr",            /* Character set */
// 			       "type"};                /* Type of field. See mysql_com.h for types */

//   MYSQL_RES   * mysql_res_ptr   = NULL;
//   MYSQL_FIELD * mysql_field_ptr = NULL;

//   CheckRhs(1,1);
//   CheckLhs(0,1);

//   GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
//   mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

//   if (mysql_res_ptr==NULL)
//     {
//       Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
//       return 0;
//     }

//   mysql_field_ptr = mysql_fetch_fields(mysql_res_ptr);

//   if (mysql_field_ptr)
//     {
//       if (mysql_field_ptr->name)
// 	{
// 	  m_name = 1; n_name = strlen(mysql_field_ptr->name);
// 	  CreateVar(Rhs+1, "c", &m_name, &n_name, &l_name);
// 	  strncpy(cstk(l_name), mysql_field_ptr->name, n_name);
// 	}
//       else
// 	{
// 	  m_name = 0; n_name = 0;
// 	  CreateVar(Rhs+1, "d", &m_name, &n_name, &l_name);
// 	}

//       if (mysql_field_ptr->org_name)
// 	{
// 	  m_org_name = 1; n_org_name = strlen(mysql_field_ptr->org_name);
// 	  CreateVar(Rhs+2, "c", &m_org_name, &n_org_name, &l_org_name);
// 	  strcpy(cstk(l_org_name),  mysql_field_ptr->org_name);
// 	}
//       else
// 	{
// 	  m_org_name = 0; n_org_name = 0;
// 	  CreateVar(Rhs+2, "d", &m_org_name, &n_org_name, &l_org_name);
// 	}

//       if (mysql_field_ptr->table)
// 	{
// 	  m_table = 1; n_table = strlen(mysql_field_ptr->table);
// 	  CreateVar(Rhs+3, "c", &m_table, &n_table, &l_table);
// 	  strcpy(cstk(l_table), mysql_field_ptr->table);
// 	}
//       else
// 	{
// 	  m_table = 0; n_table = 0;
// 	  CreateVar(Rhs+3, "d", &m_table, &n_table, &l_table);
// 	}
	 
//       if (mysql_field_ptr->org_table)
// 	{
// 	  m_org_table = 1; n_org_table = strlen(mysql_field_ptr->org_table);
// 	  CreateVar(Rhs+4, "c", &m_org_table, &n_org_table, &l_org_table);
// 	  strcpy(cstk(l_org_table), mysql_field_ptr->org_table);
// 	}
//       else
// 	{
// 	  m_org_table = 0; n_org_table = 0;
// 	  CreateVar(Rhs+4, "d", &m_org_table, &n_org_table, &l_org_table);
// 	}

//       if (mysql_field_ptr->db)
// 	{
// 	  m_db = 1; n_db = strlen(mysql_field_ptr->db);
// 	  CreateVar(Rhs+5, "c", &m_db, &n_db, &l_db);
// 	  strcpy(cstk(l_db), mysql_field_ptr->db);
// 	}
//       else
// 	{
// 	  m_db = 0; n_db = 0;
// 	  CreateVar(Rhs+5, "d", &m_db, &n_db, &l_db);
// 	}

//       if (mysql_field_ptr->catalog)
// 	{
// 	  m_catalog = 1; n_catalog = strlen(mysql_field_ptr->catalog);
// 	  CreateVar(Rhs+6, "c", &m_catalog, &n_catalog, &l_catalog);
// 	  strcpy(cstk(l_catalog), mysql_field_ptr->catalog);
// 	}
//       else
// 	{
// 	  m_catalog = 0; n_catalog = 0;
// 	  CreateVar(Rhs+6, "d", &m_catalog, &n_catalog, &l_catalog);
// 	}

//       if (mysql_field_ptr->def)
// 	{
// 	  m_def = 1; n_def = strlen(mysql_field_ptr->def);      
// 	  CreateVar(Rhs+7, "c", &m_def, &n_def, &l_def);
// 	  strcpy(cstk(l_def), mysql_field_ptr->def);
// 	}
//       else
// 	{
// 	  m_def = 0; n_def = 0;
// 	  CreateVar(Rhs+7, "d", &m_def, &n_def, &l_def);
// 	}

//       CreateVar(Rhs+8,  "i", &m_length,     &n_length,     &l_length);
//       CreateVar(Rhs+9,  "i", &m_max_length, &n_max_length, &l_max_length);
//       CreateVar(Rhs+10, "i", &m_flags,      &n_flags,      &l_flags);
//       CreateVar(Rhs+11, "i", &m_decimals,   &n_decimals,   &l_decimals);
//       CreateVar(Rhs+12, "i", &m_charsetnr,  &n_charsetnr,  &l_charsetnr);
//       CreateVar(Rhs+13, "i", &m_type,       &n_type,       &l_type);
            
//       *istk(l_length)     = mysql_field_ptr->length;
//       *istk(l_max_length) = mysql_field_ptr->max_length;
//       *istk(l_flags)      = mysql_field_ptr->flags;
//       *istk(l_decimals)   = mysql_field_ptr->decimals;
//       *istk(l_charsetnr)  = mysql_field_ptr->charsetnr;
//       *istk(l_type)       = (int)mysql_field_ptr->type;
      
//       CreateVar(Rhs+14, "m", &m_extra, &n_extra, &l_extra);
//       CreateListVarFromPtr(Rhs+14, 1, "S",  &m_list_labels, &n_list_labels, FieldNames);

//       if (mysql_field_ptr->name)
// 	{
// 	  CreateListVarFrom(Rhs+14, 2, "c", &n_name, &m_name, &l_list_name, &l_name);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 2, "d", &n_name, &m_name, &l_list_name, &l_name);
// 	}

//       if (mysql_field_ptr->org_name)
// 	{
// 	  CreateListVarFrom(Rhs+14, 3, "c", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 3, "d", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
// 	}

//       if (mysql_field_ptr->table)
// 	{
// 	  CreateListVarFrom(Rhs+14, 4, "c", &n_table, &m_table, &l_list_table, &l_table);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 4, "d", &n_table, &m_table, &l_list_table, &l_table);
// 	}

//       if (mysql_field_ptr->org_table)
// 	{
// 	  CreateListVarFrom(Rhs+14, 5, "c", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 5, "d", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
// 	}

//       if (mysql_field_ptr->db)
// 	{
// 	  CreateListVarFrom(Rhs+14, 6, "c", &n_db, &m_db, &l_list_db, &l_db);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 6, "d", &n_db, &m_db, &l_list_db, &l_db);
// 	}

//       if (mysql_field_ptr->catalog)
// 	{
// 	  CreateListVarFrom(Rhs+14, 7, "c", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 7, "d", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
// 	}
	
//       if (mysql_field_ptr->def)
// 	{
// 	  CreateListVarFrom(Rhs+14, 8, "c", &n_def, &m_def, &l_list_def, &l_def);
// 	}
//       else
// 	{
// 	  CreateListVarFrom(Rhs+14, 8, "d", &n_def, &m_def, &l_list_def, &l_def);
// 	}

//       CreateListVarFrom(Rhs+14, 9,  "i", &m_length,     &n_length,     &l_list_length,     &l_length);
//       CreateListVarFrom(Rhs+14, 10, "i", &m_max_length, &n_max_length, &l_list_max_length, &l_max_length);
//       CreateListVarFrom(Rhs+14, 11, "i", &m_flags,      &n_flags,      &l_list_flags,      &l_flags);
//       CreateListVarFrom(Rhs+14, 12, "i", &m_decimals,   &n_decimals,   &l_list_decimals,   &l_decimals);
//       CreateListVarFrom(Rhs+14, 13, "i", &m_charsetnr,  &n_charsetnr,  &l_list_charsetnr,  &l_charsetnr);
//       CreateListVarFrom(Rhs+14, 14, "i", &m_type,       &n_type,       &l_list_type,       &l_type);
      
//       LhsVar(1) = Rhs + 14;
//     }
//   else
//     {
//       m_extra = 0; n_extra = 0;

//       CreateVar(Rhs + 1, "d", &m_extra, &n_extra, &l_extra);

//       LhsVar(1) = Rhs + 1;
//     }

//   return 0;
// }

extern "C" int sci_mysql_fetch_field_direct(char * fname)
{
  int m_sql_res_pointer_in,     n_sql_res_pointer_in,     l_sql_res_pointer_in;
  int m_fieldnr,                n_fieldnr,                l_fieldnr;

  int m_extra       = 14, n_extra       = 1, l_extra;
  int m_list_labels = 14, n_list_labels = 1, l_list_list_labels, l_list_labels;
  int m_name        = 1,  n_name        = 1, l_list_name,        l_name;
  int m_org_name    = 1,  n_org_name    = 1, l_list_org_name,    l_org_name;
  int m_table       = 1,  n_table       = 1, l_list_table,       l_table;
  int m_org_table   = 1,  n_org_table   = 1, l_list_org_table,   l_org_table;
  int m_db          = 1,  n_db          = 1, l_list_db,          l_db;
  int m_catalog     = 1,  n_catalog     = 1, l_list_catalog,     l_catalog;
  int m_def         = 1,  n_def         = 1, l_list_def,         l_def;
  int m_length      = 1,  n_length      = 1, l_list_length,      l_length;
  int m_max_length  = 1,  n_max_length  = 1, l_list_max_length,  l_max_length;
  int m_flags       = 1,  n_flags       = 1, l_list_flags,       l_flags;
  int m_decimals    = 1,  n_decimals    = 1, l_list_decimals,    l_decimals;
  int m_charsetnr   = 1,  n_charsetnr   = 1, l_list_charsetnr,   l_charsetnr;
  int m_type        = 1,  n_type        = 1, l_list_type,        l_type;

  const char * FieldNames[] = {"mysql_field", 
			       "name",                 /* Name of column */
			       "org_name",             /* Original column name, if an alias */
			       "table",                /* Table of column if column was a field */
			       "org_table",            /* Org table name, if table was an alias */
			       "db",                   /* Database for table */
			       "catalog",              /* Catalog for table */
			       "def",                  /* Default value (set by mysql_list_fields) */
			       "length",               /* Width of column (create length) */
			       "max_length",           /* Max width for selected set */
			       "flags",                /* Div flags */
			       "decimals",             /* Number of decimals in field */
			       "charsetnr",            /* Character set */
			       "type"};                /* Type of field. See mysql_com.h for types */

  MYSQL_RES   * mysql_res_ptr   = NULL;
  MYSQL_FIELD * mysql_field_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_fieldnr, &n_fieldnr, &l_fieldnr);

  mysql_field_ptr = mysql_fetch_field_direct(mysql_res_ptr, *istk(l_fieldnr));  

  if (mysql_field_ptr)
    {
      if (mysql_field_ptr->name)
	{
	  m_name = 1; n_name = strlen(mysql_field_ptr->name);
	  CreateVar(Rhs+1, "c", &m_name, &n_name, &l_name);
	  strncpy(cstk(l_name), mysql_field_ptr->name, n_name);
	}
      else
	{
	  m_name = 0; n_name = 0;
	  CreateVar(Rhs+1, "d", &m_name, &n_name, &l_name);
	}

      if (mysql_field_ptr->org_name)
	{
	  m_org_name = 1; n_org_name = strlen(mysql_field_ptr->org_name);
	  CreateVar(Rhs+2, "c", &m_org_name, &n_org_name, &l_org_name);
	  strcpy(cstk(l_org_name),  mysql_field_ptr->org_name);
	}
      else
	{
	  m_org_name = 0; n_org_name = 0;
	  CreateVar(Rhs+2, "d", &m_org_name, &n_org_name, &l_org_name);
	}

      if (mysql_field_ptr->table)
	{
	  m_table = 1; n_table = strlen(mysql_field_ptr->table);
	  CreateVar(Rhs+3, "c", &m_table, &n_table, &l_table);
	  strcpy(cstk(l_table), mysql_field_ptr->table);
	}
      else
	{
	  m_table = 0; n_table = 0;
	  CreateVar(Rhs+3, "d", &m_table, &n_table, &l_table);
	}
	 
      if (mysql_field_ptr->org_table)
	{
	  m_org_table = 1; n_org_table = strlen(mysql_field_ptr->org_table);
	  CreateVar(Rhs+4, "c", &m_org_table, &n_org_table, &l_org_table);
	  strcpy(cstk(l_org_table), mysql_field_ptr->org_table);
	}
      else
	{
	  m_org_table = 0; n_org_table = 0;
	  CreateVar(Rhs+4, "d", &m_org_table, &n_org_table, &l_org_table);
	}

      if (mysql_field_ptr->db)
	{
	  m_db = 1; n_db = strlen(mysql_field_ptr->db);
	  CreateVar(Rhs+5, "c", &m_db, &n_db, &l_db);
	  strcpy(cstk(l_db), mysql_field_ptr->db);
	}
      else
	{
	  m_db = 0; n_db = 0;
	  CreateVar(Rhs+5, "d", &m_db, &n_db, &l_db);
	}

      if (mysql_field_ptr->catalog)
	{
	  m_catalog = 1; n_catalog = strlen(mysql_field_ptr->catalog);
	  CreateVar(Rhs+6, "c", &m_catalog, &n_catalog, &l_catalog);
	  strcpy(cstk(l_catalog), mysql_field_ptr->catalog);
	}
      else
	{
	  m_catalog = 0; n_catalog = 0;
	  CreateVar(Rhs+6, "d", &m_catalog, &n_catalog, &l_catalog);
	}

      if (mysql_field_ptr->def)
	{
	  m_def = 1; n_def = strlen(mysql_field_ptr->def);      
	  CreateVar(Rhs+7, "c", &m_def, &n_def, &l_def);
	  strcpy(cstk(l_def), mysql_field_ptr->def);
	}
      else
	{
	  m_def = 0; n_def = 0;
	  CreateVar(Rhs+7, "d", &m_def, &n_def, &l_def);
	}

      CreateVar(Rhs+8,  "i", &m_length,     &n_length,     &l_length);
      CreateVar(Rhs+9,  "i", &m_max_length, &n_max_length, &l_max_length);
      CreateVar(Rhs+10, "i", &m_flags,      &n_flags,      &l_flags);
      CreateVar(Rhs+11, "i", &m_decimals,   &n_decimals,   &l_decimals);
      CreateVar(Rhs+12, "i", &m_charsetnr,  &n_charsetnr,  &l_charsetnr);
      CreateVar(Rhs+13, "i", &m_type,       &n_type,       &l_type);
            
      *istk(l_length)     = mysql_field_ptr->length;
      *istk(l_max_length) = mysql_field_ptr->max_length;
      *istk(l_flags)      = mysql_field_ptr->flags;
      *istk(l_decimals)   = mysql_field_ptr->decimals;
      *istk(l_charsetnr)  = mysql_field_ptr->charsetnr;
      *istk(l_type)       = (int)mysql_field_ptr->type;
      
      CreateVar(Rhs+14, "m", &m_extra, &n_extra, &l_extra);
      CreateListVarFromPtr(Rhs+14, 1, "S",  &m_list_labels, &n_list_labels, FieldNames);

      if (mysql_field_ptr->name)
	{
	  CreateListVarFrom(Rhs+14, 2, "c", &n_name, &m_name, &l_list_name, &l_name);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 2, "d", &n_name, &m_name, &l_list_name, &l_name);
	}

      if (mysql_field_ptr->org_name)
	{
	  CreateListVarFrom(Rhs+14, 3, "c", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 3, "d", &n_org_name, &m_org_name, &l_list_org_name, &l_org_name);
	}

      if (mysql_field_ptr->table)
	{
	  CreateListVarFrom(Rhs+14, 4, "c", &n_table, &m_table, &l_list_table, &l_table);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 4, "d", &n_table, &m_table, &l_list_table, &l_table);
	}

      if (mysql_field_ptr->org_table)
	{
	  CreateListVarFrom(Rhs+14, 5, "c", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 5, "d", &n_org_table, &m_org_table, &l_list_org_table, &l_org_table);
	}

      if (mysql_field_ptr->db)
	{
	  CreateListVarFrom(Rhs+14, 6, "c", &n_db, &m_db, &l_list_db, &l_db);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 6, "d", &n_db, &m_db, &l_list_db, &l_db);
	}

      if (mysql_field_ptr->catalog)
	{
	  CreateListVarFrom(Rhs+14, 7, "c", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 7, "d", &n_catalog, &m_catalog, &l_list_catalog, &l_catalog);
	}
	
      if (mysql_field_ptr->def)
	{
	  CreateListVarFrom(Rhs+14, 8, "c", &n_def, &m_def, &l_list_def, &l_def);
	}
      else
	{
	  CreateListVarFrom(Rhs+14, 8, "d", &n_def, &m_def, &l_list_def, &l_def);
	}

      CreateListVarFrom(Rhs+14, 9,  "i", &m_length,     &n_length,     &l_list_length,     &l_length);
      CreateListVarFrom(Rhs+14, 10, "i", &m_max_length, &n_max_length, &l_list_max_length, &l_max_length);
      CreateListVarFrom(Rhs+14, 11, "i", &m_flags,      &n_flags,      &l_list_flags,      &l_flags);
      CreateListVarFrom(Rhs+14, 12, "i", &m_decimals,   &n_decimals,   &l_list_decimals,   &l_decimals);
      CreateListVarFrom(Rhs+14, 13, "i", &m_charsetnr,  &n_charsetnr,  &l_list_charsetnr,  &l_charsetnr);
      CreateListVarFrom(Rhs+14, 14, "i", &m_type,       &n_type,       &l_list_type,       &l_type);
      
      LhsVar(1) = Rhs + 14;
    }
  else
    {
      m_extra = 0; n_extra = 0;

      CreateVar(Rhs + 1, "d", &m_extra, &n_extra, &l_extra);

      LhsVar(1) = Rhs + 1;
    }

  return 0;
}

extern "C" int sci_mysql_fetch_lengths(char * fname)
{
  int m_sql_res_pointer_in, n_sql_res_pointer_in, l_sql_res_pointer_in;
  int m_lengths_out,  n_lengths_out,  l_lengths_out;
  MYSQL_RES   * mysql_res_ptr   = NULL;
  unsigned long * lengths = NULL;
  unsigned int num_fields = 0;
  int i;

  CheckRhs(1,1);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  lengths = mysql_fetch_lengths(mysql_res_ptr);
  if (lengths)
    {
      num_fields = mysql_num_fields(mysql_res_ptr);
      
      m_lengths_out = num_fields;
      n_lengths_out = 1;
      
      CreateVar(Rhs+1, "i", &m_lengths_out, &n_lengths_out, &l_lengths_out);
      for(i=0;i<num_fields;i++)
	{
	  *istk(l_lengths_out + i) = lengths[i];
	}
    }
  else
    {
      m_lengths_out = 0;
      n_lengths_out = 0;
      
      CreateVar(Rhs+1, "d", &m_lengths_out, &n_lengths_out, &l_lengths_out);
    }

  LhsVar(1) = Rhs+1;
      
  return 0;
}

extern "C" int sci_mysql_fetch_row(char * fname)
{
  int m_sql_res_pointer_in,   n_sql_res_pointer_in,   l_sql_res_pointer_in;
  int m_sql_row_pointer_out,  n_sql_row_pointer_out,  l_sql_row_pointer_out;
  MYSQL_RES * mysql_res_ptr = NULL;
  MYSQL_ROW mysql_row_ptr; // It's a char **
  char ** res_copy = NULL;
  unsigned long * lengths = NULL;
  int i;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  mysql_row_ptr = mysql_fetch_row(mysql_res_ptr);
  lengths = mysql_fetch_lengths(mysql_res_ptr);

  sciprint("DEBUG: mysql_row_ptr = %ld, lengths = %ld\n", mysql_row_ptr, lengths);

  if (mysql_row_ptr)
    {
      //YC: à vérifier !!
      n_sql_row_pointer_out = 1;
      m_sql_row_pointer_out = mysql_num_fields(mysql_res_ptr);
      
      res_copy = (char **)MALLOC(m_sql_row_pointer_out*sizeof(char *));
      sciprint("m_sql_row_pointer_out = %d\n", m_sql_row_pointer_out);
      
      for(i=0; i<m_sql_row_pointer_out;i++)
	{
	  sciprint("DEBUG: orig = !%s!, length = %d\n", mysql_row_ptr[i],lengths[i]);
	  res_copy[i] = strndup(mysql_row_ptr[i],lengths[i]);
	}
      
      CreateVarFromPtr(Rhs+1,"S",&m_sql_row_pointer_out, &n_sql_row_pointer_out, res_copy);

      freeArrayOfString(res_copy, m_sql_row_pointer_out);
    }
  else
    {
      n_sql_row_pointer_out = 0;
      m_sql_row_pointer_out = 0;

      CreateVar(Rhs+1,"d",&m_sql_row_pointer_out, &n_sql_row_pointer_out, &l_sql_row_pointer_out);
    }
  LhsVar(1) = Rhs+1;


  return 0;
}

extern "C" int sci_mysql_field_count(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_field_count_out = 1,  n_field_count_out = 1,  l_field_count_out;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_field_count_out, &n_field_count_out, &l_field_count_out);

  *istk(l_field_count_out) =  mysql_field_count(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_field_seek(char * fname)
{
  int m_sql_res_pointer_in,   n_sql_res_pointer_in,   l_sql_res_pointer_in;
  int m_offset_in,  n_offset_in,  l_offset_in;
  int m_offset_out = 1,  n_offset_out = 1,  l_offset_out;
  MYSQL_RES * mysql_res_ptr = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_offset_in, &n_offset_in, &l_offset_in);

  CreateVar(Rhs+1, "i", &m_offset_out, &n_offset_out, &l_offset_out);

  *istk(l_offset_out) =  mysql_field_seek(mysql_res_ptr, *istk(l_offset_in));

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_field_tell(char * fname)
{
  int m_sql_res_pointer_in,   n_sql_res_pointer_in,   l_sql_res_pointer_in;
  int m_offset_out = 1,  n_offset_out = 1,  l_offset_out;
  MYSQL_RES * mysql_res_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_offset_out, &n_offset_out, &l_offset_out);

  *istk(l_offset_out) =  mysql_field_tell(mysql_res_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_free_result(char * fname)
{
  int m_sql_res_pointer_in,   n_sql_res_pointer_in,   l_sql_res_pointer_in;
  MYSQL_RES * mysql_res_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  if (mysql_res_ptr) mysql_free_result(mysql_res_ptr);

  return 0;
}

extern "C" int sci_mysql_get_client_info(char * fname)
{
  int m_client_info_out, n_client_info_out, l_client_info_out;
  char * result_string = NULL;

  result_string = (char *)mysql_get_client_info();

  m_client_info_out = 1;
  n_client_info_out = strlen(result_string);

  CreateVar(1, "c", &m_client_info_out, &n_client_info_out, &l_client_info_out);
  strncpy(cstk(l_client_info_out),result_string, strlen(result_string));

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_get_client_version(char * fname)
{
  int m_client_version_out = 1, n_client_version_out = 1, l_client_version_out;

  CreateVar(1, "i", &m_client_version_out, &n_client_version_out, &l_client_version_out);

  *istk(l_client_version_out) = mysql_get_client_version();

  LhsVar(1) = 1;

  return 0;
}

extern "C" int sci_mysql_get_host_info(char * fname)
{
  int m_host_info_out, n_host_info_out, l_host_info_out;
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;
  char * result_string = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  result_string = (char *)mysql_get_host_info(mysql_ptr);

  m_host_info_out = 1;
  n_host_info_out = strlen(result_string);

  CreateVar(Rhs+1, "c", &m_host_info_out, &n_host_info_out, &l_host_info_out);

  strncpy(cstk(l_host_info_out),result_string,strlen(result_string));

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_get_proto_info(char * fname)
{
  int m_proto_info_out = 1, n_proto_info_out = 1, l_proto_info_out;
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1,"i", &m_proto_info_out, &n_proto_info_out, &l_proto_info_out);

  *istk(l_proto_info_out) =  mysql_get_proto_info(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_get_server_info(char * fname)
{
  int m_server_info_out, n_server_info_out, l_server_info_out;
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;
  char * result_string = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  result_string = (char *)mysql_get_server_info(mysql_ptr);

  m_server_info_out = 1;
  n_server_info_out = strlen(result_string);

  CreateVar(Rhs+1,"c", &m_server_info_out, &n_server_info_out, &l_server_info_out);
  strncpy(cstk(l_server_info_out),result_string, strlen(result_string));

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_get_server_version(char * fname)
{
  int m_server_version_out = 1, n_server_version_out = 1, l_server_version_out;
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_server_version_out, &n_server_version_out, &l_server_version_out);

  *istk(l_server_version_out) =  mysql_get_server_version(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_hex_string(char * fname)
{
  int m_to_out, n_to_out, l_to_out;
  int m_from_in, n_from_in, l_from_in;
  int m_length_in, n_length_in, l_length_in;
  int m_res_out = 1, n_res_out = 1, l_res_out;
  char * to = NULL;
  int res_out = 0;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "c", &m_from_in,   &n_from_in,   &l_from_in);
  GetRhsVar(2, "c", &m_length_in, &n_length_in, &l_length_in);

  to = (char *)MALLOC((2*(*istk(l_length_in))+1)*sizeof(char));

  res_out = mysql_hex_string(to, cstk(l_from_in), *istk(l_length_in));

  m_to_out = 1;
  n_to_out = strlen(to);
  CreateVar(Rhs+1, "c", &m_to_out, &n_to_out, &l_to_out);
  strncpy(cstk(l_to_out),to,strlen(to));

  FREE(to);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_info(char * fname)
{
  int m_mysql_info_out, n_mysql_info_out, l_mysql_info_out;
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  MYSQL * mysql_ptr = NULL;
  char * result_string = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  result_string = (char *)mysql_info(mysql_ptr);

  if (result_string!=NULL)
    {
      m_mysql_info_out = 1;
      n_mysql_info_out = strlen(result_string);
      
      CreateVar(Rhs+1, "c", &m_mysql_info_out, &n_mysql_info_out, &l_mysql_info_out);
      strncpy(cstk(l_mysql_info_out),result_string, strlen(result_string));
    }
  else
    {
      m_mysql_info_out = 0;
      n_mysql_info_out = 0;

      CreateVar(Rhs+1, "d", &m_mysql_info_out, &n_mysql_info_out, &l_mysql_info_out);
    }
  
  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_init(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_sql_pointer_out = 1,   n_sql_pointer_out = 1,   l_sql_pointer_out;
  MYSQL * mysql_ptr = NULL;
  MYSQL * mysql_ptr_out = NULL;

  CheckRhs(0,1);
  CheckLhs(1,1);

  if (Rhs>=1)
    {
      GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
      
      if ((m_sql_pointer_in!=0)&&(n_sql_pointer_in!=0)) mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
      else                                              mysql_ptr = NULL;
    }
  else
    {
      mysql_ptr = NULL;
    }

  mysql_ptr_out = mysql_init(mysql_ptr);

  if (mysql_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_pointer_out, &n_sql_pointer_out, mysql_ptr_out);
    }
  else
    {
      m_sql_pointer_out = 0;
      n_sql_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_pointer_out, &n_sql_pointer_out, l_sql_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_insert_id(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_id = 1, n_id = 1, l_id;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_id, &n_id, l_id);

  *istk(l_id) =  mysql_insert_id(mysql_ptr);

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_kill(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_pid, n_pid, l_pid;
  int m_out = 1, n_out = 1, l_out = 1;
  MYSQL * mysql_ptr = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_pid, &n_pid, &l_pid);

  CreateVar(Rhs+1, "i", &m_out, &n_out, l_out);

  *istk(l_out) =  mysql_kill(mysql_ptr, *istk(l_pid));

  LhsVar(1) = Rhs + 1;

  return 0;
}

// YC: ajout de parametres a mysql_library_init
extern "C" int sci_mysql_library_init(char * fname)
{
  if (mysql_library_init(0, NULL, NULL))
    {
      Scierror(999, "could not initialize MySQL library\n");
      return 0;
    }
  return 0;
}

extern "C" int sci_mysql_library_end(char * fname)
{
  mysql_library_end();

  return 0;
}

extern "C" int sci_mysql_list_dbs(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_sql_res_pointer_out = 1,   n_sql_res_pointer_out = 1,   l_sql_res_pointer_out;
  int m_wild, n_wild, l_wild;
  MYSQL * mysql_ptr_in  = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;
  char * wild = NULL;

  CheckRhs(1,2);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  if (Rhs>=2)
    {
      GetRhsVar(2, "c", &m_wild, &n_wild, &l_wild);
      wild = cstk(l_wild);
    }
  else
    {
      wild = NULL;
    }

  mysql_res_ptr_out = mysql_list_dbs(mysql_ptr_in, wild);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0;
      n_sql_res_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_list_fields(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_sql_res_pointer_out = 1,   n_sql_res_pointer_out = 1,   l_sql_res_pointer_out;
  int m_wild, n_wild, l_wild;
  int m_table, n_table, l_table;
  MYSQL * mysql_ptr = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;
  char * wild = NULL;

  CheckRhs(2,3);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_table, &n_table, &l_table);
  if (Rhs>=3)
    {
      GetRhsVar(3, "c", &m_wild, &n_wild, &l_wild);
      wild = cstk(l_wild);
    }
  else
    {
      wild = NULL;
    }

  mysql_res_ptr_out = mysql_list_fields(mysql_ptr, cstk(l_table), wild);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0;
      n_sql_res_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, &l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_list_processes(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_sql_res_pointer_out = 1,   n_sql_res_pointer_out = 1,   l_sql_res_pointer_out;
  MYSQL * mysql_ptr = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;

  CheckRhs(1,1);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  mysql_res_ptr_out = mysql_list_processes(mysql_ptr);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0;
      n_sql_res_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, &l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_list_tables(char * fname)
{
  int m_sql_pointer_in,   n_sql_pointer_in,   l_sql_pointer_in;
  int m_sql_res_pointer_out = 1,   n_sql_res_pointer_out = 1,   l_sql_res_pointer_out;
  int m_wild, n_wild, l_wild;
  MYSQL * mysql_ptr = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;
  char * wild = NULL;

  CheckRhs(1,2);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  if (Rhs>=2)
    {
      GetRhsVar(2, "c", &m_wild, &n_wild, &l_wild);
      wild = cstk(l_wild);
    }
  else
    {
      wild = NULL;
    }

  mysql_res_ptr_out = mysql_list_tables(mysql_ptr, wild);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0; n_sql_res_pointer_out = 0;
      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, &l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_num_fields(char * fname)
{
  int m_sql_res_pointer_in, n_sql_res_pointer_in, l_sql_res_pointer_in;
  int m_num_fields = 1, n_num_fields = 1, l_num_fields;
  MYSQL_RES * mysql_res_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr_in = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs + 1, "i", &m_num_fields, &n_num_fields, &l_num_fields);

  *istk(l_num_fields) =  mysql_num_fields(mysql_res_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_num_rows(char * fname)
{
  int m_sql_res_pointer_in, n_sql_res_pointer_in, l_sql_res_pointer_in;
  int m_num_rows = 1, n_num_rows = 1, l_num_rows;
  MYSQL_RES * mysql_res_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr_in = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs + 1, "i", &m_num_rows, &n_num_rows, &l_num_rows);

  *istk(l_sql_res_pointer_in) = mysql_num_rows(mysql_res_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_options(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_option, n_option, l_option;
  int m_arg, n_arg, l_arg;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(3,3);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_option, &n_option, &l_option);

  CreateVar(Rhs + 1, "i", &m_status, &n_status, &l_status);

  switch((enum mysql_option)*istk(l_option))
    {
    case MYSQL_INIT_COMMAND:
      // MYSQL_INIT_COMMAND (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_OPT_COMPRESS:
      // MYSQL_OPT_COMPRESS (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_CONNECT_TIMEOUT:
      // MYSQL_OPT_CONNECT_TIMEOUT (argument type: unsigned int *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), istk(l_arg));
      break;
    case MYSQL_OPT_GUESS_CONNECTION:
      // MYSQL_OPT_GUESS_CONNECTION (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_LOCAL_INFILE:
      // MYSQL_OPT_LOCAL_INFILE (argument type: optional pointer to unsigned int)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), istk(l_arg));
      break;
    case MYSQL_OPT_NAMED_PIPE:
      // MYSQL_OPT_NAMED_PIPE (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_PROTOCOL:
      // MYSQL_OPT_PROTOCOL (argument type: unsigned int *)
      // Type of protocol to use. Should be one of the enum values of mysql_protocol_type defined in mysql.h.

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), istk(l_arg));
      break;
    case MYSQL_OPT_READ_TIMEOUT:
      // MYSQL_OPT_READ_TIMEOUT (argument type: unsigned int *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), istk(l_arg));
      break;
    case MYSQL_OPT_RECONNECT:
      // MYSQL_OPT_RECONNECT (argument type: my_bool *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), (my_bool*)istk(l_arg));
      break;
    case MYSQL_SET_CLIENT_IP:
      // MYSQL_SET_CLIENT_IP (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_OPT_SSL_VERIFY_SERVER_CERT:
      // MYSQL_OPT_SSL_VERIFY_SERVER_CERT (argument type: my_bool *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), (my_bool *)istk(l_arg));
      break;
    case MYSQL_OPT_USE_EMBEDDED_CONNECTION:
      // MYSQL_OPT_USE_EMBEDDED_CONNECTION (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_USE_REMOTE_CONNECTION:
      // MYSQL_OPT_USE_REMOTE_CONNECTION (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_USE_RESULT:
      // MYSQL_OPT_USE_RESULT (argument: not used)

      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), NULL);
      break;
    case MYSQL_OPT_WRITE_TIMEOUT:
      // MYSQL_OPT_WRITE_TIMEOUT (argument type: unsigned int *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), istk(l_arg));
      break;
    case MYSQL_READ_DEFAULT_FILE:
      // MYSQL_READ_DEFAULT_FILE (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_READ_DEFAULT_GROUP:
      // MYSQL_READ_DEFAULT_GROUP (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_REPORT_DATA_TRUNCATION:
      // MYSQL_REPORT_DATA_TRUNCATION (argument type: my_bool *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), (my_bool *)istk(l_arg));
      break;
    case MYSQL_SECURE_AUTH:
      // MYSQL_SECURE_AUTH (argument type: my_bool *)

      GetRhsVar(3, "i", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), (my_bool *)istk(l_arg));
      break;
    case MYSQL_SET_CHARSET_DIR:
      // MYSQL_SET_CHARSET_DIR (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_SET_CHARSET_NAME:
      // MYSQL_SET_CHARSET_NAME (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    case MYSQL_SHARED_MEMORY_BASE_NAME:
      // MYSQL_SHARED_MEMORY_BASE_NAME (argument type: char *)

      GetRhsVar(3, "c", &m_arg, &n_arg, &l_arg);
      *istk(l_status) =  mysql_options(mysql_ptr_in, (enum mysql_option)*istk(l_option), cstk(l_arg));
      break;
    default:
      sciprint("%s: wrong option: %s\n", fname, cstk(l_option));
    }

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_ping(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) = mysql_ping(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_query(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_query, n_query, l_query;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_query, &n_query, &l_query);

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) = mysql_query(mysql_ptr_in, cstk(l_query));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_real_connect(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_sql_pointer_out = 1, n_sql_pointer_out = 1, l_sql_pointer_out;
  int m_host, n_host, l_host;
  int m_user, n_user, l_user;
  int m_passwd, n_passwd, l_passwd;
  int m_db, n_db, l_db;
  int m_port, n_port, l_port;
  int m_socket, n_socket, l_socket;
  int m_flag, n_flag, l_flag;
  int m_status = 1, n_status = 1, l_status;
  char * host = NULL, * user = NULL, * db = NULL, * socket = NULL, * passwd = NULL;
  unsigned int port = 0, flag = 0;
  MYSQL * mysql_ptr_in = NULL;
  MYSQL * mysql_ptr_out = NULL;

  CheckRhs(1,8);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  if (Rhs>=2) 
    {
      GetRhsVar(2, "c", &m_host, &n_host, &l_host);
      if ((m_host!=0)&&(l_host!=0)) host = cstk(l_host);
      else                          host = NULL;
    }
  else
    {
      host = NULL;
    }
  if (Rhs>=3) 
    {
      GetRhsVar(3, "c", &m_user, &n_user, &l_user);
      if ((m_user!=0)&&(n_user!=0)) user = cstk(l_user);
      else                          user = NULL;
    }
  else
    {
      user = NULL;
    }
  if (Rhs>=4)
    {
      GetRhsVar(4, "c", &m_passwd, &n_passwd, &l_passwd);
      if ((m_passwd!=0)&&(n_passwd!=0)) passwd = cstk(l_passwd);
      else                              passwd = NULL;
    }
  else
    {
      passwd = NULL;
    }
  if (Rhs>=5)
    {
      GetRhsVar(5, "c", &m_db, &n_db, &l_db);
      if ((m_db!=0)&&(n_db!=0)) db = cstk(l_db);
      else                      db = NULL;
    }
  else
    {
      db = NULL;
    }
  if (Rhs>=6)
    {
      GetRhsVar(6, "i", &m_port, &n_port, &l_port);
      if ((m_port!=0)&&(n_port!=0)) port = *istk(l_port);
      else                          port = 0;
    }
  else
    {
      port = 0;
    }
  if (Rhs>=7)
    {
      GetRhsVar(7, "c", &m_socket, &n_socket, &l_socket);
      if ((m_socket!=0)&&(n_socket!=0)) socket = cstk(l_socket);
      else                              socket = NULL;
    }
  else
    {
      socket = NULL;
    }
  if (Rhs>=8)
    {
      GetRhsVar(8,"i",&m_flag, &n_flag, &l_flag);
      if ((m_flag!=0)&&(n_flag!=0)) flag = *istk(l_flag);
      else                          flag = 0;
    }
  else
    {
      flag = 0;
    }

  mysql_ptr_out = mysql_real_connect(mysql_ptr_in, host, user, passwd, db, port, socket, flag);

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  if (mysql_ptr_out==NULL) *istk(l_status) = -1;
  else                     *istk(l_status) = 0;

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_real_escape_string(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_de, n_de, l_de;
  int m_en, n_en, l_en;
  MYSQL * mysql_ptr_in = NULL;
  char * en = NULL;
  unsigned long res_long;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_de, &n_de, &l_de);

  en = (char *)MALLOC((2*strlen(cstk(l_de))+1)*sizeof(char));

  res_long = mysql_real_escape_string(mysql_ptr_in, en, cstk(l_de), strlen(cstk(l_de)));

  m_en = 1;
  n_en = strlen(en);

  CreateVar(Rhs+1, "c", &m_en, &n_en, &l_en);
  strncpy(cstk(l_en),en,strlen(en));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_real_query(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_query, n_query, l_query;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_query, &n_query, &l_query);

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) =  mysql_real_query(mysql_ptr_in, cstk(l_query), strlen(cstk(l_query)));

  LhsVar(1) = Rhs + 1;

  return 0;
}

#ifdef USE_OLD_FUNCTIONS
extern "C" int sci_mysql_reload(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) =  mysql_reload(mysql_ptr_in);

  LhsVar(1) = Rhs+1;

  return 0;
}
#endif

extern "C" int sci_mysql_row_seek(char * fname)
{
  // MYSQL_ROW_OFFSET is an equivalent to MYSQL_ROWS *
  int m_sql_res_pointer_in, n_sql_res_pointer_in, l_sql_res_pointer_in;
  int m_sql_row_offset_pointer_in, n_sql_row_offset_pointer_in, l_sql_row_offset_pointer_in;
  int m_sql_row_offset_pointer_out = 1, n_sql_row_offset_pointer_out = 1, l_sql_row_offset_pointer_out;
  MYSQL_RES * mysql_res_ptr_in = NULL;
  MYSQL_ROW_OFFSET mysql_row_offset_ptr_in = NULL;
  MYSQL_ROW_OFFSET mysql_row_offset_ptr_out;
  
  CheckRhs(2,2);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr_in = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "p", &m_sql_row_offset_pointer_in,  &n_sql_row_offset_pointer_in,  &l_sql_row_offset_pointer_in);
  mysql_row_offset_ptr_in = (MYSQL_ROW_OFFSET)((unsigned long)*stk(l_sql_row_offset_pointer_in));

  if (mysql_row_offset_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysqlrow_offset pointer\n", fname);
      return 0;
    }

  mysql_row_offset_ptr_out = mysql_row_seek(mysql_res_ptr_in, mysql_row_offset_ptr_in);

  if (mysql_row_offset_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_row_offset_pointer_out, &n_sql_row_offset_pointer_out, mysql_row_offset_ptr_out);
    }
  else
    {
      m_sql_row_offset_pointer_out = 0;
      n_sql_row_offset_pointer_out = 0;
      CreateVar(Rhs+1, "d", &m_sql_row_offset_pointer_out, &n_sql_row_offset_pointer_out, &l_sql_row_offset_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_row_tell(char * fname)
{
  // MYSQL_ROW_OFFSET is an equivalent to MYSQL_ROWS *
  int m_sql_res_pointer_in, n_sql_res_pointer_in, l_sql_res_pointer_in;
  int m_sql_row_offset_pointer_out = 1, n_sql_row_offset_pointer_out = 1, l_sql_row_offset_pointer_out;
  MYSQL_RES * mysql_res_ptr_in = NULL;
  MYSQL_ROW_OFFSET mysql_row_offset_ptr_out;
  
  CheckRhs(1,1);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_res_pointer_in,  &n_sql_res_pointer_in,  &l_sql_res_pointer_in);
  mysql_res_ptr_in = (MYSQL_RES *)((unsigned long)*stk(l_sql_res_pointer_in));

  if (mysql_res_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql_res pointer\n", fname);
      return 0;
    }

  mysql_row_offset_ptr_out = mysql_row_tell(mysql_res_ptr_in);

  if (mysql_row_offset_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_row_offset_pointer_out, &n_sql_row_offset_pointer_out, mysql_row_offset_ptr_out);
    }
  else
    {
      m_sql_row_offset_pointer_out = 0;
      n_sql_row_offset_pointer_out = 0;
      CreateVar(Rhs+1, "d", &m_sql_row_offset_pointer_out, &n_sql_row_offset_pointer_out, &l_sql_row_offset_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_select_db(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_db, n_db, l_db;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_db, &n_db, &l_db);

  sciprint("DEBUG: db = %s\n", cstk(l_db));

  CreateVar(Rhs + 1, "i", &m_status, &n_status, &l_status);

  // YC: before using this function, we must check that we are connected ....
  //     otherwise: crash
  *istk(l_status) = mysql_select_db(mysql_ptr_in, cstk(l_db));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_set_server_option(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_option, n_option, l_option;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_option, &n_option, &l_option);

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) = mysql_set_server_option(mysql_ptr_in, (enum enum_mysql_set_option)*istk(l_option));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_shutdown(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_shutdown_level, n_shutdown_level, l_shutdown_level;
  int m_status = 1, n_status = 1, l_status;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_shutdown_level, &n_shutdown_level, &l_shutdown_level);

  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);

  *istk(l_status) = mysql_shutdown(mysql_ptr_in, (mysql_enum_shutdown_level)*istk(l_shutdown_level));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_sqlstate(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_state, n_state, l_state;
  char * state = NULL;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  state = (char *)mysql_sqlstate(mysql_ptr_in);

  m_state = 1;
  n_state = strlen(state);

  CreateVar(Rhs+1, "c", &m_state, &n_state, &l_state);
  strncpy(cstk(l_state),state,strlen(state));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_ssl_set(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_key, n_key, l_key;
  int m_cert, n_cert, l_cert;
  int m_ca, n_ca, l_ca;
  int m_capath, n_capath, l_capath;
  int m_cipher, n_cipher, l_cipher;
  int m_status = 1, n_status = 1, l_status;
  char * key = NULL, * cert = NULL, * ca = NULL, * capath = NULL, * cipher = NULL;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,6);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  if (Rhs>=2)
    {
      GetRhsVar(2,"c", &m_key, &n_key, &l_key);
      if ((m_key!=0)&&(n_key!=0)) key = cstk(l_key);
      else                        key = NULL;
    }
  else
    {
      key = NULL;
    }
  if (Rhs>=3)
    {
      GetRhsVar(3,"c", &m_cert, &n_cert, &l_cert);
      if ((m_cert!=0)&&(n_cert!=0)) cert = cstk(l_cert);
      else                          cert = NULL;
    }
  else
    {
      cert = NULL;
    }
  if (Rhs>=4)
    {
      GetRhsVar(4,"c", &m_ca, &n_ca, &l_ca);
      if ((m_ca!=0)&&(n_ca!=0)) ca = cstk(l_ca);
      else                      ca = NULL;
    }
  else
    {
      ca = NULL;
    }
  if (Rhs>=5)
    {
      GetRhsVar(5,"c", &m_capath, &n_capath, &l_capath);
      if ((m_capath!=0)&&(n_capath!=0)) capath = cstk(l_capath);
      else                              capath = NULL;
    }
  else
    {
      capath = NULL;
    }
  if (Rhs>=6)
    {
      GetRhsVar(6,"c", &m_cipher, &n_cipher, &l_cipher);
      if ((m_cipher!=0)&&(n_cipher!=0)) cipher = cstk(l_cipher);
      else                              cipher = NULL;
    }
  else
    {
      cipher = NULL;
    }
  
  CreateVar(Rhs+1, "i", &m_status, &n_status, &l_status);
	  
  *istk(l_status) = mysql_ssl_set(mysql_ptr_in, key, cert, ca, capath, cipher);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_stat(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_stat, n_stat, l_stat;
  MYSQL * mysql_ptr_in = NULL;
  char * stat = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  stat = (char *)mysql_stat(mysql_ptr_in);

  m_stat = 1;
  n_stat = strlen(stat);

  CreateVar(Rhs + 1, "c", &m_stat, &n_stat, &l_stat);
  strncpy(cstk(l_stat),stat,strlen(stat));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_store_result(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_sql_res_pointer_out = 1, n_sql_res_pointer_out = 1, l_sql_res_pointer_out;
  MYSQL * mysql_ptr_in = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;

  CheckRhs(1,1);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  mysql_res_ptr_out = mysql_store_result(mysql_ptr_in);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0;
      n_sql_res_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, &l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_thread_id(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_thread_id = 1, n_thread_id = 1, l_thread_id;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_thread_id, &n_thread_id, &l_thread_id);

  *istk(l_thread_id) = mysql_thread_id(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_use_result(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_sql_res_pointer_out = 1, n_sql_res_pointer_out = 1, l_sql_res_pointer_out;
  MYSQL * mysql_ptr_in = NULL;
  MYSQL_RES * mysql_res_ptr_out = NULL;

  CheckRhs(1,1);
  CheckLhs(1,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  mysql_res_ptr_out = mysql_use_result(mysql_ptr_in);

  if (mysql_res_ptr_out)
    {
      CreateVarFromPtr(Rhs+1, "p", &m_sql_res_pointer_out, &n_sql_res_pointer_out, mysql_res_ptr_out);
    }
  else
    {
      m_sql_res_pointer_out = 0;
      n_sql_res_pointer_out = 0;

      CreateVar(Rhs+1, "d", &m_sql_res_pointer_out, &n_sql_res_pointer_out, &l_sql_res_pointer_out);
    }

  LhsVar(1) = Rhs+1;

  return 0;
}

extern "C" int sci_mysql_warning_count(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_warning = 1, n_warning = 1, l_warning;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_warning, &n_warning, &l_warning);

  *istk(l_warning) = mysql_warning_count(mysql_ptr_in);

  LhsVar(1) = Rhs+1;
  
  return 0;
}

extern "C" int sci_mysql_commit(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_commit = 1, n_commit = 1, l_commit;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_commit, &n_commit, &l_commit);

  *istk(l_commit) = mysql_commit(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_rollback(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_rollback = 1, n_rollback = 1, l_rollback;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_rollback, &n_rollback, &l_rollback);

  *istk(l_rollback) = mysql_rollback(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_autocommit(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_mode, n_mode, l_mode;
  int m_rollback = 1, n_rollback = 1, l_rollback;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));
  
  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_mode, &n_mode, &l_mode);

  CreateVar(Rhs+1, "i", &m_rollback, &n_rollback, &l_rollback);

  *istk(l_rollback) = mysql_autocommit(mysql_ptr_in, *istk(l_mode));

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_more_results(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_more_results = 1, n_more_results = 1, l_more_results;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_more_results, &n_more_results, &l_more_results);

  *istk(l_more_results) = mysql_more_results(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_next_result(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_next_result = 1, n_next_result = 1, l_next_result;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  CreateVar(Rhs+1, "i", &m_next_result, &n_next_result, &l_next_result);

  *istk(l_next_result) = mysql_next_result(mysql_ptr_in);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_check_null(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_test = 1, n_test = 1, l_test;
  void * mysql_ptr_in = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (void *)((unsigned long)*stk(l_sql_pointer_in));

  CreateVar(Rhs+1,"i", &m_test, &n_test, &l_test);

  *istk(l_test) = (mysql_ptr_in==NULL);

  LhsVar(1) = Rhs + 1;

  return 0;
}

extern "C" int sci_mysql_get_character_set_info(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_extra       = 9, n_extra       = 1, l_extra;
  int m_list_labels = 9, n_list_labels = 1, l_list_list_labels, l_list_labels;
  int m_number,   n_number,   l_list_number,   l_number;
  int m_state,    n_state,    l_list_state,    l_state;
  int m_csname,   n_csname,   l_list_csname,   l_csname;
  int m_name,     n_name,     l_list_name,     l_name;
  int m_comment,  n_comment,  l_list_comment,  l_comment;
  int m_dir,      n_dir,      l_list_dir,      l_dir;
  int m_mbminlen, n_mbminlen, l_list_mbminlen, l_mbminlen;
  int m_mbmaxlen, n_mbmaxlen, l_list_mbmaxlen, l_mbmaxlen;

  MYSQL * mysql_ptr_in = NULL;
  MY_CHARSET_INFO cs;
  const char * FieldNames[] = {"mysql_cs", 
			       "number",      /* character set number              */
			       "state",       /* character set state               */
			       "csname",      /* collation name                    */
			       "name",        /* character set name                */
			       "comment",     /* comment                           */
			       "dir",         /* character set directory           */
			       "mbminlen",    /* min. length for multibyte strings */
			       "mbmaxlen"};   /* max. length for multibyte strings */

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  mysql_get_character_set_info(mysql_ptr_in, &cs);

  m_number = 1; n_number = 1;
  CreateVar(Rhs+1, "i", &m_number, &n_number, &l_number);
  *istk(l_number) = cs.number;

  m_state = 1; n_state = 1;
  CreateVar(Rhs+2, "i", &m_state, &n_state, &l_state);
  *istk(l_state) = cs.state;

  m_csname = 1; n_csname = strlen(cs.csname);
  CreateVar(Rhs+3, "c", &m_csname, &n_csname, &l_csname);
  strcpy(cstk(l_csname), cs.csname);

  m_name = 1; n_name = strlen(cs.name);
  CreateVar(Rhs+4, "c", &m_name, &n_name, &l_name);
  strcpy(cstk(l_name), cs.name);
  
  m_comment = 1; n_comment = strlen(cs.comment);
  CreateVar(Rhs+5, "c", &m_comment, &n_comment, &l_comment);
  strcpy(cstk(l_comment), cs.comment);

  m_dir = 1; n_dir = strlen(cs.dir);
  CreateVar(Rhs+6, "c", &m_dir, &n_dir, &l_dir);
  strcpy(cstk(l_dir), cs.dir);

  m_mbminlen = 1; n_mbminlen = 1;
  CreateVar(Rhs+7, "i", &m_mbminlen, &n_mbminlen, &l_mbminlen);
  *istk(l_mbminlen) = cs.mbminlen;

  m_mbmaxlen = 1; n_mbmaxlen = 1;
  CreateVar(Rhs+8, "i", &m_mbmaxlen, &n_mbmaxlen, &l_mbmaxlen);
  *istk(l_mbmaxlen) = cs.mbmaxlen;

  CreateVar(Rhs+9, "m", &m_extra, &n_extra, &l_extra);
  CreateListVarFromPtr(Rhs+9, 1, "S",  &m_list_labels, &n_list_labels, FieldNames);
  CreateListVarFrom(Rhs+9, 2, "i", &n_number,   &m_number,   &l_list_number,   &l_number);
  CreateListVarFrom(Rhs+9, 3, "i", &n_state,    &m_state,    &l_list_state,    &l_state);
  CreateListVarFrom(Rhs+9, 4, "c", &n_csname,   &m_csname,   &l_list_csname,   &l_csname);
  CreateListVarFrom(Rhs+9, 5, "c", &n_name,     &m_name,     &l_list_name,     &l_name);
  CreateListVarFrom(Rhs+9, 6, "c", &n_comment,  &m_comment,  &l_list_comment,  &l_comment);
  CreateListVarFrom(Rhs+9, 7, "c", &n_dir,      &m_dir,      &l_list_dir,      &l_dir);
  CreateListVarFrom(Rhs+9, 8, "i", &n_mbminlen, &m_mbminlen, &l_list_mbminlen, &l_mbminlen);
  CreateListVarFrom(Rhs+9, 9, "i", &n_mbmaxlen, &m_mbmaxlen, &l_list_mbmaxlen, &l_mbmaxlen);

  LhsVar(1) = Rhs + 9;

  return 0;
}

int sci_mysql_get_ssl_cipher(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_result = 1, n_result = 1, l_result;
  MYSQL * mysql_ptr_in = NULL;
  char * result = NULL;

  CheckRhs(1,1);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  result = (char *)mysql_get_ssl_cipher(mysql_ptr_in);

  if (result)
    {
      m_result = 1; n_result = strlen(result);
      CreateVar(Rhs+1, "c", &m_result, &n_result, &l_result);
      strcpy(cstk(l_result),result);
    }
  else
    {
      m_result = 0; n_result = 0;
      CreateVar(Rhs+1, "d", &m_result, &n_result, &l_result);
    }

  LhsVar(1) = 0;

  return 0;
}

int sci_mysql_refresh(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_option, n_option, l_option;
  int m_result = 1, n_result = 1, l_result;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "i", &m_option, &n_option, &l_option);

  CreateVar(Rhs+1, "i", &m_result, &n_result, &l_result);

  *istk(l_result) = mysql_refresh(mysql_ptr_in, *istk(l_option));
  
  LhsVar(1) = Rhs + 1;

  return 0;
}

int sci_mysql_set_character_set(char * fname)
{
  int m_sql_pointer_in, n_sql_pointer_in, l_sql_pointer_in;
  int m_charset, n_charset, l_charset;
  int m_result = 1, n_result = 1, l_result;
  MYSQL * mysql_ptr_in = NULL;

  CheckRhs(2,2);
  CheckLhs(0,1);

  GetRhsVar(1, "p", &m_sql_pointer_in,  &n_sql_pointer_in,  &l_sql_pointer_in);
  mysql_ptr_in = (MYSQL *)((unsigned long)*stk(l_sql_pointer_in));

  if (mysql_ptr_in==NULL)
    {
      Scierror(999,"%s: problem with the mysql pointer\n", fname);
      return 0;
    }

  GetRhsVar(2, "c", &m_charset, &n_charset, &l_charset);

  CreateVar(Rhs + 1, "i", &m_result, &n_result, &l_result);

  *istk(l_result) =  mysql_set_character_set(mysql_ptr_in, cstk(l_charset));

  return 0;
}
