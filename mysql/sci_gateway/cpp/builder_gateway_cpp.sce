// ====================================================================
// Yann COLLETTE
// Copyright 2010-2012
// This file is released into the public domain
// ====================================================================

UseOldFunctions = %f; // Force the use of some old deprecated functions

list_add_inter = [];

if UseOldFunctions then
  list_add_inter = ['mysql_connect',            'sci_mysql_connect'; ...
                    'mysql_create_db',          'sci_mysql_create_db'; ...
		    'mysql_drop_db',            'sci_mysql_drop_db'; ...
		    'mysql_reload',             'sci_mysql_reload'];
end

// TO BE DONE:
//		  'mysql_fetch_fields',       'sci_mysql_fetch_fields'; ...

list_add_inter = [list_add_inter;
                  'mysql_affected_rows',          'sci_mysql_affected_rows'; ...
                  'mysql_change_user',            'sci_mysql_change_user'; ...
		  'mysql_character_set_name',     'sci_mysql_character_set_name'; ...
		  'mysql_close',                  'sci_mysql_close'; ...
		  'mysql_data_seek',              'sci_mysql_data_seek'; ...
		  'mysql_debug',                  'sci_mysql_debug'; ...
		  'mysql_dump_debug_info',        'sci_mysql_dump_debug_info'; ...
		  'mysql_eof',                    'sci_mysql_eof'; ...
		  'mysql_errno',                  'sci_mysql_errno'; ...
		  'mysql_error',                  'sci_mysql_error'; ...
		  'mysql_fetch_field',            'sci_mysql_fetch_field'; ...
		  'mysql_fetch_field_direct',     'sci_mysql_fetch_field_direct'; ...
		  'mysql_fetch_lengths',          'sci_mysql_fetch_lengths'; ...
		  'mysql_fetch_row',              'sci_mysql_fetch_row'; ...
		  'mysql_field_count',            'sci_mysql_field_count'; ...
		  'mysql_field_seek',             'sci_mysql_field_seek'; ...
		  'mysql_field_tell',             'sci_mysql_field_tell'; ...
		  'mysql_free_result',            'sci_mysql_free_result'; ...
		  'mysql_get_client_info',        'sci_mysql_get_client_info'; ...
		  'mysql_get_client_version',     'sci_mysql_get_client_version'; ...
		  'mysql_get_host_info',          'sci_mysql_get_host_info'; ...
		  'mysql_get_proto_info',         'sci_mysql_get_proto_info'; ...
		  'mysql_get_server_info',        'sci_mysql_get_server_info'; ...
		  'mysql_get_server_version',     'sci_mysql_get_server_version'; ...
		  'mysql_hex_string',             'sci_mysql_hex_string'; ...
		  'mysql_info',                   'sci_mysql_info'; ...
		  'mysql_init',                   'sci_mysql_init'; ...
		  'mysql_insert_id',              'sci_mysql_insert_id'; ...
		  'mysql_kill',                   'sci_mysql_kill'; ...
		  'mysql_list_dbs',               'sci_mysql_list_dbs'; ...
		  'mysql_list_fields',            'sci_mysql_list_fields'; ...
		  'mysql_list_processes',         'sci_mysql_list_processes'; ...
		  'mysql_list_tables',            'sci_mysql_list_tables'; ...
		  'mysql_num_fields',             'sci_mysql_num_fields'; ...
		  'mysql_num_rows',               'sci_mysql_num_rows'; ...
		  'mysql_options',                'sci_mysql_options'; ...
		  'mysql_ping',                   'sci_mysql_ping'; ...
		  'mysql_query',                  'sci_mysql_query'; ...
		  'mysql_real_connect',           'sci_mysql_real_connect'; ...
		  'mysql_real_escape_string',     'sci_mysql_real_escape_string'; ...
		  'mysql_real_query',             'sci_mysql_real_query'; ...
		  'mysql_row_seek',               'sci_mysql_row_seek'; ...
		  'mysql_row_tell',               'sci_mysql_row_tell'; ...
		  'mysql_select_db',              'sci_mysql_select_db'; ...
		  'mysql_set_server_option',      'sci_mysql_set_server_option'; ...
		  'mysql_shutdown',               'sci_mysql_shutdown'; ...
		  'mysql_sqlstate',               'sci_mysql_sqlstate'; ...		    
		  'mysql_ssl_set',                'sci_mysql_ssl_set'; ...
		  'mysql_stat',                   'sci_mysql_stat'; ...
		  'mysql_store_result',           'sci_mysql_store_result'; ...
		  'mysql_thread_id',              'sci_mysql_thread_id'; ...
		  'mysql_use_result',             'sci_mysql_use_result'; ...
		  'mysql_warning_count',          'sci_mysql_warning_count'; ...
		  'mysql_commit',                 'sci_mysql_commit'; ...
		  'mysql_rollback',               'sci_mysql_rollback'; ...
		  'mysql_autocommit',             'sci_mysql_autocommit'; ...
		  'mysql_more_results',           'sci_mysql_more_results'; ...
		  'mysql_next_result',            'sci_mysql_next_result'; ...
		  'my_init',                      'sci_my_init'; ...
		  'mysql_library_init',           'sci_mysql_library_init'; ...
		  'mysql_library_end',            'sci_mysql_library_end'; ...
		  'mysql_check_null',             'sci_mysql_check_null'; ...
		  'mysql_get_character_set_info', 'sci_mysql_get_character_set_info'; ...
		  'mysql_get_ssl_cypher',         'sci_mysql_get_ssl_cypher'; ...
		  'mysql_refresh',                'sci_mysql_refresh'; ...
		  'mysql_set_character_set',      'sci_mysql_set_character_set' ];

files_to_compile    = ['sci_mysql.cpp'];

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');

if MSDOS then
  include_mysql       = '""c:/Program Files/MySQL/MySQL Server 5.1/include/""';
  mysql_lib           = '""c:/Program Files/MySQL/MySQL Server 5.1/lib/opt/libmysql.lib""';

  link_options        = '';//'/NODEFAULTLIB:LIBCMT'; // /NODEFAULTLIB:MSVCPRT /NODEFAULTLIB:MSVCRT 
  cflags              = '/D__USE_DEPRECATED_STACK_FUNCTIONS__ /I ' + path_builder + ' /I ' + include_mysql;
  ldflags             = mysql_lib + ' ' + link_options;
else
  include_mysql       = '/usr/include/mysql';

  mysql_lib           = '-L/usr/lib64 -lmysqlclient -lmysqlclient_r';

  cflags              = '-ggdb -D__USE_DEPRECATED_STACK_FUNCTIONS__ -I' + path_builder + ' -I' + include_mysql;
  ldflags             = mysql_lib;
end

tbx_build_gateway('mysql_cpp', list_add_inter, files_to_compile, get_absolute_file_path('builder_gateway_cpp.sce'), [], ldflags, cflags);

clear tbx_build_gateway;
