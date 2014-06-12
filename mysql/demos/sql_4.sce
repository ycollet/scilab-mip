// Don't forget to start mysql using /etc/init.d/mysqld start or via service mysqld start
// Create a data base using mysqladmin create database
// Fill it. The SQL query is based on a "table1" table and a "dvar" variable

// Some database examples (retrieved from the glpk project) are shipped in this directory
// To add this databases, proceed the following steps:
//
// If you wanted to test mysql interface on sudoku database: 
// mysql < sudoku.sql
// and set UseSudoku to %T and UseTransp to %F
// If you wanted to test mysql interface on transportation database:
// mysql < transp.sql
// and set UseSudoku to %F and UseTransp to %T
// These scripts will create a database glpk.
// They will create a user glpk with password gnu also.
// Once you have finished playing with mysql, you can remove the test database using erase.sql:
// mysql < erase.sql

/////////////////////
// File sudoku.sql //
/////////////////////
// Database: glpk
// Table sudoku
// Columns: ID (INT), COL (INT), LIN (INT), VAL (INT)
// Table: sudoku_solution
// Columns: ID (INT), COL (INT), LIN (INT), VAL (INT)

////////////////////
// File trans.sql //
////////////////////
// Database: glpk
// Table transp_capa (production capacity)
// Columns: PLANT (TEXT(127)), CAPA (REAL)
// Table: transp_demand (demand)
// Columns:  MARKET (TEXT(127)), DEMAND (REAL)
// Table trans_dist (distance)
// Columns: LOC1 (TEXT(127)), LOC2 (TEXT(127)), DIST (REAL)
// Table transp_result (result)
// Columns: LOC1 (TEXT(127)), LOC2 (TEXT(127)), QUANTITY (REAL)

username = 'glpk'; // Put your username
password = 'gnu';  // Put your password
database = 'glpk'; 
port     = 3306;   // use netstat -a | grep mysql to locate the mysql port
                   // or ps -elf | grep mysql and locate --port
myhost   = 'localhost'; // localhost most of the time

sql_ptr = mysql_init();
status  = mysql_real_connect(sql_ptr, myhost, username, password, database, port);

timer();

// Exemple 4
if (status) then
  printf('Error message: %s, error number: %d\n', mysql_error(sql_ptr), mysql_errno(sql_ptr));
end

// For sudoku
step       = 1.0;
start_step = 0.0;
end_step   = 7.0;
res        = 0.0;
row        = [];

for var=start_step:step:end_step do
  sql = sprintf("select val from sudoku where val>=%f and val<%f;\n", var, end_step);

  printf('query: ' + sql + ' -> ');
  status = mysql_real_query(sql_ptr, sql);
  disp(mysql_error(sql_ptr));
  rs   = mysql_store_result(sql_ptr);
  _row = mysql_fetch_row(rs);
  for i=1:length(_row)
    printf('%s ', _row(i));
  end
  printf('\n');
  res = res + sum(evstr(_row));
  mysql_free_result(rs);
end

disp(timer());

printf('Cumulated sum for the selected items: %f\n',res);

printf('Number of affected rows: %d\n', mysql_affected_rows(sql_ptr));

mysql_close(sql_ptr);
