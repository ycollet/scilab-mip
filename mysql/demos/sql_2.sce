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

// Example 2
status  = mysql_select_db(sql_ptr,'glpk');

printf('Number of affected rows: %d\n', mysql_affected_rows(sql_ptr));

mysql_close(sql_ptr);
