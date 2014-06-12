#!/bin/sh
# This file can be used to create database glpk in MySQL.
USERNAME=$1
echo MySQL is called for user $USERNAME. 
mysql -f -u $USERNAME -p < scripts/sudoku.sql
#echo MySQL is called for user $USERNAME.
#mysql -f -u $USERNAME -p < scripts/transp.sql
