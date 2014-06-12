#!/bin/sh
# This file can be used to create database glpk in MySQL.
USERNAME=$1
echo MySQL cleaning is called for user $USERNAME. 
mysql -f -u $USERNAME -p < scripts/erase.sql
