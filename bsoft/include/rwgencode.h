/**
@file	rwgencode.h
@brief	Header file for reading a genetic code from a file
@author Bernard Heymann
@date	Created: 20030316
@date	Modified: 20250601
**/

#include <map>
#include "string_util.h"

/* Function prototypes */
map<string,char>	get_genetic_code(string& filename);
int 		write_genetic_code(string& filename, map<string,char>& gc);


