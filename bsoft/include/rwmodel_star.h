/**
@file	rwmodel_star.h
@brief	Header file for reading and writing STAR model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20210219
**/

#include "rwmodel.h"

/* Function prototypes */
map<string, int>	comptype_tags();
map<string, int>	linktype_tags();
map<string, int>	angletype_tags();
Bmodel*		read_model_star(vector<string> file_list);
int			write_model_star(string& filename, Bmodel* model, int splt);


