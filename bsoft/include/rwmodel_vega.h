/**
@file	rwmodel_vega.h
@brief	Header file for reading and writing Vega model parameters
@author 	Bernard Heymann
@date	Created: 20060919
@date	Modified: 20221115
**/

#include "rwmodel.h"

/* Function prototypes */
Bmodel*		read_model_vega(vector<string> file_list);
int			write_model_vega(string& filename, Bmodel* model, int splt);


