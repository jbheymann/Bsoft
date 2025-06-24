/**
@file	rwmodel.h
@brief	Header file for reading and writing atomic model files
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20230623
**/

#include "Bmodel.h"
#include "json.h"

/* Function prototypes */
Bmodel*		read_model(string filename);
Bmodel*		read_model(vector<string> file_list);
Bmodel*		read_model(string filename, string paramfile);
Bmodel*		read_model(vector<string> file_list, string paramfile);
Bmodel*		read_model(string filename, int type_select);
Bmodel*		read_model(vector<string> file_list, int type_select);
Bmodel*		read_model(string filename, string paramfile, int type_select);
Bmodel*		read_model(vector<string> file_list, string paramfile, int type_select);
int			write_model(string filename, Bmodel* model);
int			write_model(string filename, Bmodel* model, int split);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2, double length, double radius);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2);
int			model_check(Bmodel* model, string path);
JSvalue		model_types(Bmodel *model);
int			models_assign_types(Bmodel* model, JSvalue& types);
Bmodel*		model_list_copy(Bmodel* model);
int			comp_associated_links_kill(Bcomponent* comp, Blink** link);


