/**
@file	rwmodel.h
@brief	Header file for reading and writing atomic model files
@author Bernard Heymann
@date	Created: 20060919
@date	Modified: 20230621
**/

#include "Bmodel.h"

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
Bmodel*		model_add(Bmodel** model, string id);
Bcomponent*	component_add(Bcomponent** comp, string& id);
Bcomponent*	component_add(Bcomponent** comp, unsigned long number);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2, double length, double radius);
Blink*		link_add(Blink** link, Bcomponent* comp1, Bcomponent* comp2);
int			model_check(Bmodel* model, string path);
Bmodel*		model_list_copy(Bmodel* model);
int			component_list_kill(Bcomponent* comp);
int			comp_type_list_kill(Bcomptype* type);
int			model_link_list_kill(Bmodel* model);
int			link_kill(Blink** link_list, Bcomponent* comp, int i);
int			link_kill(Blink** link_list, Bcomponent* comp, Bcomponent* comp2);
int			poly_list_kill(Bpolygon* poly);
int			comp_associated_links_kill(Bcomponent* comp, Blink** link);
int 		model_kill(Bmodel* model);


