/**
@file	model_links.h
@brief	Library routines used for model link processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250426
**/

#include "rwmodel.h"

// Function prototypes
long		models_generate_links(Bmodel* waters, int distance_type);
long		models_generate_angles(Bmodel* model);
long		models_setup_links(Bmodel* model);
long		models_link_list_generate(Bmodel* model, double maxlength);
long		models_link_list_generate(Bmodel* model, double maxlength,
				string type1, string type2, int flag);
long		models_set_link_length(Bmodel* model, double linklength);
long		models_set_link_radius(Bmodel* model, double linkrad);
long		model_reduce_linked(Bmodel* model, string submodname, int flags);
long		models_links_minimum_valency(Bmodel* model, long valency);

