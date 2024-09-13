/**
@file	model_links.h
@brief	Library routines used for model link processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20230717
**/

#include "rwmodel.h"

// Function prototypes
long		model_generate_links(Bmodel* waters, int distance_type);
long		model_generate_angles(Bmodel* model);
long		model_setup_links(Bmodel* model);
long		model_link_list_generate(Bmodel* model, double maxlength);
long		model_link_list_generate(Bmodel* model, double maxlength,
				string type1, string type2, int flag);
long		model_set_link_length(Bmodel* model, double linklength);
long		model_set_link_radius(Bmodel* model, double linkrad);
long		model_reduce_linked(Bmodel* model, string submodname, int flags);
long		model_links_minimum_valency(Bmodel* model, long valency);

