/**
@file	model_views.h
@brief	Library routines used for analysing model component views
@author 	Bernard Heymann
@date	Created: 20081120
@date	Modified: 20230512
**/

#include "Bmodel.h"
#include "Bstring.h"

// Function prototypes
vector<View2<double>>	views_from_model(Bmodel* model);
vector<View2<double>>	views_from_models(Bmodel* model);
long		model_set_views(Bmodel* model, View2<double> view);
long		model_invert_views(Bmodel* model);
long		model_find_views(Bmodel* model, Bstring& reffile, Bstring& paramfile);
long		model_calculate_views(Bmodel* model, Bstring& mode);
long		model_calculate_local_views(Bmodel* model);
long		model_view_directions(Bmodel* model, int bin_width, int ref_flag);
int			component_hand(Bstring s);

