/**
@file	model_compare.h
@brief	Functions to compare models and components
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250615
**/

#include "Bmodel.h"
#include "Matrix.h"
#include "rwresprop.h"

// Function prototypes
long		model_component_number_difference(Bmodel* model1, Bmodel* model2);
long		model_maxnum_components(Bmodel* model);
double		model_compare(Bmodel* model1, Bmodel* model2);
double		model_compare_by_distance(Bmodel* model1, Bmodel* model2);
double		model_compare_corresponding(Bmodel* model1, Bmodel* model2);
vector<double>	models_compare_corresponding(Bmodel* model1, Bmodel* model2);
Matrix		models_interfaces(Bmodel* model1, Bmodel* model2, double dcut);
long		model_interface(Bmodel* model1, Bmodel* model2, double dcut);
long		model_interface(Bmodel* model1, Bmodel* model2, map<string,Bresidue_type> res_prop);
Matrix		model_distance_matrix(Bmodel* model, int view_flag);
Matrix		model_distance_matrix(Bmodel* m1, Bmodel* m2);
Matrix		model_adjacency_matrix(Bmodel* model);
long		model_consolidate(Bmodel* model, double distance);
Bmodel*		models_consensus(Bmodel* model, double distance);
int			model_fit(Bmodel* model, Bmodel* refmod, string id);

