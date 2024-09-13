/**
@file	model_shell.h
@brief	Library routines used for shell model processing
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20161013
**/

#include "Bmodel.h"

// Function prototypes
long		model_add_shell(Bmodel* model, double add_distance, string new_type);
int			model_adjust_shell_to_guide(Bmodel* model, Bmodel* gmod, double fraction, int curv_flag);
Bmodel*		model_components_to_shells(Bmodel* model, double distance, string nutype, int twod);
double		model_sphericity(Bmodel* model);
double		model_ellipsoidicity(Bmodel* model);
int			model_curvature(Bmodel* model);
double		model_inside_outside(Vector3<double> vec, Bmodel* model, int curv_flag, int fast);

