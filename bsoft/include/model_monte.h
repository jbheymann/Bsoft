/**
@file	model_monte.h
@brief	Headers of functions using a monte carlo metroplis algorithm to energy minimize models.
@author 	Bernard Heymann
@date	Created: 20041230
@date	Modified: 20250620
**/

#include "rwimg.h"
#include "rwmodel.h"
#include "rwmodel_param.h"
#include "rwmd.h"


// Function prototypes 
Bmodel*	monte_carlo_metropolis(Bmodel* model, Bmodparam& md, Bimage* map, 
				double beta, double max_angle, double max_shift, long max_iter, 
				double (Efunc)(Bmodel*, Bimage*, Bmodparam&),
				int (Tfunc)(Bmodel*, double, double, int));
long		model_test_if_within_box(Bmodel* model, Vector3<double> min, Vector3<double> max);
double		monte_rigid_body_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md);
double		monte_component_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md);
double		monte_link_fit_energy(Bmodel* model, Bimage* map, Bmodparam& md);
double		models_component_overlap(Bmodel* model, Bmodparam& md);
double		models_clashes(Bmodel* model, Bmodparam& md);
int			model_rigid_body_transform(Bmodel* model, double max_angle, double shift_std, int rigid);
int			model_move_components_down_energy(Bmodel* model, double max_angle, double max_shift, int rigid);

