/**
@file	model_water.h
@brief	Generating and managing water
@author 	Bernard Heymann
@date	Created: 20001014
@date	Modified: 20230704
**/

#include "Bmodel.h"
#include "Vector3.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Function prototypes
Bmodel*		model_generate_one_water(string& watername, Vector3<double> Ocoord);
Bmodel*		model_generate_regular_water(Vector3<double> size, int type);
Bmodel*		model_generate_random_water(Vector3<double> size);
//Bangle*		water_angle_list(Bmolgroup* molgroup);
int			model_calc_water_rdf(Bmodel* waters, double interval, double cutoff);

