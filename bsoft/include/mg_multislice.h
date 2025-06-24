/**
@file	mg_multislice.h
@brief	Multislice simulation 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20241211
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "rwmolecule.h"

// Function prototypes
Bimage*		img_calc_wave_propagator(Vector3<long> size, Vector3<double> sam,
				double thickness, double volt);
Bimage*		img_calc_multi_slice(Bimage* pgrate, double thickness, double volt, bool norm);
Bimage*		img_project_multi_slice(Bimage* pgrate, double thickness, double volt);
Bimage*		img_calc_potential(Bmolgroup* molgroup, Vector3<long> size, Vector3<double> origin,
				Vector3<double> sam, double thickness, double resolution, double Bfactor, 
				Bstring& paramfile, int type);
int			img_calc_phase_grating(Bimage* p, double volt, bool norm);
