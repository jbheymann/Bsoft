/**
@file	mg_simulate.h
@brief	Generates projects and simulate images 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20250407
**/

#include "mg_processing.h"
#include "rwimg.h"
#include "rwmolecule.h"

// Function prototypes
//int			img_apply_complex_CTF(Bimage* p, CTFparam& cp);
Bproject*	project_generate(int nfield, int nmg, int npart,
				Vector3<double>  pixel_size, double img_origin,
				CTFparam& cp, double def_min, double def_max, double dose,
				double tsigma, Bstring& fieldbase, Bstring& mgbase, Bstring& partbase,
				int fieldnumber, int mgnumber, int partnumber);
Bproject*	project_generate_asu(string& symmetry_string,
				Vector3<double>  pixel_size, double img_origin,
				double theta_step, double phi_step,
				CTFparam& cp, double defocus, double dose, Bstring& mgbase, Bstring& partbase);
int			project_generate_potential(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution, 
				double Bfactor, int pottype, Bstring& paramfile);
int			project_generate_projections(Bmolgroup* molgroup, Bmolgroup* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int type, Bstring& paramfile);
int			project_generate_projections(Bmodel* model, Bmodel* water, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				Vector3<long> size, double thickness, double resolution,
				double Bfactor, int type, int ab_flag, Bstring& paramfile);
int			project_generate_projections(Bimage* map, Bproject* project, 
				Bstring& fieldname, Bstring& mgname, int partselect,
				double resolution, int ew_flag, int ab_flag);
int			project_generate_image(Bproject* project, double thickness, double resolution);
int			project_apply_distortions(Bproject* project, int poisson, double gauss, double kmtf);
