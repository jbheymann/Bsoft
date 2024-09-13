/**
@file	bpot.cpp
@brief	Calculating potential from atomic models
@author Bernard Heymann
@date	Created: 20230409
@date 	Modified: 20230614
**/

#include "model_map.h"
#include "rwmodel.h"
#include "rwimg.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_mol.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bpot [options] in.pdb out.mrc",
"------------------------------------",
"Calculates a potential map from atomic models.",
"A coordinate file is required.",
"Note: The occupancy for each atom is used - make sure it is set properly.",
" ",
"Actions:",
"-center                  Center coordinates before calculations.",
"-View 0.3,-0.5,0.8,33    View to rotate the molecule to.",
"-realspace               Transform the output back to real space.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-size 100,80,70          Size (default automatic voxels).",
"-origin 0,0,0            Origin placement within image (default 0,0,0).",
"-sampling 2.5,2.5,2.5    Sampling (default 1 angstrom/voxel).",
"-resolution 3            High resolution limit (default 0.1 angstrom).",
"-density 1.27g           Material density (g=g/cm3, d=Da/A3, n=#/A3).",
" ",
"Input:",
"-parameters scat.star    Parameter file with scattering coefficients.",
" ",
"Output:",
"-rps file.txt            File for radial power spectrum.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    /* Initialize variables */
	int				center(0);				// Flag to center coordinates
	View2<double>	view;					// View to generate
	int 			set_backtransform(0);	// Flag for back transformation
	Vector3<double>	origin;					// Coordinate origin placement
	int				set_origin(0);			// Flag to set origin
	Vector3<long>	size;					// New map size
	Vector3<double>	sam(1,1,1);    			// Sampling in angstrom/voxel side
	double			hires(0.1);				// High resolution limit (angstrom)
	double			density(0);				// Material density
	DensityUnit		density_units(G_CM3);	// Density units (default g/cm3)
	int				ps_flags(0);			// Radial power spectrum flags, bits: 1=zero origin, 2=average, 8=log
	Bstring			paramfile;				// Use default parameter file
	Bstring			rpsfile;				// File name for radial power spectrum

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "View" )
			view = curropt->view();
		if ( curropt->tag == "realspace" )
			set_backtransform = 1;
		if ( curropt->tag == "size" )
			size = curropt->size();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "resolution" )
			if ( ( hires = curropt->value.real() ) < 0.001 )
				cerr << "-resolution: A high resolution limit must be specified!" << endl;
		if ( curropt->tag == "density" )
			density = curropt->density(density_units);
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "rps" )
			rpsfile = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();
	
	string			coorfile(argv[optind++]);
	Bmodel*			model = NULL;
	Vector3<double> box;
	vector<string>	file_list;
	
	if ( coorfile.length() ) {
		file_list = split(coorfile, ',');
		model = read_model(file_list, paramfile.str(), 1);
		
		if ( !model ) {
			cerr << "Error: Problem with coordinate file " << coorfile << ", exiting!" << endl;
			bexit(-1);
		}
	}
	
	if ( !model ) {
		cerr << "No coordinates read!" << endl;
		bexit(-1);
	}

	model_selection_stats(model);

	long			nimg(1);
	string			imgfile;
	if ( optind < argc ) imgfile = argv[optind++];
	
	if ( center )
		models_shift(model, -models_center_of_coordinates(model));

	if ( view[2] < 1 )
		model_rotate(model, view);

	vector<Vector3<double>>	bounds = models_calculate_bounds(model);

	if ( size.volume() < 100 ) {
		model->calculate_bounds();
		size = sam * (bounds[1] - bounds[0]);
	}
	
	if ( size[2] < 2 ) sam[2] = 1;

	if ( density < 1e-3 ) {
		density = RHO;
		density_units = DA_A3;
	}

	Bimage*			p = new Bimage(Float, TComplex, size, nimg);
	p->file_name(imgfile);
	p->sampling(sam);
	if ( set_origin != 1 ) origin = p->size()/2;
	p->fourier_type(Standard);
	p->image->view(view);

	Vector3<double>		start = -origin*sam;
	Vector3<double>		end = start + sam*size;
	long				nsel = models_select_within_bounds(model, start, end);
	if ( verbose ) {
		cout << "Model bounds:                   " << bounds[0] << tab << bounds[1] << endl;
		cout << "Volume bounds:                  " << start << tab << end << endl;
		cout << "Components selected:            " << nsel << endl;
	}

	model_show_selection(model);

//	img_potential_from_model(model, p, paramfile, hires, density, density_units);
	double		mip = img_potential_from_model_structure_factors(model, -1, p, paramfile.str(), hires);
	mip *= p->real_size().volume()*density/model_mass(model);

	if ( verbose ) {
		cout << "Model mass:                     " << model_mass(model) << " Da" << endl;
		cout << "Density:                        " << density << " Da/A3" << endl;
		cout << "Mean inner potential:           " << mip << " V" << endl << endl;
	}

	if ( rpsfile.length() ) {
		double		sampling_ratio(1);
		Bimage*		prad = p->fspace_radial_power(hires, sampling_ratio);
//		if ( ps_flags & 2 ) p->average_images();
		if ( ps_flags & 8 ) prad->logarithm();
		Bplot*		plot = prad->plot_radial_powerspectrum(hires, ps_flags);
		ps_plot(rpsfile, plot);
		delete prad;
	}

	if ( set_backtransform ) {
		if ( verbose )
			cout << "Back transforming" << endl;
		p->phase_shift(origin);
		p->fft(FFTW_BACKWARD, 0, Real);
	}

	// Write output file
	if ( imgfile.length() )
		write_img(imgfile, p, 0);

 	delete model;
	if ( p ) delete p;

	
		timer_report(ti);
	
	bexit(0);
}

