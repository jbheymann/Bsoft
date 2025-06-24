/**
@file	bmonte.cpp
@brief	Program to use a monte carlo metroplis algorithm to energy minimize models.
@author Bernard Heymann
@date	Created: 20041230
@date 	Modified: 20250618
**/

#include "model_monte.h"
#include "model_links.h"
#include "model_util.h"
#include "rwmodel.h"
#include "rwimg.h"
#include "rwmd.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmonte [options] in1.pdb [in2.cif ...]",
"---------------------------------------------",
"Minimizes the energy of a set of models using a Monte Carlo Metropolis algorithm.",
" ",
"Actions:",
"-rigid comp              Rigidity: all/whole (default), model, component.",
"-iterations 60           Maximum number of iterations (default 0).",
"-SS 2.05                 Add disulphides with this reference link length.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-componentradius 0.5     Component display radius.",
"-linkradius 0.5          Link display radius.",
" ",
"Parameters for mechanics:",
"-bbox 0,82.5,50,50,80,50 Bounding box center and size (default from coordinates or map).",
"-wrap                    Wrap around (periodic boundaries).",
"-Klink 150               Bond strength (default 0).",
"-Kangle 4                Angle strength (default 0).",
"-Kdistance 0.1           Distance (Van der Waals) strength (default 0).",
"-Kelectrostatic 0.4      Electrostatic strength (default 0).",
"-Kseparation 0.02        Separation energy constant (default 0).",
"-Kmap 1.5                Map energy constant (default 1 if -Map option is used).",
"-separation 4.5          Separation distance for overlap calculation (default 4 A).",
"-cutoff 7.8              Distance cutoff for non-linked forces (default 5 A).",
"-beta 12                 Inverse of mean energy per atom and degree of freedom (default 10).",
"-angle 1.5               Maximum angular increment per iteration (default 1 degree).",
"-shift 3.1               Maximum shift per iteration (default 1 angstrom).",
"-linksteps 5             Number of steps along a link for map fitting (default none).",
"-location 12,5.3,6       Location of harmonic force to apply to center-of-mass.",
"-Klocation 16.8,0.2      Magnitude of harmonic force and decay constant.",
" ",
"Input:",
"-parameters md.star      Molecular dynamics parameter file.",
"-Map file.map            Map to use as an additional restraint.",
//"-Mask mask.mrc           Mask to limit grid-searches (must be byte data type).",
" ",
"Output:",
"-output file.cif         Output model file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
    // Initialize variables
	Vector3<double>	sam;    				// Map sampling
	double			compradius(0);			// Component display radius
	double			linkradius(0);			// Link display radius
	Vector3<double>	bbox_center;			// Bounding box center
	Vector3<double>	bbox_size;				// Bounding box size
	long			max_iter(0);			// Maximum number of iterations
	double			beta(10);				// Equivalent of 1/kT
	double			max_angle(M_PI/180.0);	// Maximum angular deviation
	double			max_shift(1);			// Maximum allowed shift
	int				link_steps(0);			// Number of steps along link
	Vector3<double>	location;				// Point for harmonic force
	double			ss(0);					// SS link length
	int				type_select(0);			// Selection
	string			mapfile;				// Input map file
	string			maskfile;				// Input map file
	string			paramfile;				// Parameter file
 	string			outfile;				// Output model file
    
	random_seed();

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	Bmodparam	md;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
	}
	
	if ( paramfile.length() )
		md = read_dynamics_parameters(paramfile);

	md.rigid = 0;			// All models rigid
	md.distancetype = 2;	// Soft potential

	for ( curropt = option; curropt; curropt = curropt->next ) {
 		if ( curropt->tag == "sampling" )
        	sam = curropt->scale();
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 0.001 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 0.001 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "rigid" ) {
			if ( curropt->value[0] == 'm' ) md.rigid = 1;
			if ( curropt->value[0] == 'c' ) md.rigid = 2;
		}
		if ( curropt->tag == "bbox" )
			if ( curropt->box(bbox_center, bbox_size) < 6 )
				cerr << "-bbox: All 6 values must be specified!" << endl;
		if ( curropt->tag == "wrap" )
			md.wrap = 1;
		if ( curropt->tag == "iterations" )
			if ( ( max_iter = curropt->value.integer() ) < 1 )
				cerr << "-iterations: A number must be specified!" << endl;
		if ( curropt->tag == "Klink" )
			if ( ( md.Klink = curropt->value.real() ) < 1e-30 )
				cerr << "-Klink: The link strength must be specified!" << endl;
		if ( curropt->tag == "Kangle" )
			if ( ( md.Kangle = curropt->value.real() ) < 1e-30 )
				cerr << "-Kangle: The angle strength must be specified!" << endl;
		if ( curropt->tag == "Kelectrostatic" )
			if ( ( md.Kelec = curropt->value.real() ) < 1e-30 )
				cerr << "-Kelectrostatic: The electrostatic strength must be specified!" << endl;
		if ( curropt->tag == "Kdistance" )
			if ( ( md.Kdistance = curropt->value.real() ) < 1e-30 )
				cerr << "-Kdistance: The distance force constant must be specified!" << endl;
		if ( curropt->tag == "Kseparation" )
			if ( ( md.Ksep = curropt->value.real() ) < 1e-30 )
				cerr << "-Kseparation: The separation energy constant strength must be specified!" << endl;
		if ( curropt->tag == "Kmap" )
			if ( ( md.Kmap = curropt->value.real() ) < 1e-30 )
				cerr << "-Kmap: The map energy constant strength must be specified!" << endl;
		if ( curropt->tag == "separation" )
			if ( ( md.sepdist = curropt->value.real() ) < 1e-30 )
				cerr << "-separation: The separation distance must be specified!" << endl;
		if ( curropt->tag == "cutoff" )
			if ( ( md.cutoff = curropt->value.real() ) < 1e-30 )
				cerr << "-cutoff: The cutoff distance must be specified!" << endl;
		if ( curropt->tag == "beta" )
			if ( ( beta = curropt->value.real() ) < 1e-30 )
				cerr << "-beta: A value must be specified!" << endl;
		if ( curropt->tag == "angle" ) {
			if ( ( max_angle = curropt->value.real() ) < 1e-30 )
				cerr << "-angle: An angle must be specified!" << endl;
			else
				max_angle = angle_set_negPI_to_PI(max_angle*M_PI/180.0);
		}
		if ( curropt->tag == "shift" )
			if ( ( max_shift = curropt->value.real() ) < 1e-30 )
				cerr << "-shift: A distance must be specified!" << endl;
		if ( curropt->tag == "linksteps" )
			if ( ( link_steps = curropt->value.integer() ) < 1 )
				cerr << "-linksteps: The link step length must be specified!" << endl;
		if ( curropt->tag == "location" )
			md.point = curropt->vector3();
		if ( curropt->tag == "Klocation" )
			if ( curropt->values(md.Kpoint, md.pointdecay) < 1 )
				cerr << "-Klocation: The location force constant must be specified!" << endl;
		if ( curropt->tag == "SS" )
			if ( ( ss = curropt->value.real() ) < 1e-30 )
				cerr << "-SS: A link length must be specified!" << endl;
		if ( curropt->tag == "Map" )
			mapfile = curropt->filename().str();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
    }
	option_kill(option);

	double		ti = timer_start();
	
	// Read all the model files
	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( file_list.size() < 1 ) {
		cerr << "Error: No model files specified!" << endl;
		bexit(-1);
	}

	Bmodel*			model = read_model(file_list, paramfile, type_select);
	if ( !model ) {
		cerr << "Error: No models read!" << endl;
		bexit(-1);
	}

	if ( compradius > 0 ) models_set_component_radius(model, compradius);

	if ( linkradius > 0 ) models_set_link_radius(model, linkradius);
	
	if ( ss ) md.add_linktype("S", "S", ss);

	if ( md.rigid > 1 )
		models_link_list_generate(model, 2);
	
	model_update_reference_parameters(model, md);
	
	vector<Vector3<double>>	bounds = models_calculate_bounds(model);

	cout << bbox_center << tab << bbox_size << endl;
	if ( bbox_size.volume() < 1 ) {
		md.min = bounds[0];
		md.max = bounds[1];
	} else {
		md.min = bbox_center - (bbox_size * 0.5);
		md.max = bbox_center + (bbox_size * 0.5);
	}
	cout << md.min << tab << md.max << endl;
	
	// Read the map file
	long		natom_in(0);
	Bimage*		map = NULL;
	Bimage*		pmask = NULL;
	
	if ( mapfile.length() ) {
		map = read_img(mapfile, 1, 0);
		if ( sam.volume() ) map->sampling(sam);
		if ( bbox_size.volume() < 1 ) {
			bbox_size = map->sampling(0)*map->size();
			md.min = -map->image->origin()*map->sampling(0);
			md.max = md.min + bbox_size;
		}
		natom_in = model_test_if_within_box(model, md.min, md.max);
		if ( natom_in < 2 ) {
			cerr << "Error: Model not within map boundaries!" << endl;
			bexit(-1);
		}
		if ( md.Kmap <= 0 ) md.Kmap = 1;
	}

	if ( maskfile.length() ) {
		pmask = read_img(maskfile, 1, 0);
		if ( sam.volume() ) pmask->sampling(sam);
	}

	if ( max_iter ) {
		if ( md.rigid < 2 ) {
			model = monte_carlo_metropolis(model, md, map, beta, max_angle, max_shift,
					max_iter, monte_rigid_body_fit_energy, model_rigid_body_transform);
		} else {
			if ( link_steps > 0 )
				model = monte_carlo_metropolis(model, md, map, beta, max_angle, max_shift, 
						max_iter, monte_link_fit_energy, model_move_components_down_energy);
			else
				model = monte_carlo_metropolis(model, md, map, beta, max_angle, max_shift, 
						max_iter, monte_component_fit_energy, model_move_components_down_energy);
		}
	}

	if ( outfile.length() )
		write_model(outfile, model);
	
	delete map;
	delete pmask;
	delete model;
		
	timer_report(ti);
	
	bexit(0);
}

