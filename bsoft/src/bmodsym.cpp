/**
@file	bmodsym.cpp
@brief	Deals with model and component symmetries.
@author Bernard Heymann
@date	Created: 20060908
@date 	Modified: 20250618
**/

#include "rwmodel.h"
#include "model_util.h"
#include "model_transform.h"
#include "model_symmetry.h"
#include "model_map.h"
#include "model_select.h"
#include "model_links.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodsym [options] in1.star [in2.star...]",
"-----------------------------------------------",
"Manipulates models.",
" ",
"Actions for preparation:",
"-all                     Reset selection to all models and components before other selections.",
" ",
"Selections:",
"-select #232@14          Select models and components.",
" ",
"Actions:",
"-center                  Center before all other operations.",
"-axes                    Show symmetry axes.",
"-setasu D8               Set components to within an asymmetric unit.",
"-apply T                 Apply point group symmetry.",
"-check I                 Check point group symmetry.",
"-find C5                 Find standard orientation for this point group symmetry.",
"-Bfactor D5              Calculate B factors for this point group symmetry.",
"-consolidate C5          Average symmetry-related poisitions.",
"-lattice 2,5,3           Generate a lattice with the number of unit cells in each direction",
"-cyclic 3,12             Find component cyclic symmetry in the given order range.",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components before other selections.",
"-merge                   Merge models before writing.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
"-origin 0,22.5,30        Set the symmetry origin.",
"-view 0,1.5,-0.2,35      View relative to standard view (default 0,0,1,0).",
"-map image.pif,2         Map and image number associated with model.",
" ",
"Parameters for generating a lattice:",
"-unitcell 50,50,50,90,90,90 Unit cell parameters (angstrom & degrees)",
" ",
"Parameters for finding symmetry:",
"-annuli 3,12             Annular range to find cyclic symmetry (pixels).",
"                         If not specified, taken from component radius.",
"-width 10                Annular width to find cyclic symmetry (pixels).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-split 3                 Split models into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": model ID's are used as file names.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	bool 			all(0);						// Keep selection as read from file
	bool			reset(0);					// Keep selection as ouput
	bool			merge(0);					// Flag to merge models
	string			mod_select;					// Model and component selection
	bool			center(0);					// Flag to center the structure
	bool			show_axes(0);				// Flag to show axes
	string			map_name;					// Density map reference
	int				img_num(0);					// Image number in density map file
	Vector3<double>	origin;						// Translate
	long			minorder(0), maxorder(0);	// Order range for finding cyclic symmetry
	long			ann_min(0), ann_max(0);		// Annular range in pixels for finding cyclic symmetry
	long			ann_width(0);				// Annular width in pixels for finding cyclic symmetry
	string			action;						// Symmetry operation: asu, app, chk, fnd, Bfc, con
	Bsymmetry		sym;						// Symmetry for various options
	Vector3<long>	lattice;					// Crystal lattice size
	UnitCell		uc;							// Unit cell parameters for lattice
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
	View2<double>	view;						// Reference view
    string    		atom_select("all");
	string			paramfile;					// Input parameter file name
	string			outfile;					// Output parameter file name
	int				splt(0);					// Sets output of multiple single-model files
	string			astr;
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value.str();
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "map" ) {
			vector<string>	vs = split(curropt->value.str(), ',');
			map_name = vs[0];
			if ( vs.size() > 1 )
			img_num = to_integer(vs[1]);
			if ( map_name.length() < 1 )
				cerr << "-map: A file name must be specified!" << endl;
		}
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
		if ( curropt->tag == "axes" ) show_axes = 1;
		if ( curropt->tag == "setasu" ) {
			sym = curropt->symmetry();
			action = "asu";
		}
		if ( curropt->tag == "apply" ) {
			sym = curropt->symmetry();
			action = "app";
		}
		if ( curropt->tag == "check" ) {
			sym = curropt->symmetry();
			action = "chk";
		}
		if ( curropt->tag == "find" ) {
			sym = curropt->symmetry();
			action = "fnd";
		}
		if ( curropt->tag == "consolidate" ) {
			sym = curropt->symmetry();
			action = "con";
		}
		if ( curropt->tag == "Bfactor" ) {
			sym = curropt->symmetry();
			action = "Bfc";
		}
		if ( curropt->tag == "lattice" ) {
			lattice = curropt->vector3();
			if ( lattice.volume() < 1 )
				cerr << "-lattice: Three values must be specified" << endl;
		}
		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
		if ( curropt->tag == "cyclic" )
			if ( curropt->values(minorder, maxorder) < 2 )
				cerr << "-cyclic: Minimum and maximum orders must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 1 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "annuli" )
			if ( curropt->values(ann_min, ann_max) < 1 )
				cerr << "-annuli: An annular range must be specified!" << endl;
		if ( curropt->tag == "width" )
			if ( ( ann_width = curropt->value.integer() ) < 1 )
				cerr << "-width: An annular width must be specified!" << endl;
		if ( curropt->tag == "view" )
			view = curropt->view();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) splt = 9;
			else if ( ( splt = curropt->value.integer() ) < 1 )
				cerr << "-splt: An integer must be specified!" << endl;
			else
				if ( splt > 6 ) splt = 6;
		}
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	// Read all the parameter files
	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( file_list.size() < 1 ) {
		cerr << "Error: No model files specified!" << endl;
		bexit(-1);
	}

	Bmodel*			model = read_model(file_list, paramfile);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( all ) models_select_all(model);

	if ( mod_select.length() ) models_select(model, mod_select);
	
	if ( map_name.length() ) {
		model->mapfile(map_name);
		model->image_number(img_num);
	}

	if ( comprad > 0 ) models_set_component_radius(model, comprad);

	if ( linkrad > 0 ) models_set_link_radius(model, linkrad);

	if ( center ) models_center(model);

	if ( show_axes ) model_symmetry_axes(model);
	
	if ( action == "asu" ) models_find_asymmetric_unit(model, sym);
	
	if ( action == "chk" )
		model_symmetry_RMSD(model, sym, view);
	else if ( action == "con" )
		models_symmetrize(model, sym);
	else if ( action == "app" )
		models_apply_point_group(model, sym, origin, view);
	else if ( action == "fnd" ) {
		model_orient_to_standard_view(model, sym, view);
	} else if ( action == "Bfc" )
		model_symmetry_B(model, sym, view);

	if ( lattice.volume() > 1 )
		models_generate_lattice(model, uc, lattice);
		
	if ( maxorder )
		model_component_symmetry(model, 360, ann_min, ann_max, ann_width, 0, 0, 0, minorder, maxorder);
		
	if ( reset ) models_select_all(model);

	if ( merge ) model_merge(model);
	
	models_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( model && ( outfile.length() || splt == 9 ) )
		write_model(outfile, model, splt);

	delete model;
		
	timer_report(ti);
	
	bexit(0);
}

