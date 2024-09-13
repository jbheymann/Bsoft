/**
@file	bmodsym.cpp
@brief	Deals with model and component symmetries.
@author Bernard Heymann
@date	Created: 20060908
@date 	Modified: 20221115
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
"-map image.pif,2         Map and image number associated with model.",
"-setasu D8               Set components to within an asymmetric unit.",
"-apply T                 Apply point group symmetry.",
"-symmetrize C5           Symmetrize for a point group symmetry.",
"-lattice 2,5,3           Generate a lattice with the number of unit cells in each direction",
"-separate D7             Generate separate symmetry-related models.",
"-find 3,12               Find component cyclic symmetry in the given order range.",
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
	int 			all(0);						// Keep selection as read from file
	int 			reset(0);					// Keep selection as ouput
	int				merge(0);					// Flag to merge models
	Bstring			mod_select;					// Model and component selection
	int				center(0);					// Flag to center the structure
	Bstring			map_name;					// Density map reference
	int				img_num(0);					// Image number in density map file
	Vector3<double>	origin;						// Translate
	long			minorder(0), maxorder(0);	// Order range for finding cyclic symmetry
	long			ann_min(0), ann_max(0);		// Annular range in pixels for finding cyclic symmetry
	long			ann_width(0);				// Annular width in pixels for finding cyclic symmetry
	string			asu_sym;					// Set coordinates within the ASU
	string			symmetry_apply_string;		// Point group string
	string			symmetrize_string;			// Point group string
	Vector3<long>	lattice;					// Crystal lattice size
	UnitCell		uc;							// Unit cell parameters for lattice
	string			separate_symmetrize;		// Point group string
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
	View2<double>	ref_view;					// Reference view
    Bstring    		atom_select("all");
	Bstring			paramfile;					// Input parameter file name
	Bstring			outfile;					// Output parameter file name
	int				split(0);					// Sets output of multiple single-model files
	Bstring			astr;
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "map" ) {
			map_name = curropt->value;
			img_num = (map_name.post(',')).integer();
			map_name = map_name.pre(',');
			if ( map_name.length() < 1 )
				cerr << "-map: A file name must be specified!" << endl;
		}
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "setasu" )
			asu_sym = curropt->symmetry_string();
		if ( curropt->tag == "apply" )
			symmetry_apply_string = curropt->symmetry_string();
		if ( curropt->tag == "symmetrize" )
			symmetrize_string = curropt->symmetry_string();
		if ( curropt->tag == "lattice" ) {
			lattice = curropt->vector3();
			if ( lattice.volume() < 1 )
				cerr << "-lattice: Three values must be specified" << endl;
		}
		if ( curropt->tag == "unitcell" )
			uc = curropt->unit_cell();
		if ( curropt->tag == "separate" )
			separate_symmetrize = curropt->symmetry_string();
		if ( curropt->tag == "find" )
			if ( curropt->values(minorder, maxorder) < 2 )
				cerr << "-find: Minimum and maximum orders must be specified!" << endl;
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
		if ( curropt->tag == "reference" )
			ref_view = curropt->view();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) split = 9;
			else if ( ( split = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( split > 6 ) split = 6;
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

	Bmodel*			model = read_model(file_list, paramfile.str());		

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( all ) models_process(model, model_reset_selection);

	if ( mod_select.length() ) model_select(model, mod_select);
	
	if ( map_name.length() ) {
		model->mapfile(map_name.str());
		model->image_number(img_num);
	}

	if ( comprad > 0 ) models_process(model, comprad, model_set_component_radius);

	if ( linkrad > 0 ) models_process(model, linkrad, model_set_link_radius);

	if ( center ) models_process(model, model_center);

	if ( asu_sym.size() ) models_process(model, asu_sym, model_find_asymmetric_unit);
	
	if ( symmetry_apply_string.size() )
		model_apply_point_group(model, symmetry_apply_string, origin, ref_view);
	
	if ( symmetrize_string.size() )
		models_process(model, symmetrize_string, model_symmetrize);

	if ( lattice.volume() > 1 )
		model_generate_lattice(model, uc, lattice);
		
	if ( separate_symmetrize.size() )
		model_symmetry_related(model, separate_symmetrize);
	
	if ( maxorder )
		model_component_symmetry(model, 360, ann_min, ann_max, ann_width, 0, 0, 0, minorder, maxorder);
		
	if ( reset ) models_process(model, model_reset_selection);

	if ( merge ) model_merge(model);
	
	model_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( model && ( outfile.length() || split == 9 ) )
		write_model(outfile.str(), model, split);

	model_kill(model);
		
	
		timer_report(ti);
	
	bexit(0);
}

