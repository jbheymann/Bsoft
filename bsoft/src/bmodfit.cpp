/**
@file	bmodfit.cpp
@brief	A tool to fit models.
@author Bernard Heymann
@date	Created: 20220223
@date	Modified: 20250516
**/

#include "rwmodel.h"
#include "rwimg.h"
#include "model_compare.h"
#include "model_mol.h"
#include "model_mechanics.h"
#include "model_select.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodfit [options] in.star",
"--------------------------------",
"Calculate fits of model to a reference model or map.",
" ",
"Actions:",
"-fit mod_id              Fit a model with the given id to a reference.",
"-aligned                 Align and fit a model to a reference.",
"-scale 0.04,0.001        Fit a map: maximum scale adjustment and increment.",
" ",
"Selections:",
"-select #Mod1@14         Select models and components.",
"-type 2                  Type for fit and transformation: 0,1,2.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-gap 10,0.2              Gap opening and extending penalties for alignment (default 20,0.2).",
" ",
"Input:",
"-reference ref.pdb       Input reference model file.",
"-map ref.pdb             Input map file.",
"-Similarity sim.star     Residue similarity matrix (default blosum62.star).",
" ",
"Output:",
"-output newmod.star      Output model file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	string			fit_id;						// Model ID to fit
	bool			aligned(0);					// Flag to alig model sequences
	double			scale_max(0), scale_inc(0);	// Maximu scale for map fit
	Bstring			mod_select;					// Model and component selection
	int				fit_type(0);				// Flag to select type of fit and transformation
	double			gapopen(20), gapextend(0.2);	// Gap opening and extending penalties
	string			paramfile;					// Input parameter file name
	string			reffile;					// Input reference model file name
	string			mapfile;					// Input map file name
	string			simfile;					// Use default similarity matrix file
	string			outfile;					// Output model file name
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "fit" ) fit_id = curropt->value.str();
		if ( curropt->tag == "aligned" ) aligned = 1;
		if ( curropt->tag == "scale" )
        	if ( curropt->values(scale_max, scale_inc) < 2 )
				cerr << "-scale: A maximum scale adjustment must be specified!" << endl;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "type" ) fit_type = curropt->value.integer();
		if ( curropt->tag == "gap" )
			if ( curropt->values(gapopen, gapextend) < 1 )
				cerr << "-gap: At least one gap penalty must be specified!" << endl;
		if ( curropt->tag == "reference" )
			reffile = curropt->filename().str();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename().str();
		if ( curropt->tag == "Similarity" )
			simfile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
    }
	option_kill(option);
	
	double			ti = timer_start();

	// Read all the parameter files
	Bmodel*			model = NULL;
	Bmodel*			refmod = NULL;
	Bimage*			map = NULL;

	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( file_list.size() )
		model = read_model(file_list, paramfile);
	
	if ( verbose )
		cout << "Models: " << model->count() << endl;
	
	if ( reffile.length() ) {
		refmod = read_model(reffile, paramfile);
		if ( !refmod ) {
			error_show("Error: No reference model read!", __FILE__, __LINE__);
			bexit(-1);
		}
	}

	if ( mapfile.length() ) {
		map = read_img(mapfile, 1, 0);
//		if ( sam.volume() ) map->sampling(sam);
	}

	Bresidue_matrix	simat;
	
	if ( aligned )
		simat = get_residue_matrix(simfile);	

	if ( mod_select.length() ) models_select(model, mod_select);

	if ( refmod ) {
//		if ( aligned )
//			models_select_aligned_CA(model, refmod, gapopen, gapextend, simat);
//		else
//			models_select_corresponding_CA(model, refmod);
		if ( fit_id.length() )
			model_fit(model, refmod, fit_id);
		else
			models_fit_CA(model, refmod, fit_type);
	}

	if ( map && scale_max > 0 ) model_find_map_scale(model, map, scale_max, scale_inc);

	if ( verbose )
		cout << "Models: " << model->count() << endl;
	
	// Write an output model file if a name is given
    if ( outfile.length() && model )
		write_model(outfile, model);

	delete model;
	delete refmod;
	delete map;
	
	timer_report(ti);
	
	bexit(0);
}


