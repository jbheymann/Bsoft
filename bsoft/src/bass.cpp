/**
@file	bass.cpp
@brief	Program to assemble molecular components.
@author Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250318
**/

#include "rwmodel.h"
#include "model_assembly.h"
#include "model_select.h"
#include "model_views.h"
#include "model_links.h"
#include "model_util.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bass [options] in1.star [in2.star ...]",
"---------------------------------------------",
"Manipulate an assembly of models.",
" ",
"Actions for preparation:",
"-create                  Generate an assembly model from a set of models.",
"-assemble                Generate a full model from an assembly model.",
"-all                     Reset selection to all components before other selections.",
" ",
"Selections:",
"-select #232@14          Select models and components.",
"-first 22                Select the first number of components.",
" ",
"Actions:",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-associate TRS,trs.pdb   Associate a component type with a file name.",
"-untangle 3.5,0.2        Eliminate overlaps by moving molecules apart (sampling and damping factor).",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-separate                Each model is defined as a separate molecule group.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components before other selections.",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-coordinates all.pdb     Output coordinate files.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	bool			create(0);					// Flag to generate an assembly model
	bool			assemble(0);				// Flag to generate a full model
	bool 			all(0);						// Keep selection as read from file
	bool 			reset(0);					// Keep selection as ouput
	Bstring			mod_select;					// Model and component selection
	int				first(0);					// First number of components to select
	Bstring			calc_views;					// Mode to calculate component views
	string			associate_type;				// Component type
	string			associate_file;				// Component file name
	double			untangle(0);				// Untangling grid sampling
	double			lambda(0.1);				// Untangling damping factor
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
	string			paramfile;					// Input parameter file name
	string			outfile;					// Output parameter file name
    
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" ) create = 1;
		if ( curropt->tag == "assemble" ) assemble = 1;
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "first" )
			if ( ( first = curropt->value.integer() ) < 1 )
				cerr << "-first: An integer must be specified!" << endl;
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower();
		if ( curropt->tag == "associate" ) {
			vector<string>	vs = split(curropt->value.str(), ',');
			associate_type = vs[0];
			associate_file = vs[1];
		}
		if ( curropt->tag == "untangle" )
			if ( curropt->values(untangle, lambda) < 1 )
				cerr << "-untangle: Grid sampling in angstrom must be specified!" << endl;
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 1 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
	}
	option_kill(option);
	
	double			ti = timer_start();
	
	Bmodel*			model = NULL;

	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( file_list.size() < 1 ) {
		cerr << "Error: No model files specified!" << endl;
		bexit(-1);
	}

	if ( create )
		model = model_generate_assembly(file_list, paramfile);
	else
		model = read_model(file_list, paramfile);
	
	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( all ) models_select_all(model);

	if ( mod_select.length() ) models_select(model, mod_select);
	
	if ( first ) models_select_first(model, first);
	
	if ( associate_file.length() )
		model_associate(model, associate_type, associate_file);
	
	if ( comprad > 0 ) models_set_component_radius(model, comprad);

	if ( linkrad > 0 ) models_set_link_radius(model, linkrad);

	if ( calc_views.length() ) model_calculate_views(model, calc_views);
	
	if ( assemble ) {
		Bmodel*		numod = model_assemble(model, paramfile);
		delete model;
		model = numod;
	}

	if ( reset ) models_select_all(model);

	models_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile, model);
	}

	delete model;
		
	timer_report(ti);
	
	bexit(0);
}

