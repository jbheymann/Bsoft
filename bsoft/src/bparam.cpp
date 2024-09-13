/**
@file	bparam.cpp
@brief	A tool to extract parameters from coordinate files
@author Bernard Heymann
@date	Created: 20050304
@date	Modified: 20230706
**/

#include "rwmodel.h"
#include "rwmodel_param.h"
#include "model_links.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bparam [options] param.star",
"----------------------------------",
"Manipulate molecular model parameter files.",
" ",
"Actions:",
//"-select CA               Atom selection (default all).",
//"-links                   List all link lengths.",
//"-angles                  List all angles.",
"-length 1.85             Length to define links.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-type elements           Model type: 0=custom(default), 1=elements, 2=atoms.",
" ",
"Input:",
"-from input.pdb          Extract parameters from an atomic coordinate file.",
"-parameters parm.star    Molecular parameter file (default atom_prop.star)",
" ",
"Output:",
"-output param.star       Output parameter file.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
    int				type_select(0);			// Model type: 0=custom, 1=elements, 2=atoms
//	string    		atom_select("ALL");
//	int				show_links(0);
//	int				show_angles(0);
	double			length_cutoff(0);
	string			paramfile;					// Use default parameter file
	string			modelfile;					// Atom coordinate file
	string			outfile;					// Output parameter file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'e' ) type_select = 1;
			else if ( curropt->value[0] == 'a' ) type_select = 2;
			else type_select = curropt->value.integer();
			if ( type_select < 0 ) type_select = 0;
			if ( type_select > 2 ) type_select = 0;
		}
//		if ( curropt->tag == "select" ) {
//			atom_select = curropt->value;
//			if ( atom_select.length() < 1 )
//				cerr << "-select: A selection must be specified!" << endl;
//		}
//		if ( curropt->tag == "links" ) show_links = 1;
//		if ( curropt->tag == "angles" ) show_angles = 1;
		if ( curropt->tag == "length" )
			length_cutoff = curropt->value.real();
//		if ( curropt->tag == "elements" ) elements = 1;
		if ( curropt->tag == "from" )
			modelfile = curropt->filename().str();
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	Bmodparam	md;
	string		filename;
	
	if ( optind < argc ) {
		filename = argv[optind];
		md = read_dynamics_parameters(filename);
		if ( md.comptype.size() < 1 ) {
			cerr << "Error: Parameter file " << filename << " not read!" << endl;
			bexit(-1);
		}
	}

 	if ( modelfile.length() ) {
		Bmodel*		model = read_model(modelfile, paramfile, type_select);
		if ( !model ) {
			cerr << "Error: Model file " << modelfile << " not read!" << endl;
			bexit(-1);
		}
		if ( length_cutoff ) {
			model->clear_links();
			model_link_list_generate(model, length_cutoff);
		}
//		md = md_calculate_parameters(molgroup, elements, show_links, show_angles);
		md = model_param_generate(model);
	}
	
	if ( outfile.length() )
		write_dynamics_parameters(outfile, md);
	
	timer_report(ti);
	
	bexit(0);
}

