/**
@file	bmodcomp.cpp
@brief	A tool to compare polyhedra.
@author Bernard Heymann
@date	Created: 20080102
@date 	Modified: 20250623
**/

#include "rwmodel.h"
#include "model_select.h"
#include "model_poly.h"
#include "model_compare.h"
#include "model_mol.h"
#include "model_util.h"
#include "rwimg.h"
#include "rwresprop.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodcomp [options] in.star",
"---------------------------------",
"Compares models internally or to reference models.",
" ",
"Actions:",
"-all                     Select all models for comparison.",
"-select #Mod1@14         Select models and components.",
"-distance                Calculate a distance matrix.",
"-corresponding           Calculate the RMSDs between corresponding models.",
"-interface 2.5           List interface components based on cutoff distance.",
"-orientations            Calculate relative orientations of all models.",
"-unknown                 Select unknown models for comparison.",
"-closed order,3          Select based on valency (valency,<n>) or polygon order (order,<n>).",
"-consensus 12.5          Calculate a consensus model: components are considered the same if within the given distance.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
" ",
"Input:",
"-parameters parm.star    Molecular parameter file (default atom_prop.star).",
"-reference ref.star      File with reference models.",
" ",
"Output:",
"-output new.star         Output model file.",
"-writereference ref.star Reference output model file.",
"-matrix file.mat         Distance matrix.",
"-image file.map          Distance matrix as an image.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
 	int				all(0);					// Flag to select all models
	Bstring			mod_select;				// Model and component selection
	bool			corresponding;			// Falg to calculate corresponding RMSDs
 	bool			calc_dist(0);			// Flag to calculate a distance matrix
 	bool			unknown(0);				// Flag to select unknown models
 	bool			orientations(0);		// Flag to calculate relative orientations
	int				closure_rule(0);		// Closure rule: 1=valency, 2=order
	int				val_order(0);			// Valency or order - depending on rule
	double			cons_dist(0);			// Cutoff distance for consensus model
	double			interface_dist(-1);		// Cutoff distance for interfaces
	Bstring    		atom_select("ALL");
	string			paramfile;				// Use default parameter file
	string			reffile;				// Reference model file
	string			outfile;				// Output model file
	string			refoutfile;				// Reference output model file
	string			matfile;				// Output matrix file
	string			imgfile;				// Output image file
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "corresponding" ) corresponding = 1;
		if ( curropt->tag == "distance" ) calc_dist = 1;
		if ( curropt->tag == "unknown" ) unknown = 1;
		if ( curropt->tag == "orientations" ) orientations = 1;
		if ( curropt->tag == "closed" ) {
			if ( curropt->value[0] == 'v' ) closure_rule = 1;
			if ( curropt->value[0] == 'o' ) closure_rule = 2;
			if ( curropt->value.contains(",") )
				val_order = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "consensus" )
			if ( ( cons_dist = curropt->real() ) < 1 )
				cerr << "-consensus: A cutoff distance must be specified!" << endl;
		if ( curropt->tag == "interface" )
			if ( ( interface_dist = curropt->real() ) < 0 )
				cerr << "-interface: A cutoff distance must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
		if ( curropt->tag == "reference" )
			reffile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
		if ( curropt->tag == "writereference" )
			refoutfile = curropt->filename().str();
		if ( curropt->tag == "matrix" )
			matfile = curropt->filename().str();
		if ( curropt->tag == "image" )
			imgfile = curropt->filename().str();
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

	Bmodel*		model = read_model(file_list, paramfile);		

	if ( all ) models_select_all(model);

	if ( mod_select.length() ) models_select(model, mod_select);
	
	Bmodel*		refmodel = NULL;		
	if ( reffile.length() ) {
		refmodel = read_model(reffile, paramfile);
		if ( all ) models_select_all(refmodel);
		if ( mod_select.length() ) models_select(refmodel, mod_select);
	}
	
//	if ( !model->poly ) model_poly_generate(model);

	Matrix		mat;

	if ( calc_dist ) {
		if ( refmodel )
			mat = model_distance_matrix(model, refmodel);
		else
			mat = model_distance_matrix(model, 0);
//		for ( mp = model; mp; mp = mp->next ) {
//			mat = model_distance_matrix(mp, 0);
//			cout << mp->identifier() << endl << mat << endl;
//			cout << mp->identifier() << tab << mp->mapfile() << tab << mat[0][1] << endl;
//		}
	}
	
	if ( unknown ) {
		models_select_unknowns(model);
		refmodel = model;
	}
	
	if ( closure_rule ) models_select_closed(model, closure_rule, val_order);

	map<string,Bresidue_type>	rtv = get_residue_properties_code3(paramfile);
	
	if ( refmodel ) {
		if ( model->poly )
			model_poly_compare(model, refmodel);
		else if ( interface_dist > 0 )
			model_interface(model, refmodel, interface_dist);
		else if ( interface_dist == 0 )
			model_interface(model, refmodel, rtv);
		else if ( corresponding )
			models_compare_corresponding(model, refmodel);
		else
			model_compare_by_distance(model, refmodel);
	} else if ( interface_dist > 0 ) {
		mat = models_interfaces(model, model, interface_dist);
	} else if ( orientations ) {
		models_compare_orientations(model);
	}
	
	if ( imgfile.length() && mat.rows() ) {
		Bimage*		pimg = new Bimage(mat, 1);
//		pimg->change_type(nudatatype);
		write_img(imgfile, pimg, 0);
		delete pimg;
	}

	Bmodel*		numod = NULL;
	if ( cons_dist > 0 ) {
		numod = models_consensus(model, cons_dist);
		delete model;
		model = numod;
	}
	
	if ( outfile.length() ) {
		write_model(outfile, model);
	}

	if ( refoutfile.length() ) {
		write_model(refoutfile, refmodel);
	}

	if ( model != refmodel ) delete model;
	delete refmodel;
	
	timer_report(ti);
	
	bexit(0);
}

