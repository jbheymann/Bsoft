/**
@file	bmodel.cpp
@brief	Manipulates models.
@author Bernard Heymann
@date	Created: 20060908
@date 	Modified: 20250617
**/

#include "rwmodel.h"
#include "rwsequence.h"
#include "model_util.h"
#include "model_shell.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_views.h"
#include "model_mol.h"
#include "model_links.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bmodel [options] in1.star [in2.star...]",
"----------------------------------------------",
"Manipulates models.",
" ",
"Actions for preparation:",
"-addnumber               Add numbers to the model id's.",
"-merge                   Merge models in different files rather than concatenate.",
"-all                     Reset selection to all models and components before other selections.",
"-type 1                  Component type designation: 0=from file, 1=elements, 2=atoms.",
" ",
"Selections:",
"-select #232@14          Select models and components.",
"-number 24,36            Select models based on the number of components.",
"-closed order,3          Select models based on valency (valency,<n>) or polygon order (order,<n>).",
"-fullerene               Select fullerene type models.",
"-fom 0.25                Select: FOM cutoff.",
" ",
"Actions:",
"-listmodels              List model information in table form.",
"-listcomponents          List component counts in table form.",
"-listsequences           List molecular sequences in table form.",
"-showmass                Show model masses.",
"-axes                    Calculate principal axes.",
"-first B                 Reorder models, setting the indicated model first.",
"-reorder B,C,A,Z         Copy and reorder selected models.",
"-delete HO               Delete a component type.",
"-setfom 0.8              Set the FOM of selected components to a specific value.",
"-center                  Center before all other operations.",
"-average 4               Average sequential components and delete all but one in each set.",
"-histogram 0.2           FOM histogram with the given step size.",
"-id model1               Set model identifier.",
"-map image.pif,2         Map and image number associated with model.",
"-path ../map             Path to map.",
"-views local             Calculate views for components. Modes: origin: current origin,",
"                         com: center-of-mass origin, map: map origin, local: neigbor plane.",
"-invert                  Invert views for selected components.",
"-rdf 1.2                 Calculate the radial distribution function with the given sampling (angstrom).",
"-curvature               Calculate the curvature associated with each component.",
"-settype VER             Set all component types.",
"-associate TRS,trs.pdb   Associate a component type with a file name.",
"-mass TRS,288.7          Associate a component type with a mass.",
"-translate 12.3,95,-10.5 Shift along a vector (x,y,z).",
"-rotate 0.5,0,5.8,45     Rotate around a vector (x,y,z) by an angle.",
"-toview 0.3,0.5,0.8,33   Rotate to view: vector {xyz} and angle (default 0,0,1,0).",
"-fromview 0.3,0.5,0.8,33 Rotate from view: vector {xyz} and angle (default 0,0,1,0).",
"-scale 1.5,0.8,1         Scale around the origin.",
"-reflect 0.4,-.22,0.5    Reflect through a mirror plane defined by the given normal.",
"-linklength 56.3,CA,OH   Generate links below a maximum length between indicated types (all if omitted).",
"-thickness               Calculate effective thickness from z coordinates.",
" ",
"Actions for finishing:",
"-reset                   Reset selection to all components before other selections.",
" ",
"Parameters:",
"-verbose 7               Verbose output.",
"-componentradius 8.4     Set display radius for all components.",
"-linkradius 5.1          Set display radius for all links.",
"-origin 0,22.5,30        Set the origin for rotation and scaling.",
"-separate                Each model is defined as a separate molecule group.",
"-reference 0,1.5,-0.2,35 Reference view (default 0,0,1,0).",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
"-newtypes new.json       New model types to assign (JSON).",
//"-reference file.pdb      Reference coordinate file to set views.",
" ",
"Output:",
"-output file.star        Output model parameter file.",
"-split 3                 Split models into individual files:",
"                         Argument: 1-6: number of digits inserted before extension",
"                         Argument: \"id\": model ID's are used as file names.",
"-types all.json          Output model types file.",
"-sequences all.fasta     Output sequence file.",
" ",
NULL
};

int 	main(int argc, char **argv)
{
	/* Initialize variables */
	int 			num_id(0);					// Flag to add numbers to model id's
	int 			all(0);						// Keep selection as read from file
	int 			reset(0);					// Keep selection as ouput
	int				type_select(0);				// Component type designation
	Bstring			mod_select;					// Model and component selection
	bool			list_models(0);				// List model information in table form
	bool			list_comps(0);				// List component information in table form
	bool			list_seq(0);				// List molecular sequences in table form
	int				show_mass(0);				// Show model masses
	int				show_axes(0);				// Calculate principal axes
	string			id_first;					// ID of first model in reordering
	string			reorder;					// String for copying and reordering
	int				merge(0);					// Flag to merge models rather than concatenate
	long			ncomp_min(0), ncomp_max(0);	// Minimum and maximum number of components to select for
	int				closure_rule(0);			// Closure rule: 1=valency, 2=order
	int				val_order(0);				// Valency or order - depending on rule
	int				fullerene(0);				// Flag to select fullerenes
	Bstring			comp_delete;				// Component type selection for deletion
	double			setfom = -1;				// Value to set FOM to
	long			average(0);					// Number of sequential components to average
	double			fom_cutoff(0);				// No selection based on FOM
	double			fom_step(0);				// Histogram step size
	string			model_id;					// Model identifier
	string			map_name;					// Density map reference
	long			img_num(0);					// Image number in density map file
	string			map_path;					// Density map path
	string			calc_views;					// Mode to calculate component views
	int				inv_views(0);				// Flag to invert views for selected components
	double			rdf_interval(0);			// Calculate RDF at this sampling (0=not)
	int				curvature(0);				// Flag to calculate curvature
	Bstring			associate_type;				// Component type
	Bstring			associate_file;				// Component type file name
	double			associate_mass(0);			// Component type mass
	Bstring			set_type;					// Component type to set
	double			comprad(0);					// Component display radius
	double			linkrad(0);					// Link display radius
	double			linklength(0);				// Link length for generating links
	Bstring			link_type1, link_type2;		// Component type IDs for specifying links
	int				link_closest(0);			// Only generate one link per first component
	bool			center(0);					// Flag to center model
    Transform		t;							// Transformation details
    View2<double>	view;						// View for rotations
	Vector3<double>	shift;						// Translate
	Vector3<double>	scale;						// Scale around origin
	Vector3<double>	reflect;					// Mirror plane normal
	View2<double>	ref_view;					// Reference view
	bool			eff_thick(0);				// Flag for effective thickness calculation
	string			paramfile;					// Input parameter file name
	string			nutypefile;					// Input new model types file (JSON)
	string			outfile;					// Output parameter file name
	string			typefile;					// Output model types file (JSON)
	string			seqfile;					// Output sequence file name
	int				splt(0);					// Sets output of multiple single-model files
	Bstring			astr;
    
	int				i, optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "addnumber" ) num_id = 1;
		if ( curropt->tag == "all" ) all = 1;
//		if ( curropt->tag == "type" )
//			if ( ( type_select = curropt->value.integer() ) < 1 )
//				cerr << "-type: A type designation must be specified!" << endl;
		if ( curropt->tag == "type" ) {
			if ( curropt->value[0] == 'e' ) type_select = 1;
			else if ( curropt->value[0] == 'a' ) type_select = 2;
			else type_select = curropt->value.integer();
			if ( type_select < 0 ) type_select = 0;
			if ( type_select > 2 ) type_select = 0;
		}
		if ( curropt->tag == "merge" ) merge = 1;
		if ( curropt->tag == "reset" ) reset = 1;
		if ( curropt->tag == "select" )
			mod_select = curropt->value;
		if ( curropt->tag == "listmodels" ) list_models = 1;
		if ( curropt->tag == "listcomponents" ) list_comps = 1;
		if ( curropt->tag == "listsequences" ) list_seq = 1;
		if ( curropt->tag == "showmass" ) show_mass = 1;
		if ( curropt->tag == "axes" ) show_axes = 1;
		if ( curropt->tag == "first" ) {
			id_first = curropt->value.str();
			if ( id_first.length() < 1 )
				cerr << "-first: An identifier must be specified!" << endl;
		}
		if ( curropt->tag == "reorder" ) {
			reorder = curropt->value.str();
			if ( reorder.length() < 1 )
				cerr << "-reorder: At least one identifier must be specified!" << endl;
		}
		if ( curropt->tag == "number" )
			if ( curropt->values(ncomp_min, ncomp_max) < 1 )
				cerr << "-number: A number of components must be specified!" << endl;
		if ( curropt->tag == "closed" ) {
			if ( curropt->value[0] == 'v' ) closure_rule = 1;
			if ( curropt->value[0] == 'o' ) closure_rule = 2;
			if ( curropt->value.contains(",") )
				val_order = curropt->value.post(',').integer();
		}
		if ( curropt->tag == "fullerene" ) fullerene = 1;
		if ( curropt->tag == "delete" )
			comp_delete = curropt->value;
		if ( curropt->tag == "setfom" )
			if ( ( setfom = curropt->value.real() ) < 0.001 )
				cerr << "-setfom: A FOM value must be specified!" << endl;
		if ( curropt->tag == "average" )
			if ( ( average = curropt->value.integer() ) < 1 )
				cerr << "-average: A number of components must be specified!" << endl;
		if ( curropt->tag == "fom" )
			if ( ( fom_cutoff = curropt->value.real() ) < 0.001 )
				cerr << "-fom: A FOM cutoff must be specified!" << endl;
		if ( curropt->tag == "histogram" )
			if ( ( fom_step = curropt->value.real() ) < 0.001 )
				cerr << "-histogram: A FOM step size must be specified!" << endl;
		if ( curropt->tag == "id" ) {
			model_id = curropt->value.str();
			if ( model_id.length() < 1 )
				cerr << "-id: An identifier must be specified!" << endl;
		}
		if ( curropt->tag == "map" ) {
			map_name = curropt->value.str();
//			img_num = (map_name.post(',')).integer();
//			map_name = map_name.pre(',');
			img_num = to_integer(map_name.substr(map_name.find(',')));
			map_name = map_name.substr(0,map_name.find(','));
			if ( map_name.length() < 1 )
				cerr << "-map: A file name must be specified!" << endl;
		}
		if ( curropt->tag == "path" ) map_path = curropt->value.str();
		if ( curropt->tag == "views" ) calc_views = curropt->value.lower().str();
		if ( curropt->tag == "invert" ) inv_views = 1;
		if ( curropt->tag == "rdf" )
			if ( ( rdf_interval = curropt->value.real() ) < 0.1 )
				cerr << "-rdf: A sampling interval must be specified!" << endl;
		if ( curropt->tag == "curvature" ) curvature = 1;
		if ( curropt->tag == "associate" ) {
			astr = curropt->value;
			associate_type = astr.pre(',');
			associate_file = astr.post(',');
			astr = 0;
		}
		if ( curropt->tag == "mass" ) {
			astr = curropt->value;
			associate_type = astr.pre(',');
			astr = astr.post(',');
//			associate_mass = get_option_mass(astr);
			associate_mass = curropt->real_unit(astr);
			astr = 0;
		}
		if ( curropt->tag == "settype" )
			set_type = curropt->value;
		if ( curropt->tag == "origin" )
			t.origin = curropt->origin();
		if ( curropt->tag == "reference" )
			ref_view = curropt->view();
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "translate" ) {
			shift = curropt->vector3();
        	if ( shift.length() < 0.1 )
				cerr << "-translate: Three values must be specified!" << endl;
		}
		if ( curropt->tag == "rotate" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			if ( d.size() < 4 )
				cerr << "-rotate: All three vector elements and an angle must be specified!" << endl;
			else {
				t.axis = Vector3<double>(d[0], d[1], d[2]);
				t.angle = d[3] * M_PI/180.0;
			}
		}
		if ( curropt->tag == "toview" ) {
			view = curropt->view();
			t.axis = view.rotation_axis();
			t.angle = view.rotation_angle();
		}
		if ( curropt->tag == "fromview" ) {
			view = curropt->view();
			t.axis = view.rotation_axis();
			t.angle = -view.rotation_angle();
		}
		if ( curropt->tag == "scale" ) {
			scale = curropt->vector3();
			for ( i=0; i<3 && scale[i]; i++ ) ;
			if ( i < 1 )
				cerr << "-scale: At least one scale value must be specified!" << endl;
			else {
				if ( i == 1 ) scale[2] = scale[1] = scale[0];
				else if ( i == 2 ) scale[2] = 1;
			}
		}
		if ( curropt->tag == "reflect" ) {
			reflect = curropt->vector3();
			if ( reflect.length() < 0 )
				cerr << "-reflect: A mirror plane normal must be specified!" << endl;
		}
		if ( curropt->tag == "linklength" ) {
			if ( curropt->value.contains(",") ) {
				Bstring*	sv = curropt->value.split(",");
				linklength = sv->real();
				if ( sv->next ) {
					link_type1 = *(sv->next);
					link_type1.next = NULL;
					if ( sv->next->next ) {
						link_type2 = *(sv->next->next);
						if ( sv->next->next->next ) link_closest = 1;
					}
				}
				string_kill(sv);
			} else {
				if ( ( linklength = curropt->value.real() ) < 1 )
					cerr << "-linklength: A link length must be specified!" << endl;
			}
		}
		if ( curropt->tag == "componentradius" )
			if ( ( comprad = curropt->value.real() ) < 1 )
				cerr << "-componentradius: A display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkrad = curropt->value.real() ) < 1 )
				cerr << "-linkradius: A radius must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename().str();
		if ( curropt->tag == "thickness" ) eff_thick = 1;
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
		if ( curropt->tag == "split" ) {
			if ( curropt->value.contains("id") || curropt->value.contains("ID") ) splt = 9;
			else if ( ( splt = curropt->value.integer() ) < 1 )
				cerr << "-split: An integer must be specified!" << endl;
			else
				if ( splt > 6 ) splt = 6;
		}
		if ( curropt->tag == "newtypes" )
			nutypefile = curropt->filename().str();
		if ( curropt->tag == "types" )
			typefile = curropt->filename().str();
		if ( curropt->tag == "sequences" ) {
			seqfile = curropt->filename().str();
			list_seq = 1;
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

	// Make sure masses are assigned
	if ( show_mass && type_select < 1 ) type_select = 1;
	
	Bmodel*			model = read_model(file_list, paramfile, type_select);

	if ( !model ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}

//	for ( Bcomponent* comp = model->comp; comp; comp = comp->next )
//		cout << comp->view() << endl;

	if ( num_id ) model_number_ids(model);
	
	if ( all ) models_select_all(model);

	if ( merge ) model_merge(model);
	
	if ( id_first.length() ) models_reorder_circular(&model, id_first);
	
	if ( reorder.length() ) models_copy(&model, reorder);
	
	if ( mod_select.length() ) models_select(model, mod_select);
	
	if ( list_models ) models_list(model);

	if ( list_comps ) models_list_comp(model);
	
	if ( typefile.length() ) {
		JSvalue 	types = model_types(model);
		types.write(typefile);
	}

	if ( nutypefile.length() ) {
		JSvalue 	types = JSparser(nutypefile).parse();
		models_assign_types(model, types);
	}

	if ( list_seq ) {
		vector<Bsequence>	seqs = models_sequence(model);
		if ( seqfile.length() )
			write_sequence(seqfile, seqs);
	}
	
	if ( model_id.length() ) model->identifier(model_id);
	
	if ( map_name.length() ) {
		model->mapfile(map_name);
		model->image_number(img_num);
		model_check(model, map_path);
	} else if ( map_path.length() ) {
		model_check(model, map_path);
	}

	if ( set_type.length() ) models_set_type(model, set_type.str());
	
	if ( associate_file.length() )
		model_associate(model, associate_type.str(), associate_file.str());
	
	if ( associate_mass > 0 )
		model_associate_mass(model, associate_type.str(), associate_mass);
	
	if ( ncomp_min ) models_select_number_of_components(model, ncomp_min, ncomp_max);

	if ( closure_rule ) models_select_closed(model, closure_rule, val_order);

	if ( fullerene ) models_select_fullerene(model);
	
	if ( comp_delete.length() ) models_delete_comp_type(model, comp_delete.str());
	
	if ( setfom >= 0 ) model->set_component_fom(setfom);
	
	if ( linklength > 0 ) models_link_list_generate(model, linklength, link_type1.str(), link_type2.str(), link_closest);
//	if ( linklength > 0 ) models_link_list_generate(model, linklength);
	
	if ( comprad > 0 ) models_set_component_radius(model, comprad);

	if ( linkrad > 0 ) models_set_link_radius(model, linkrad);

	if ( center ) models_center(model);

	if ( fom_cutoff > 0 ) model_fom_deselect(model, fom_cutoff);

	if ( average > 0 ) model_average_components(model, average);
		
	if ( shift.length() ) models_shift(model, shift);

	if ( t.angle ) models_rotate(model, t);

	if ( scale.volume() > 0 ) models_scale(model, scale, t.origin);
	
	if ( reflect.length() > 0 ) models_reflect(model, reflect, t.origin);
	
	if ( fom_step > 0 ) model_fom_histogram(model, fom_step);

	if ( calc_views.length() ) model_calculate_views(model, calc_views);
	
	if ( rdf_interval > 0 ) model_radial_distribution(model, rdf_interval);

	if ( inv_views ) model_invert_views(model);
	
	if ( curvature ) model_curvature(model);

	if ( show_mass ) model_mass_all(model);

	if ( show_axes ) model_principal_axes(model);
	
	if ( eff_thick )
		cout << "Effective thickness:            " << model_effective_thickness(model) << endl << endl;
	
	if ( reset ) models_select_all(model);

	models_selection_stats(model);
	
	// Write an output parameter format file if a name is given
    if ( model && ( outfile.length() || splt == 9 ) )
		write_model(outfile, model, splt);

	delete model;
	
	timer_report(ti);
	
	bexit(0);
}

