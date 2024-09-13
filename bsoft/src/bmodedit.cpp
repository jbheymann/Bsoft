/**
@file	bmodedit.cpp
@brief	A tool to create and edit models.
@author Bernard Heymann
@date	Created: 20090714
@date	Modified: 20220210
**/

#include "rwmodel.h"
#include "model_create.h"
#include "model_map.h"
#include "model_multifit.h"
#include "model_transform.h"
#include "model_compare.h"
#include "model_color.h"
#include "model_select.h"
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
"Usage: bmodedit [options] in.star",
"---------------------------------",
"Creates models.",
"The size of a new model is defined by the -radius option (except -gaussian).",
"If input files are specified, the new model is added or the new components replace old ones.",
" ",
"Actions for creation:",
"-create mod1             New model with the given ID.",
"-symmetry O-3            New model type based on symmetry (T, O and I).",
//"-type tetrahedron        New model type: tetrahedron, cube, octahedron, dodecahedron,",
//"                         icosahedron, helix.",
"-box 120,80,50           Generate model within a box starting at {0,0,0}.",
"-shell 134.2             Generate a spherical shell with the given radius.",
"-fibonacci 244           Generate a fibonacci sphere with the given radius.",
"-icosahedron 238.5,3,1   Generate an icosahedral shell with the given radius and divisions,",
"                         and a flag to make it spherical.",
"-circle 527,243          Generate a circle with the given radius and z offset.",
"-ellipse 620.5,314.9,5.2 Generate an ellipse with the given x and y semi-axes lengths and z offset.",
"-ellipsoid 24.3,51.2,68  Generate an ellipsoid with the given x,y,z semi-axes.",
"-cylinder 88.3,155       Generate a cylinder with the given radius and length.",
"-spindle 88.3,9.5        Generate a spindle with the radius and packing separation.",
"-lattice 3,5,2,1         Generate a lattice with the number of unit cells and a flag for hexagonal packing.",
"-randomfill 56           Random spherical model with the given radius (use with -points & -componentradius).",
"-randomshell 87          Random shell model with the given radius (use with -points & -componentradius).",
"-gaussian 30.7           Random model with gaussian distributed components with the given standard deviation.",
"-peaks cc.map            Generate a model from peaks in a cross correlation map.",
" ",
"Actions:",
"-all                     Reset selection to all models and components before other selections.",
"-replace                 Replace existing components with newly created ones (default add).",
"-reduce submod.star      Replace each set of linked components with one new component and write sub-model file.",
"-id model1               Set model identifier.",
"-settype VER             Set all component types.",
"-changetype VER,KRN      Change component type name.",
"-add -56,3.45,-123       Add a component (define model with -id and type with -settype, can be repeated).",
"-associate TRS,trs.pdb   Associate a component type with a file name.",
"-linklength 56.3         Generate links with this maximum length.",
"-rdf 1.2                 Calculate the radial distribution function with the given sampling (angstrom).",
"-consolidate 85.4        Consolidate models: components are considered the same if within the given distance.",
"-translate 12.3,95,-10.5 Shift along a vector (x,y,z).",
"-trim 200,230,180        Trim to new box size (use with -translate).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-helix 31.6,43.2,35      Helix with given rise, rotation angle and number of components.",
"-direction 0.2,-0.5,0.3  Direction vector (for -cylinder and -spindle).",
"-componentradius 0.5     Component display radius (angstrom).",
"-linkradius 0.5          Link display radius (angstrom).",
"-color 0.1,0.5,0.2,1     Color new models.",
" ",
"Parameters for generating a model:",
"-radius 54.7             Model radius (default 100 angstrom).",
"-points 152              Number of points.",
"-separation 2.5          Separation distance between points.",
"-density 0.3             Density in number per unit volume.",
" ",
"Input:",
"-parameters param.star   Input parameter file.",
"-map file.map            Input map file to associate with shell model.",
" ",
"Output:",
"-output newmod.star      New model file.",
"-catenate catmap.pif     New file name for catenated model maps.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
    /* Initialize variables */
	Bstring			create(0);					// ID of new model
	string			symmetry_string;			// Point group for new model type
	Vector3<double>	box;						// Size of box
	double			shell_radius(-1);			// Radius of sphere
	long			fibo_radius(0);				// Number of components for Fibonacci sphere
	int				ico_divisions(0);			// Number of divisions for icosahedral sphere
	int				sphere_flag(0);				// Flag to make the icosahedral shell spherical
	double			circle_radius(0);			// Circle radius
	double			circle_offset(0);			// Offset of circle in z
	Vector3<double>	ellipse_semi_z;				// Ellipse semi-axes and offset in z
	Vector3<double> ellipsoid_semi;				// Ellipsoid semi-axes
	Vector3<double>	dir;						// Cylinder or spindle direction
	double			cyl_rad(0), cyl_len(0);		// Cylinder radius and length
	double			spindle_rad(0), packing(0);	// Spindle radius and packing separation
	Vector3<long>	lattice;					// Size of lattice
	int				lattice_type(0);			// Lattice type: 0=cubic, 1=hexagonal
	double			density(0);					// Density of components in number per unit volume
	double			gauss_stdev(0);				// Standard deviation for random gaussian model
	double			random_sphere_radius(0);	// Radius for a random sphere model
	double			random_shell_radius(0);		// Radius for a random shell model
	int 			all(0);						// Keep selection as read from file
	int				replace(0);					// Default action is to add models
	Bstring			reduce;						// Reduce linked components and write sub-model
	Bstring			model_id;					// Model identifier
	Bstring			set_type;					// Component type to set
	Bstring			change_type;				// Component type to change and new type name
	Bstring			associate_type;				// Component type
	Bstring			associate_file;				// Component type file name
	Vector3<double>	shift;						// Translate
	Vector3<double>	trim;						// New enclosing box
	Bstring			peakmap;					// Map with peaks to generate a new model
	vector<Vector3<double>>	newloc;				// New component locations
	double			linklength(0);				// Link length for generating links
	double			rdf_interval(0);			// Calculate RDF at this sampling (0=not)
	double			consolidate(0);				// Cutoff to consider components the same
	double 			model_radius(100);			// Model radius
	double			helix_rise(0);				// Helical rise
	double			helix_angle(0);				// Helical rotation angle
	long			helix_comp(0);				// Number of helical components
	double			compradius(0);				// Component display radius
	double			linkradius(0);				// Link display radius
	RGBA<double>	color(1,1,1,1);				// Color specification
	int				set_color(0);				// Flag to set color
	long			points(0);					// Number of points on the sphere
	double			separation(0);				// Distance between points
	Bstring			paramfile;					// Input parameter file name
	Bstring			mapfile;					// Input map file name
	Bstring			outfile;					// Output parameter file name
	Bstring			catname;					// Output catenated map file name
	Bstring			astr;
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "create" ) create = curropt->value;
		if ( curropt->tag == "symmetry" )
			symmetry_string = curropt->symmetry_string();
		if ( curropt->tag == "box" )
			box = curropt->vector3();
		if ( curropt->tag == "shell" )
			if ( ( shell_radius = curropt->value.real() ) < 1 )
				cerr << "-shell: The radius must be specified!" << endl;
		if ( curropt->tag == "fibonacci" )
			if ( ( fibo_radius = curropt->value.real() ) < 1 )
				cerr << "-fibonacci: A radius must be specified!" << endl;
		if ( curropt->tag == "icosahedron" ) {
			if ( curropt->values(shell_radius, ico_divisions, sphere_flag) < 1 )
				cerr << "-icosahedron: The radius and divisions must be specified!" << endl;
			else if ( ico_divisions < 1 ) ico_divisions = 1;
		}
		if ( curropt->tag == "circle" )
			if ( curropt->values(circle_radius, circle_offset) < 1 )
				cerr << "-circle: The circle radius must be specified!" << endl;
		if ( curropt->tag == "ellipse" )
			if ( curropt->values(ellipse_semi_z[0], ellipse_semi_z[1],
					ellipse_semi_z[2]) < 2 )
				cerr << "-ellipse: The major and minor semiaxes must be specified!" << endl;
		if ( curropt->tag == "ellipsoid" )
			if ( curropt->values(ellipsoid_semi[0], ellipsoid_semi[1],
					ellipsoid_semi[2]) < 3 )
				cerr << "-ellipsoid: The 3 semiaxes must be specified!" << endl;
		if ( curropt->tag == "cylinder" )
			if ( curropt->values(cyl_rad, cyl_len) < 2 )
				cerr << "-cylinder: Both radius and length must be specified!" << endl;
		if ( curropt->tag == "spindle" )
			if ( curropt->values(spindle_rad, packing) < 2 )
				cerr << "-spindle: Both radius and packing separation must be specified!" << endl;
		if ( curropt->tag == "lattice" )
			if ( curropt->values(lattice[0], lattice[1], lattice[2], lattice_type) < 3 )
				cerr << "-lattice: All three dimensions must be specified!" << endl;
		if ( curropt->tag == "randomfill" )
			if ( ( random_sphere_radius = curropt->value.real() ) < 1 )
				cerr << "-randomfill: A sphere radius must be specified!" << endl;
		if ( curropt->tag == "randomshell" )
			if ( ( random_shell_radius = curropt->value.real() ) < 1 )
				cerr << "-randomshell: A shell radius must be specified!" << endl;
		if ( curropt->tag == "gaussian" )
			if ( ( gauss_stdev = curropt->value.real() ) < 1 )
				cerr << "-gaussian: A standard deviation must be specified!" << endl;
		if ( curropt->tag == "all" ) all = 1;
		if ( curropt->tag == "replace" ) replace = 1;
		if ( curropt->tag == "reduce" )
			reduce = curropt->filename();
		if ( curropt->tag == "id" ) {
			model_id = curropt->value;
			if ( model_id.length() < 1 )
				cerr << "-id: An identifier must be specified!" << endl;
		}
		if ( curropt->tag == "settype" )
			set_type = curropt->value;
		if ( curropt->tag == "changetype" )
			change_type = curropt->value;
		if ( curropt->tag == "associate" ) {
			astr = curropt->value;
			associate_type = astr.pre(',');
			associate_file = astr.post(',');
			astr = 0;
		}
		if ( curropt->tag == "peaks" )
			peakmap = curropt->filename();
		if ( curropt->tag == "add" )
			newloc.push_back(curropt->vector3());
		if ( curropt->tag == "translate" ) {
			shift = curropt->vector3();
        	if ( shift.length() < 0.1 )
				cerr << "-translate: Three values must be specified!" << endl;
		}
		if ( curropt->tag == "trim" )
			trim = curropt->size();
		if ( curropt->tag == "linklength" )
			if ( ( linklength = curropt->value.real() ) < 0.1 )
				cerr << "-linklength: A link length must be specified!" << endl;
		if ( curropt->tag == "rdf" )
			if ( ( rdf_interval = curropt->value.real() ) < 0.001 )
				cerr << "-rdf: A sampling interval must be specified!" << endl;
		if ( curropt->tag == "consolidate" )
			if ( ( consolidate = curropt->real() ) < 1 )
				cerr << "-consolidate: A cutoff distance must be specified!" << endl;
		if ( curropt->tag == "radius" )
			if ( ( model_radius = curropt->value.real() ) < 1 )
				cerr << "-radius: A value must be specified!" << endl;
		if ( curropt->tag == "helix" ) {
			if ( curropt->values(helix_rise, helix_angle, helix_comp) < 3 )
				cerr << "-helix: The rise, angle and number must be specified!" << endl;
			else
				helix_angle *= M_PI/180.0;
		}
		if ( curropt->tag == "direction" )
			dir = curropt->vector3();
		if ( curropt->tag == "componentradius" )
			if ( ( compradius = curropt->value.real() ) < 0.1 )
				cerr << "-componentradius: The component display radius must be specified!" << endl;
		if ( curropt->tag == "linkradius" )
			if ( ( linkradius = curropt->value.real() ) < 0.1 )
				cerr << "-linkradius: The link display radius must be specified!" << endl;
		if ( curropt->tag == "color" ) {
			if ( curropt->values(
					color[0], color[1], color[2], color[3]) < 3 ) {
				cerr << "-color: At least three color values must be specified!" << endl;
			} else {
				set_color = 1;
			}
		}
		if ( curropt->tag == "points" )
			if ( ( points = curropt->value.integer() ) < 1 )
				cerr << "-points: The number of points must be specified!" << endl;
		if ( curropt->tag == "density" )
			if ( ( density = curropt->value.real() ) < 1 )
				cerr << "-density: The number per unit volume be specified!" << endl;
		if ( curropt->tag == "separation" )
			if ( ( separation = curropt->value.real() ) < 0.1 )
				cerr << "-separation: The separation distance must be specified!" << endl;
		if ( curropt->tag == "parameters" )
			paramfile = curropt->filename();
		if ( curropt->tag == "map" )
			mapfile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
		if ( curropt->tag == "catenate" )
			catname = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();

	// Read all the parameter files
	Bmodel*			model = NULL;
	Bmodel*			mp = NULL;
	Bmodel*			newmod = NULL;
	Bsymmetry		sym;
	
	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( file_list.size() )
		model = read_model(file_list, paramfile.str());
	
	if ( create.length() ) newmod = new Bmodel(create.str());

	if ( peakmap.length() ) {
		Bimage*			map = read_img(peakmap, 1, 0);
		if ( !map ) bexit(-1);
		newmod = model_from_peaks(map, 0, 0);
		delete map;
	}

	if ( symmetry_string.length() ) {
		if ( symmetry_string[0] == 'H' ) {
			newmod = model_helix(model_radius, helix_rise, helix_angle, helix_comp);
		} else {
			sym = Bsymmetry(symmetry_string);
			newmod = model_platonic(sym, model_radius);
		}
//	} else if ( model_type ) {
//		if ( model_type < 6 ) newmod = model_platonic(model_type, model_radius);
//		else newmod = model_helix(model_radius, helix_rise, helix_angle, helix_comp);
	} else if ( box.volume() > 0 ) {
		if ( density > 0 ) points = density*box.volume();
		Vector3<double>	start;
		newmod = model_random(points, compradius, start, box);
	} else if ( random_sphere_radius > 0 ) {
		newmod = model_random(points, compradius, random_sphere_radius);
	} else if ( random_shell_radius > 0 ) {
		if ( separation ) newmod = model_random_shell(points, random_shell_radius, separation);
		else newmod = model_random_shell(points, random_shell_radius);
	} else if ( gauss_stdev > 0 ) {
		newmod = model_random_gaussian(points, gauss_stdev);
	} else if ( fibo_radius > 0 ) {
		newmod = model_create_fibonacci_sphere(points, fibo_radius);
	} else if ( shell_radius >= 0 ) {
		if ( ico_divisions > 0 )
			newmod = model_create_icosahedron(shell_radius, ico_divisions, sphere_flag);
		else
			newmod = model_create_shell(points, shell_radius, separation);
	} else if ( circle_radius > 0) {
		newmod = model_create_circle(circle_radius, circle_offset, separation);
	} else if ( ellipse_semi_z[0] > 0 && ellipse_semi_z[1] > 0) {
		newmod = model_create_ellipse(ellipse_semi_z, separation);
	} else if ( ellipsoid_semi.volume() > 0) {
		newmod = model_create_ellipsoid(ellipsoid_semi, separation);
	} else if ( cyl_len > 0) {
		newmod = model_create_cylinder(dir, cyl_rad, cyl_len, separation);
	} else if ( spindle_rad > 0 && packing > 0) {
		newmod = model_create_spindle(dir, spindle_rad, separation, packing);
	} else if ( lattice.volume() ) {
		if ( lattice_type ) newmod = model_create_hexagonal_lattice(lattice, separation);
		else newmod = model_create_cubic_lattice(lattice, separation);
	}
	
	if ( model_id.length() ) newmod->identifier(model_id.str());
	else model_id = "create";
	
	if ( set_type.length() ) model_set_type(newmod, set_type.str());
	
	if ( associate_file.length() )
		model_associate(newmod, associate_type.str(), associate_file.str());
	
	if ( linklength > 0 ) model_link_list_generate(newmod, linklength);

	if ( mapfile.length() ) newmod->mapfile(mapfile.str());

	if ( replace && model && newmod ) {
		model_kill(model);
		model = newmod;
		newmod = NULL;
	} else {
		if ( model )
			for ( mp = model; mp->next; mp = mp->next ) ;	// Pointer to last model
		if ( mp ) mp->next = newmod;
		else model = newmod;
		mp = newmod;
	}
	
	if ( all ) models_process(model, model_reset_selection);
	
	if ( newloc.size() ) model_add_components(model, model_id.str(), set_type.str(), newloc);

	if ( change_type.length() ) model_change_type(model, change_type.str());
	
	if ( reduce.length() ) model_reduce_linked(model, reduce.str(), 1);

	if ( consolidate ) model_consolidate(model, consolidate);

//	if ( !model->poly ) model_poly_generate(model);
	
	if ( rdf_interval > 0 ) model_radial_distribution(model, rdf_interval);

	if ( compradius ) models_process(model, compradius, model_set_component_radius);

	if ( linkradius ) models_process(model, linkradius, model_set_link_radius);

	if ( set_color ) model_color_uniformly(model, color);

	if ( catname.length() ) model_catenate_maps(model, catname.str());
	
	if ( shift.length() ) models_shift(model, shift);

	if ( trim.volume() ) models_trim(model, trim);

	model_selection_stats(model);

	// Write an output parameter format file if a name is given
    if ( outfile.length() && model ) {
		write_model(outfile.str(), model);
	}

	model_kill(model);
	
	
		timer_report(ti);
	
	bexit(0);
}

