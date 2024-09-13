/**
@file	bmgrec.cpp
@brief	Micrograph reconstration from oriented particle projections
@author	Bernard Heymann
@date	Created: 20240406
@date	Modified: 20240408
**/

#include "mg_reconstruct.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bmgrec [options] input.star [input2.star ...]",
"----------------------------------------------------",
"Calculates projections of a reference map based on particle orientations and.",
"reintegrate them into the framework of the original micrograph.",
" ",
//"Actions:",
//"-removemarkers 14        Mask out markers with this radius (pixels) before transformation.",
//" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
//"-datatype u              Force writing of a new data type.",
"-resolution 20           Resolution limit (angstrom, default Nyquist).",
"-Volt 200k               Acceleration voltage (default 0).",
" ",
"Input:",
"-reference map.grd       Input 3D map for generating projections.",
" ",
"Output:",
"-output file.star        Output STAR file name.",
" ",
NULL
};

int 		main(int argc, char **argv)
{
//	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	double 			resolution(0); 				// Must be set > 0 to limit resolution
	double			volt(0);					// Acceleration voltage
	Bstring			reffile;					// Input reference map
	Bstring			outfile;					// Output STAR file

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "resolution" )
			if ( ( resolution = curropt->value.real() ) < 0.01 )
				cerr << "-resolution: A resolution limit must be specified!" << endl;
		if ( curropt->tag == "Volt" )
			if ( ( volt = curropt->real_units() ) < 1 )
				cerr << "-Volt: The acceleration voltage must be specified!" << endl;
		if ( curropt->tag == "reference" )
			reffile = curropt->filename();
		if ( curropt->tag == "output" )
			outfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();

	// Read all the parameter files
	Bstring*			file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);
	string_kill(file_list);

	if ( project == NULL )  {
		cerr << "Error: No input file read!" << endl;
		bexit(-1);
	}
	
	Bimage*			pref = NULL;
	if ( reffile.length() )
		pref = read_img(reffile, 1, 0);
	
	if ( pref ) {
		project_micrograph_reconstruct(project, pref, resolution, volt);
		delete pref;
	}
	
	if ( project && outfile.length() ) {
		write_project(outfile, project, 0, 0);
	}
	
	project_kill(project);

	timer_report(ti);

	bexit(0);
}

