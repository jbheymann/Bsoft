/**
@file	bcat.cpp
@brief	Program to catenate image files
@author Bernard Heymann
@date	Created: 19991221
@date	Modified: 20230524
**/

#include "rwimg.h"
#include "img_combine.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen


// Usage assistance
const char* use[] = {
" ",
"Usage: bcat [options] input.images",
"----------------------------------",
"Catenates a list of images into a multi-image file.",
"Any number of input images may be given and may include the wild card \"*\".",
"(However, VMS does not support Unix style usage of wild cards).",
"All images must have the same data type.",
" ",
"Actions:",
"-rescale -0.1,5.2        Rescale input images to average and standard deviation.",
" ",
"Parameters:",
"-verbose 3               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-reverse                 Reverse the order of catenation.",
"-slices                  Pack 2D images as z-slices (default pack as separate images).",
"-size 10,50,8            New image size (pixels).",
"-fill 125.5              Fill value for padding (default average).",
"-raw d=f#x=120,120,1     Formatting to reinterpret image files.",
" ",
"Output:",
"-output output.img       Output image (default temp.miff).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	// Initialize variables
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	string			catfile("temp.pif");
	string			rawstring;
	int 			setZslices(0);             	// Pack slices as separate 2D images
	Vector3<long> 	nusize;
	bool			reverse_list(0);
	double			nuavg(0), nustd(0);
	int 			fill_type(FILL_AVERAGE);	// Fill type for smaller images
	double			fill(0);
	
	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "datatype" )
			nudatatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "origin" ) {
			if ( curropt->value[0] == 'c' ) {
				set_origin = 2;
			} else {
				origin = curropt->origin();
				set_origin = 1;
			}
		}
		if ( curropt->tag == "reverse" )
        	reverse_list = 1;
		if ( curropt->tag == "slices" )
        	setZslices = 1;
		if ( curropt->tag == "size" )
			nusize = curropt->size();
		if ( curropt->tag == "rescale" ) {
			if ( curropt->values(nuavg, nustd) < 2 )
				cerr << "-rescale: Both averge and standard deviation must be specified!" << endl;
			if ( nustd < 0 ) nustd = 0;
		}
		if ( curropt->tag == "fill" )
			fill = curropt->fill(fill_type);
		if ( curropt->tag == "raw" ) {
			rawstring = curropt->value.str();
			if ( rawstring.length() < 1 )
				cerr << "-raw: A valid format reinterpretation string must be specified!" << endl;
		}
		if ( curropt->tag == "output" ) {
			catfile = curropt->value.str();
			if ( catfile.length() < 1 )
				cerr << "-output: An output file must be specified!" << endl;
		}
    }
	option_kill(option);

	if ( verbose && argc < 3 ) {
		cerr << "Error: No input files given!" << endl;
		bexit(-1);
	}
	
	double		ti = timer_start();

	// Read all the file names
	vector<string>	file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);
	if ( !file_list.size() ) {
		cerr << "Error: No image files specified!" << endl;
		bexit(-1);
	}

	if ( reverse_list ) reverse(file_list.begin(), file_list.end());
	
	Bimage*			pcat = img_catenate(file_list, rawstring, nudatatype, 
						nusize, setZslices, fill_type, fill, nuavg, nustd);

	if ( sam.volume() > 0 ) pcat->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) pcat->origin(pcat->default_origin());
		else pcat->origin(origin);
	}
	
	pcat->calculate_background();
	write_img(catfile, pcat, 0);
	
	delete pcat;
	
	
		timer_report(ti);
	
	bexit(0);
}

