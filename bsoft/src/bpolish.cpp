/**
@file	bpolish.cpp
@brief	Program to align and sum particle frames from micograph frames
@author	Bernard Heymann
@date	Created: 20040407
@date	Modified: 20240215
**/

#include "mg_processing.h"
#include "mg_align.h"
#include "mg_extract.h"
//#include "mg_ctf.h"
#include "rwmg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bpolish [options] input.star [input.star]",
"------------------------------------------------",
"Aligns particle frames extracted from micrograph frames.",
" ",
"Actions:",
"-align 15,3,2            Align to reference image (first=0), moving sum and interval (default 1,1).",
//"-snr 5                   Estimate SSNR over a summation window.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
//"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-mgpath dir/subdir       Set micrograph file paths.",
"-framepath dir/subdir    Set micrograph frames file paths.",
"-partpath dir/subdir     Set the particle file paths.",
" ",
"Parameters for alignment:",
"-counts                  Flag to rescale images based on their counts histogram.",
"-resolution 900,300      High and low resolution limits for cross-correlation (default 0.1,1000 angstrom).",
"-shiftlimit 3.5          Limit on origin shift relative to nominal center (default 10% of box edge size).",
"-subset 2-8,12           Subset of micrographs to align and average.",
" ",
"Input:",
"-Gainreference gr.tif    Gain reference to correct the input micrographs.",
" ",
"Output:",
"-output file.star        Output parameter file, if -average, also output average parameter file.",
" ",
NULL
};


int			main(int argc, char** argv)
{
	// Initializing variables
//	DataType 		datatype(Unknown_Type);		// Conversion to new type
	int				ref_img(-1);				// Reference image for alignment, <0 means don't align
	int				window(1), step(1);			// Moving sum window for alignment
	int				flags(0);					// Flags: 1=rescale based on histogram; 2=weigh by dose; 4=write aligned frames; 8=write frame sum; 16=local
//	long			shift_window(1);			// Frame window for shift envelope calculation
//	double			shift_res(10);				// Resolution for shift envelope calculation
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	double			hi_res(0), lo_res(1e10);	// Default resolution range
	double			shift_limit(-1);			// Maximum shift from nominal image origin
	Vector3<double>	origin;						// Tilt axis origin
	double			edge_width(0), gauss_width(0);	// Edge parameters
	long		 	bin(1);						// Binning before alignment and analysis
	Bstring			mgpath;						// Micrograph file path
	Bstring			framepath;					// Micrograph frames file path
	Bstring			partpath;					// Particle file path
	Bstring			subset;						// Subset of micrographs to average
	Bstring			paramfile;					// Output parameter file
    Bstring			grfile;						// Gain reference
    Bstring			maskfile;					// Mask to use for cross-correlation

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;
	for ( curropt = option; curropt; curropt = curropt->next ) {
//		if ( curropt->tag == "datatype" )
//			datatype = curropt->datatype();
		if ( curropt->tag == "sampling" )
			sam = curropt->scale();
		if ( curropt->tag == "align" )
			if ( curropt->values(ref_img, window, step) < 1 )
 				cerr << "-align: A reference image number must be specified." << endl;
//		if ( curropt->tag == "snr" ) {
//			if ( ( snr_window = curropt->value.integer() ) < 2 )
// 				cerr << "-snr: A window larger than 1 must be specified." << endl;
//			if ( snr_prog ) snr_prog = 2;
//		}
		if ( curropt->tag == "counts" )
			flags |= 1;
		if ( curropt->tag == "resolution" ) {
    	    if ( curropt->values(hi_res, lo_res) < 1 )
				cerr << "-resolution: Resolution limits must be specified." << endl;
			else if ( hi_res > lo_res )
				swap(hi_res, lo_res);
        }
		if ( curropt->tag == "shiftlimit" )
			if ( ( shift_limit = curropt->value.real() ) < 1 )
				cerr << "-shiftlimit: A maximum shift in pixels must be specified!" << endl;
		if ( curropt->tag == "origin" )
			origin = curropt->origin();
		if ( curropt->tag == "edge" )
    	    if ( curropt->values(edge_width, gauss_width) < 1 )
				cerr << "-edge: An edge width must be specified." << endl;
		if ( curropt->tag == "bin" )
			if ( ( bin = curropt->value.integer() ) < 1 )
				cerr << "-bin: An ineteger greater than zero must be specified!" << endl;
		if ( curropt->tag == "subset" )
			subset = curropt->value;
		if ( curropt->tag == "mgpath" ) {
			mgpath = curropt->value;
			if ( mgpath.length() < 1 )
				cerr << "-mgpath: The micrograph file path must be specified!" << endl;
			else
				if ( mgpath[-1] != '/' ) mgpath += "/";
		}
		if ( curropt->tag == "framepath" ) {
			framepath = curropt->value;
			if ( framepath.length() < 1 )
				cerr << "-framepath: The micrograph frames file path must be specified!" << endl;
			else
				if ( framepath[-1] != '/' ) framepath += "/";
		}
		if ( curropt->tag == "partpath" ) {
			partpath = curropt->value;
			if ( partpath.length() < 1 )
				cerr << "-partpath: The particle file path must be specified!" << endl;
			else
				if ( partpath[-1] != '/' ) partpath += "/";
		}
		if ( curropt->tag == "Gainreference" )
			grfile = curropt->filename();
		if ( curropt->tag == "Mask" )
			maskfile = curropt->filename();
		if ( curropt->tag == "output" )
			paramfile = curropt->filename();
    }
	option_kill(option);

	double		ti = timer_start();

	// Read all the parameter files
	Bstring*		file_list = NULL;
	while ( optind < argc ) string_add(&file_list, argv[optind++]);
	if ( !file_list ) {
		cerr << "Error: No parameter or image files specified!" << endl;
		bexit(-1);
	}

	Bproject*		project = read_project(file_list);		
	string_kill(file_list);

	if ( !project ) {
		cerr << "Error: Input file not read!" << endl;
		bexit(-1);
	}
	
	if ( sam[0] > 0 )
		project_set_mg_pixel_size(project, sam);

	if ( mgpath.length() )
		project_set_micrograph_path(project, mgpath);
	if ( framepath.length() )
		project_set_frame_path(project, framepath);
	if ( partpath.length() )
		project_set_particle_path(project, partpath);

	Bimage*			pgr = NULL;
	if ( grfile.length() )
		pgr = read_img(grfile, 1, 0);
	
	Bimage*			pmask = NULL;
	if ( maskfile.length() )
		pmask = read_img(maskfile, 1, 0);
		
	if ( ref_img > -1 )
		project_align_particle_frames(project, ref_img, window, step, pgr, pmask, origin, hi_res, lo_res,
				shift_limit, edge_width, gauss_width, bin, subset, flags);
	
//	if ( snr_window > 1 )
//		project_frames_snr(project, hi_res, snr_window, subset, sampling_ratio, (flags&1));
	
    if ( paramfile.length() )
		write_project(paramfile, project, 0, 0);
	
	project_kill(project);
	if ( pgr ) delete pgr;
	delete pmask;

	timer_report(ti);

	bexit(0);
}

