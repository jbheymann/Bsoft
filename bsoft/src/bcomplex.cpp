/**
@file	bcomplex.cpp
@brief	Program for handing complex images.
@author Bernard Heymann
@date	Created: 19990321
@date	Modified: 20240905
**/

#include "rwimg.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Usage assistance
const char* use[] = {
" ",
"Usage: bcomplex [options] input.img output.img",
"----------------------------------------------",
"Converts complex image forms.",
" ",
"Actions:",
"-size 120,102,200        Resize the image as a Fourier transform.",
"-tocomplex               Convert a phase image to a complex image (phases in radians).",
"-from2                   Convert two sub-images to a complex image.",
"-conjugate               Convert a complex image to its conjugate.",
"-invert                  Invert a complex image.",
"-logarithm               Calculate the logarithm of the image.",
"-convert real            Convert the complex transform: real, imag, amp, int, phase (default not).",
//"-signed                  Convert a complex image to amplitudes, with sign based on phase.",
"-positive                Convert a complex image to be postive definite.",
"-square                  Convert a complex image to its square.",
"-addcomplex 2.3,-4.1     Add a constant complex value to a complex image.",
"-addphase 2.5            Add a constant to the phase (can be pi or pi/2).",
"-friedel                 Check and apply Friedel symmetry.",
"-ewald opp               Convert an Ewald projection to the opposite or combined.",
"-shift -3,4.5,-0.3       Shift phases (pixels, can be half).",
"-center                  Shift image to center origin.",
"-correlate other.sup     Calculate correlation.",
"-color 1.5               Generate a phase-coloured power spectrum: amplitude scaling.",
"-partcolor 0.3,-0.4,0    Generate an image with part intensities and part phase-coloured.",
"-sideband -0.2,0.8,0.1   Retain a single side band image with division by plane normal.",
"-plot 2478               Polar plot up to a maximum amplitude.",
"-phasegrating 100k       Convert the image to a phase grating for the given voltage.",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-datatype u              Force writing of a new data type.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
" ",
NULL
};

int 		main(int argc, char **argv)
{
	DataType 		nudatatype(Unknown_Type);	// Conversion to new type
	Vector3<long>	nusize;						// New size if > 0
	Vector3<double>	sam;    					// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	bool			tocomplex(0);				// Flag to convert a phase image to complex
	bool			from2(0);					// Flag to convert a phase image to complex
	bool			conjugate(0);				// Flag to calculate conjugate
	bool			invert(0);					// Flag to invert
	bool 			setlogarithm(0);			// Flag for logarithmic power spectrum
    ComplexConversion	conv(NoConversion);		// Conversion from complex transform
	bool			positive(0);				// Flag to convert to positive definite
	bool			square(0);			  		// Flag to convert to its square
	Complex<double>	addcomplex;					// Constant complex to add
	double			addphase(0);				// Phase to add
	bool			friedel(0);					// Flag for applying Fiedel symmetry
	int				ewald(0);					// Flag for converting an Ewald projection
	bool			half(0);					// Flag to shift half of the image size
	Vector3<double>	shift;						// Shift vector
	bool			center(0);					// Center image
	double	 		colour_phase_scale(0);		// Scale for coloured phases
	bool			int_part_phase(0);			// Flag to calculate part phase colored
	Vector3<double> plane_normal;				// Plane normal for part phase image
	Vector3<double> sideband_normal;			// Plane normal for side band image
	double			plot_amp(0);				// Polar plot maximum amplitude
	double			phase_grating_volt(0);		// Voltage for phase grating
	Bstring			file_for_corr;				// Complex image to calculate difference
	
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
		if ( curropt->tag == "size" )
			nusize = curropt->size();
		if ( curropt->tag == "tocomplex" )
        	tocomplex = 1;
		if ( curropt->tag == "from2" )
        	from2 = 1;
		if ( curropt->tag == "conjugate" )
        	conjugate = 1;
		if ( curropt->tag == "invert" )
        	invert = 1;
		if ( curropt->tag == "logarithm" )
			setlogarithm = 1;
		if ( curropt->tag == "convert" )
			conv = curropt->complex_conversion();
		if ( curropt->tag == "positive" )
        	positive = 1;
		if ( curropt->tag == "square" )
        	square = 1;
 		if ( curropt->tag == "addcomplex" )
			if ( curropt->values(addcomplex[0], addcomplex[1]) < 2 )
				cerr << "-addcomplex: A complex value must be specified" << endl;
		if ( curropt->tag == "addphase" ) {
			if ( curropt->value[0] == 'p' ) {
				addphase = M_PI;
				if ( curropt->value.contains("/2") ) addphase /= 2 ;
			} else addphase = curropt->value.real();
		}
		if ( curropt->tag == "friedel" )
        	friedel = 1;
		if ( curropt->tag == "ewald" ) {
			if ( curropt->value[0] == 'o' ) ewald = 1;
			else if ( curropt->value[0] == 'c' ) ewald = 2;
			else
				cerr << "-ewald: other or combined must be specified" << endl;
		}
		if ( curropt->tag == "shift" ) {
			if ( curropt->value[0] == 'h' ) half = 1;
			else shift = curropt->vector3();
		}
 		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "color" )
			if ( ( colour_phase_scale = curropt->value.real() ) < 0.001 )
				cerr << "-color: A scale for the phase colours must be specified" << endl;
		if ( curropt->tag == "partcolor" ) {
		    plane_normal = curropt->vector3();
		    int_part_phase = 1;
		}
		if ( curropt->tag == "sideband" )
		    sideband_normal = curropt->vector3();
		if ( curropt->tag == "plot" )
			if ( ( plot_amp = curropt->value.real() ) < 0.001 )
				cerr << "-plot: A maximum amplitude must be specified" << endl;
		if ( curropt->tag == "phasegrating" )
			if ( ( phase_grating_volt = curropt->real_units() ) < 0.001 )
				cerr << "-phasegrating: An acceleration voltage must be specified" << endl;
		if ( curropt->tag == "correlate" )
			file_for_corr = curropt->filename();
    }
	option_kill(option);
	
 	double		ti = timer_start();
	
	int 		dataflag(friedel);
	if ( optind < argc - 1 ) dataflag = 1;
	Bimage*		p = read_img(argv[optind++], dataflag, -1);
	if ( !dataflag ) bexit(0);
	if ( p == NULL ) bexit(-1);
	
	if ( nudatatype == Unknown_Type )
		nudatatype = p->data_type();
	else if ( nudatatype > p->data_type() )
		p->change_type(nudatatype);
	
	if ( sam.volume() > 0 ) p->sampling(sam);
	
	if ( set_origin ) {
		if ( set_origin == 2 ) p->origin(p->default_origin());
		else p->origin(origin);
	}
	
	if ( tocomplex ) p->phase_to_complex();
	
	if ( from2 ) p->two_to_complex();
	
	if ( nusize.volume() > 0 ) p->change_transform_size(nusize);
	
	if ( conjugate ) p->complex_conjugate();
	
	if ( invert ) p->complex_invert();
	
	if ( half ) shift = p->size()/2;
	if ( shift.length() ) {
		p->phase_shift(shift);
	}
	
	if ( center )
		p->center_wrap();
	
	if ( friedel ) {
		p->friedel_check();
		p->friedel_difference();
//		p->friedel_apply();
	}
	
	if ( ewald == 1 ) {
		if ( p->sizeZ() < 2 ) p->opposite_ewald();
		else p->complex_geometric_invert();
	} else if ( ewald == 2 ) {
		if ( p->sizeZ() < 2 ) p->combine_ewald();
		else {
			Bimage*	pc = p->copy();
			pc->complex_geometric_invert();
			p->add(pc);
			delete pc;
		}
	}
	
	if ( phase_grating_volt )
		p->phase_grating(phase_grating_volt);
	
	if ( file_for_corr.length() ) {
		Bimage*		pref = read_img(file_for_corr, 1, -1);
		p->complex_conjugate_product(pref, 1);
		delete pref;
	}
	
	if ( positive )
		p->fspace_positive();
	
	if ( addcomplex.real() || addcomplex.imag() )
		p->complex_add(addcomplex);
	
	if ( addphase )
		p->complex_add_phase(addphase);
	
	if ( square )
		p->complex_square();
	
	if ( sideband_normal.length() )
		p->side_band(sideband_normal);
		
	if ( conv && p->compound_type() == TComplex ) {
		p->complex_convert(conv);
		p->statistics();
	} else if ( colour_phase_scale > 0 ) {
		Bimage*		pnu = p->intensities_phase_colored(colour_phase_scale, !half);
		if ( pnu ) {
			delete p;
			p = pnu;
			nudatatype = p->data_type();
		}
		Bstring		filename("ps_color_wheel.png");
		pnu = new Bimage(UCharacter, TRGB, 512, 512, 1, 1);
		pnu->phase_colour_wheel();
		write_img(filename, pnu, 0);
		delete pnu;
	} else if ( int_part_phase ) {
		Bimage*		pnu = p->intensities_with_part_phase_colored(plane_normal);
			delete p;
			p = pnu;
			nudatatype = p->data_type();
	}
	
	if ( setlogarithm ) p->logarithm();
	
	if ( plot_amp ) {
		Bimage*		ppol = p->polar_plot(plot_amp);
		delete p;
		p = ppol;
	}
	
    // Write an output file if a file name is given
    if ( optind < argc ) {
		p->change_type(nudatatype);
    	write_img(argv[optind], p, 0);
	}

	delete p;
	
	timer_report(ti);
	
	bexit(0);
}
