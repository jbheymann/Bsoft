/**
@file	bfocal.cpp
@brief	Processing focal series
@author Bernard Heymann
@date	Created: 20220808
@date	Modified: 20240506
**/

#include "rwimg.h"
#include "mg_ctf_fit.h"
#include "mg_ctf_focal.h"
#include "utilities.h"
#include "options.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

Bimage*		img_focal_aberration_phase_difference(Bimage* p, Bimage* pref, CTFparam cp, vector<double>& dfocus);
Bimage*		img_ctf_focal_series_tilted(CTFparam& cp, vector<double>& dfocus,
				Vector3<long> size, Vector3<double> sam, double lores, double hires,
				double beam_tilt, double beam_dir);

/* Usage assistance */
const char* use[] = {
" ",
"Usage: bfocal [options] img.mrc out.pif",
"---------------------------------------",
"Processes micrographs taken as a focal series.",
" ",
"Actions (mutually exclusive):",
"-series 20               Calculate a focal series function for the number of images.",
"-fit 1000,20             Fit a focal series CTF: maximum iterations and truncation maximum (default 20).",
//"-extractsphere           Extract sphere from 3D frequency space.",
//"-weighsphere             Weigh a 3D transform with a focal coherent sphere.",
"-apply                   Apply the CTF to the input image.",
"-reconstruct             Reconstruct focal series.",
"-exitwave 100,1e-5,1     Reconstruct an exit wave from a focal series; maximum iterations, tolerance, flag.",
" ",
"Actions:",
"-zft                     Back transform z columns.",
"-center                  Center the simulated image.",
"-convert real            Convert the complex transform: real, imag, amp, int, phase (default not).",
" ",
"Parameters:",
"-verbose 7               Verbosity of output.",
"-sampling 1.5,1.5,1.5    Sampling (A/pixel; a single value can be given).",
"-origin 0.8,-10,15.7     Set the origin (default from input image).",
"-size 200,300,50         Size of focal series (must be specified for focal series).",
"-resolution 2.3,120      Resolution limits for fitting (A).",
"-focusstep 37,2,1        Focus increment between images, step to refine, flag for linear search.",
"-ctf 3                   CTF type: 0=complex, 1=sine, 2=envelope.",
"-Bfactor 44              B-factor (default 0 A^2).",
"-beamtilt 5mr,45d        Beam tilt and direction (radians, mr=millirad, d=degrees).",
"-Defocus 1.2,1.0,47      Defocus average & deviation, and astigmatism angle (default 2 um, 0, 0).",
"-Astigmatism 0.3,-34     Set defocus deviation and astigmatism angle.",
" ",
#include "use_ctf.inc"
" ",
"Input:",
"-jsin file.json          Input JSON file with CTF parameters.",
" ",
"Output:",
"-jsout file.json         Output CTF parameters to a JSON file.",
//"-envelope env.mrc        The Ewald sphere envelope map.",
" ",
NULL
};

void		ctf_apply_beam_tilt(CTFparam& cp, double tilt, double dir)
{
	cp.show_aberration();
	
	if ( verbose ) {
		cout << "Beam tilt angle:                " << tilt << " rad (" << tilt*180.0/M_PI << " degrees)" << endl;
		cout << "Beam direction angle:           " << dir << " rad (" << dir*180.0/M_PI << " degrees)" << endl;
	}
	
	long		i;
	double		tamp(sin(tilt)), v;
	double		tx(tamp*cos(dir)), ty(tamp*sin(dir));
	
	double		a = cp.aberration_coefficient(0,0);
//	cout << "lambda=" << cp.lambda() << endl;
//	cout << "amp=" << a << endl;
	for ( i=1, v=tamp; i<5; ++i, v*=tamp ) a += v/i;
//	cout << "amp=" << a << endl;
	cp.aberration_coefficient(0,0,a);
	
	double		cs = cp.aberration_coefficient(4,0);
	
	cp.beam_tilt(tx, ty);
	
	double		d = cp.aberration_coefficient(2,0);
	d += 2*cs*tamp*tamp;
	cp.aberration_coefficient(2,0,d);
	
	double		cx = cp.aberration_coefficient(2,-2);
	double		cy = cp.aberration_coefficient(2,2);
	cp.aberration_coefficient(2, -2, cy + 2*cs*tx*ty);
	cp.aberration_coefficient(2, 2, cx + cs*(tx*tx-ty*ty));
	
}

int 		main(int argc, char **argv)
{
    /* Initialize variables */
    bool			series(0);					// Flag to calculate a focal series CTF
 	long			fit(0);						// Flag to fit the CTF and number of iterations
 	double			tmax(20);					// Truncation maximum for power spectrum preparation
	bool			center(0);					// Flag to center image
//	bool			extract_sphere(0);			// Flag to extract a sphere based on the voltage
//	bool			weigh_sphere(0);			// Flag to weigh a 3D transform with a sphere based on the voltage
	bool			apply(0);					// Flag to apply the CTF
	bool			zft(0);						// Flag to do z transform
    ComplexConversion	conv(NoConversion);		// Conversion from complex transform
	Vector3<double>	sam;						// Units for the three axes (A/pixel)
	Vector3<double>	origin;						// New image origin
	int				set_origin(0);				// Flag to set origin
	Vector3<long>	size;						// Size for new focal series image
	double			hires(0), lores(0);			// High and low resolution limits
	int				ctf_flag(0);				// Flag for CTF form: 0=complex, 1=envelope, 2=sine
	bool			rec_flag(0);				// Flag to reconstruct
	long			exit_iter(0);				// Number of iterations for exit wave reconstruction
	double			tolerance(1e-6);			// Tolerance for exit wave reconstruction
	int				exit_flag(0);				// Flag for exit wave reconstruction: 0=intensity, 1=sqrt(int), 2=real
	double			Bfactor(0); 				// B-factor
	double			beam_tilt(0), beam_dir(0);	// Beam tilt and direction
	double			def_avg(0);					// In angstrom
	double			def_dev(0);					// In angstrom
	double			ast_angle(0);	 			// Used to limit astigmatism
/*	int				basetype(1);				// Baseline type: 1=poly, 2=double_gauss, 3=EMAN
	int				setbase(0);
	vector<double>	base = {1,0,0,0,0};			// Baseline coefficients
	int				envtype(1);					// Envelope type: 1=gauss, 2=gauss+, 3=double_gauss, 4=double_gauss+
	int				setenv(0);
	vector<double>	env = {1,-1,0,0,0};			// Envelope coefficients
	*/
 	double			focus_step(0), focus_inc(0);	// Focus step and search increment (angstrom)
 	double			focus_linear(0);			// Linear focus search extent (±angstrom)
 	long			nfoc(1);					// Number of images in focal series
	Bstring			jsin, jsout;				// JSON files

	double			v;
	JSvalue			jsctf(JSobject);

	int				optind;
	Boption*		option = get_option_list(use, argc, argv, optind);
	Boption*		curropt;

	// If the JSON file with CTF parameters is specified, read it first
	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "jsin" ) {
			jsin = curropt->filename();
			if ( jsin.length() ) jsctf = JSparser(jsin.c_str()).parse();
			if ( !jsctf.exists("Volt") ) {
				cerr << "CTF parameters not found in file " << jsin << endl;
				bexit(-1);
			}
		}
	}

	for ( curropt = option; curropt; curropt = curropt->next ) {
		if ( curropt->tag == "series" ) {
			if ( ( nfoc = curropt->integer() ) < 1 )
				cerr << "-series: A number of images must be specified!" << endl;
			else
				series = 1;
		}
		if ( curropt->tag == "fit" )
			if ( curropt->values(fit, tmax) < 1 )
				cerr << "-fit: A number of iterations must be specified!" << endl;
//		if ( curropt->tag == "extractsphere" ) extract_sphere = 1;
//		if ( curropt->tag == "weighsphere" ) weigh_sphere = 1;
		if ( curropt->tag == "apply" ) apply = 1;
		if ( curropt->tag == "zft" ) zft = 4;
		if ( curropt->tag == "center" ) center = 1;
		if ( curropt->tag == "convert" )
			conv = curropt->complex_conversion();
		if ( curropt->tag == "reconstruct" ) rec_flag = 1;
		if ( curropt->tag == "exitwave" )
			if ( curropt->values(exit_iter, tolerance, exit_flag) < 1 )
				cerr << "-exitwave: A number of iterations must be specified!" << endl;
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
			size = curropt->size();
 		if ( curropt->tag == "resolution" )
			if ( curropt->values(hires, lores) < 1 )
				cerr << "-resolution: A high resolution limit must be specified" << endl;
		if ( curropt->tag == "focusstep" )
			if ( curropt->values(focus_step, focus_inc, focus_linear) < 1 )
				cerr << "-focusstep: A focus step size and increment must be specified!" << endl;
		if ( curropt->tag == "ctf" )
			if ( ( ctf_flag = curropt->integer() ) < 0 )
				cerr << "-ctf: A CTF flag must be specified!" << endl;
		if ( curropt->tag == "Bfactor" )
			if ( ( Bfactor = curropt->real() ) < 0.1 )
				cerr << "-Bfactor: A B-factor must be specified!" << endl;
		if ( curropt->tag == "beamtilt" ) {
			if ( curropt->real_units(beam_tilt, beam_dir) < 2 )
				cerr << "-beamtilt: A beam tilt angle and direction angle must be specified!" << endl;
			else {
				if ( fabs(beam_tilt) > M_PI ) beam_tilt *= M_PI/180;	// Assume degrees
				if ( fabs(beam_tilt) > M_PI ) beam_dir *= M_PI/180;		// Assume degrees
			}
		}
		if ( curropt->tag == "Defocus" ) {
			if ( curropt->real_units(def_avg, def_dev, ast_angle) < 1 )
				cerr << "-Defocus: At least the defocus average must be specified!" << endl;
			else {
				ast_angle *= M_PI/180;				// Assume degrees
			}
		}
		if ( curropt->tag == "Astigmatism" ) {
			if ( curropt->real_units(def_dev, ast_angle) < 1 )
				cerr << "-Astigmatism: A defocus value must be specified!" << endl;
			else {
				ast_angle *= M_PI/180.0;						// Assume degrees
			}
		}
/*		if ( curropt->tag == "basetype" ) {
			basetype = curropt->value.integer();
			if ( basetype < 1 || basetype > 6 ) {
				basetype = 1;
				cerr << "Warning: The baseline type must be 1, 2 or 3. Reset to 1." << endl;
			} else
				setbase = 1;
		}
		if ( curropt->tag == "baseline" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			for ( size_t i=0; i<d.size(); i++ ) base[i] = d[i];
			if ( d.size() < 1 )
				cerr << "-baseline: At least one coefficient must be specified!" << endl;
			else
				setbase = 2;
		}
		if ( curropt->tag == "envtype" ) {
			envtype = curropt->value.integer();
			if ( envtype < 1 || envtype > 4 ) {
				envtype = 4;
				cerr << "Warning: The envelope type must be 1, 2, 3 or 4. Reset to 4." << endl;
			} else
				setenv = 1;
		}
		if ( curropt->tag == "envelope" ) {
			vector<double>	d = curropt->value.split_into_doubles(",");
			for ( size_t i=0; i<d.size(); i++ ) env[i] = d[i];
			if ( d.size() < 1 )
				cerr << "-envelope: At least an envelope amplitude must be specified!" << endl;
			else
				setenv = 2;
		}*/
#include "ctf.inc"
		if ( curropt->tag == "jsout" )
			jsout = curropt->filename();
    }
	option_kill(option);
	
	double			ti = timer_start();
	
//	jsctf.write("t.json");

    vector<CTFparam>	cpa;			// Array of CTF parameters
	CTFparam			cp;
	if ( jsctf.type() == JSarray ) {
		cpa = ctf_from_json_array(jsctf);
		cp = cpa[0];
	} else {
		cp = ctf_from_json(jsctf);
		cpa.push_back(cp);
	}

//	cp.beam_tilt(sin(beam_tilt)*cos(beam_dir), sin(beam_tilt)*sin(beam_dir));
	if ( def_avg ) cp.defocus_average(def_avg);
	if ( def_dev ) cp.astigmatism(def_dev, ast_angle);
	if ( focus_step ) cp.focus_step(focus_step);
	else focus_step = cp.focus_step();
//	if ( beam_tilt ) ctf_apply_beam_tilt(cp, beam_tilt, beam_dir);
//	if ( setbase ) cp.baseline(basetype, base);
//	if ( setenv ) cp.envelope(envtype, env);

//	if ( def_ser_inc ) cpa = ctf_series(cp, def_ser_start, def_ser_end, def_ser_inc);
	
//	if ( verbose & VERB_FULL )
		cp.show();
		
    Bimage* 		p = NULL;
    
	if ( size.volume() < 1 && optind < argc )
		p = read_img(argv[optind++], 1, -1);
    
	if ( !p && size.volume() < 1 ) {
		cerr << "Error: No input file read or size specified!" << endl;
		bexit(-1);
	}

	Bstring			outfile;
	if ( optind < argc ) outfile = argv[optind++];

	if ( p ) {
		if ( sam.volume() ) p->sampling(sam);
		else sam = p->image->sampling();
		size = p->size();
		if ( set_origin ) {
			if ( set_origin == 2 ) p->origin(p->default_origin());
			else p->origin(origin);
		}
		if ( nfoc != p->images()*p->sizeZ() ) {
			nfoc = p->images()*p->sizeZ();
			cerr << "Warning: Setting the number of focus values to " << nfoc << endl;
		}
		if ( p->images()*p->sizeZ() > 1 && fabs(focus_step) < 1 ) {
			focus_step = p->sampling(0)[2];
			cerr << "Warning: Setting the focus increment to " << focus_step << endl;
		}
	}
	
	if ( sam.volume() < 1e-3 ) sam = Vector3<double>(1,1,1);
	
	vector<double>		dfocus;
	for ( long i=0; i<nfoc; ++i )
		dfocus.push_back(focus_step*(i-nfoc/2));
//		dfocus.push_back(focus_step*i);

	if ( fit ) {		// Fitting power spectra
		Bimage*	pfit = img_ctf_focal_fit(p, cp, hires, lores, tmax, Bfactor, fit);
    	if ( outfile.length() )
	    	write_img(outfile, pfit, 0);
		bexit(0);
	}
	
	if ( focus_inc ) {		// Estimating the focus step size
		img_ctf_refine_focus_step(p, cp.lambda(), focus_step, focus_inc, hires, lores, focus_linear);
		cp.focus_step(focus_step);
		for ( long i=0; i<nfoc; ++i )
			dfocus[i] = focus_step*(i-nfoc/2);
	}

	Bimage*		pctf = NULL;
	if ( ( series || apply || rec_flag ) && size.volume() ) {
		if ( beam_tilt ) pctf = img_ctf_focal_series_tilted(cp, dfocus, size, sam, hires, lores, beam_tilt, beam_dir);
		else pctf = img_ctf_focal_series(cp, nfoc, size, sam, hires, lores, ctf_flag);
		if ( !p ) p = pctf;
	}

	if ( apply ) {
		if ( p->compound_type() != TComplex )
			p->fftxy();
		p->complex_product(pctf);
/*	} else if ( extract_sphere ) {
		if ( verbose )
			for ( auto v: dfocus ) cout << v << endl;
		Bimage*		ps = img_fspace_extract_sphere(p, cp, dfocus);
		delete p;
		p = ps;
	} else if ( weigh_sphere ) {
		img_fspace_weigh_sphere(p, cp.volt());
	} else if ( phi_diff ) {
		Bimage*		pphi = img_focal_aberration_phase_difference(p, ps, cp, dfocus);
		delete p;
		p = pphi;*/
	} else if ( rec_flag ) {
		if ( p->compound_type() != TComplex )
			p->fftxy();
		p->complex_conjugate_product(pctf);
		Bimage*		pr = p->fspace_sum(0);
		Bplot*		plot = pr->fspace_ssnr(nfoc, hires, 1);
		pr->multiply(1.0/nfoc);
		if ( outfile.length() ) {
			Bstring			psfile = outfile.base() + "_ssnr.ps";
			ps_plot(psfile, plot);
		}
		delete plot;
		delete p;
		p = pr;
	} else if ( exit_iter ) {
		Bimage*		pew = img_ctf_focal_exit_wave_reconstruct(p, cp, hires, exit_iter, tolerance, ctf_flag, exit_flag);
		delete p;
		p = pew;
	}

	if ( jsout.length() ) {
		jsctf = ctf_to_json(cp);
		jsctf.write(jsout.str());
	}

	if ( Bfactor != 0 ) p->fspace_weigh_B_factor(Bfactor, hires);

	if ( zft ) p->fftz(FFTW_BACKWARD, 1);

	if ( center ) p->center_wrap();

	p->complex_convert(conv);
	

    // Write an output file if a file name is given
    if ( outfile.length() ) {
//		p->change_type(nudatatype);
    	write_img(outfile, p, 0);
	}
	
	delete p;
	
	
		timer_report(ti);
	
	bexit(0);
}

Bimage*		img_focal_aberration_phase_difference(Bimage* p, Bimage* pref, CTFparam cp, vector<double>& dfocus)
{
	p->slices_to_images();
	
	if ( dfocus.size() != p->images() ) {
		cerr << "Error: The number of focal values (" << dfocus.size() <<
			") must equal the number of images (" << p->images() << ")" << endl;
		return NULL;
	}
	
	if ( p->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the images" << endl;
		p->fft();
	}

	if ( pref->fourier_type() != Standard ) {
		if ( verbose )
			cout << "Transforming the reference map" << endl;
		pref->fft();
	}

	long			nn;
	Bimage*			p1 = NULL;
	
	if ( verbose ) {
		cout << "Calculating phase difference for a focal series:" << endl;
		cout << "Focus step size:                " << dfocus[1] - dfocus[0] << " A" << endl;
	}
		
	Bimage*			psum = new Bimage(Float, 6, pref->sizeX(), pref->sizeY(), 1, 1);
	psum->sampling(pref->sampling(0));
	
	double			def(cp.defocus_average());
	
	for ( nn=0; nn<p->images(); ++nn ) {
		cp.defocus_average(def + dfocus[nn]);
		if ( verbose )
			cout << "Image " << nn << ": " << cp.defocus_average() << " A" << endl;
		p1 = p->extract(nn);
		img_add_aberration_terms(psum, pref, p1, cp);
		delete p1;
	}

	Bimage*		pphi = img_calculate_phase_differences(psum);

	delete psum;

	return pphi;
}

Bimage*		img_ctf_focal_series_tilted(CTFparam& cp, vector<double>& dfocus,
				Vector3<long> size, Vector3<double> sam, double lores, double hires,
				double beam_tilt, double beam_dir)
{
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;
	
	long			nimg(dfocus.size());
	if ( size[2] > 1 ) {
		nimg = 1;
		size[2] = dfocus.size();
	}
	
	double			shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	double			shi2(shi*shi), slo2(slo*slo);
	
	double			tx(sin(beam_tilt)*cos(beam_dir)/cp.lambda());
	double			ty(sin(beam_tilt)*sin(beam_dir)/cp.lambda());
//	double			t2(tx*tx+ty*ty);
	
	Bimage*			p = new Bimage(Float, TComplex, size, nimg);
	if ( sam.volume() > 0 ) p->sampling(sam);
	p->origin(0,0,size[2]/2);
	p->fourier_type(Standard);
	
	long 			i, x, y, z, n;
	double			sx, sy, s, s2, a, dphi;
	Complex<double>	cv(1,0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	
	double			def(cp.defocus_average());

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF focal series:" << endl;
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Defocus start, end, increment:  " << def+dfocus[0] << " - " << def+dfocus.back() << " ∆ " << dfocus[1] - dfocus[0] << endl;
		cout << "First sinc node:                " << sqrt(fabs(dfocus.back()-dfocus[0])*cp.lambda()/2) << " A" << endl;
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		cout << "Beam tilt vector:               " << tx << tab << ty << endl;
		cp.show_envelope();
		cout << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[1]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[1]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}

	for ( i=n=0; n<nimg; ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			cp.defocus_average(def+dfocus[n+z]);
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				sy += ty;
				for ( x=0; x<p->sizeX(); ++x, ++i ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					sx += tx;
					s2 = sx*sx + sy*sy;
					if ( s2 >= slo2 && s2 <= shi2 ) {
						s = sqrt(s2);
						a = atan2(sy,sx);
						dphi = cp.calculate_aberration_even(s2, a);
						cv = cp.aberration_odd_complex(s2, a);
						cv = cv.conj() * (sinl(dphi) * cp.calc_envelope(s));
						p->add(i, cv);
					}
				}
			}
		}
	}
	
//	p->multiply(1.0L/n);

	cp.defocus_average(def);

	return p;
}
