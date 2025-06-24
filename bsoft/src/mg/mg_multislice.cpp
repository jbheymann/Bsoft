/**
@file	mg_multislice.cpp
@brief	Multislice simulation 
@author	Bernard Heymann
@date	Created: 20030805
@date	Modified: 20241211
**/

#include "mg_processing.h"
#include "mg_ctf.h"
#include "molecule_to_map.h"
#include "mol_transform.h"
#include "mol_edit.h"
#include "Complex.h"
#include "Vector3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the wave propagation function between slices.
@param 	size		size of projection image (z = 1).
@param 	sam			pixel size in x and y, slice thickness in z.
@param 	thickness	slice thickness (in angstrom).
@param 	volt		acceleration voltage (volt).
@return Bimage*	 	wave propagation function image.
**/
Bimage*		img_calc_wave_propagator(Vector3<long> size, Vector3<double> sam,
				double thickness, double volt)
{
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->sampling(sam);
	
	long			i, h, k, x, y;
	double			wavelength = electron_wavelength(volt);
	double			arg, sx2, sy2;
	double			afactor = M_PI*wavelength*thickness;
//	double			nfactor = 1.0/(p->sizeX()*p->sizeY());	// Normalization for passing through FFT's
	double			xscale2 = 1.0/(p->sizeX()*p->sizeX()*sam[0]*sam[0]);
	double			yscale2 = 1.0/(p->sizeY()*p->sizeY()*sam[1]*sam[1]);
	
	for ( i=y=0; y<p->sizeY(); ++y ) {
		k = y;
		if ( k > (p->sizeY()-1)/2 ) k -= p->sizeY();
		sy2 = k*k*yscale2;
		for ( x=0; x<p->sizeX(); ++x, ++i ) {
			h = x;
			if ( h > (p->sizeX()-1)/2 ) h -= p->sizeX();
			sx2 = h*h*xscale2;
			arg = afactor*(sx2 + sy2);
//			p->set(i, Complex<double>(cos(arg), sin(arg)) * nfactor);
			p->set(i, Complex<double>(cos(arg), sin(arg)));
		}
	}
	
//	write_img("propagator.map", p);
	
	return p;
}

/**
@brief 	Simulates the electron imaging process using a multi-slice approach.
@param 	*pgrate			phase grating multi-image.
@param 	thickness		slice thickness (in angstrom).
@param 	volt	 		acceleration voltage (volt).
@param	norm			flag to normalize accross slices.
@return Bimage*	 		simulated projection image transform.

	The passage of the electron beam is simulated as the interaction of
	a planar wave with successive planar phase gratings spaced at regular 
	intervals, with the wave propagated between the 2D gratings. The
	phase gratings are derived from slabs of the atomic potential 
	calculated from the atomic structure using scattering profiles.
	Note: The final product is a 2D transform of the exit wave.
	Reference: Cowley, J. M. (1995) Diffraction Physics. 3rd Rev. Ed. 
		Elsevier Science, Amsterdam.

**/
Bimage*		img_calc_multi_slice(Bimage* pgrate, double thickness, double volt, bool norm)
{
	Vector3<long>	size(pgrate->size());
	
	// wavelength & size - loop one slice -> propagator (one 2D image)
	Bimage*			prop = img_calc_wave_propagator(size, pgrate->sampling(0), thickness, volt);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Multislice simulation:" << endl;
		cout << "Acceleration voltage:           " << volt << " V" << endl;
		cout << "Slice thickness:                " << thickness << " A" << endl;
//		cout << "Normalization:                  " << norm << endl << endl;
	}
	
	// For all slices: Multiply previous exitwave with phasegrating - FFT - 
	// multiply with propagator - FFT-1 - next slice -> exitwave (one 2D image)
	long   			n;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->label(pgrate->label());
	p->fourier_type(NoTransform);
	p->origin(pgrate->image->origin());
	p->sampling(pgrate->sampling(0));
	
	bool			flag(1);
	long			slice_size(size[0]*size[1]);
	Complex<double>	cv(1,0);
	for ( n=0; n<slice_size; ++n ) p->set(n, cv);
	
	Bimage* 		ptemp;
	Bimage*			pp = NULL;
	if ( flag ) pp = pgrate->copy();
	
	for ( n=0; n<pgrate->images(); ++n ) {
//		cout << n << endl;
		ptemp = pgrate->extract(n);			// Get phase grating slice
		p->complex_product(ptemp);			// Multiply with current wave
		delete ptemp;
//		p->fft(FFTW_FORWARD, 0);			// Normalization in propagator
		p->fft();
		if ( flag ) pp->replace(n, p);
		if ( n < pgrate->images()-1 ) {
			p->complex_product(prop);		// Multiply with propagator
//			p->fft(FFTW_BACKWARD, 0);		// Normalization in propagator
			p->fft(FFTW_BACKWARD, 1);
		}
	}
	
	delete prop;
	
//	if ( norm ) p->multiply(1.0/n);

	if ( flag ) {
		write_img("ms.grd", pp, 0);
		delete pp;
	}
	
	return p;
}

Bimage*		img_project_multi_slice(Bimage* pgrate, double thickness, double volt)
{
	Vector3<long>	size(pgrate->size());
	
	Bimage*			prop = img_calc_wave_propagator(size, pgrate->sampling(0), thickness, volt);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Projection multislice simulation:" << endl;
		cout << "Acceleration voltage:           " << volt << " V" << endl;
		cout << "Slice thickness:                " << thickness << " A" << endl;
	}
	
	bool			flag(1);
	long   			n;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, 1);
	p->label(pgrate->label());
	p->fourier_type(Standard);
	p->origin(pgrate->image->origin());
	p->sampling(pgrate->sampling(0));
	
	Bimage*			prop2 = prop->copy();
	Bimage* 		ptemp;
	Bimage*			pp = NULL;
	if ( flag ) pp = pgrate->copy();

	for ( n=0; n<pgrate->images(); ++n ) {
		ptemp = pgrate->extract(n);			// Get phase grating slice
		ptemp->fft();
		if ( n ) {
			ptemp->complex_product(prop2);	// Multiply with propagator
			prop2->multiply(prop); 
		}
		p->add(ptemp);						// Add to current wave
		if ( flag ) pp->replace(n, p);
		delete ptemp;
	}
	
	delete prop;
	delete prop2;
	
	p->multiply(1.0/n);

//	p->fft(FFTW_BACKWARD, 1);

	if ( flag ) {
		write_img("mp.grd", pp, 0);
		delete pp;
	}

	return p;
}

/**
@brief 	Calculates the atomic potential.
@param 	*molgroup	set of molecules.
@param 	size		size of projection image (z = 1).
@param 	origin		origin in x and y.
@param 	sam			voxel size.
@param 	thickness	slice thickness (in angstrom).
@param 	resolution	resolution limit (angstrom).
@param 	Bfactor		overall temperature factor.
@param 	&paramfile	parameter file for atomic scattering coefficients.
@param 	type		type of potential calculation: 0=reciprocal space, 1=real space, 2=gaussian
@return Bimage*	 	complex potential image.
**/
Bimage*		img_calc_potential(Bmolgroup* molgroup, Vector3<long> size, Vector3<double> origin,
				Vector3<double> sam, double thickness, double resolution, double Bfactor,
				Bstring& paramfile, int type)
{
	// Atoms & scat cross-sections - Loop slices - calc 2D reciprocal space structure factors
	long			i;
	int				nslices(0);
	Bmolgroup**		slice_molgroup = molgroup_split_into_slices(molgroup, thickness, nslices);
	
	if ( verbose )
		cout << "Calculating the atomic potential" << endl;
		
	int				spacegroup(0);
	UnitCell		unit_cell(sam[0]*size[0],sam[1]*size[1],sam[2]*size[2],M_PI_2,M_PI_2,M_PI_2);
	Bimage*			ptemp = NULL;
	Bimage*			p = new Bimage(Float, TComplex, size[0], size[1], 1, nslices);
	p->label(molgroup->comment.str());
	p->unit_cell(unit_cell);
	p->sampling(sam[0], sam[1], thickness);
	p->space_group(spacegroup);
	
	for ( i=0; i<nslices; i++ ) {
		if ( type ) {
			origin[2] = -slice_molgroup[i]->min[2]/sam[2] - 0.5;
			size[2] = (int) (thickness/sam[2] + 1);
			ptemp = img_from_molecule(slice_molgroup[i], origin, size, sam, 
					resolution, 0.001, 0, 2 - type, spacegroup, unit_cell);
			ptemp->project('z', 1);
			ptemp->simple_to_complex();
		} else {
//			origin[2] = -slice_molgroup[i]->min[2]/sam[2] - 0.5;
//			origin[2] = -i*thickness;
//			origin[2] = 0;
			size[2] = 1;
			sam[2] = thickness;
			ptemp = img_sf_from_molecule(slice_molgroup[i], origin, size, sam, 
					resolution, spacegroup, unit_cell, 2, Bfactor, paramfile);
			ptemp->fft(FFTW_BACKWARD, 0);	// No normalization - get actual potential at the beam
		}
		p->replace(i, ptemp);
		delete ptemp;
		p->origin(i, origin);
		if ( verbose )
			cout << "Slice " << i+1 << " done" << endl << endl;
	}

	p->fourier_type(NoTransform);
	
	for ( i=0; i<nslices; i++ )
		molgroup_kill(slice_molgroup[i]);
	
	delete[] slice_molgroup;
	
	return p;
}



/**
@brief 	Calculates the phase grating approximation from the atomic potential.
@param 	*p				atomic potential image (modified).
@param 	volt			acceleration voltage (volt).
@param	norm			flag to normalize accross slices.
@return int 				0.

	All calculations are complex.

**/
int			img_calc_phase_grating(Bimage* p, double volt, bool norm)
{
	p->simple_to_complex();
	
	long   			i, j, n;
	long   			imgsize = p->sizeX()*p->sizeY();
	double			arg;

//	double			wavelength = electron_wavelength(volt);
//	double			sigma = 0.0208886 * wavelength *
//						sqrt(1 + 5.8886579e-4/(wavelength*wavelength));	// According to Roar Kilaas
	double			sigma = 1e-10*TWOPI*ECHARGE/(PLANCK*LIGHTSPEED*beta(volt));
		
//	if ( verbose & VERB_FULL )
	if ( verbose )
		cout << "Sigma (interaction factor):     " << sigma << endl;

	if ( norm ) sigma /= p->images();

	for ( i=n=0; n<p->images(); n++ ) {
		for ( j=0; j<imgsize; i++, j++ ) {
			arg = -sigma*(p->complex(i)).real();
			p->set(i, Complex<double>(cos(arg), sin(arg)));
		}
	}
	
	return 0;
}

