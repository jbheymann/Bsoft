/**
@file	scatter.cpp
@brief	Functions for calculating electron scattering profiles
@author Bernard Heymann
@date	Created: 20190521
@date	Modified: 20230507
**/

#include "scatter.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/*
	The intensity is composed of the direct and scattering components.
	The average for the scattering component (mu) is modified
	with aperture-based fraction k for the elastic part
	and the use of an energy filter for the inelastic part.
*/
double			signal_integrated(double emu, double imu,
		double ke, double ki, int ef)
{
	double			sig(1), f(1), mu(ke*emu), pmu(1);
	if ( !ef ) mu += ki*imu;
	
	for ( long n=1; n<=10; ++n ) {
		f *= n;			// Factorial
		pmu *= mu;		// Powers of average
		sig += pmu/f;
//		cout << n << tab << f << tab << pmu << tab << sig << endl;
	}

	sig *= exp(-emu-imu);
	
	return sig;
}

double			signal_integrated_new(double emu, double imu, double k, int ef)
{
	double			sig(0), f(1), pmu(1);
	
	// Elastic scattering
	for ( long n=1; n<=10; ++n ) {
		f *= n;			// Factorial
		pmu *= emu;		// Powers of average
		sig += pmu/f;
	}

	if ( ef )
		sig = (1 + k*sig)*exp(-emu-imu);
	else
		sig = 1 - (1 - k)*sig*exp(-emu-imu);
	
	return sig;
}

auto string_is_empty = []( const std::string &s )
{
	return s.size() < 1;
};

vector<string>	all_elements(map<string,Bcomptype>& types)
{
	long			mx(0);
	
	for ( auto ct: types )
		if ( mx < ct.second.index() ) mx = ct.second.index();
	
//	cout << "Maximum Z = " << mx << endl;

	vector<string>	el(mx);
	
	for ( auto ct: types )
		el[ct.second.index()-1] = ct.second.identifier();

	el.erase( remove_if( el.begin(), el.end(), string_is_empty ), el.end() );

	return el;
}

/**
@brief 	Calculates elastic electron scattering profiles for a subset of component types.
@param 	&types			component types.
@param 	ds 				frequency space sampling increment.
@param 	scut 			maximum spacial frequency.
@return map<string, vector<double>>	array with scattering profile.
**/
map<string, vector<double>>		calculate_elastic_scattering_curves(
				map<string,Bcomptype>& types, double ds, double scut)
{
	map<string, vector<double>>	scat;
	
	for ( auto ct: types ) {
		if ( verbose & VERB_FULL )
			cout << "Calculating scattering curve for " << ct.first << endl;
		vector<double>	scat1 = ct.second.elastic_scattering_curve(ds, scut);
		scat[ct.first] = scat1;
	}
		
	return scat;
}

/**
@brief 	Calculates the atomic potential for elastic electron scattering for a component type.
@param 	&ct				component type.
@param 	dr 				radius increment.
@param 	rmax 			maximum radius.
@return vector<double>	array with potential profile.
**/
vector<double>		calculate_potential_curve(Bcomptype& ct, double dr, double rmax)
{
	if ( dr < 1e-10 ) dr = 0.001;
	if ( rmax < 2 ) rmax = 2;
	
	long			i, j, k, imax(rmax/dr + 0.5);
	double			r2, dr2(dr*dr);
    
    // Calculate the atomic potential profile
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG calculate_potential_curve: z = " << ct.index() << endl;
		
	vector<double>	pot;
	vector<double>	coef(ct.coefficients());

	for ( j=0, k=5; j<5; ++j, ++k ) {
		coef[j] *= sqrt(M_PI/coef[k]);
		coef[k] = M_PI*M_PI/coef[k];
	}

	for ( i=0; i<imax; ++i ) {
		r2 = i*i*dr2;
		pot.push_back(coef[10]);
		for ( j=0, k=5; j<5; ++j, ++k )
			pot[i] += coef[j]*exp(-coef[k]*r2);
	}
		
	return pot;
}

/**
@brief 	Calculates potential profiles for a subset of component types.
@param 	&types			component type.
@param 	dr 				radius increment.
@param 	rmax 			maximum radius.
@return map<string, vector<double>>	array with potential profile.
**/
map<string, vector<double>>	calculate_potential_curves(map<string,Bcomptype>& types,
								double dr, double rmax)
{
	map<string, vector<double>>	pot;
	
	for ( auto ct: types ) {
		if ( verbose & VERB_FULL )
			cout << "Calculating potential curve for " << ct.first << endl;
		vector<double>	pot1 = calculate_potential_curve(ct.second, dr, rmax);
		pot[ct.first] = pot1;
	}
		
	return pot;
}

/**
@brief 	Integrates the elastic scattering curve up to the aperture cutoff.
@param 	&ct				component type.
@param	ctf				CTF and microscope parameters.
@return double			integral.
**/
double		elastic_cross_section_integrated(Bcomptype& ct, CTFparam& ctf)
{
	long			i(0);
	double			ds(0.001), f, cs(0);
	double			wlr = ctf.lambda()*lorentz(ctf.volt());
	double			scut(ctf.frequency_cutoff());
	double			scale(8*M_PI*ds*ds*wlr*wlr);

//	vector<double>	scat = calculate_elastic_scattering_curve(ct, ds, scut);
	vector<double>	scat = ct.elastic_scattering_curve(ds, scut);

	for ( auto it = scat.begin(); it != scat.end(); ++it, ++i ) {
		f = *it;
		cs += f*f*i;
//		cout << i*ds << tab << cs << endl;
	}
	
	cs *= scale;
	
	return cs;
}

/**
@brief 	Calculates the combined elastic cross section for a defined composition.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		elastic_cross_section(map<string,Bcomptype>& types, double volt)
{
	double			cs(0);

	for ( auto ct: types )
		cs += ct.second.elastic_cross_section(volt) * ct.second.component_count();

	return cs;
}


/**
@brief 	Calculates the combined elastic cross section for a defined composition.
@param 	&types			component types.
@param	ctf				CTF and microscope parameters.
@return double			cross section.
**/
double		elastic_cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double			cs(0);
	double			scut(ctf.frequency_cutoff());

	for ( auto ct: types )
		cs += ct.second.elastic_cross_section(ctf.volt(), 0, scut) * ct.second.component_count();

	return cs;
}

/**
@brief 	Calculates the combined excluded elastic cross section for a defined composition.
@param 	&types			component types.
@param	ctf				CTF and microscope parameters.
@return double			excluded cross section.
**/
double		elastic_excluded_cross_section(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double			cs(0);
	double			scut(ctf.frequency_cutoff());

	for ( auto ct: types )
		cs += ct.second.elastic_cross_section(ctf.volt(), scut, 1e10) * ct.second.component_count();

	return cs;
}

double		elastic_cross_section_lenz(long Z, double volts)
{
	double			f = 1e10*PLANCK/(EMASS*LIGHTSPEED);
	double			b2 = beta2(volts);

	double			cs = (f*f/(M_PI*b2))*pow(Z, 4.0/3.0);
	
	return cs;
}

double		elastic_cross_section_langmore(long Z, double volts)
{
	double		b2 = beta2(volts), b = sqrt(b2);

	double		cs = 1.4e-4*(pow(Z, 1.5)/b2)*(1-0.26*Z/(137*b));
	
	return cs;
}


/**
@brief 	Calculates the full inelastic cross section for an component type.
@param 	Z				atomic number.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		inelastic_cross_section_langmore(long Z, double volt)
{
	double		b2 = beta2(volt);

	double		cs = 1.5e-4*
		(sqrt(Z)/b2)*log(2*b2*
		(ECHARGE*volt+EMASS*LIGHTSPEED*LIGHTSPEED)/(20*ECHARGE));
	
	if ( Z == 1 ) cs *= 0.154;	// Correction for hydrogen
	
	return cs;
}

/**
@brief 	Calculates the effective inelastic cross section for an component type.
@param 	Z				atomic number.
@param	&ctf			microscope parameters.
@return double			integral.
**/
double		inelastic_cross_section_langmore(long Z, CTFparam& ctf)
{
	if ( ctf.slit_width() ) return 0;

	double		cs = inelastic_cross_section_langmore(Z, ctf.volt());

	double		scut(ctf.frequency_cutoff());
	
	cs  *= scut/(scut + 0.054);
	
	return cs;
}

/**
@brief 	Calculates the excluded inelastic cross section for an component type.
@param 	Z				atomic number.
@param	&ctf			microscope parameters.
@return double			integral.
**/
double		inelastic_excluded_cross_section_langmore(long Z, CTFparam& ctf)
{
//	if ( ctf.slit_width() ) return 0;

	double		cs = inelastic_cross_section_langmore(Z, ctf.volt());

	if ( ctf.slit_width() ) return cs;

	double		scut(ctf.frequency_cutoff());
	
	cs  *= 0.054/(scut + 0.054);
	
	return cs;
}

/**
@brief 	Calculates the combined full inelastic cross section for a defined composition.
@param 	&types			component types.
@param	volt			acceleration voltage.
@return double			integral.
**/
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, double volt)
{
	double		cs(0);

	for ( auto ct: types )
		cs += ct.second.inelastic_cross_section_langmore(volt) * ct.second.component_count();
	
	return cs;
}

/**
@brief 	Calculates the combined effective inelastic cross section for a defined composition.
@param 	&types			component types.
@param	&ctf			microscope parameters.
@return double			integral.
**/
double		inelastic_cross_section_langmore(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double		cs(0);
	double		scut(ctf.frequency_cutoff());

	for ( auto ct: types )
		cs += ct.second.inelastic_cross_section_langmore(ctf.volt(), 0, scut) * ct.second.component_count();
	
	return cs;
}

/**
@brief 	Calculates the combined excluded inelastic cross section for a defined composition.
@param 	&types			component types.
@param	&ctf			microscope parameters.
@return double			integral.
**/
double		inelastic_excluded_cross_section_langmore(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double		cs(0);
	double		scut(ctf.frequency_cutoff());

	if ( ctf.slit_width() ) scut = 0;

	for ( auto ct: types )
		cs += ct.second.inelastic_cross_section_langmore(ctf.volt(), scut, 1e10) * ct.second.component_count();
	
	return cs;
}

/**
@brief 	Calculates the combined total cross section for a defined composition.
@param 	&types			component types.
@param	&ctf			CTF and microscope parameters.
@return double			integral.

	If the slit width is specified, the energy filter is assumed to be used
	and only the elastic cross section is returned.
**/
double		cross_section_integrated(map<string,Bcomptype>& types, CTFparam& ctf)
{
	double			ecs = elastic_cross_section_integrated(types, ctf);
	double			ics = inelastic_cross_section_langmore(types, ctf);

	if ( verbose ) {
		cout << "Elastic cross section:          " << ecs << " A2" << endl;
		cout << "Inelastic cross section:        " << ics << " A2" << endl;
	}
	
	return ecs + ics;
}

/**
@brief 	Calculates the half-maximal frequency for a component type.
@param 	&ct				component type.
@return double			half-maximal frequency.
**/
double		cross_section_half_maximal_frequency(Bcomptype& ct)
{
	long			i(0);
	double			ds(0.001), f, csi(0);
	double			cs = ct.elastic_scattering_integral(0,1e10);
//	vector<double>	scat = calculate_elastic_scattering_curve(ct, ds, 1);
	vector<double>	scat = ct.elastic_scattering_curve(ds, 1);
	
	if ( verbose & VERB_FULL )
		cout << "Full scattering integral = " << cs << endl;
	
	for ( auto it = scat.begin(); it != scat.end() && csi < 0.5*cs; ++it, ++i ) {
		f = *it * ds;
		csi += 2*f*f*i;
//		cout << i*ds << tab << *it << tab << f << tab << csi << endl;
	}
	
	return ds*i;
}

/**
@brief 	Calculates the half-maximal frequency for an component type.
@param 	&types			component types.
@return double			half-maximal frequency.
**/
double		cross_section_half_maximal_frequency(map<string,Bcomptype>& types)
{
	double			c, n(0), k(0);

	for ( auto ct: types ) {
		c = ct.second.component_count();
		k += c * cross_section_half_maximal_frequency(ct.second);
//		cout << ct.second.identifier() << tab << ct.second.component_count() << tab << k << endl;
		n += c;
	}

	if ( n ) k /= n;
	else k = 1;
	
	return k;
}

/**
@brief 	Calculates and shows the half-maximal frequencies for a material.
@param 	&material		material with an elemental composition.
@return int				0.

	Calculates the half-maximal cross section to estimate the mid-cutoff
	frequency coefficient for each element to model the aperture effect.
**/
int			material_cross_section_half_maximal_frequencies(Bmaterial& material)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	double			cs, sm;

	cout << "Material:                       " << material.identifier() << endl;
	cout << "Elastic cross section half-maximal frequencies:" << endl;
	cout << "Element\tZ\tIntensity\tFrequency(/A)" << endl;

	for ( auto el1: comp ) {
		Bcomptype&		ct = el1.second;
		cs = ct.elastic_scattering_integral(0,1e10);
		sm = cross_section_half_maximal_frequency(ct);
		cout << setprecision(5) << ct.identifier() << tab << ct.index() << tab
			<< cs << tab << sm << endl;
	}
	
	return 0;
}

/**
@brief 	Calculates the elastic scattering cross sections for a material.
@param 	&material		material with an elemental composition.
@param 	&ctf			CTF parameters.
@return double			scattering cross section (A^2).

**/
double		material_elastic_cross_section(Bmaterial& material, CTFparam& ctf)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	double			csel(0);
	
	for ( auto el: comp ) {
		Bcomptype&	ct = el.second;
		csel += ct.component_count() * elastic_cross_section_integrated(ct, ctf);
	}
	
	return csel;
}

/**
@brief 	Calculates the scattering cross sections for a material.
@param 	&material		material with an elemental composition.
@param 	&ctf			CTF parameters.
@return double			scattering cross section (A^2).

**/
double		material_cross_section(Bmaterial& material, CTFparam& ctf)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	double			csel(0), csin(0);
	
	for ( auto el: comp ) {
		Bcomptype&	ct = el.second;
		csel += ct.component_count() * elastic_cross_section_integrated(ct, ctf);
		if ( ctf.slit_width() < 0.1 )
			csin += ct.component_count() * inelastic_cross_section_langmore(ct.index(), ctf);
	}
	
	return csel + csin;
}

/**
@brief 	Calculates and shows the scattering cross sections for a material.
@param 	&material		material with an elemental composition.
@param 	&ctf			CTF parameters.
@return double			scattering cross section (A^2).

**/
double		material_show_cross_section(Bmaterial& material, CTFparam& ctf)
{
	double			csel, csin, csel_sum(0), csin_sum(0), cnt(0);
	double			lambda = electron_wavelength(ctf.volt());
	double			scut(ctf.frequency_cutoff());

	map<int,Bcomptype>		cn = material.composition_numbered();

//	if ( verbose ) {
		cout << "Cross sections:" << endl;
		cout << "Material:                       " << material.identifier() << endl;
		cout << "Density:                        " << material.dalton_per_angstrom3() << " Da/A3" << endl;
		cout << "Density:                        " << material.number_per_angstrom3() << " /A3" << endl;
		cout << "Acceleration voltage:           " << ctf.volt()*1e-3 << " kV" << endl;
		cout << "Wavelength:                     " << lambda << endl;
		cout << "Lorentz factor:                 " << lorentz(ctf.volt()) << endl;
		cout << "Beta:                           " << beta(ctf.volt()) << endl;
		cout << "Aperture diameter:              " << ctf.objective_aperture()*1e-4 << " µm" << endl;
		cout << "Focal length:                   " << ctf.focal_length()*1e-7 << " mm" << endl;
		cout << "Cutoff frequency:               " << scut << " /Å (" << 1.0/scut << " A)" << endl;
		if ( ctf.slit_width() )
			cout << "Slit width:                     " << ctf.slit_width() << " eV" << endl;
		cout << "Volume:                         " << material.volume() << " A3" << endl;
		cout << "Z\tElement\tCount\tMass(Da)\tcsel\tcsin\tcs" << endl;
//	}

	for ( auto el: cn ) {
		Bcomptype&	ct = el.second;
		csel = ct.component_count() * elastic_cross_section_integrated(ct, ctf);
		if ( ctf.slit_width() < 0.1 )
			csin = ct.component_count() * inelastic_cross_section_langmore(ct.index(), ctf);
		csel_sum += csel;
		csin_sum += csin;
		cnt += ct.component_count();
		cout << ct.index() << tab << ct.identifier() << tab << ct.component_count() 
			<< tab << setprecision(5) << ct.mass()*ct.component_count() 
			<< tab << setprecision(6) << csel << tab << csin << tab << csel+csin << endl;
//		cout << ct.identifier() << tab << ct.component_count() << tab
//				<< ct.mass()*ct.component_count() << tab << ct.elastic_cross_section(ctf.volt(),0,scut) << tab
//				<< ct.inelastic_cross_section_langmore(ctf.volt(),0,scut) << tab << csel+csin << endl;
	}
	
	cout << "Total" << tab << tab << setprecision(1) << cnt 
		<< tab << setprecision(5) << material.mass() 
		<< tab << setprecision(6) << csel_sum << tab << csin_sum << tab << csel_sum+csin_sum << endl << endl;
	
	return csel_sum + csin_sum;
}

/**
@brief 	Calculates the excluded scattering cross section for a material.
@param 	&material		material with an elemental composition.
@param 	&ctf			CTF parameters.
@return double			scattering cross section (A^2).

**/
double		material_excluded_cross_section(Bmaterial& material, CTFparam& ctf)
{
	map<string,Bcomptype>&	comp = material.composition();

	double		csel = elastic_excluded_cross_section(comp, ctf);
	double		csin = inelastic_excluded_cross_section_langmore(comp, ctf);
	
	return csel + csin;
}

/**
@brief 	Calculates the average effective mean free path.
@param 	&material		material with types.
@param	&ctf			microscope parameters.
@return double			EMFP average.
**/
double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf)
{
	double			thick_step(100), thick_max(5000), t(thick_step), emfp_avg(0);

	vector<double>	I = signal_intensity(material, thick_step, thick_max, ctf);

	for ( long i=0; i<I.size(); ++i, t+=thick_step )
		emfp_avg += -t/log(I[i]);
		
	emfp_avg /= I.size();

	return emfp_avg;
}


/*
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&material		material with types.
@param	&ctf			microscope parameters.
@return double			effective mean free path.
**/
/*double		effective_mean_free_path(Bmaterial& material, CTFparam& ctf)
{
	return effective_mean_free_path(material.composition(), ctf);
	map<string,Bcomptype>& comp = material.composition();
	
	double			d(material.number_per_angstrom3());
	bool			ef(ctf.slit_width()>0);	// Energy filter
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
//	double			ki = sc/(sc+0.054);
	double			ki = sc/(sc+0.2);
//	double			ki = sc*sc/(sc*sc+0.45*0.45);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	double			ics = inelastic_cross_section_langmore(comp, ctf.volt());
//	cout << d << tab << ecs << tab << ics << tab << ke << tab << ki << endl;
//	double			emfp = (ef)? 1/(d*ke*ecs): 1/(d*(ke*ecs+ki*ics));
	double			emfp = 1/(d*ke*ecs);
	if ( !ef ) emfp += 1/(d*ki*ics);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Material:                       " << material.identifier() << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
	}
	
	return emfp;
}*/

/**
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&materials		list of materials with types.
@param 	&fractions		fractional contributions.
@param	&ctf			microscope parameters.
@return double			effective mean free path.
**/
double		effective_mean_free_path(vector<Bmaterial>& materials, vector<double> fractions, CTFparam& ctf)
{
	Bmaterial		m = material_combine(materials, fractions);
	
	return effective_mean_free_path(m, ctf);
}

/**
@brief 	Calculates the expected intensity given a thickness of vitreous ice.
@param 	&material		material with types.
@param 	thickness		specimen thickness.
@param	ctf				microscope parameters.
@return double			intensity.
**/
double		signal_intensity(Bmaterial& material, double thickness, CTFparam& ctf)
{
	if ( thickness < 1 ) thickness = 1000;
	
	material.show();

	map<string,Bcomptype>& comp = material.composition();
	
	double			d(material.number_per_angstrom3());
	
	int				ef(ctf.slit_width()>0);	// Energy filter
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
	double			ki = sc/(sc+0.054);
//	double			ki = sc/(sc+0.2);
//	double			ki = sc*sc/(sc*sc+0.45*0.45);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	double			ics = inelastic_cross_section_langmore(comp, ctf);
	double			mfp = 1/(d*(ecs+ics));
	
	double			emu = ecs*d*thickness;
	double			imu = ics*d*thickness;
	double			i = signal_integrated(emu, imu, ke, ki, ef);
//	cout << i << endl;
	double			emfp = -thickness/log(i);
//	double			emfp = 1/(d*(ke*ecs+ki*ics));
	
	if ( verbose ) {
		cout << "Mean free path:                 " << mfp << " A" << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
		cout << "Half_maximal frequency:         " << sm << " 1/Å" << endl;
		cout << "Frequency cutoff:               " << sc << " 1/Å" << endl;
		cout << "Aperture effect:                " << ke << tab << ki << endl;
		cout << "Thickness:                      " << thickness << " A" << endl;
		cout << "Intensity:                      " << i << endl;
		cout << endl;
	}
	
	return i;
}

/**
@brief 	Calculates the expected intensity vs thickness of vitreous ice.
@param 	&material		material with types.
@param 	thick_step		specimen tickness step size (angstrom).
@param 	thick_max		maximum specimen thickness (angstrom).
@param	&ctf			microscope parameters.
@return vector<double>	array of intensities.
**/
vector<double>	signal_intensity(Bmaterial& material, double thick_step, double thick_max, CTFparam& ctf)
{
	if ( thick_step < 1 ) thick_step = 100;
	if ( thick_max < 1 ) thick_max = 5000;

	map<string,Bcomptype>& comp = material.composition();

	int				ef(ctf.slit_width()>0);	// Energy filter
	double			d(material.number_per_angstrom3());
	double			sm = cross_section_half_maximal_frequency(comp);
	double			sc = ctf.frequency_cutoff();
	double			ke = sc*sc/(sc*sc+sm*sm);
	double			ki = sc/(sc+0.054);

	double			ecs = elastic_cross_section(comp, ctf.volt());
	
	double			ics = inelastic_cross_section_langmore(comp, ctf);
	
	double			t, emu, imu, i, emfp;
	vector<double>	I;

//	cout << "Aperture effect:                " << k << endl;

	if ( verbose & VERB_FULL )
		cout << "t(A)\tI\tEMFP" << endl;
	for ( t=thick_step; t<=thick_max; t+=thick_step ) {
		emu = ecs*d*t;
		imu = ics*d*t;
		i = signal_integrated(emu, imu, ke, ki, ef);
		emfp = -t/log(i);
		I.push_back(i);
		if ( verbose & VERB_FULL )
			cout << t << tab << i << tab << emfp << endl;
	}
	
	return I;
}

/**
@brief 	Calculates and writes atomic scattering profiles to a file.
@param 	paramfile	file with scattering coefficients.
@param 	outfile		file to write curves to.
@param 	selection	element selection comma-delimited string.
@param 	resolution	resolution limit.
@return int			0.

	The scattering coefficients are obtained from an input parameter file.

**/
int			write_scattering_curves(string paramfile, string outfile, string selection, double resolution)
{
 	map<string,Bcomptype> 	types = read_atom_properties(paramfile);
 	
	Bmaterial				material = material_from_elements(selection);
	
	material.update_parameters(types);
	
	return write_scattering_curves(outfile, material.composition(), resolution);
}

/**
@brief 	Calculates and writes atomic scattering profiles to a file.
@param 	outfile		file to write curves to.
@param 	&types		element types.
@param 	resolution	resolution limit.
@return long			number of curves.

	The scattering coefficients are obtained from an input parameter file.

**/
long		write_scattering_curves(string outfile, map<string,Bcomptype>& types, double resolution)
{
	if ( resolution < 0.01 || resolution > 2 ) resolution = 2;
	
	long			i, j, nscat(100);
	double			scut(1.0/resolution), ds(scut/nscat);
    
	map<string, vector<double>>	scat = calculate_elastic_scattering_curves(types, ds, scut);

	ofstream		fd(outfile.c_str());
	if ( fd.fail() ) return -1;

	if ( verbose )
		cout << "Writing scattering curves to:   " << outfile << endl;
	
	fd << "Scattering coefficients:" << endl;
	for ( auto ct: types )
		fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( j=0; j<5; j++ ) {
		fd << "a" << j+1;
		for ( auto ct: types )
			fd << tab << ct.second.coefficients()[j];
		fd << endl;
	}
	for ( j=0; j<5; j++ ) {
		fd << "b" << j+1;
		for ( auto ct: types )
			fd << tab << ct.second.coefficients()[j+5];
		fd << endl;
	}
	fd << "c";
	for ( auto ct: types )
		fd << tab << ct.second.coefficients()[10];
	fd << endl;
	
	fd << "Scattering curves:" << endl << "s";
	for ( auto ct: types )
		fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( i=0; i<nscat; i++ ) {
		fd << ds*i;
		for ( auto ct: types )
			fd << tab << scat[ct.first][i];
		fd << endl;
	}
	fd << endl;
	
	fd.close();
	
	return scat.size();
}

/**
@brief 	Calculates and writes atomic scattering profiles to a file.
@param 	outfile		file to write curves to.
@param 	&scat		scattering curves.
@param 	ds			spatial frequency increments (1/Å).
@return long			number of curves.

	The scattering coefficients are obtained from an input parameter file.

**/
long		write_scattering_curves(string outfile, map<string, vector<double>>& scat, double ds)
{
	long			i, nscat(scat.begin()->second.size());
    
	ofstream		fd(outfile.c_str());
	if ( fd.fail() ) return -1;

	if ( verbose )
		cout << "Writing scattering curves to:   " << outfile << endl;
		
	fd << "Scattering curves:" << endl << "s";
	for ( auto ct: scat )
		fd << tab << ct.first;
	fd << endl;
	for ( i=0; i<nscat; i++ ) {
		fd << ds*i;
		for ( auto ct: scat )
			fd << tab << ct.second[i];
		fd << endl;
	}
	fd << endl;
	
	fd.close();
	
	return scat.size();
}

/**
@brief 	Calculates and writes atomic potential profiles to a file.
@param 	paramfile		file with scattering coefficients.
@param 	outfile			file to write curves to.
@param 	selection		element selection.
@param 	radius			maximum radius.
@return int				0.

	The scattering coefficients are obtained from an input parameter file.

**/
int			write_potential_curves(string paramfile, string outfile, string selection, double radius)
{
 	map<string,Bcomptype> 	types = read_atom_properties(paramfile);
 	
	Bmaterial				material = material_from_elements(selection);
	
	material.update_parameters(types);
	
	return write_potential_curves(outfile, material.composition(), radius);
}

/**
@brief 	Calculates and writes atomic potential profiles to a file.
@param 	outfile			file to write curves to.
@param 	types			component types.
@param 	radius			maximum radius.
@return int				0.

	The scattering coefficients are obtained from an input parameter file.

**/
int			write_potential_curves(string outfile, map<string,Bcomptype>& types, double radius)
{
	if ( radius < 0.1 || radius > 10 ) radius = 10;
	
	long			i, j, nrad(100);
	double			dr(radius/nrad);
    
	map<string, vector<double>>	pot = calculate_potential_curves(types, dr, radius);

	ofstream		fd(outfile.c_str());
	if ( fd.fail() ) return -1;

	if ( verbose )
		cout << "Writing potential curves to:    " << outfile << endl;
	
	fd << "Scattering coefficients:" << endl;
	for ( auto ct: types )
		fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( j=0; j<5; j++ ) {
		fd << "a" << j+1;
		for ( auto ct: types )
			fd << tab << ct.second.coefficients()[j];
		fd << endl;
	}
	for ( j=0; j<5; j++ ) {
		fd << "b" << j+1;
		for ( auto ct: types )
			fd << tab << ct.second.coefficients()[j+5];
		fd << endl;
	}
	fd << "c";
	for ( auto ct: types )
		fd << tab << ct.second.coefficients()[10];
	fd << endl;
	
	fd << "Potential curves:" << endl << "s";
	for ( auto ct: types )
		fd << tab << ct.second.identifier() << "(" << ct.second.index() << ")";
	fd << endl;
	for ( i=0; i<nrad; i++ ) {
		fd << dr*i;
		for ( auto ct: types )
			fd << tab << pot[ct.first][i];
		fd << endl;
	}
	fd << endl;
	
	fd.close();
	
	return 0;
}

