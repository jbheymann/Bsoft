/**
@file	bscatter.cpp
@brief	Program to calculate scattering cross sections
@author	Bernard Heymann
@date	Created: 20190521
@date	Modified: 20230804
**/

#include "rwmodel.h"
#include "rwresprop.h"
#include "ctf.h"
#include "scatter.h"
#include "json.h"
#include "file_util.h"
#include "ps_plot.h"
#include "options.h"
#include "utilities.h"
#include "timer.h"


// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

// Constants
#define MAXSCAT 	1000			// Maximum number of atomic scattering data points

double		particle_signal(map<string,Bmaterial>& material, CTFparam& cp, Vector3<double> vol);
double		micrograph_recorded_intensity(Bmaterial& material, Bmaterial ice, CTFparam& cp, double area, double thickness);
int			micrograph_recorded_intensity(Bmaterial& material, Bmaterial& ice, CTFparam& cp, double area);
int			aperture_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> apser);
int			collection_angle_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> angser);
double		particle_snr(Bmaterial& material, double mass, double radius, double thickness, CTFparam& ctf);
Bplot*		particle_spectral_signal(Bmaterial& material, double mass,
				double radius, double thickness, CTFparam& ctf);

// Usage assistance
const char* use[] = {
" ",
"Usage: bscatter [options] input.pdb",
"-----------------------------------",
"Calculating properties from electron scattering cross sections.",
"Materials are specified either by input, from a materials ",
"property file, or by providing elements and their counts.",
" ",
"Actions:",
"-list                    List available materials from reference file.",
"-crosssections           List cross sections for a material.",
"-halfmaximum             Calculate the half-maximal frequency for each element.",
"-mfp                     Calculate the true mean free path.",
"-emfp                    Calculate the effective mean free path.",
"-series 30,50,70         Calculate EMFP using series of aperture sizes (um).",
"-angles 2.6,4.6,11.8     Calculate EMFP using series of collection angles (mrad).",
" ",
"Material parameters:",
"-material \"Vitreous ice\" Material name to read from properties file.",
"-elements H,C,O          Create material from elements.",
"-counts 10,3,2           Element counts.",
"-density 1.27g           Set density for all materials (g=g/cm3, d=Da/A3, n=#/A3).",
"-mass 234m               Set molecular weight for all materials (k=kDa, m=MDa).",
"-protein 65k             Add a protein with this mass (k=kDa, m=MDa).",
"-addhydrogens            Add hydrogens to a protein.",
"-ice                     Fill volume with ice.",
" ",
"Spatial parameters:",
"-area 120.5,45.7         Area of imaging (angstrom).",
"-thickness 1500          Specimen thickness (angstrom).",
"-radius 284              Particle radius (angstrom).",
" ",
"Parameters:",
"-verbose 1               Verbosity of output.",
"-resolution 1.5          High resolution for calculating curves (A).",
"-dose 25                 Electron dose (electrons/pixel).",
"-defocus 1.2,3.5,0.1     Defocus range an steps (µm).",
" ",
#include "use_ctf.inc"
" ",
"Input:",
"-coordinates prot.pdb    Input atomic coordinate model.",
"-atoms atomprop.star     Input atom properties file.",
"-residues resprop.star   Input residue properties file.",
"-properties matter.star  Input material density and composition file.",
" ",
"Output:",
"-output material.star    File with a list of materials.",
"-curves outfile.txt      File for scattering curve output (use with -elements option).",
"-rps powerspec.ps        Posrtscript file for radial power spectrum.",
" ",
NULL
};


int		main(int argc, char** argv)
{
	// Initialize variables
	bool			list_materials(0);	// Flag to display available materials
	bool			cross(0);			// Flag to list cross sections
	bool			halfmax(0);			// Flag to calculate half-maximal frequency
	bool			calc_mfp(0);		// Flag to calculate the mean free path
	bool			calc_emfp(0);		// Flag to calculate the effective mean free path
	double			add_protein(0);		// Mass for added protein
	bool			make_ice(0);		// Flag to fill volume with ice
	vector<double>	apser;				// Series of apertures
	vector<double>	angser;				// Series of collection angels
	string			material_str;		// Material string (name or coordinate file)
	string			elements;			// Elements
	string			count_str;			// String with element counts
	double			density(0);			// Material density
	DensityUnit		density_units(G_CM3);	// Density units (default g/cm3)
	Vector3<double>	volume;				// Imaging area and specimen thickness
	double			mass(0);			// Molecular weight
	bool			add_hydrogens(0);	// Flag to add hydrogens
	double			radius(0);			// Particle radius
	double			dose(0);			// Fluence
	double			def_min(0), def_max(0), def_step(0.1); // Defocus range
	double			hires(0);			// High resolution limit
	string			atom_select("all");	// Selection
	string			coorfile;			// Atomic coordinate file
	string 			atompropfile;		// Atom properties file
	string 			respropfile;		// Residue properties file
	string 			propfile("material.star");	// Material density and composition file
	string			paramfile;
	string			outfile;			// Output material file
	string			curvefile;			// Output file for scattering curves
	Bstring			rpsfile;			// Output file for scattering curves
	Bstring			jsin;				// JSON file

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
		if ( curropt->tag == "list" ) list_materials = 1;
		if ( curropt->tag == "crosssections" ) cross = 1;
		if ( curropt->tag == "halfmaximum" ) halfmax = 1;
		if ( curropt->tag == "mfp" ) calc_mfp = 1;
		if ( curropt->tag == "emfp" ) calc_emfp = 1;
		if ( curropt->tag == "series" )
			apser = curropt->value.split_into_doubles(",");
		if ( curropt->tag == "angles" )
			angser = curropt->value.split_into_doubles(",");
		if ( curropt->tag == "material" )
			material_str = curropt->value.str();
		if ( curropt->tag == "elements" )
			elements = curropt->value.str();
		if ( curropt->tag == "counts" )
			count_str = curropt->value.str();
		if ( curropt->tag == "density" )
			density = curropt->density(density_units);
/*		if ( curropt->tag == "density" ) {
			if ( ( density = curropt->value.real() ) < 0.01 ) {
				cerr << "-density: A density must be specified!" << endl;
			} else {
				if ( curropt->value.contains("g") ) density_units = G_CM3;
				if ( curropt->value.contains("d") ) density_units = DA_A3;
				if ( curropt->value.contains("n") ) density_units = NUM_A3;
			}
		}*/
		if ( curropt->tag == "mass" )
			mass = curropt->real_units();
		if ( curropt->tag == "protein" )
			add_protein = curropt->real_units();
		if ( curropt->tag == "addhydrogens" ) add_hydrogens = 1;
		if ( curropt->tag == "ice" ) make_ice = 1;
#include "ctf.inc"
		if ( curropt->tag == "area" )
			if ( curropt->values(volume[0], volume[1]) < 2 )
				cerr << "-area: An imaging area x and y coordinates must be specified!" << endl;
		if ( curropt->tag == "thickness" )
			if ( ( volume[2] = curropt->value.real() ) < 1 )
				cerr << "-thickness: A specimen thickness must be specified!" << endl;
		if ( curropt->tag == "radius" )
			if ( ( radius = curropt->value.real() ) < 1 )
				cerr << "-radius: A particle radius must be specified!" << endl;
		if ( curropt->tag == "dose" )
			if ( ( dose = curropt->value.real() ) < 1 )
				cerr << "-dose: A total electron dose per pixel must be specified!" << endl;
		if ( curropt->tag == "defocus" ) {
			if ( curropt->values(def_min, def_max, def_step) < 1 )
				cerr << "-defocus: At least one defocus value must be specified!" << endl;
			else {
				if ( def_min < 20 ) def_min *= 1e4;	// Assume µm
				if ( def_max < 20 ) def_max *= 1e4;	// Assume µm
				if ( def_max < def_min ) def_max = def_min;
				if ( def_step < 10 ) def_step *= 1e4;
			}
		}
		if ( curropt->tag == "resolution" )
			if ( ( hires = curropt->value.real() ) < 0.01 )
				cerr << "-resolution: A high resolution limit must be specified!" << endl;
		if ( curropt->tag == "coordinates" )
			coorfile = curropt->filename().str();
		if ( curropt->tag == "atoms" )
			atompropfile = curropt->filename().str();
		if ( curropt->tag == "residues" )
			respropfile = curropt->filename().str();
		if ( curropt->tag == "properties" )
			propfile = curropt->filename().str();
		if ( curropt->tag == "output" )
			outfile = curropt->filename().str();
		if ( curropt->tag == "curves" )
			curvefile = curropt->filename().str();
		if ( curropt->tag == "rps" )
			rpsfile = curropt->filename();
    }
	option_kill(option);
	
	double		ti = timer_start();
	
	map<string,Bmaterial>	mprop = read_material_properties(propfile);

	if ( list_materials ) {
		materials_list(mprop);
		bexit(0);
	}
	
	CTFparam			cp = ctf_from_json(jsctf);
	cp.defocus_average(def_min);
//	cp.defocus_deviation(def_dev);
//	cp.astigmatism_angle(ast_angle);

	if ( verbose & VERB_FULL )
		cp.show();

	vector<string>			file_list;
	while ( optind < argc ) file_list.push_back(argv[optind++]);

	map<string,Bmaterial> 	material = read_material_properties(file_list, paramfile, add_hydrogens);

	if ( density )
		for ( auto& m: material )
			m.second.density(density, density_units);

	map<string,Bcomptype>	atompar = read_atom_properties(atompropfile);
	
	if ( material_str.length() ) {
		if ( material_str[0] != '"' )
			material_str = "\"" + material_str + "\"";
		if ( verbose )
			cout << "Requested material: " << material_str << endl;
		if ( mprop.find(material_str) != mprop.end() ) {
			Bmaterial	m = mprop[material_str];
			if ( density ) m.density(density, density_units);
			if ( mass ) m.mass(mass);
			m.update_parameters(atompar);
			material[material_str] = m;
		}
	}
	
	if ( elements.length() ) {
		Bmaterial			m("Elements");
		map<string,Bcomptype>&	comp = m.composition();
		if ( elements == "all" ) {
			for ( auto ct: atompar ) {
				comp[ct.first] = ct.second;
				comp[ct.first].component_count(1);
			}
		} else {
			vector<string>		el = split(elements,',');
			vector<long>		cnt(el.size(),1);
			if ( count_str.length() )
				cnt = parse_integer_vector(count_str);
			for ( long i=0; i<el.size(); ++ i ) {
				comp[el[i]].identifier(el[i]);
				if ( i < cnt.size() ) comp[el[i]].component_count(cnt[i]);
				else comp[el[i]].component_count(1);
//				cout << comp[el[i]].identifier() << tab << comp[el[i]].component_count() << endl;
			}
			m.update_parameters(atompar);
		}
		if ( density ) m.density(density, density_units);
		if ( mass ) m.mass(mass);
		material[elements] = m;
	}

	if ( add_protein ) {
		Bmaterial			m = material_protein(1000);
		m.update_parameters(atompar);
		m.mass(add_protein);
		material["Protein"] = m;
	}

	if ( make_ice && volume.volume() ) {
		double		v(0);
		for ( auto& m: material ) v += m.second.volume();
		if ( v < volume.volume() ) {
			Bmaterial			m = material_ice(1);
			m.update_parameters(atompar);
			m.volume(volume.volume() - v);
			material["Vitreous ice"] = m;
		} else {
			cerr << "Error: The specimen volume is less than the combined material volume!" << endl;
			cerr << "Specimen volume:              " << volume.volume() << " A3" << endl;
			cerr << "Combined material volume:     " << v << " A3" << endl << endl;
		}
	}

	if ( verbose )
		for ( auto m: material )
			m.second.show();

	vector<double>	numbers;
	Bmaterial		mcomb = material_combine(material, numbers);
	
	if ( volume.volume() < 1 ) volume = mcomb.volume();
	
	if ( verbose )
		mcomb.show();

	if ( curvefile.length() )
		write_scattering_curves(curvefile, mcomb.composition(), hires);

	if ( halfmax )
		material_cross_section_half_maximal_frequencies(mcomb);
	
	if ( cross )
		material_show_cross_section(mcomb, cp);
	
	if ( calc_mfp )
		material_mean_free_path(material, cp.volt());

	if ( calc_emfp )
		material_effective_mean_free_path(material, cp.volt(), cp.frequency_cutoff(), cp.slit_width());
	
	if ( rpsfile.length() ) {
		Bplot* 	plot = material_radial_power_spectrum(material, cp.volt(), cp.frequency_cutoff(), cp.slit_width());
		ps_plot(rpsfile, plot);
		delete plot;
	}
	
	if ( material.find("Vitreous ice") != material.end() )
		particle_signal(material, cp, volume);

//		double		cs = material_show_cross_section(material, cp);
//		cout << "Signal:                         " << cs/area << endl << endl;
//		double		cs_ice = material_cross_section(ice, cp);
//		cout << "Signal(ice):                      " << cs_ice/area << endl << endl;
//	}
	
	// Calculate radius from mass, assuming a globular/spherical protein
//	if ( radius < 1 )
//		radius = pow((3*mass)/(4*M_PI*mcomb.density(DA_A3)), 1.0/3.0);
	
//	if ( thickness < 1 ) thickness = 2*radius;
	
//	if ( area < 1 ) area = 4*radius*radius;
	
//	if ( thickness < 1 ) thickness = 2*radius;
/*
	if ( verbose ) {
		cout << "Specimen dimensions:" << endl;
		cout << "Area:                           " << area << " A2" << endl;
		cout << "Thickness:                      " << thickness << " A" << endl;
		cout << "Volume:                         " << area*thickness << " A3" << endl;
		cout << "Mass:                           " << material.mass() << " Da" << endl;
		cout << "Density:                        " << material.dalton_per_angstrom3() << " Da/A3" << endl;
		cout << "Mean inner potential:           " << material.mean_inner_potential() << " V" << endl;
		cout << "Elastic cross section:          " << material.elastic_cross_section(cp.volt()) << " A2" << endl << endl;
	}
*/
//	Bmaterial	ice = material_ice(1);
//	ice.update_parameters(atompar);
	
//	if ( area && thickness ) ice.volume(area*thickness - mcomb.volume());
/*
	if ( cross ) {
		double		cs = material_show_cross_section(material, cp);
		cout << "Signal:                         " << cs/area << endl << endl;
//		double		cs_ice = material_cross_section(ice, cp);
//		cout << "Signal(ice):                      " << cs_ice/area << endl << endl;
	}
*/
//	material = material_combine(material, ice, 1, 1);
/*
	if ( thickness ) micrograph_recorded_intensity(material, ice, cp, area, thickness);
	else micrograph_recorded_intensity(material, ice, cp, area);
	

	if ( outfile )
		write_material_properties(outfile, material);
*/
/*	bexit(0);

	if ( verbose ) {
		cout << "General parameters:" << endl;
		cout << "Beta:                           " << beta(cp.volt()) << endl;
//		cout << "Frequency cutoff:               " << scut << " /A (" <<
//			1/scut << " A)" << endl;
		cout << "Specimen density:               " << material.density(DA_A3) << " Da/A3" << endl;
		cout << "Specimen thickness:             " << thickness << " A" << endl;
		cout << "Dose/fluence:                   " << dose << " e/px" << endl;
		if ( def_min ) {
			cout << "Defocus range:                  " << def_min*1e-4 << " - " << def_max*1e-4 << " µm" << endl;
			cout << "Defocus step:                   " << def_step*1e-4 << " µm" << endl;
		}
		cout << endl;
	}

	double		def_fac = defocus_factor(cp, def_min, def_max, def_step);
	double		intensity = signal_intensity(ice, thickness, cp);
	double		emfp = effective_mean_free_path(ice, cp);
	double		snre = particle_snr(material, mass, radius, thickness, cp);
	
	if ( dose ) {
		cout << "Total intensity:                " << dose << " e/px" << endl;
		cout << "Ice EMFP:                       " << emfp << endl;
		cout << "Direct beam attenuation:        " << intensity << endl;
		cout << "Defocus factor:                 " << def_fac << endl;
		cout << "Particle SNR:                   " << dose*def_fac*snre*intensity << endl << endl;
	}

	if ( apser.size() )
		aperture_series(thickness, ice, cp, apser);

	if ( angser.size() )
		collection_angle_series(thickness, ice, cp, angser);

//	Bplot*		plot = particle_spectral_signal(protcomp, mass,
//		radius, thickness, atompar, ctf);
*/

	if ( outfile.length() )
		write_material_properties(outfile, material);

	timer_report(ti);

	return 0;
}

/*
	The material is assumed to fill the volume
*/
double		material_scattering_average(map<string,Bmaterial>& material, CTFparam& cp, double area, double thickness)
{
	if ( area <  1) {
		cerr << "Error: The area must be specified!" << endl;
		bexit(-1);
	}
	
	double			scat_avg(0);
	double			scut(cp.frequency_cutoff());
	
	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
		scat_avg += m1.elastic_cross_section(cp.volt(), 0, scut);
		if ( cp.slit_width() < 0.01 )
			scat_avg += m1.inelastic_cross_section(cp.volt(), 0, scut);
	}

	return scat_avg/area;
}

double		micrograph_emfp(Bmaterial& material, Bmaterial ice, CTFparam& cp, double area, double thickness)
{
	if ( area <  1) {
		cerr << "Error: The area must be specified!" << endl;
		bexit(-1);
	}

	ice.volume(area*thickness - material.volume());
	
	Bmaterial	mat_ice = material_combine(material, ice, 1, 1);

	double		cs = material_excluded_cross_section(mat_ice, cp);
	
	double		emfp = area*thickness/cs;
	
	return emfp;
}

double		micrograph_signal(Bmaterial& material, Bmaterial ice, CTFparam& cp, double area, double thickness)
{
	if ( area <  1) {
		cerr << "Error: The area must be specified!" << endl;
		bexit(-1);
	}

	ice.volume(area*thickness - material.volume());
	
	Bmaterial	mat_ice = material_combine(material, ice, 1, 1);

	double		cs = material_cross_section(mat_ice, cp);

	double		csel = material_elastic_cross_section(mat_ice, cp);

	double		signal = (csel/area)*exp(-cs/area);
	
	return signal;
}

double		particle_signal(map<string,Bmaterial>& material, CTFparam& cp, Vector3<double> vol)
{
	if ( material.find("Vitreous ice") == material.end() ) {
		cerr << "Error: Vitreous ice parameters must be specified!" << endl;
		return 0;
	}
	
	double			cs(0), csf(0), cse, cses(0), mu, muf, mus, muss(0), area(vol[0]*vol[1]), I, snr(0);

	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
		cs += m1.cross_section(cp.volt());
		csf += m1.elastic_cross_section(cp.volt(), 0, cp.frequency_cutoff());
		if ( cp.slit_width() < 0.01 )
			csf += m1.inelastic_cross_section(cp.volt(), 0, cp.frequency_cutoff());
	}
	mu = cs/area;
	muf = csf/area;
	I = exp(muf-mu);
	
	if ( verbose ) {
		cout << "Calculating the signal for each material:" << endl;
		cout << "Area:                           " << area << " A2" << endl;
		cout << "Thickness:                      " << vol[2] << " A" << endl;
		cout << "Intensity recorded:             " << I << endl;
		cout << "Scattering average:             " << mu << endl;
		cout << "Material\tCross(A2))\tSignal\tSNR" << endl;
	}

	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
//		cse = m1.cross_section(cp.volt(), 0, cp.frequency_cutoff(), cp.slit_width());
		cse = m1.elastic_cross_section(cp.volt(), 0, cp.frequency_cutoff());
//		if ( cp.slit_width() < 0.01 )
//			cse += m1.inelastic_cross_section(cp.volt(), 0, cp.frequency_cutoff());
		cses += cse;
		mus = (cse/area) * exp(-mu);
		snr += mus/I;
		muss += mus;
		if ( verbose )
			cout << setw(15) << m1.identifier() << tab << cse << tab << mus << tab << mus/I << endl;
	}
	
	if ( verbose )
		cout << setw(15) << "Total" << tab << cses << tab << muss << tab << snr << endl << endl;
	
	return mu;
}

double		particle_signal(Bmaterial& material, Bmaterial ice, CTFparam& cp, double area, double thickness)
{
	if ( area <  1) {
		cerr << "Error: The area must be specified!" << endl;
		bexit(-1);
	}

	ice.volume(area*thickness - material.volume());
	
	Bmaterial	mat_ice = material_combine(material, ice, 1, 1);

	double		cs = material_cross_section(mat_ice, cp);

	double		csel = material_elastic_cross_section(material, cp);

	double		signal = (csel/area)*exp(-cs/area);
	
	return signal;
}

double		micrograph_recorded_intensity(Bmaterial& material, Bmaterial ice, CTFparam& cp, double area, double thickness)
{
	if ( area <  1) {
		cerr << "Error: The area must be specified!" << endl;
		bexit(-1);
	}

	double		emfp = micrograph_emfp(material, ice, cp, area, thickness);

	double		signal = micrograph_signal(material, ice, cp, area, thickness);
	double		part_signal = particle_signal(material, ice, cp, area, thickness);

	double		intensity = exp(-thickness/emfp);
	
	cout << "Material mass:                  " << material.mass() << " Da" << endl;
	cout << "Ice mass:                       " << ice.mass() << " Da" << endl;
	cout << "Area:                           " << area << " A2" << endl;
	cout << "Thickness:                      " << thickness << " A" << endl;
	cout << "Volume:                         " << area*thickness << " A3" << endl;
	cout << "Effective mean free path:       " << emfp << " A" << endl;
	cout << "Relative intensity:             " << intensity << endl;
	cout << "Micrograph signal:              " << signal << endl;
	cout << "Particle signal:                " << part_signal << endl << endl;

	return intensity;
}

int			micrograph_recorded_intensity(Bmaterial& material, Bmaterial& ice, CTFparam& cp, double area)
{
	double		t, emfp, intensity, signal, part_signal;
	
	double		radius = pow((3*material.mass())/(4*M_PI*material.density(DA_A3)), 1.0/3.0);
	if ( area < radius*radius ) area = 4*radius*radius;
	
	double		t_min(100*long(radius/50+1)), t_max(10*t_min);

	cout << "Material mass:                  " << material.mass() << " Da" << endl;
	cout << "Area:                           " << area << " A2" << endl;

	cout << "T(A)\tEMFP\tInt\tS(mg)\tS(part)" << setprecision(5) << endl;
	for ( t = t_min; t <= t_max; t += t_min ) {
		emfp = micrograph_emfp(material, ice, cp, area, t);
		intensity = exp(-t/emfp);
		signal = micrograph_signal(material, ice, cp, area, t);
		part_signal = particle_signal(material, ice, cp, area, t);
		cout << t << tab << emfp << tab << intensity << tab << signal << tab << part_signal << endl;
	}
	cout << endl;
	
	return 0;
}


int			aperture_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> apser)
{
	double		oa, ang, Iunf, Ifil, lam_unf, lam_fil;
	
	cout << "Aper\ta(mr)\tLunf\tLfil\tIunf\tIfil\tln(Iu)\tln(If)\tIf/Iu\tln(If/Iu)" << endl;
	for ( auto it = apser.begin(); it != apser.end(); ++it ) {
		oa = *it;
		if ( oa < 1e4 ) oa *= 1e4;
		ctf.objective_aperture(oa);
		ang = 1e3*oa/(2*ctf.focal_length());
		ctf.slit_width(0);
		Iunf = signal_intensity(material, thickness, ctf);
		lam_unf = -thickness/log(Iunf);
		ctf.slit_width(20);
		Ifil = signal_intensity(material, thickness, ctf);
		lam_fil = -thickness/log(Ifil);
		cout << *it << tab << setprecision(4) << ang << tab
			<< lam_unf << tab << lam_fil << tab
			<< Iunf << tab << Ifil << tab
			<< -log(Iunf) << tab << -log(Ifil) << tab
			<< Ifil/Iunf << tab << -log(Ifil/Iunf) << endl;
	}

	return 0;
}

int			collection_angle_series(double thickness, Bmaterial& material, CTFparam& ctf, vector<double> angser)
{
//	map<string,Bcomptype>& 	comp = material.composition();

	double		oa, ang, Iunf, Ifil, lam_unf, lam_fil;
	
	cout << "Aper\ta(mr)\tLunf\tLfil\tIunf\tIfil\tln(Iu)\tln(If)\tIf/Iu\tln(If/Iu)" << endl;
	for ( auto it = angser.begin(); it != angser.end(); ++it ) {
		ang = *it;
		if ( ang > 1 ) ang *= 1e-3;
		oa = 2*ang*ctf.focal_length();
		ctf.objective_aperture(oa);
		ctf.slit_width(0);
		Iunf = signal_intensity(material, thickness, ctf);
		lam_unf = -thickness/log(Iunf);
		ctf.slit_width(20);
		Ifil = signal_intensity(material, thickness, ctf);
		lam_fil = -thickness/log(Ifil);
		cout << oa*1e-4 << tab << setprecision(4) << ang*1e3 << tab
			<< lam_unf << tab << lam_fil << tab
			<< Iunf << tab << Ifil << tab
			<< -log(Iunf) << tab << -log(Ifil) << tab
			<< Ifil/Iunf << tab << -log(Ifil/Iunf) << endl;
	}

	return 0;
}


double		particle_snr(Bmaterial& material, double mass, double radius, double thickness, CTFparam& ctf)
{
	map<string,Bcomptype>& comp = material.composition();

	if ( mass < 1 ) mass = material.mass();

	// Calculate radius from mass, assuming a globular/spherical protein
	if ( radius < 1 ) radius = 50;
		radius = pow((3*mass)/(4*M_PI*material.density(DA_A3)), 1.0/3.0);

	if ( thickness < 1 ) thickness = 2*radius;

	double			cs = elastic_cross_section_integrated(comp, ctf);
	double			ics = inelastic_cross_section_langmore(comp, ctf);
	double			emfp = effective_mean_free_path(material, ctf);
	double			snre = 0.5*cs/(M_PI*radius*radius);

	if ( verbose ) {
		cout << "Particle:" << endl;
		cout << "Molecular weight:               " << mass << " Da" << endl;
		cout << "Particle radius:                " << radius << " A" << endl;
		cout << "Density:                        " << material.density(G_CM3) << " g/cm3" << endl;
		cout << "Density:                        " << material.density(DA_A3) << " Da/A3" << endl;
		cout << "Density:                        " << material.density(NUM_A3) << " #/A3" << endl;
		cout << "Occupancy:                      " << mass/(M_PI*radius*radius*thickness) << " Da/A3" << endl;
		cout << "Elastic cross section:          " << cs << " A2" << endl;
		cout << "Inelastic cross section:        " << ics << " A2" << endl;
		cout << "Particle SNR per electron:      " << snre << " /e" << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
		cout << endl;
	}
	
	return snre;
}

Bplot*		particle_spectral_signal(Bmaterial& material, double mass,
	double radius, double thickness, CTFparam& ctf)
{
	map<string,Bcomptype>& comp = material.composition();

	double			ds(0.5/radius);
	double			scut(ctf.frequency_cutoff());

	map<string, vector<double>>		sc = calculate_elastic_scattering_curves(comp, ds, scut);

	string			e;
	long			i, ns(sc.begin()->second.size());
	double			pf(ctf.lambda()/sqrt(1-beta2(ctf.volt())));
	double			w, area(M_PI * radius * radius);
	Bplot*			plot = new Bplot(1,ns,2);
	
	cout << "Prefactor: " << pf << endl;
	
	for ( i=0; i<ns; ++i ) {
		(*plot)[i] = ds*i;
		(*plot)[i+ns] = 0;
	}

	for ( auto it: comp ) {
		e = it.first;
		w = pf * it.second.component_count();
		cout << e << tab << w << endl;
		i = ns;
		for ( auto is = sc[e].begin(); is != sc[e].end(); ++is )
			(*plot)[i++] += *is * w;
	}
	
	vector<double>	dp = defocus_range_profile(ctf, ds, 5e3, 2e4, 1e3);
	
	for ( i=0; i<ns; ++i ) {
		(*plot)[i+ns] *= (*plot)[i+ns];
		(*plot)[i+ns] *= dp[i]/area;
	}

	cout << "s\tAmp\t" << endl;
	for ( i=0; i<ns; ++i )
		cout << (*plot)[i] << tab << (*plot)[i+ns] << tab << dp[i] << endl;

	return plot;
}

