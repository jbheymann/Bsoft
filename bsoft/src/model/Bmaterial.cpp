/**
@file	Bmaterial.cpp
@brief	Functions for material roperties and calculating electron scattering profiles
@author Bernard Heymann
@date	Created: 20190521
@date	Modified: 20230720
**/

#include "Bmaterial.h"
#include "string_util.h"
#include "Color.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Prints a list of materials.
@param 	&mlist				list of materials.
@return long					number of materials.
**/
long		materials_list(map<string,Bmaterial>& mlist)
{
	cout << "Material" << tab << "Density(g/cm3)" << tab << "Density(Da/A3)" << tab << "Density(#/A3)" << endl;
	for ( auto m: mlist )
		cout << m.first << tab << m.second.gram_per_cm3() << tab << m.second.dalton_per_angstrom3() << tab << m.second.number_per_angstrom3() << endl;
	
	return mlist.size();
}

/**
@brief 	Prints a list of materials.
@param 	elements			comma-separated list of elements.
@return Bmaterial				new material.
**/
Bmaterial	material_from_elements(string elements)
{
	Bmaterial				material("FromElements: " + elements);
	map<string,Bcomptype>&	comp = material.composition();
	map<string,Bcomptype>	types;
	
	vector<string>			vs = split(elements,',');
	
	for ( auto el: vs )
		comp[el] = 1;

	return material;
}

/**
@brief 	Retruns a combined composition.
@param 	&material			list of materials.
@return map<string,Bcomptype>	new composition.
**/
map<string,Bcomptype>	material_combined_composition(map<string,Bmaterial>& material)
{
	map<string,Bcomptype>	nucomp;
	
	for ( auto m: material ) {
		Bmaterial&		m1 = m.second;
		for ( auto ct: m1.composition() ) {
			if ( nucomp.find(ct.first) == nucomp.end() ) {
				nucomp[ct.first] = ct.second;
				nucomp[ct.first].component_count(0);
			}
			nucomp[ct.first].component_count_add(ct.second.component_count());
		}
	}
	
	return nucomp;
}

/**
@brief 	Assembles a component composition from a set of materials with fractional contributions.
@param 	&mlist				list of materials.
@param 	&numbers			array of copies of each material.
@return Bmaterial				new combined material.
**/
Bmaterial	material_combine(vector<Bmaterial>& mlist, vector<double>& numbers)
{
	long					i(0);
	double					sum(0), d(0);
	Bmaterial				material("Combined");
//	material.identifier("Combined");
	map<string,Bcomptype>&	nucomp = material.composition();
	map<string,Bcomptype>	types;
	
	for ( auto n: numbers ) sum += n;
	
	for ( auto m: mlist ) {
		d += numbers[i] * m.dalton_per_angstrom3();
		map<string,Bcomptype>	comp = m.composition();
		for ( auto ct: comp ) {
			if ( nucomp.find(ct.first) == nucomp.end() ) {
				nucomp[ct.first] = ct.second;
				nucomp[ct.first].component_count(0);
			}
			nucomp[ct.first].component_count_add(numbers[i] * ct.second.component_count());
		}
		i++;
	}
	
	material.density(d/sum, DA_A3);
	
	return material;
}

/**
@brief 	Assembles a component composition from a set of materials with fractional contributions.
@param 	&mlist				list of materials.
@param 	&numbers			array of copies of each material.
@return Bmaterial				new combined material.
**/
Bmaterial	material_combine(map<string,Bmaterial>& mlist, vector<double>& numbers)
{
	long					i(0);
	double					sum(0), vol(0), mass(0);
	Bmaterial				material("Combined");
	map<string,Bcomptype>&	nucomp = material.composition();
	
	if ( numbers.size() < mlist.size() ) numbers.resize(mlist.size(),1);
	for ( auto n: numbers ) sum += n;
	
	if  ( verbose ) {
		cout << "Combining materials:" << endl;
		cout << "Material\tFraction" << endl;
	}
	
	for ( auto m: mlist ) {
		Bmaterial&		m1 = m.second;
		if ( verbose )
			cout << m1.identifier() << tab << numbers[i]/sum << endl;
//		m1.show();
		mass += numbers[i] * m1.mass();
		vol += numbers[i] * m1.volume();
		for ( auto ct: m1.composition() ) {
			if ( nucomp.find(ct.first) == nucomp.end() ) {
				nucomp[ct.first] = ct.second;
				nucomp[ct.first].component_count(0);
			}
			nucomp[ct.first].component_count_add(numbers[i] * ct.second.component_count());
		}
		i++;
	}
	
	if ( verbose ) cout << endl;
	
	material.density(mass/vol, DA_A3);
	
	return material;
}

/**
@brief 	Assembles a component composition from two materials with fractional contributions.
@param 	&m1				first materials.
@param 	&m2				second materials.
@param 	n1				copies of first material.
@param 	n2				copies of second material.
@return Bmaterial			new combined material.
**/
Bmaterial	material_combine(Bmaterial& m1, Bmaterial& m2, double n1, double n2)
{
	vector<Bmaterial>	mlist;
	vector<double>		numbers;
	mlist.push_back(m1);
	mlist.push_back(m2);
	numbers.push_back(n1);
	numbers.push_back(n2);
	return material_combine(mlist, numbers);
}

/**
@brief 	Default protein composition.
@param	natom			number of atoms.
@return Bmaterial			material.
**/
Bmaterial	material_protein(long natom)
{
	Bmaterial			material;
	map<string,Bcomptype>&	comp = material.composition();
	
	material.identifier("Protein");
	material.density(RHO, DA_A3);

	comp["H"] = Bcomptype("H");
	comp["C"] = Bcomptype("C");
	comp["N"] = Bcomptype("N");
	comp["O"] = Bcomptype("O");
	comp["S"] = Bcomptype("S");

	comp["H"].component_count(0.498*natom);
	comp["C"].component_count(0.320*natom);
	comp["N"].component_count(0.085*natom);
	comp["O"].component_count(0.095*natom);
	comp["S"].component_count(0.002*natom);
	
	return material;
}

/**
@brief 	Adds hydrogens to a protein composition.
@return long			number of hydrigens.
**/
long		material_protein_add_hydrogens(Bmaterial& material)
{
	map<string,Bcomptype>&	comp = material.composition();
	
	comp["H"].component_count(comp["C"].component_count()*0.498/0.320);
	
	return comp["H"].component_count();
}

/**
@brief 	Vitreous ice material.
@param	nmol			number of molecules.
@return Bmaterial			material.
**/
Bmaterial	material_ice(long nmol)
{
	Bmaterial		ice;
	map<string,Bcomptype>&	comp = ice.composition();

	ice.identifier("Vitreous ice");
	ice.density(0.93, G_CM3);	// Density in g/cm3
	ice.unit(DA_A3);			// Density in Da/A3
	
	comp["H"] = Bcomptype("H");
	comp["O"] = Bcomptype("O");

//	comp["H"].mass(1);
//	comp["O"].mass(16);

	comp["H"].component_count(2*nmol);
	comp["O"].component_count(nmol);

	return ice;
}

/**
@brief 	Calculating the true mean free path of a material.
@param 	material		list of materials.
@param	volt			acceleration voltage (V).
@return double			mean free path.
**/
double		material_mean_free_path(map<string,Bmaterial>& material, double volt)
{
	double			cs(0), cs_comb(0), vol(0);
	
	if ( verbose ) {
		cout << "Calculating the mean free path at " << volt << " V" << endl;
		cout << "Material\tCross(A2)\tVolume" << endl;
	}
	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
		cs = m1.cross_section(volt);
		cs_comb += cs;
		vol += m1.volume();
		if ( verbose )
			cout << setw(15) << m1.identifier() << tab << cs << tab << m1.volume() << endl;
	}
	
	double			mfp = vol/cs_comb;
	
	if ( verbose ) {
		cout << setw(15) << "Total" << tab << cs_comb << tab << vol << endl;
		cout << "Mean free path:                 " << mfp << " A" << endl << endl;
	}

	return mfp;
}

/**
@brief 	Calculating the true mean free path of a material.
@param	material		list of materials.
@param	volt			acceleration voltage (V).
@param	scut			frequency cutoff (/A).
@param	slit			energy filter slit width (eV).
@return double			mean free path.
**/
double		material_effective_mean_free_path(map<string,Bmaterial>& material, double volt, double scut, double slit)
{
	double			cs(0), cs_comb(0), vol(0);
	
	if ( verbose ) {
		cout << "Calculating the effective mean free path at " << volt << " V" << endl;
		cout << "Material\tCross(A2)\tVolume(A3)" << endl;
	}
	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
//		cs = m1.cross_section(volt, scut, 1e10, slit);
		cs = m1.elastic_cross_section(volt, scut, 1e10);
		if ( slit ) cs += m1.inelastic_cross_section(volt);
		else cs += m1.inelastic_cross_section(volt, scut, 1e10);
		cs_comb += cs;
		vol += m1.volume();
		if ( verbose )
			cout << setw(15) << m1.identifier() << tab << cs << tab << m1.volume() << endl;
	}
	
	double			emfp = vol/cs_comb;
	
	if ( verbose ) {
		cout << "Total" << tab << cs_comb << tab << vol << endl << endl;
		cout << "Effective mean free path:       " << emfp << " A" << endl;
		cout << "Aperture cutoff frequency:      " << scut << " /A" << endl;
		if ( slit )
			cout << "Energy filter slit width:       " << slit << " eV" << endl;
		cout << endl;
	}

	return emfp;
}

Bplot*		material_radial_power_spectrum(map<string,Bmaterial>& material, double volt, double scut, double slit)
{
	long			i, j, k;
	double			ds(0.01), vol(0);
	long			n(scut/ds+1);
	vector<double>	rps(n,0);

	for ( auto m: material ) {
		Bmaterial& 		m1 = m.second;
		vector<double>	cs = m1.elastic_cross_section_curve(volt, ds, scut);
		for ( i=0; i<n; ++i ) rps[i] += cs[i];
		vol += m1.volume();
	}

	Bstring			title("Radial power spectrum");
	long			nc(2);
	RGB<float>		color;

	Bplot*			plot = new Bplot(1, n, nc);
	plot->title(title);
	plot->page(0).title(title);
	plot->page(0).columns(nc);
	for ( i=0; i<nc; i++ ) plot->page(0).column(i).number(i);
	plot->page(0).column(0).label("SpatialFrequency(1/A)");
	plot->page(0).column(0).axis(1);
	plot->page(0).axis(0).min(0);
	plot->page(0).axis(0).max(scut);
	for ( i=1; i<nc; ++i ) {
		plot->page(0).column(i).type(2);
		plot->page(0).column(i).label(Bstring(i, "%d"));
		plot->page(0).column(i).axis(3);
		if ( nc > 2 ) color.spectrum(i,1,nc-1);
		plot->page(0).column(i).color(color.r(),color.g(),color.b());
	}
//	plot->page(0).axis(3).flags(2);
	plot->page(0).axis(3).label("Power");
	for ( i=0; i<n; ++i ) (*plot)[i] = ds*i;
	for ( i=0, k=n; i<n; ++i, ++k ) (*plot)[k] = rps[i]/vol;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Spatial Frequency (1/A)";
		for ( k=1; k<nc; ++k ) cout << tab << k;
		cout << scientific << endl;
		for ( i=0; i<n; ++i ) {
			cout << (*plot)[i];
			for ( j=n+i, k=1; k<nc; ++k, j+=n ) cout << tab << (*plot)[j];
			cout << endl;
		}
	}
		
	return plot;
}
