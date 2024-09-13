/**
@file	Bmaterial.h
@brief	Header to read and write model dynamics parameters in STAR format
@author Bernard Heymann
@date	Created: 20100305
@date	Modified: 20230808
**/

#include "ctf.h"
#include "Bplot.h"
#include "utilities.h"

#ifndef _Bmaterial_
#define _Bmaterial_
#define	MATERIAL_NAME		"material.name"
#define	MATERIAL_DENSITY	"material.density"
#define	MATERIAL_DENS_UNIT	"material.density_unit"

enum DensityUnit {
	G_CM3 = 0,
	DA_A3 = 1,
	NUM_A3 = 2
} ;

/************************************************************************
@Object: class Bcomptype
@Description:
	Model component type parameter structure.
@Features:
	Reference to a model file, either atomic coordinates or a 3D map.
	Number of components of this type.
	Figure-of-merit.
*************************************************************************/
class Bcomptype {
private:
	void	initialize() {
		next = NULL;
		id = "?";
		ind = 0;
		fmod = "?";
		num = 0;
		cnt = 0;
		mas = 1;
		chrg = 0;
		lrad = 1;
		rad = 1;
		fom = 0;
		sel = 1;
		for ( int i=0; i<11; ++i ) coef[i] = 0;
	}
public:
	Bcomptype*		next;			// Next component type in list
private:
	string			id;				// Component type identifier
	int				ind;			// Type index for reference in a distance matrix, element number
	string			fmod;			// Model file name (molecule or map)
	int				num;			// Image number in model or map file
	long			cnt;			// Number of components
	float			mas;			// Component mass
	float			chrg;			// Component charge
	float			lrad;			// Link radius (bond radius; half link length)
	float			rad;			// Characteristic radius (VdW radius)
	float			coef[11];		// Cross section coefficients
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Bcomptype() { initialize(); }
	Bcomptype(string s) { initialize(); id = s; }
	Bcomptype(long i) { initialize(); id = to_string(i); }
	Bcomptype(string fn, long img_num) { initialize(); id = "1"; fmod = fn; num = img_num; }
	Bcomptype(string s, string fn, long img_num) { initialize(); id = s; fmod = fn; num = img_num; }
	Bcomptype(const Bcomptype& ct) {
		id = ct.id;
		ind = ct.ind;
		fmod = ct.fmod;
		num = ct.num;
		cnt = ct.cnt;
		mas = ct.mas;
		chrg = ct.chrg;
		lrad = ct.lrad;
		rad = ct.rad;
		fom = ct.fom;
		sel = ct.sel;
		coefficients(ct.coefficients());
//		cout << "In Bcomptype:" << endl;
//		show();
	}
	Bcomptype(Bcomptype* ct) {
		next = NULL;
		id = ct->id;
		ind = ct->ind;
		fmod = ct->fmod;
		num = ct->num;
		cnt = ct->cnt;
		mas = ct->mas;
		chrg = ct->chrg;
		lrad = ct->lrad;
		rad = ct->rad;
		fom = ct->fom;
		sel = ct->sel;
		coefficients(ct->coefficients());
//		cout << "In Bcomptype:" << endl;
//		show();
	}
	void			identifier(string s) { id = s; }
//	void			identifier(Bstring s) { id = s.str(); }
	string&			identifier() { return id; }
	void			index(long i) { ind = i; }
	long			index() { return ind; }
	long			index(string s) {
		Bcomptype*		c = find(s);
		if ( c ) return c->ind;
		else return -1;
	}
	long			maximum_index() {
		long		 mx(0);
		for ( Bcomptype* c=this; c; c=c->next )
			if ( mx < c->index() ) mx = c->index();
		return mx;
	}
	void			file_name(string s) { fmod = s; }
	string&			file_name() { return fmod; }
	void			image_number(long i) { num = i; }
	long			image_number() { return num; }
	void			mass(double m) { mas = m; }
	double			mass() { return mas; }
	void			charge(double c) { chrg = c; }
	double			charge() { return chrg; }
	void			link_radius(double r) { lrad = r; }
	double			link_radius() { return lrad; }
	void			radius(double r) { rad = r; }
	double			radius() { return rad; }
	template <typename T>
	void			coefficients(vector<T> c) {
		for ( int i=0; i<11 && i<c.size(); ++i ) coef[i] = c[i];
	}
	vector<double>	coefficients() const {
		vector<double>	c(11);
		for ( int i=0; i<11; ++i ) c[i] = coef[i];
		return c;
	}
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	double			FOM_add(double d) { fom += d; return fom; }
	double			FOM_multiply(double d) { fom *= d; return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
	long			select_increment() { return sel++; }
	void			component_count(long n) { cnt = n; }
	long			component_count() { return cnt; }
	void			component_count_increment() { cnt++; }
	void			component_count_decrement() { cnt--; }
	void			component_count_add(long n) { cnt += n; }
	Bcomptype*		add(string s) {
		Bcomptype* 		c(this);
		while ( c->id != s && c->next ) c = c->next;
		if ( c->next ) return c;
		return c->next = new Bcomptype(s);
	}
	Bcomptype*		add(long i) {
		Bcomptype* 		c(this);
		string			s(to_string(i));
		while ( c->id != s && c->next ) c = c->next;
		if ( c->next ) return c;
		return c->next = new Bcomptype(s);
	}
	Bcomptype*		add(string fn, long img_num) {
		Bcomptype* 		c(this);
		while ( !( c->fmod == fn && c->num == img_num ) && c->next ) c = c->next;
		if ( c->fmod == fn && c->num == img_num ) return c;
		return c->next = new Bcomptype(fn, img_num);
	}
	Bcomptype*		add(string s, string fn, long img_num) {
		Bcomptype* 		c(this);
		while ( !( c->id == s && c->fmod == fn && c->num == img_num ) && c->next ) c = c->next;
		if ( c->fmod == fn && c->num == img_num ) return c;
		return c->next = new Bcomptype(s, fn, img_num);
	}
	Bcomptype*		add(Bcomptype* ct) {
		Bcomptype* 		c(this);
		while ( c->next ) c = c->next;
		return c->next = new Bcomptype(ct);
	}
	Bcomptype*		find(string s) {
		Bcomptype* 		c(this);
		while ( c && c->id != s ) c = c->next;
		return c;
	}
	long			count() {
		long			n(0);
		for ( Bcomptype* c=this; c; c=c->next ) n++;
		return n;
	}
	long			count_selected() {
		long			n(0);
		for ( Bcomptype* c=this; c; c=c->next ) if ( sel ) n++;
		return n;
	}
	void		update(Bcomptype& ct) {
		id = ct.id;
		ind = ct.ind;
		fmod = ct.fmod;
		num = ct.num;
		cnt = ct.cnt;
		mas = ct.mas;
		chrg = ct.chrg;
		lrad = ct.lrad;
		rad = ct.rad;
		fom = ct.fom;
		sel = ct.sel;
		coefficients(ct.coefficients());
	}
	bool			check() {
		if ( id.empty() ) id = "?";
		if ( fmod.empty() ) fmod = "?";
		if ( mas < 1 ) mas = 1;
		return 1;
	}
	double			elastic_scattering_amplitude(double s) {
		long			j, k;
		double			s2(s*s), amp(coef[10]);
		for ( j=0, k=5; j<5; ++j, ++k )
			amp += coef[j]*exp(-coef[k]*s2);
		return amp;
	}
	vector<double>	elastic_scattering_curve(double ds, double scut) {
		if ( ds < 1e-10 ) ds = 0.01;
		long			i, imax(scut/ds + 1);
		vector<double>	scat;
		for ( i=0; i<imax; ++i )
			scat.push_back(elastic_scattering_amplitude(ds*i));
		return scat;
	}
	// Burge & Smith (1962) formula
	double			elastic_scattering_integral(double slo, double shi) {
		long			i, j, k, l;
		double			a, b, cs(0);
		for ( i=0, k=5; i<5; ++i, ++k ) {
			for ( j=0, l=5; j<5; ++j, ++l ) {
				b = coef[k] + coef[l];
				a = coef[i]*coef[j]/b;
				cs += a*(exp(-b*slo*slo) - exp(-b*shi*shi));
			}
		}
		return cs;
	}
	double			elastic_cross_section(double volt) {
		return elastic_cross_section(volt, 0, 1e10);
	}
	double			elastic_cross_section(double volt, double slo, double shi) {
		double		wlr = electron_wavelength_relativistic(volt);
		return 4*M_PI*wlr*wlr*elastic_scattering_integral(slo, shi);
	}
	double			inelastic_cross_section_langmore(double volt) {
		return inelastic_cross_section_langmore(volt, 0, 1e10);
	}
	double			inelastic_cross_section_langmore(double volt, double slo, double shi) {
		double		wlr = electron_wavelength_relativistic(volt);
		double		g = lorentz(volt);
		double		cs = 0.25*wlr*wlr*sqrt(ind)*log(2*(1+1/g)*volt/20);
		if ( ind == 1 ) cs *= 0.154;	// Correction for hydrogen
		cs *= shi/(shi + 0.054) - slo/(slo + 0.054);
		return cs;
	}
	void			show() {
		cout << "Identity:                       " << id << endl;
		cout << "Index:                          " << ind << endl;
		cout << "File name:                      " << fmod << endl;
		cout << "Number:                         " << num << endl;
		cout << "Count:                          " << cnt << endl;
		cout << "Mass:                           " << mas << endl;
		cout << "Charge:                         " << chrg << endl;
		cout << "Link radius:                    " << lrad << endl;
		cout << "Radius:                         " << rad << endl;
		cout << "FOM:                            " << fom << endl;
		cout << "Select:                         " << sel << endl;
		cout << "Coefficients:" << endl;
		for ( int i=0; i<11; ++i ) cout << i+1 << tab << coef[i] << endl;
	}
} ;

/************************************************************************
@Object: class Blinktype
@Description:
	A structure for a link type.
@Features:
	A link is defined between two components.
*************************************************************************/
class Blinktype {
private:
	void	initialize() {
		len = dist = 1;
		Klen = Kdist = 1;
		fom = 0;
		sel = 1;
	}
	string		id[2];		// Component type identifiers
	float		len;		// Covalent bond length
	float		dist;		// Van der Waals interaction distance
	float		Klen;		// Covalent bond strenght
	float		Kdist;		// Van der Waals interaction strength
	float		fom;		// Figure-of-merit
	int			sel;		// Selection flag
public:
	Blinktype() { initialize(); }
	Blinktype(string& id1, string& id2) {
		initialize();
		id[0] = id1;
		id[1] = id2;
	}
	Blinktype(string& id1, string& id2, double d) {
		initialize();
		id[0] = id1;
		id[1] = id2;
		len = d;
	}
	void			identifier(string s, int i) { id[i] = s; }
	string&			identifier(int i) { return id[i]; }
	void			length(double d) { len = d; }
	double			length() { return len; }
	void			distance(double d) { dist = d; }
	double			distance() { return dist; }
	void			Klength(double d) { Klen = d; }
	double			Klength() { return Klen; }
	void			Kdistance(double d) { Kdist = d; }
	double			Kdistance() { return Kdist; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
} ;

/************************************************************************
@Object: class Bangletype
@Description:
	Model angle parameter structure.
@Features:
	An angle is defined for three components.
*************************************************************************/
class Bangletype {
private:
	void	initialize() {
		ang = 0;
		Kang = 1;
		fom = 0;
		sel = 1;
	}
	string			id[3];			// List of component type identifiers
	float			ang;			// Reference angle
	float			Kang;			// Angle interaction strength
	float			fom;			// Figure-of-merit
	int				sel;			// Selection flag
public:
	Bangletype() { initialize(); }
	Bangletype(string& id1, string& id2, string& id3) {
		initialize();
		id[0] = id1;
		id[1] = id2;
		id[2] = id3;
	}
	void			identifier(string s, int i) { id[i] = s; }
	string&			identifier(int i) { return id[i]; }
	void			angle(double d) { ang = d; }
	double			angle() { return ang; }
	void			Kangle(double d) { Kang = d; }
	double			Kangle() { return Kang; }
	void			FOM(double d) { fom = d; }
	double			FOM() { return fom; }
	void			select(long i) { sel = i; }
	long			select() { return sel; }
} ;

/************************************************************************
@Object: class Bmaterial
@Description:
	Material parameter structure.
@Features:
	A material has an elemental composition and density.
*************************************************************************/
class Bmaterial {
private:
	string				id;		// Material name
	double				den;	// Density in g/cm3
	DensityUnit			unt;	// Density units: 0 = g/cm3, 1 = Da/Å3, 2 = number/Å3
	long				sel;	// Selection
	map<string,Bcomptype>	comp;	// Elemental composition
	void		initialize() {
		den = 0;
		unt = G_CM3;
		sel = 1;
	}
public:
	Bmaterial() { initialize(); }
	Bmaterial(string s) { initialize(); id = s; }
	Bmaterial(string s, double d, DensityUnit u) {
		initialize(); id = s; den = d; unt = u;
	}
	void		identifier(string s) { id = s; }
	string		identifier() { return id; }
	void		select(long i) { sel = i; }
	long		select() { return sel; }
	double		density(double d, DensityUnit u) { den = d; unt = u; return den; }
	double		density(DensityUnit u) {
		if ( u == G_CM3 ) return gram_per_cm3();
		else if ( u == DA_A3 ) return dalton_per_angstrom3();
		else return number_per_angstrom3();
	}
	double		gram_per_cm3() {
		if ( unt == DA_A3 ) return den*1.0e24/AVOGADRO;
		else if ( unt == NUM_A3 ) return den*mass()*1.0e24/AVOGADRO;
		return den;
	}
	double		dalton_per_angstrom3() {
		if ( unt == G_CM3 ) return AVOGADRO*den/1.0e24;
		else if ( unt == NUM_A3 ) return den*mass();
		return den;
	}
	double		number_per_angstrom3() {
		if ( unt == G_CM3 ) return AVOGADRO*den/(1.0e24*mass());
		else if ( unt == DA_A3 ) return den/mass();
		return den;
	}
	double		unit(DensityUnit u) {
		if ( u == unt ) return den;
		if ( u == G_CM3 ) den = gram_per_cm3();
		else if ( u == DA_A3 ) den = dalton_per_angstrom3();
		else den = number_per_angstrom3();
		unt = u;
		return den;
	}
	DensityUnit	unit() { return unt; }
	map<string,Bcomptype>&	composition() { return comp; }
	map<int,Bcomptype>		composition_numbered() {
		map<int,Bcomptype>	cn;
		for ( auto ct: comp )
			cn[ct.second.index()] = ct.second;
		return cn;
	}
	Bcomptype&	operator[](string s) { return comp[s]; }
	long		type_count() { return comp.size(); }
	long		number() {
		long		n(0);
		for ( auto i: comp ) n += i.second.component_count();
		return n;
	}
	double		mass() {
		double		m(0);
		for ( auto ct: comp ) m += ct.second.mass() * ct.second.component_count();
		return m;
	}
	void		mass(double m) {
		double		r(m/mass());
		for ( auto &ct: comp ) ct.second.component_count(ct.second.component_count()*r);
	}
	double		volume() {
		if ( den < 1e-10 )
			cerr << "Error: The density for " << id << " is not specified!" << endl;
		return mass()/dalton_per_angstrom3();
	}
	void		volume(double v) {
		double		r(v*number_per_angstrom3());
		for ( auto &ct: comp ) ct.second.component_count(ct.second.component_count()*r);
	}
	double		mean_inner_potential() {
		double			scale(1e20*PLANCK*PLANCK/(TWOPI*ECHARGE*EMASS));
		double			mip(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			mip += ct.component_count() * ct.elastic_scattering_amplitude(0);
		}
		mip *= scale*dalton_per_angstrom3()/mass();
		return mip;
	}
	double		elastic_cross_section(double volt) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * ct.elastic_cross_section(volt);
		}
		return cs;
	}
	double		elastic_cross_section(double volt, double slo, double shi) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * ct.elastic_cross_section(volt, slo, shi);
		}
		return cs;
	}
	double		inelastic_cross_section(double volt) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * ct.inelastic_cross_section_langmore(volt);
		}
		return cs;
	}
	double		inelastic_cross_section(double volt, double slo, double shi) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * ct.inelastic_cross_section_langmore(volt, slo, shi);
		}
		return cs;
	}
	double		cross_section(double volt) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * (ct.elastic_cross_section(volt) + ct.inelastic_cross_section_langmore(volt));
		}
		return cs;
	}
/*	double		cross_section(double volt, double slo, double shi, double slit=0) {
		double			cs(0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			cs += ct.component_count() * ct.elastic_cross_section(volt, slo, shi);
			if ( slit < 0.001 ) cs += ct.component_count() * ct.inelastic_cross_section_langmore(volt, slo, shi);
			else cs += ct.component_count() * ct.inelastic_cross_section_langmore(volt);
		}
		return cs;
	}*/
	vector<double>	elastic_cross_section_curve(double volt, double ds, double scut) {
		long			n(scut/ds+1);
		vector<double>	cs(n,0);
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			long			cnt = ct.component_count();
			vector<double>	cs1 = ct.elastic_scattering_curve(ds, scut);
			for ( long i=0; i<n; ++i )
				cs[i] += cnt * cs1[i] * cs1[i];
		}
		double		wlr = electron_wavelength_relativistic(volt);
		for ( long i=0; i<n; ++i )
			cs[i] *= 8*M_PI*wlr*wlr;
		return cs;
	}
	void		update_parameters(map<string,Bcomptype>& types) {
		for ( auto &it: comp ) {
			Bcomptype&		ct = it.second;
			if ( types.find(it.first) != types.end() ) {
				Bcomptype&	at = types[it.first];
				ct.identifier(at.identifier());
				ct.index(at.index());
				ct.mass(at.mass());
				ct.charge(at.charge());
				ct.coefficients(at.coefficients());
			}
		}
	}
	void		show() {
		long		cnt(0);
		cout << "Material:                       " << id << endl;
		cout << "Density:                        " << fixed << den;
		if ( unt == G_CM3 ) cout << " g/cm3" << endl;
		else if ( unt == DA_A3 ) cout << " Da/A3" << endl;
		else cout << " #/A3" << endl;
		cout << "Composition:" << endl;
		cout << "Element\tCount\tMass(Da)" << endl;
		for ( auto ct: comp ) {
			cnt += ct.second.component_count();
			cout << ct.second.identifier() << tab << ct.second.component_count()
				<< tab << ct.second.mass()*ct.second.component_count() << endl;
		}
		cout << "Total" << tab << cnt << tab << mass() << endl;
		cout << "Volume:                         " << volume() << " A3" << endl;
		cout << "Mean inner potential:           " << mean_inner_potential() << " V" << endl << endl;
	}
} ;

#endif

// Function prototypes
long		materials_list(map<string,Bmaterial>& mlist);
Bmaterial	material_from_elements(string elements);
map<string,Bcomptype>	material_combined_composition(map<string,Bmaterial>& material);
Bmaterial	material_combine(vector<Bmaterial>& mlist, vector<double>& numbers);
Bmaterial	material_combine(map<string,Bmaterial>& mlist, vector<double>& numbers);
Bmaterial	material_combine(Bmaterial& m1, Bmaterial& m2, double n1, double n2);
Bmaterial	material_protein(long natom);
long		material_protein_add_hydrogens(Bmaterial& material);
Bmaterial	material_ice(long nmol);
Bmaterial	material_ice_from_volume(double vol);
double		material_mean_free_path(map<string,Bmaterial>& material, double volt);
double		material_effective_mean_free_path(map<string,Bmaterial>& material, double volt, double scut, double slit);
Bplot*		material_radial_power_spectrum(map<string,Bmaterial>& material, double volt, double scut, double slit);

