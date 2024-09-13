/**
@file	rwmodel_param.h
@brief	Header to read and write model dynamics parameters in STAR format
@author 	Bernard Heymann
@date	Created: 20100305
@date	Modified: 20230717
**/

//#include "Bmaterial.h"
#include "rwmodel.h"

#ifndef _Bmodparam_
#define _Bmodparam_
#define	DISTMAT_COMPTYPE1	"distance_matrix.component_type1"
#define	DISTMAT_COMPTYPE2	"distance_matrix.component_type2"
#define	DISTMAT_LLEN		"distance_matrix.link_length"
#define	DISTMAT_DIST		"distance_matrix.distance"
#define	DISTMAT_KL			"distance_matrix.link_constant"
#define	DISTMAT_KD			"distance_matrix.distance_constant"
#define	DISTMAT_POTTYPE		"distance_matrix.potential_type"

#define ANGLEMAT_COMPTYPE1	"angle_matrix.component_type1"
#define ANGLEMAT_COMPTYPE2	"angle_matrix.component_type2"
#define ANGLEMAT_COMPTYPE3	"angle_matrix.component_type3"
#define ANGLEMAT_ANGLE		"angle_matrix.angle"
#define ANGLEMAT_KA			"angle_matrix.angle_constant"

/************************************************************************
@Object: class Bmodparam
@Description:
	A structure used for model mechanics.
@Features:
	This defines variables and constants for doing model dynamics or
	other model mechanics.
*************************************************************************/
class Bmodparam {
public:
	string		comment;
	double		timestep;		// Integration time step for MD
	double		Kfriction;		// Friction coefficient
	double		Kdistance;		// Distance/Van der Waals force constant
	double		Kelec;			// Electrostatic force constant
	double		Klink;			// Link/bond force constant
	double		Kangle;			// Angular force constant
	double		Kpolyangle;		// Angular force constant for polygons
	double		Kpolygon;		// Polygon regularity force constant
	double		Kpolyplane;		// Polygon planarity force constant
	double		Kpoint;			// Point force constant
	double		Kradial;		// Radial force constant
	double		Kplane;			// Neighbor plane force constant
	double		Kguide;			// Polyhedron guide force constant
	double		Kmap;			// Density map force constant
	double		sepdist;		// Separation distance grid sampling
	double		cutoff;			// Distance cutoff for non-bonded calculations
	double		pointdecay;		// Decay constant for the point force
	double		radius;			// Radius for radial force
	Vector3<double>	point;		// Center of point or radial force
	Bmodel*		guide;			// Polyhedron guide model
	int			distancetype;	// 1=harmonic, 2=soft, 3=lennard-jones, 4=morse
	int			linksteps;		// Number of sampling intervals along a link
	Vector3<double>	min, max;	// Boundary of box
	bool		wrap;			// Flag to turn periodic boundaries on
	double		sigma;			// Gaussian decay for density fitting
	double		Edistance;		// Non-linked distance energy
	double		Eelec;			// Electrostatic energy
	double		Elink;			// Link energy
	double		Eangle;			// Angle energy
	double		Epolyangle;		// Angle energy for polygons
	double		Epolygon;		// Polygon regularity energy
	double		Epolyplane;		// Polygon planarity energy
	double		Epoint;			// Point force energy
	double		Eradial;		// Radial force energy
	double		Eplane;			// Neighbor plane force energy
	double		Eguide;			// Polyhedron guide force energy
	double		Emap;			// Density map associated energy
	double		Ekin;			// Kinetic energy
	double		Epot;			// Potential energy
	map<string,Bcomptype>		comptype;
	vector<vector<Blinktype>>	linktype;
	vector<vector<vector<Bangletype>>>	angletype;
private:
	void	initialize() {
		comment = "?";
		timestep = 1;
		Kfriction = 1;
		Kdistance = 0;
		Kelec = 0;
		Klink = 0;
		Kangle = 0;
		Kpolyangle = 0;
		Kpolygon = 0;
		Kpolyplane = 0;
		Kpoint = 0;
		Kradial = 0;
		Kplane = 0;
		Kguide = 0;
		Kmap = 0;
		sepdist = 10;
		cutoff = 10;
		pointdecay = 0.01;
		radius = 0;
		distancetype = 1;
		linksteps = 3;
		wrap = 0;
		sigma = 0;
	}
public:
	Bmodparam()	{ initialize(); }
	void			minimum(Vector3<double> v) { min = v; }
	Vector3<double>	minimum() { return min; }
	void			maximum(Vector3<double> v) { max = v; }
	Vector3<double>	maximum() { return max; }
	Vector3<double>	box() { return max - min; }
	void		show() {
//		cout << "Kfriction:                      " << Kfriction << endl;
		cout << "Kdistance:                      " << Kdistance << " (" << distancetype << ")" << endl;
		cout << "Klink:                          " << Klink << endl;
		cout << "Kangle:                         " << Kangle << endl;
		cout << "Kpolyangle:                     " << Kpolyangle << endl;
		cout << "Kpolygon:                       " << Kpolygon << endl;
		cout << "Kpolyplane:                     " << Kpolyplane << endl;
		cout << "Kpoint and decay:               " << Kpoint << " (" << pointdecay << ")" << endl;
		cout << "Kradial and radius:             " << Kradial << " (" << radius << ")" << endl;
		cout << "Kplane:                         " << Kplane << endl;
		cout << "Box mimimum:                    " << min << endl;
		cout << "Box maximum:                    " << max << endl;
	}
	void		show_types() {
		cout << "Component types:" << endl << "#\tName\tMass" << endl;;
		for ( auto c: comptype )
			cout << c.second.index() << tab << c.first << tab << c.second.mass() << endl;
	}
	void		show_links() {
		long		i(0), j(0);
		cout << "Link types:" << endl << "#1\t#2\tLength\tDistance" << endl;
		for ( auto lr: linktype ) {
			j = 0;
			for ( auto l: lr ) {
				cout << i << tab << j << tab << l.length() << tab << l.distance() << endl;
				j++;
			}
			i++;
		}
	}
	void		show_angles() {
		long		i(0), j(0), k(0);
		cout << "Angle types:" << endl << "#1\t#2\t#3\tAngle" << endl;
		for ( auto am: angletype ) {
			j = 0;
			for ( auto ar: am ) {
				k = 0;
				for ( auto a: ar ) {
					cout << i << tab << j << tab << k << tab << a.angle()*180.0/M_PI << endl;
					k++;
				}
				j++;
			}
			i++;
		}
	}
} ;
#endif

// Function prototypes
Bmaterial	material_from_model(Bmodel* model);
Bmodparam	model_param_generate(Bmodel* model);
int			model_param_generate(Bmodparam& md, Bmodel* model);
int			model_param_set_type_indices(Bmodel* model, Bmodparam& mp);
int			model_update_reference_parameters(Bmodel* model, Bmodparam& md);
map<string,Bcomptype> 	read_atom_properties(string& filename);
map<string,Bmaterial> 	read_material_properties(string& filename);
map<string,Bmaterial> 	read_material_properties(vector<string> file_list, string paramfile, int flags);
int			write_material_properties(string filename, map<string,Bmaterial> material);
Bmodparam 	read_dynamics_parameters(string filename);
int			update_dynamics_parameters(Bmodparam& md, string filename);
int		 	write_dynamics_parameters(string filename, Bmodparam& md);

