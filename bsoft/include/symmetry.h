/**
@file	symmetry.h
@brief	Header file for general symmetry functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20230623
**/

#include "View2.h"
#include "string_util.h"

#ifndef _Bsymmetry_
/************************************************************************
@Object: class Bsymop
@Description:
	General symmetry operator.
@Features:
	Both point groups and space groups are covered.
	A symmetry operator describes the complete series of rotations around 
		an axis, with the number of operations given by the order field
		and the rotation angle given by the angle field. The shift field
		is used in gliding operations, such as needed for helical symmetry.
*************************************************************************/
class Bsymop {
private:
	Vector3<double>	ax;		// Symmetry axis
	int				n;		// Number of times the operation is applied, 1=identity
	double			ang;	// Rotation angle, derived from order when 0
	double			sh;		// Shift up symmetry axis for helix
	void	initialize() {
		n = 0;
		ang = 0;
		sh = 0;
	}
public:
	Bsymop() { initialize(); }
	Bsymop(double x, double y, double z, int nn, double a) :
		ax(Vector3<double>(x,y,z)), n(nn), ang(a), sh(0) {}
	void			order(int i) { n = i; }
	int				order() { return n; }
	void			angle(double a) { ang = a; }
	double			angle() { return ang; }
	void			shift(double s) { sh = s; }
	double			shift() { return sh; }
	void			axis(Vector3<double> a) { ax = a; }
	void			axis(double x, double y, double z) { ax = Vector3<double>(x,y,z); }
	Vector3<double>&	axis() { return ax; }
	double			normalize() { return ax.normalize(); }
	Matrix3			matrix() {
		return Matrix3(ax, M_PI*2.0/n);
	}
} ;

/************************************************************************
@Object: class Bsymmetry
@Description:
	General symmetry structure.
@Features:
	Both point groups and space groups are covered.
*************************************************************************/
struct Bsymmetry {
private:
	int 			pnt;		// Point group (< 1 if a crystal)
	int 			sp;			// Space group (< 1 if not a crystal)
	string			lbl;		// Symmetry label string
	vector<Bsymop>	op; 		// Symmetry operators
private:
	void	initialize() {
		pnt = 0;
		sp = 0;
		lbl = "C1";
	}
	string		clean_symstring(string& sym) {
		int					j(0);
		lbl = "C1";
		if ( sym.length() < 1 ) return lbl;
		// Remove leading blanks and convert to upper case
		lbl = sym;
		remove_spaces(lbl);
		// Alternate nomenclature for point groups
		if ( isdigit(lbl[0]) ) {
			j = to_integer(lbl);
			if ( lbl.find("532") != string::npos && lbl.find("90") != string::npos )
				lbl = "I90";
			else if ( lbl.find("532") != string::npos )
				lbl = "I";
			else if ( lbl.find("432") != string::npos )
				lbl = "O";
			else if ( lbl.find("23") != string::npos )
				lbl = "T";
			else if ( j > 1 && lbl[1] == 2 )
				lbl = "D" + lbl.substr(0,1);
			else if ( j > 0 )
				lbl = "C" + to_string(j);
			else
				lbl = "C1";
		}
		lbl[0] = toupper(lbl[0]);
		return lbl;
	}
public:
	Bsymmetry() { initialize(); }
	Bsymmetry(string sym);
//	~Bsymmetry() { op.clear(); }
	void		point(int i) { pnt = i; }
	int			point() { return pnt; }
	int			point() const { return pnt; }
	void		space(int i) { sp = i; }
	int			space() { return sp; }
	void		label(string& s) { lbl = s; }
	string&		label() { return lbl; }
	int		order() {
		int		ns(1);
		for ( auto it = op.begin(); it != op.end(); ++it ) ns *= it->order();
		return ns;
	}
	int			operations() { return op.size(); }
	Bsymop&		operator[](int i) { return op[i]; }
	vector<Matrix3>	matrices() {
		long				ns(1), k;
		Matrix3				mt;
		vector<Matrix3>		mat;
		mat.push_back(Matrix3(1));
		for ( auto it = op.begin(); it != op.end(); ++it ) {
			for ( int i = 1; i < it->order(); ++i ) {
				mt = Matrix3(it->axis(), i*M_PI*2.0L/it->order());
				for ( k=0; k<ns; ++k )
					mat.push_back(mt*mat[k]);
			}
			ns *= it->order();
		}
//		cout << "Number of matrices = " << mat.size() << endl;
		return mat;
	}
	void		transform(Matrix3& mat) {
		for ( auto it = op.begin(); it != op.end(); ++it )
			it->axis(mat * it->axis());
	}
	vector<Vector3<double>>	get_axes();
	View2<double>	reference_symmetry_view();
	vector<View2<double>>	reference_asymmetric_unit_views() {
		View2<double>	ref = reference_symmetry_view();
		return get_all_views(ref);
	}
	int 			asymmetric_unit_index(View2<double> theview);
	Matrix3			rotate_to_axis(long axis, long axis_flag);
	vector<View2<double>>	get_all_views(View2<double> asu_view);
	vector<View2<double>>	asymmetric_unit_views(double theta_step, double phi_step, int flag);
	vector<View2<double>>	asymmetric_unit_views(double theta_step, double phi_step, double alpha_step, int flag);
	View2<double>	find_asymmetric_unit_view(View2<double> theview);
	View2<double>	find_closest_symmetric_view(View2<double> view_ref, View2<double> view);
	vector<View2<double>>	side_views(double side_ang, double theta_step, double phi_step);
	vector<View2<double>>	side_views(double side_ang, double theta_step, double phi_step, double alpha_step);
	void			change_views_to_asymmetric_unit(vector<View2<double>>& views) {
		if ( pnt < 102 ) return;
		for ( auto& v: views )
			v = find_asymmetric_unit_view(v);
	}
	View2<double>	random_view(View2<double>& asu_view) {
		vector<View2<double>>	views = get_all_views(asu_view);
		long			i(order()*random()*0.99999L/get_rand_max());
		return views[i];
	}
	int				show_matrices();
	int				show_operational_matrices();
	int				show_pdb_matrices();
} ;
#define _Bsymmetry_
#endif

// Function prototypes
string		symmetry_helical_label(double helix_rise, double helix_angle,
				int dyad_axis, int cyclic, double seam_shift);

