/**
@file	View.h
@brief	View object
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20230524
**/

#ifndef _View2_
#define _View2_

#include "Vector3.h"
#include "Matrix3.h"
#include "Quaternion.h"
#include "random_numbers.h"
#include "json.h"
#include "utilities.h"

/************************************************************************
@Object: class View
@Description:
	View/orientation parameter structure.
@Features:
	Pointer for a linked list.
	Composed of a view vector and rotation angle around the vector.
*************************************************************************/
template <typename Type>
class View2 {
private:
	Vector3<Type>	v;
	Type			a;
	void	check() {	// Sets small numbers to zero and adjust the angle to [-PI,PI]
		for ( size_t i = 0; i<3; i++ ) if ( fabs(v[i]) < TRIGPRECISION ) v[i] = 0;
		a = angle_set_negPI_to_PI(a);
		if ( fabs(a) < TRIGPRECISION ) a = 0;
	}
public:
	View2() : v(0,0,1), a(0) {  }	// Initialize to {0,0,1,0} with a NULL pointer
	View2(const View2& view) : v(view.v), a(view.a) { normalize(); }	// Copies from another view, not the pointer
	template <typename T>
	View2(Vector3<T>& vec) : v(vec), a(0) { normalize(); }
	View2(vector<double>& vec) : v(vec), a(0) {
		normalize();
		if ( vec.size() > 3 ) a = vec[3];
	}
	View2(const double x, const double y, const double z, const double angle)
		: v(x,y,z), a(angle) {		// Constructs from 4 values
		normalize();
	}
	// View from another form
	View2(Quaternion& q) {	// Convert from a quaternion
		double		ang = atan2(q[3], q[0]);
		double		ca = cos(ang), sa = sin(ang);
		double		q2 = q[0]*q[0] + q[3]*q[3];
		double		f2 = 2*sqrt(q2);
		if ( q2 < SMALLFLOAT ) ang = atan2(q[1], q[2]);
		v[0] = f2*(q[2]*ca + q[1]*sa);
		v[1] = f2*(q[2]*sa - q[1]*ca);
		v[2] = 2*q2 - 1;
		a = 2*ang;
		normalize();
	}
	View2(Matrix3 m) {		// Convert from a matrix
		Quaternion	q = m.quaternion();
		*this = View2(q);
	}
	View2(double tilt, double axis) :
			v(-sin(axis)*sin(tilt), cos(axis)*sin(tilt), cos(tilt)), a(0) {
		normalize();
	}
	View2(double angle, Vector3<double> axis) {
		Matrix3		mat = Matrix3(axis, angle);
		*this = View2(mat);
	}
	// Operators
	template <typename T2> operator View2<T2>() const {
		return View2<T2>((T2)v[0], (T2)v[1], (T2)v[2], (T2)v[3]);
	}
	View2	operator=(const View2& view) {	// Copies the contents, not the pointer
		v = view.v; a = view.a;
		return *this;
	}
	View2	operator=(const vector<double>& vec) {
		size_t	i;
		for ( i=0; i<3 && i<vec.size(); ++i ) v[i] = vec[i];
		if ( i < vec.size() ) a = vec[i]; 
		return *this; 	
	}
	View2	operator-() {			// Returns the negative of the view
		View2	vn(-v[0], -v[1], -v[2], -a);
		return vn;
	}
	bool	operator==(const View2& view) {
		return ( v == view.v ) && ( a == view.a );
	}
	Type&	operator[](size_t i) { if ( i < 3 ) return v[i]; return a; }
	double	x() { return v[0]; }
	double	y() { return v[1]; }
	double	z() { return v[2]; }
	double	angle() { return a; }
	void	x(const double d) { v[0] = d; }
	void	y(const double d) { v[1] = d; }
	void	z(const double d) { v[2] = d; }
	void	angle(const double d) { a = d; }
	vector<double>	array() {
		return {v[0], v[1], v[2], a};
	}
	double	vector_size() {	// Returns the size of the vector
		return v.length();
	}
	double	normalize() {
		double		size = v.length();
		if ( size < SMALLFLOAT ) { v[0] = v[1] = 0; v[2] = 1; size = 1; }
		else { v[0] /= size; v[1] /= size; v[2] /= size; }
		check();
		return size;
	}
	// View conversions
	void	negate() { v = -v; a = -a; }	// Negates the view
	View2	backward() {			// Returns the backwards/inverse form of the view
		double		ca = cos(a);
		double		sa = sin(a);
		View2		backview;
		backview.v[0] = -v[0]*ca - v[1]*sa;
		backview.v[1] =  v[0]*sa - v[1]*ca;
		backview.v[2] =  v[2];
		backview.a =  angle_set_negPI_to_PI(-a);
		backview.check();
		return backview;
	}
	Vector3<double> vector3() { return v; }
	Quaternion	quaternion() {		// Converts to a quaternion
		normalize();
		double		ang = a/2.0L;
		double		ca = cos(ang), sa = sin(ang);
		double		z1 = v[2] + 1, f1, f2;
		Quaternion	q(ca, 0, 0, sa);			// z = 1
		if ( v[2] < 1.0L - SMALLFLOAT ) {
			if ( v[2] <= SMALLFLOAT - 1.0L ) {	// z = -1
				q[0] = 0;
				q[1] = sa;
				q[2] = ca;
				q[3] = 0;
			} else {							// -1 < z < 1
				f1 = sqrt(z1/2.0);
				f2 = sqrt(1.0/(2*z1));
				q[0] = f1*ca;
				q[1] = f2*(v[0]*sa - v[1]*ca);
				q[2] = f2*(v[0]*ca + v[1]*sa);
				q[3] = f1*sa;
			}
		}
		return q;
	}
	Matrix3	matrix() {		// Converts to a matrix
		Quaternion	q = quaternion();
		return Matrix3(q);
	}
	void	random_view() {
		double			irm = 2.0/get_rand_max();
		v[0] = random()*irm - 1;
		v[1] = random()*irm - 1;
		v[2] = random()*irm - 1;
		a = M_PI*(random()*irm - 1);
		normalize();
	}
// View comparisons
	double	distance(const View2& view) {	// Returns the distance between view vectors
		double		d = a - view.a;
		return sqrt((v - view.v).length2() + d*d);
	}
	double	angle(const View2& view) {	// Returns the angle between view vectors
		double		d(v.scalar(view.v));
		if ( d < 1 && d > -1 ) return acos(d);
		else if ( d < 0 ) return M_PI;
		else return 0;
	}
	double	residual(const View2& view) {	// Returns a residual difference between views
		double		d = 1 - cos(a - view.a);
		return sqrt((v - view.v).length2() + d*d)/2;
	}
	long		find_closest(vector<View2<double>>& views) {
		long		i, closest(0);
		double		a, amin(M_PI);
		for ( i=0; i<views.size(); ++i ) {
			a = angle(views[i]);
			if ( amin > a ) {
				amin = a;
				closest = i;
			}
		}
		return closest;
	}
} ;

//template <typename Type>
//ostream& operator<<(ostream& output, View2<Type> v);
ostream& operator<<(ostream& output, View2<float> v);
ostream& operator<<(ostream& output, View2<double> v);


#endif

// Function prototypes
int			show_views(vector<View2<float>>& views);
int			show_views(vector<View2<double>>& views);
View2<double>			random_view();
View2<double>			view_random_reslice();
vector<View2<double>>	random_views(long nviews);
double		random_view_error(View2<double>& v, double std);
vector<View2<double>>	tilt_views(double ang_min, double ang_max, double ang_step, double axis);
vector<View2<double>>	views_within_limits(View2<double> theview, double theta_step, double phi_step,
					double alpha_step, double view_angle_limit, double alpha_angle_limit);
vector<View2<double>>	view_list_expand_angles(vector<View2<double>> views, double amin, double amax, double astep);
vector<View2<double>>	view_list_subset(vector<View2<double>>& view_list, int start, int size);
vector<View2<double>>	views_for_refinement(View2<double> theview, double alpha_step);
JSvalue		js_views(vector<View2<double>>& views);
//View*		view_array(View* views, long& n);


