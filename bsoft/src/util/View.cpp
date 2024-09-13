/**
@file	View.cpp
@brief	View functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20230831
**/

#include "View2.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, View2<float> v) {
	output.setf(ios::fixed, ios::floatfield);
	output << v.vector3() << tab << setw(10) << v.angle()*180.0/M_PI;
	return output;
}

ostream& operator<<(ostream& output, View2<double> v) {
	output.setf(ios::fixed, ios::floatfield);
	output << v.vector3() << tab << setw(10) << v.angle()*180.0/M_PI;
	return output;
}

/**
@brief 	Displays a linked list of views.
@param 	&views			the list of views.
@return int				number of views.
**/
int			show_views(vector<View2<float>>& views)
{
	for ( auto& v: views )
		cout << v << endl;
	
	return views.size();
}

int			show_views(vector<View2<double>>& views)
{
	for ( auto& v: views )
		cout << v << endl;
	
	return views.size();
}

/**
@brief 	Calculates a random view.
@return View2<double>	random view.

	A random seed should already have been generated.

**/
View2<double>		random_view()
{
	double			irm = 2.0/get_rand_max();
	View2<double>	v;
	
	v[0] = random()*irm - 1;
	v[1] = random()*irm - 1;
	v[2] = random()*irm - 1;
	v[3] = M_PI*(random()*irm - 1);
	v.normalize();
	
	return v;
}

/**
@brief 	Generates a random reslicing 3x3 rotation matrix.
@return View2<double> 	new view.

	The view represents any one or more 90 degree rotations,
	randomly chosen.

**/
View2<double>	view_random_reslice()
{
	double			irm = 1.0/get_rand_max();
	View2<double>	v;
	
	int				rv = (int) (5.999*irm*random());
	int				ra = (int) (3.999*irm*random());
	switch ( rv ) {
		case 0: v = View2<double>(1,0,0,ra*M_PI_2); break;
		case 1: v = View2<double>(0,1,0,ra*M_PI_2); break;
		case 2: v = View2<double>(0,0,1,ra*M_PI_2); break;
		case 3: v = View2<double>(-1,0,0,ra*M_PI_2); break;
		case 4: v = View2<double>(0,-1,0,ra*M_PI_2); break;
		case 5: v = View2<double>(0,0,-1,ra*M_PI_2); break;
	}
	
	return v;
}

/**
@brief 	Calculates a set of random views.
@param 	nviews				number of views.
@return vector<View2<double>>	list of random views.
**/
vector<View2<double>>	random_views(long nviews)
{
	vector<View2<double>>	views;
	
	if ( nviews < 1 ) {
		error_show("Error in random_views: No views defined!", __FILE__, __LINE__);
		return views;
	}
	
	random_seed();
	
	for ( long i=0; i<nviews; ++i )
		views.push_back(random_view());
	
	return views;
}


/**
@brief 	Calculates a new view with a random error from the given view.
@param 	&v				given view (modified).
@param 	std				standard deviation of gaussian distribution.
@return double			random angle.

	The rotation between the given and new views is defined by the 
	gaussian distributed random angle and a random axis.
	A random seed should already have been generated.

**/
double		random_view_error(View2<double>& v, double std)
{
	double			angle = random_gaussian(0.0, std);

	Vector3<double>	axis = vector3_random(-1.0, 1.0);
	axis.normalize();

	Quaternion		q(v.quaternion()*Quaternion(axis, angle));
	
	v = View2<double>(q);
	
	return angle;
}

/**
@brief 	Initializes a set of views tilted around the y axis.
@param 	ang_min			starting angle (radians).
@param 	ang_max			ending angle (radians).
@param 	ang_step		angular step size (radians).
@param 	axis			tilt axis angle (radians).
@return vector<View2<double>>	a set of 4-value views.

	A set of views is calculated corresponding to tilted views imaged
	during tomography. The tilt axis angle is taken as a counter-clockwise
	rotation from the x-axis.

**/
vector<View2<double>>	tilt_views(double ang_min, double ang_max, double ang_step, double axis)
{
	if ( ang_max < ang_min ) swap(ang_max, ang_min);
	if ( ang_step < 0 ) ang_step = -ang_step;
	
	double					ang;
	int						n(0);
	vector<View2<double>>	views;
	
	for ( n=0, ang=ang_min; ang<ang_max+0.5*ang_step; ang+=ang_step, n++ )
		views.push_back(View2<double>(ang, axis));
	
	if ( verbose & VERB_PROCESS )
		cout << "Getting " << n << " views for tilting from " << ang_min*180.0/M_PI <<
			" to " << ang_max*180.0/M_PI << " by step " << ang_step*180.0/M_PI <<
			" around axis " << axis*180.0/M_PI << endl;
	
	return views;
}


/**
@brief 	Generates a list of views within an angular distance from the input view.
@param 	theview				the input view.
@param 	theta_step			theta step size (radians).
@param 	phi_step			phi step size (radians).
@param 	alpha_step			alpha step size (radians).
@param 	view_angle_limit	angular distance limit from view vector (radians).
@param 	alpha_angle_limit	angular distance limit from view rotation angle (radians).
@return vector<View2<double>>	a list of views.

	The list of views forms a 3D search grid in orientation space.

**/
vector<View2<double>>	views_within_limits(View2<double> theview, double theta_step, double phi_step,
					double alpha_step, double view_angle_limit, double alpha_angle_limit)
{
	double			theta, phi, alpha;
	double			theta_start = (theta_step)? -theta_step*floor(view_angle_limit/theta_step): 0;
	double			phi_start = (phi_step)? -phi_step*floor(view_angle_limit/phi_step): 0;
	double			alpha_start = theview.angle();
	double			alpha_end = alpha_start;
	
	if ( alpha_step > 0 && alpha_angle_limit > 0 ) {
		alpha = alpha_step*floor(alpha_angle_limit/alpha_step);
		alpha_start -= alpha;
		alpha_end += alpha;
	}
	
	if ( alpha_step <= 0 ) alpha_step = 1;
	
	Vector3<double>	ref_axis(0,0,1);
		
	theview.normalize();
	
	Vector3<double>	vv(theview[0], theview[1], theview[2]);
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	vector<View2<double>>	v;
	Quaternion    		q1, q;
	Quaternion    		qv(0, theview[0], theview[1], theview[2]);
	
	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_within_limits: The View: " << theview << endl;
		cout << "DEBUG views_within_limits: theta_step: "<< theta_step*180/M_PI
			<< ", phi_step: " << phi_step*180/M_PI
			<< ", alpha_step: " << alpha_step*180/M_PI << endl;
		cout << "DEBUG views_within_limits: Axis1: " << axis1 << endl;
		cout << "DEBUG views_within_limits: Axis2: " << axis2 << endl;
	}
	
	if ( theta_step < 1e-6 ) theta_step = 1e6;
	if ( phi_step < 1e-6 ) phi_step = 1e6;
	if ( alpha_step < 1e-6 ) alpha_step = 1e6;
	for ( theta = theta_start; theta <= view_angle_limit; theta += theta_step ) {
		q1 = Quaternion(axis1, theta);
		for ( phi = phi_start; phi <= view_angle_limit; phi += phi_step ) {
			q = Quaternion(axis2, phi);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha = alpha_start; alpha <= alpha_end; alpha += alpha_step ) {
				v.push_back(View2<float>(q[1], q[2], q[3], alpha));
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_within_limits: " << v.back() << " - " << theview.angle(v.back())*180/M_PI << endl;
			}
		}
	}
	
	return v;
}


/**
@brief 	Generates a 3x3 grid of views around an input view.
@param 	theview				the input view.
@param 	alpha_step			angular step size around view vector (radians).
@return vector<View2<double>>	a list of views.
**/
vector<View2<double>>	views_for_refinement(View2<double> theview, double alpha_step)
{
	theview.normalize();

	vector<View2<double>>	views;

	if ( alpha_step <= 0 ) {
		views.push_back(theview);
		return views;
	}
	
	double			alpha1, alpha2, alpha3;
	Vector3<double>	ref_axis(0,0,1);
			
	Vector3<double>	vv(theview.vector3());
	
	Vector3<double>	axis1 = vv.cross(ref_axis);
	if ( axis1.length() < 0.001 ) axis1[0] = 1;
	axis1.normalize();
	
	Vector3<double>	axis2 = vv.cross(axis1);
	if ( axis2.length() < 0.001 ) axis2[1] = 1;
	axis2.normalize();
	
	Quaternion    	q1, q;
	Quaternion    	qv(0, theview[0], theview[1], theview[2]);
	View2<double>	v;

	if ( verbose & VERB_FULL )
		cout << qv << endl;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "DEBUG views_for_refinement: The View: " << theview << endl;
		cout << "DEBUG views_for_refinement: alpha_step: "<< alpha_step*180/M_PI << endl;
		cout << "DEBUG views_for_refinement: Axis1: " << axis1 << endl;
		cout << "DEBUG views_for_refinement: Axis2: " << axis2 << endl;
	}
	
	for ( alpha1 = -alpha_step; alpha1 <= alpha_step; alpha1 += alpha_step ) {
		q1 = Quaternion(axis1, alpha1);
		for ( alpha2 = -alpha_step; alpha2 <= alpha_step; alpha2 += alpha_step ) {
			q = Quaternion(axis2, alpha2);
			q *= q1;
			q = q.rotate(qv);
			for ( alpha3 = -alpha_step; alpha3 <= alpha_step; alpha3 += alpha_step ) {
				v = View2<double>(q[1], q[2], q[3], theview.angle() + alpha3);
				views.push_back(v);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG views_for_refinement: " << v << " - " << theview.angle(v)*180/M_PI << endl;
			}
		}
	}
	
	return views;
}


/**
@brief 	Expands each view to several views with different rotation angles.
@param 	views			view list.
@param 	amin			minimum angle.
@param 	amax			maximum angle.
@param 	astep			angular step.
@return vector<View2<double>>	new view list.

	The new angles are added to the existing angles of the view.

**/
vector<View2<double>>	view_list_expand_angles(vector<View2<double>> views, double amin, double amax, double astep)
{
	double					a;
	vector<View2<double>>	nuviews;
	View2<double>			vn;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG view_list_expand_angles: amin=" << amin << " amax=" << amax << " astep=" << astep << endl;
	
	for ( auto& v: views ) {
		if ( verbose & VERB_FULL )
			cout << "Old view: " << v << endl;
		for ( a = amin; a <= amax; a += astep ) {
			vn = v;
			vn[3] += a;
			nuviews.push_back(vn);
			if ( verbose & VERB_FULL )
				cout << vn << endl;
		}
	}
	
	if ( verbose & VERB_FULL )
		show_views(nuviews);
	
	return nuviews;
}

/**
@brief 	Replaces a view list with a subset of it.
@param 	&view_list		view list.
@param 	start			offset of first view of subset.
@param 	size			number of views in subset.
@return vector<View2<double>>	new list ov views.

**/
vector<View2<double>>	view_list_subset(vector<View2<double>>& view_list, int start, int size)
{
	vector<View2<double>>	nu_list(view_list.begin()+start, view_list.begin()+start+size);
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of views selected:       " << nu_list.size() << endl;
	
	return nu_list;
}

/**
@brief 	Generates a JSON list of views.
@param 	&views			view list.
@return int				number of views selected, <0 on error.

**/
JSvalue		js_views(vector<View2<double>>& views)
{
	vector<double>	vv;
	JSvalue			js(JSarray);
	
	for ( auto v: views ) {
		vv = v.array();
		vv[3] *= 180.0/M_PI;
		js.push_back(JSvalue(vv));
	}
		
	return js;
}

