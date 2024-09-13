/**
@file	symmetry.cpp
@brief	General symmetry functions
@author Bernard Heymann
@date	Created: 20010420
@date	Modified: 20230524
**/

#include "symmetry.h"
#include "random_numbers.h"
#include "Euler.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Symmetry constructor from a symmetry identifier.
@param 	sym			string containing point group identifier.

	The point group symmetries are identified by the following strings:
		C<n>		cyclic point group of order n.
		D<n>			dihedral point group of order n.
		T			tetrahedral point group.
		O			octahedral/cubic point group.
		I			icosahedral/dodecahedral point group.
		H<r>,<a>,<d>	helical symmetry with rise r, rise angle a and dyad d (1/2).
	For the higher symmetries the following adjustments are available:
		T-3			no three-fold operator.
		O-2			no two-fold operator.
		O-3			no three-fold operator.
		O-4			no four-fold operator.
		I-2			no two-fold operator.
		I-3			no three-fold operator.
		I-5			no five-fold operator.
		I90			90 degrees rotated around z-axis.
		I90-3			90 degrees rotated around z-axis and no three-fold operator.
	If the point group string is empty, the default is C1 (asymmetric).

**/
Bsymmetry::Bsymmetry(string sym)
{
	int 			theorder, helix_dyad(1);
	double			helix_rise = 1, helix_angle = M_PI;
	Matrix3 		mat(0,-1,0,1,0,0,0,0,1); // 90 degree rotation for I90
	
	clean_symstring(sym);
	
 	if ( verbose & VERB_FULL )
		cout << endl << "Getting symmetry operators for " << lbl << endl;
	
	if ( lbl[0] == 'C' ) {									// Cyclic point groups
		if ( lbl == "Cs" ) theorder = 1;					// Reflection
		else theorder = to_integer(lbl.substr(1,10));
		pnt = 100 + theorder;
		op.push_back(Bsymop(0,0,1,theorder,M_PI*2.0L/theorder));
	} else if ( lbl[0] == 'D' ) {							// Dihedral point groups
		theorder = to_integer(lbl.substr(1,10));
		pnt = 200 + theorder;
		op.push_back(Bsymop(0,0,1,theorder,M_PI*2.0L/theorder));
		op.push_back(Bsymop(1,0,0,2,M_PI));					// Rotate around 2-fold x-axis
	} else if ( lbl[0] == 'T' ) {							// Tetrahedral point group
		pnt = 320;
		if ( lbl.find("-3") == string::npos )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		op.push_back(Bsymop(0,0,1,2,M_PI));					// Rotate around 2-fold z-axis
		op.push_back(Bsymop(1,0,0,2,M_PI));					// Rotate around 2-fold x-axis
	} else if ( lbl[0] == 'O' ) {							// Octahedral point group
		pnt = 432;
		if ( lbl.find("-3") == string::npos )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		if ( lbl.find("-4") != string::npos ) {
			op.push_back(Bsymop(1,-1,0,2,M_PI));			// Rotate around 2-fold axis
		} else if ( lbl.find("-2") != string::npos ) {
			op.push_back(Bsymop(0,0,1,4,M_PI_2));			// Rotate around 4-fold z-axis
		} else {
			op.push_back(Bsymop(0,0,1,4,M_PI_2));			// Rotate around 4-fold z-axis
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Rotate around 2-fold x-axis
		}
	} else if ( lbl[0] == 'I' ) {							// Icosahedral point group
		pnt = 532;
		if ( lbl.find("-3") == string::npos )
			op.push_back(Bsymop(1,1,1,3,M_PI*2.0L/3.0));	// Rotate around 3-fold (1,1,1)
		if ( lbl.find("-5") != string::npos ) {
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Rotate around 2-fold x-axis
		} else if ( lbl.find("-2") != string::npos ) {
			op.push_back(Bsymop(1,-1,0,2,M_PI));			// Rotate around 2-fold axis
		} else {
			op.push_back(Bsymop(1,1.0L/GOLDEN,GOLDEN,2,M_PI));	// Rotate around 2-fold
		}
		if ( lbl.find("-5") == string::npos ) {
			op.push_back(Bsymop(1.0L/GOLDEN,1,0,5,M_PI*2.0L/5.0));	// Rotate around 5-fold
		}
		if ( lbl.find("-2") == string::npos ) {
			op.push_back(Bsymop(0,0,1,2,M_PI));				// Rotate around 2-fold z-axis
		}
		if ( lbl.find("90") != string::npos )
			transform(mat);
	} else if ( lbl[0] == 'H' ) {							// Helical symmetry
		sscanf(lbl.c_str(), "H%lf,%lf,%d", &helix_rise, &helix_angle, &helix_dyad);
		pnt = 600;
		op.push_back(Bsymop(0,0,1,1,helix_angle*M_PI/180.0L));	// Helical axis
		op[0].shift(helix_rise);							// Rise in angstrom
		if ( helix_dyad == 2 )
			op.push_back(Bsymop(1,0,0,2,M_PI));				// Dyad axis
	} else {
		error_show("Error in Bsymmetry::Bsymmetry", __FILE__, __LINE__);
		cerr << "This is not a valid symmetry designation: " << lbl << endl;
		return;
	}
	
	for ( auto it = op.begin(); it != op.end(); ++it )
		it->normalize();
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bsymmetry::Bsymmetry: lbl=" << lbl << endl;
}

/**
@brief 	Returns the label for helical symmetry.
@param 	helix_rise		helical rise (angstroms).
@param	helix_angle		helical rotation angle (radians).
@param	dyad_axis		presence of dyad axis (1/2).
@param	cyclic			cyclic symmetry.
@param	seam_shift		fractional shift for seamed helices.
@return string			label.

	Thge symmetry order is defined as the product of all the individual
	orders of the symmetry operations, or alternatively, the number of views.

**/
string		symmetry_helical_label(double helix_rise, double helix_angle,
								   int dyad_axis, int cyclic, double seam_shift)
{
	string			label = "H" + to_string(helix_rise);
	label += "," + to_string(helix_angle*180.0/M_PI);
	label += "," + to_string(dyad_axis);
	label += "," + to_string(cyclic);
	label += "," + to_string(seam_shift);

	return label;
}

/**
@brief 	Returns an asymmetric unit reference point.
@return View2<double>		reference view.
**/
View2<double>	Bsymmetry::reference_symmetry_view()
{
	View2<double>	ref;
 
	if ( pnt < 200 ) {
		ref[0] = 1;
		ref[2] = 0;
	} else if ( pnt < 300 ) {
		ref[0] = sin(M_PI/4);
		ref[2] = cos(M_PI/4);
	} else if ( pnt == 320 ) {
//		ref[0] = sin(M_PI/4);
//		ref[2] = cos(M_PI/4);
		ref[0] = ref[1] = cos(M_PI/8.0);
		ref[2] = sin(M_PI/8.0);
	} else if ( pnt == 432 ) {
//		ref[0] = sin(M_PI/8);
//		ref[2] = cos(M_PI/8);
		ref[0] = cos(M_PI/3.0);
		ref[2] = sin(M_PI/3.0);
	} else if ( pnt == 532 ) {
//		ref[0] = sin(M_PI/18);
//		ref[2] = cos(M_PI/18);
		if ( lbl.find("I90") ) {
			ref[1] = 1;
			ref[2] = GOLDEN + 1;
		} else {
			ref[0] = 1;
			ref[2] = GOLDEN + 1;
		}
	} else if ( pnt == 600 ) {	// What is it for helical symmetry?
	}

	ref.normalize();

	return ref;
}

int 		Bsymmetry::asymmetric_unit_index(View2<double> theview)
{
	int				i(0), n(0);
	double			a, amin(TWOPI);

	vector<View2<double>>	rv = reference_asymmetric_unit_views();
	
	for ( i=0; i<rv.size(); ++i ) {
		a = theview.angle(rv[i]);
		if ( amin > a ) {
			amin = a;
			n = i;
		}
	}
		
	return n;
}

/**
@brief 	Rotation matrix to orient a symmetry axis on the z axis.
@param 	axis		desired symmetry axis order.
@param 	axis_flag	view modifier.
@return Matrix3		new rotation matrix.
**/
Matrix3		Bsymmetry::rotate_to_axis(long axis, long axis_flag)
{
	double			sqrt2(sqrt(2));
	Matrix3			mat(1), mat2(1);

	if ( axis_flag ) mat2 = Matrix3(Vector3<double>(0,0,1), M_PI_2);
	
	if ( pnt < 200 ) {							// Cyclic
	} else if ( pnt < 300 ) {					// Dihedral
		if ( axis_flag ) mat2 = Matrix3(Vector3<double>(0,0,1), M_PI/op[0].order());
		if ( axis > 2 ) mat = Matrix3(Vector3<double>(0,1,0), M_PI_2);
	} else if ( pnt < 400 ) {					// Tetrahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(-1/sqrt2,1/sqrt2,0), atan(sqrt2));
	} else if ( pnt < 500 ) {					// Octahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(-1/sqrt2,1/sqrt2,0), atan(sqrt2));
		if ( axis == 2 ) mat = Matrix3(Vector3<double>(0,1,0), M_PI/4);
	} else {											// Icosahedral
		if ( axis == 3 ) mat = Matrix3(Vector3<double>(0,1,0), atan(1/(GOLDEN*sqrt(3))));
		if ( axis == 5 ) mat = Matrix3(Vector3<double>(1,0,0), atan(1/GOLDEN));
	}
		
	mat = mat2*mat;

	return mat;
}

/**
@brief 	Get all symmetry axes.
@return vector<Vector3<double>>	array of axes.
**/
vector<Vector3<double>>	Bsymmetry::get_axes()
{
	long			i, j, k, m, ns;
	Vector3<double>	axis1;
	Matrix3			mat;

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bsymmetry::get_axes: label=" << lbl << endl;
	
	if ( pnt < 300 ) {					// Cyclic and dihedral
		axis1[2] = 1;
	} else if ( pnt < 400 ) {					// Tetrahedral
		axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
	} else if ( pnt < 500 ) {					// Octahedral
		if ( lbl.find("-4") != string::npos ) {
			axis1[2] = 1;
		} else if ( lbl.find("-3") != string::npos ) {
			axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
		} else {
			cerr << "Error: The symmetry designation must be either O-3 or O-4!" << endl;
			bexit(-1);
		}
	} else if ( pnt < 600 ) {						// Icosahedral
		if ( lbl.find("-5") != string::npos ) {
			axis1[1] = 1/sqrt(2+GOLDEN);
			axis1[2] = 1/sqrt(3-GOLDEN);
		} else if ( lbl.find("-3") != string::npos ) {
			axis1[0] = axis1[1] = axis1[2] = 1/sqrt(3.0);
		} else {
			cerr << "Error: The symmetry designation must be either I-3 or I-5!" << endl;
			bexit(-1);
		}
	}
	
	vector<Vector3<double>>	axis;
	
	axis.push_back(axis1);

	for ( i=0, k=1, ns=1; i<op.size(); i++ ) {
		for ( j=1; j<op[i].order(); j++ ) {
			mat = Matrix3(op[i].axis(), j*TWOPI/op[i].order());
			for ( m=0; m<ns; m++, k++ ) axis.push_back(mat * axis[m]);
		}
		ns = k;
	}
	
	if ( verbose & VERB_FULL )
		for ( i=0; i<axis.size(); i++ )
			cout << "Axis[" << i+1 << "]: " << axis[i] << endl;
	
	return axis;
}

/**
@brief 	Get all symmetry-related views of one given view.
@param 	asu_view			asymmetric unit vector and rotation angle.
@return vector<View2<double>>	list of views.

	The number of views generated for a point group symmetry is
	calculated as the product of the order fields in the symmetry
	structure.

**/
vector<View2<double>>	Bsymmetry::get_all_views(View2<double> asu_view)
{
	int				i, j, k;
	double			angle;
	Quaternion		q, qv;
	
	asu_view.normalize();
	
	vector<View2<double>>	view;
	view.push_back(asu_view);
	
	if ( verbose & VERB_FULL ) {
		cout << "Getting all the symmetric views:" << endl;
		cout << "Symmetry:                       " << lbl << endl;
		cout << "View:                           " << asu_view << endl;
	}
	
	int				nview(1);
	for ( i=0; i<op.size(); i++ ) {
		for ( j=1; j<op[i].order(); j++ ) {
			angle = j*TWOPI*1.0L/op[i].order();
			q = Quaternion(op[i].axis(), angle);
			for ( k=0; k<nview; ++k ) {
				qv = view[k].quaternion();
				qv = q * qv;
				view.push_back(qv);
				if ( verbose & VERB_DEBUG )
					cout << "DEBUG symmetry_get_all_views: " << view.back() << endl;
			}
		}
		nview *= op[i].order();
	}
	
	if ( verbose & VERB_FULL )
		for ( i=0; i<view.size(); ++i )
			cout << "View " << i+1 << ": " << tab << view[i] << endl;
	
	return view;
}

/**
@brief 	Initializes a well-distributed set of views in an asymmetric unit.
@param 	theta_step			angular step size from primary symmetry axis (radians).
@param 	phi_step			angular step size around primary symmetry axis (radians).
@param 	flag				flag: 0=half, 1=full, 2=no in-plane.
@return vector<View2<double>>	list of views.

	A set of views is calculated with tesselation within each asymmetric
	unit such that the views are well-distributed.
 	Flag bits:
		1: both halves of the asymmetric unit are covered.
 		2: no in-plane rotations are applied.

**/
vector<View2<double>>	Bsymmetry::asymmetric_unit_views(double theta_step, double phi_step, int flag)
{
	if ( theta_step < M_PI/1800 ) theta_step = M_PI/180;
	if ( phi_step < M_PI/1800 ) phi_step = M_PI/180;

	char			full(flag & 1);
	char			inplanerot = !(flag & 2);
	int 			i, j, n(0), ntheta, nphi, nrphi, istart;
	double			theta, phi, max_theta, theta_start, phi_start, phi_end, ang(0);
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG asymmetric_unit_views: flag=" << flag << " full=" << full << endl;

	max_theta = M_PI_2 + 1e-6;			// Most symmetries limited to upper half
	ntheta = (int) (max_theta/theta_step + 0.5);
	nphi = (int) (M_PI*2.0/(op[0].order()*phi_step) + 0.5);
	
	if ( pnt < 200 ) {
		if ( full ) {
			max_theta = M_PI + 1e-6;
			ntheta += ntheta;
		}
	} else if ( pnt == 432 ) {
		max_theta = M_PI_4 + 1e-6;
		ntheta = (int) (max_theta/theta_step + 0.5);
	} else if ( pnt == 532 ) {
		max_theta = atan(1/(GOLDEN+1)) + 1e-6;
		ntheta = (int) (max_theta/theta_step + 2.5);
		nphi = (int) (atan(1/GOLDEN)/phi_step + 1.5);
	} else if ( pnt == 600 ) {
		ntheta = 1;
		nphi = (int) (M_PI*2.0/phi_step + 1);
	}
	
	if ( ntheta < 1 ) ntheta = 1;
	if ( nphi < 1 ) nphi = 1;
		
	if ( verbose & (VERB_PROCESS | VERB_LABEL) ) {
		if ( full )
			cout << "Getting all the asymmetric unit views for symmetry " << lbl << endl;
		else
			cout << "Getting half the asymmetric unit views for symmetry " << lbl << endl;
		cout << "Theta and phi step sizes:       " <<
				theta_step*180/M_PI << " " << phi_step*180/M_PI << " degrees" << endl;
		cout << "Theta and phi steps:            " << ntheta << " " << nphi << endl;
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG asymmetric_unit_views: Views allocated = " << 2*ntheta*nphi << endl;
	
	vector<View2<double>>	views;
	
	if ( pnt < 500 ) {						// Top view for all but icosahedral and helical symmetry
		views.push_back(View2<double>(0,0,1,0));
		n = 1;
		if ( verbose & VERB_DEBUG )
			cout << "DEBUG asymmetric_unit_views: Adding the top view" << endl;
	}
	
	if ( pnt > 100 && pnt < 200 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = istart = 0;
			if ( nrphi ) {
				phi = M_PI*1.0/(op[0].order()*nrphi);
				istart = -nrphi;
				if ( pnt < 102 ) nrphi -= 1;
			}
			for ( j=istart; j<=nrphi; j++ ) {
				if ( inplanerot ) ang = j*phi;
				views.push_back(View2<float>(sin(theta)*cos(j*phi), sin(theta)*sin(j*phi), cos(theta), ang));
				n++;
			}
		}
		if ( full && views.back()[2] > -0.999999 ) {
			views.push_back(View2<double>(0,0,-1,0));
			n++;
		}
	} else if ( pnt > 200 && pnt < 300 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			nrphi = (int) (nphi*sin(theta)/2 + 0.5);	// Number of views at this radius and theta
			phi = 0;
			if ( nrphi ) phi = M_PI*1.0/(op[0].order()*nrphi);
			istart = 0;
			if ( full ) istart = -nrphi;
			for ( j=istart; j<=nrphi; j++ ) {
				if ( inplanerot ) ang = j*phi;
				views.push_back(View2<double>(sin(theta)*cos(j*phi), sin(theta)*sin(j*phi), cos(theta), ang));
				n++;
			}
		}
	} else if ( pnt == 320 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			phi_end = theta;
			if ( theta > M_PI_4 ) phi_end = M_PI_2 - theta;
			if ( full ) phi_start = -theta;
			if ( full && theta > M_PI_4 ) phi_start = theta - M_PI_2;
			for ( phi=phi_start; phi<=phi_end; phi+=phi_step ) {
				views.push_back(View2<double>(sin(theta), sin(phi), cos(theta), 0));
				n++;
			}
		}
	} else if ( pnt == 432 ) {
		for ( theta=theta_step; theta<=max_theta; theta+=theta_step ) {
			phi_start = 0;
			if ( full ) phi_start = -theta;
			for ( phi=phi_start; phi<=theta; phi+=phi_step ) {
				views.push_back(View2<double>(sin(theta), sin(phi), cos(theta), 0));
				n++;
			}
		}
	} else if ( pnt == 532 ) {
		if ( lbl.find("I90") != string::npos ) {
			theta_start = 0;
			if ( full ) theta_start = -ntheta*theta_step;
			for ( theta=theta_start; theta<=max_theta; theta+=theta_step ) {
				for ( phi=0; phi<=GOLDEN*(max_theta - fabs(theta)); phi+=phi_step ) {
					views.push_back(View2<double>(tan(phi), tan(theta), cos(theta), 0));
					n++;
				}
			}
		} else {
			for ( theta=0; theta<=max_theta; theta+=theta_step ) {
				nrphi = (int) (GOLDEN*(max_theta - theta)/phi_step);
				phi_start = 0;
				if ( full ) phi_start = -nrphi*phi_step;
				for ( phi=phi_start; phi<=GOLDEN*(max_theta - theta); phi+=phi_step ) {
					views.push_back(View2<double>(tan(theta), tan(phi), cos(theta), 0));
					n++;
				}
			}
		}
	} else if ( pnt == 600 ) {					// Helical symmetry
		for ( phi = 0; phi < TWOPI-0.001; phi += phi_step ) {
			if ( inplanerot ) ang = phi;
			views.push_back(View2<double>(cos(phi), sin(phi), 0, ang));
			n++;
		}
	} else {
		cerr << "Warning: Symmetry type " << pnt << " not supported!" << endl;
	}

	if ( verbose & ( VERB_PROCESS | VERB_LABEL ) )
		cout << "Views generated:                " << n << endl << endl;

	for ( auto& v: views ) v.normalize();
	
	if ( verbose & VERB_FULL ) {
		cout << "View\tx\ty\tz\ta" << endl;
		for ( i=0; i<views.size(); ++i )
			cout << i+1 << tab << views[i] << endl;
		cout << endl;
	}
	
	return views;
}

vector<View2<double>>	Bsymmetry::asymmetric_unit_views(double theta_step, double phi_step, double alpha_step, int flag)
{
	vector<View2<double>>	views = asymmetric_unit_views(theta_step, phi_step, flag);
	
	vector<View2<double>>	views2 = view_list_expand_angles(views, -M_PI, M_PI - alpha_step/2, alpha_step);
	
	return views2;
}

/**
@brief 	Initializes a set of views around the z-axis for helical projection.
@param 	side_ang			starting angle (radians).
@param 	theta_step			angular step size perpendicular to equator (radians).
@param 	phi_step			angular step size around equator (radians).
@return vector<View2<double>>	list of views.

	A set of views is calculated corresponding to views around the z-axis
	including some tilting to account for oblique views.

**/
vector<View2<double>>	Bsymmetry::side_views(double side_ang, double theta_step, double phi_step)
{
	if ( side_ang < 0 ) side_ang = -side_ang;
	if ( theta_step <= 0 ) {
		side_ang = 0;
		theta_step = M_PI/180;
	} else {
		side_ang = theta_step*floor(side_ang/theta_step);
	}
	if ( phi_step <= 0 ) phi_step = M_PI_2;
	
	int				order = op[0].order();
	if ( order < 1 ) order = 1;
	
	double			phi, theta;
	int				n;
	Euler			euler;
	vector<View2<double>>	views;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Getting all the asymmetric unit side views for symmetry " << lbl << endl;
		cout << "Step size around equator:       " << phi_step*180.0/M_PI << " degrees" << endl;
		if ( side_ang ) {
			cout << "Deviation from equator:         " << side_ang*180.0/M_PI << " degrees" << endl;
			cout << "Step size from equator:         " << theta_step*180.0/M_PI << " degrees" << endl;
		}
	}
	
	for ( n = 0, phi = -M_PI*1.0L/order; phi < M_PI*1.0L/order; phi += phi_step ) {
		for ( theta = M_PI_2 - side_ang; theta <= M_PI_2 + side_ang; theta += theta_step ) {
			euler = Euler(0, theta, phi);
			views.push_back(euler.view());
			n++;
//			cout << "theta=" << theta*180.0/M_PI << " phi=" << phi*180.0/M_PI << endl;
		}
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Number of views:                " << n << endl << endl;
	
	return views;
}

vector<View2<double>>	Bsymmetry::side_views(double side_ang, double theta_step, double phi_step, double alpha_step)
{
	vector<View2<double>>	views = side_views(side_ang, theta_step, phi_step);
	
	vector<View2<double>>	views2 = view_list_expand_angles(views, -M_PI, M_PI - alpha_step/2, alpha_step);
	
	return views2;
}

/**
@brief 	Finds the corresponding view in the asymmetric unit.
@param 	theview			view.
@return View2<double>	the asymmetric unit view.

	The asymmetric unit view is found and the the new view with the 
	link from the old view is returned.

**/
View2<double>	Bsymmetry::find_asymmetric_unit_view(View2<double> theview)
{
	if ( pnt < 102 ) return theview;
	if ( pnt >= 600 ) return theview;
	
	View2<double>	v, bv(theview), tv;
	vector<View2<double>>	views = get_all_views(theview);
	
	int				found(0);
	double			tol(1e-10);
	double			lim = tan(M_PI*1.0L/op[0].order());
	
	for ( auto& v: views ) {
		if ( found ) break;
		if ( ( v[0] - tol >= 0 ) || ( v[0] + tol >= 0 && v[1] >= 0 ) ) {
			tv = View2<double>(0,0,0,0);
			if ( pnt < 200 ) {
				if ( op[0].order() == 2 ) {
					tv = v;
				} else {
					if ( fabs(v[1]) <= v[0]*lim + tol ) tv = v;
				}
			} else if ( pnt > 200 && pnt < 300 ) {
				if ( v[2] + tol >= 0 ) {
					if ( op[0].order() == 2 ) {
						tv = v;
					} else {
						if ( fabs(v[1]) <= v[0]*lim + tol ) tv = v;
					}
				}
			} else if ( pnt == 320 ) {
				if ( fabs(v[1]) <= v[0] + tol && fabs(v[1]) <= v[2] + tol )
					tv = v;
			} else if ( pnt == 432 ) {
				if ( fabs(v[1]) <= v[0] + tol && v[0] <= v[2] + tol )
					tv = v;
			} else if ( pnt == 532 ) {
				if ( lbl.find("I90") != string::npos ) {
					if ( v[2] + tol >= 0 && fabs(v[1]) <= (v[2]/GOLDEN - v[0])/GOLDEN + tol )
						tv = v;
				} else {
					if ( v[2] + tol >= 0 && fabs(v[1]) <= v[2]/GOLDEN - v[0]*GOLDEN + tol )
						tv = v;
				}
			} else if ( pnt >= 600 ) {	// What is it for helical symmetry?
				tv = theview;
			}
			if ( tv.vector_size() > 0.9 ) {
				bv = tv;
				found++;
			}
		}
	}
	
	if ( !found ) {
		cerr << "ASU view not found: " << theview << endl;
	} else {
		if ( verbose & VERB_STATS ) {
			cout << "Old view:\t" << theview << endl;
			cout << "New view:\t" << bv << endl;
		}
	}

	return bv;
}

/**
@brief 	Finds the closest symmetric match between two views.
@param 	view_ref		reference view.
@param 	view			test view.
@return View2<double>	matched symmetric version of test view.

	A list of symmetry-related views of the test view is searched
	for the closest to the reference view.
	The matched view is returned.

**/
View2<double>	Bsymmetry::find_closest_symmetric_view(View2<double> view_ref, View2<double> view)
{
	if ( pnt < 102 ) return view;
	if ( pnt >= 600 ) return view;
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Bsymmetry::find_closest_symmetric_view: " << lbl << endl;
	
	double					r, rb(1e30);
	View2<double>			bv(view);
	vector<View2<double>>	viewsym = get_all_views(view);
	
	for ( auto& v: viewsym ) {
//		r = v.residual(view_ref);
		r = v.angle(view_ref);
		if ( rb > r ) {
			rb = r;
			bv = v;
		}
	}
	
	return bv;
}

/**
@brief 	Show symmetry matrices.
@return int 			number of symmetry matrices.
**/
int			Bsymmetry::show_matrices()
{
	vector<Matrix3>	mat = matrices();
	int				n(mat.size());

	cout << "\nSymmetry matrices (" << n << "):" << endl;
	for ( int i=0; i<n; ++i ) {
		cout << "Matrix " << i+1 << ":" << endl;
		matrix3_show_hp(mat[i]);
	}
	cout << endl;
	
	return n;
}

/**
@brief 	Show symmetry matrices associated with each symmetry operator.
@return int 			number of symmetry matrices.
**/
int			Bsymmetry::show_operational_matrices()
{
	int			i;
	Matrix3		mat;

	for ( i=0; i<op.size(); ++i ) {
		mat = op[i].matrix();
		cout << "Operation " << i+1 << ":" << endl;
		cout << "Axis and angle: " << op[i].axis() << " " << 360.0/op[i].order() << endl;
		cout << mat << endl;
	}
	cout << endl;
	
	return op.size();
}

/**
@brief 	Show PDB format symmetry matrices associated with each symmetry operator.
@return int 			number of symmetry matrices.
**/
int			Bsymmetry::show_pdb_matrices()
{
	int			i, j, k;
	double		t(0);

	vector<Matrix3>	mat = matrices();
	long			nsym = mat.size();

	for ( i=0; i<nsym; i++ ) {
		for ( j=0; j<3; j++ ) {
			cout << "BIOMT" << j+1 << setw(4) << i+1 << setprecision(6);
			for ( k=0; k<3; k++ )
					cout << setw(10) << mat[i][j][k];
			cout << setprecision(6) << setw(15) << t << endl;
		}
	}
	
	return nsym;
}


