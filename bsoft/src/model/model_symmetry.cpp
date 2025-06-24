/**
@file	model_symmetry.cpp
@brief	Library routines used for model symmetry operations
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250618
**/

#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_compare.h"
#include "model_mol.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Set model component locations within the asymmetric unit.
@param 	*model				model parameters.
@param 	&sym				point group symmetry.
@return long					number of components (<0 means failure).

	Only the first model is processed.

**/
long		model_find_asymmetric_unit(Bmodel* model, Bsymmetry& sym)
{
	if ( sym.point() < 102 ) return -1;

	if ( verbose )
		cout << "Setting asymmetric unit for " << model->identifier() << " to symmetry " << sym.label() << endl;

	long			ncomp(0);
	Bcomponent*		comp;
	View2<double>	view, view_asu;

	for ( comp = model->comp; comp; comp = comp->next, ncomp++ ) {
		view = View2<double>(comp->location()[0], comp->location()[1], comp->location()[2], 0);
		view_asu = sym.find_asymmetric_unit_view(view);
		comp->view(View2<float>(view_asu[0],view_asu[1],view_asu[2],view_asu[3]));
		comp->location(view_asu.vector3() * comp->location().length());
	}
	
	return ncomp;
}

long		models_find_asymmetric_unit(Bmodel* model, Bsymmetry& sym)
{
	if ( sym.point() < 102 ) return -1;

	long			ncomp(0);
	Bmodel*			mp;

	for ( mp = model; mp; mp = mp->next )
		ncomp += model_find_asymmetric_unit(model, sym);
	
	return ncomp;
}


/**
@brief 	Applying symmetry to model components.
@param 	*model				model parameters.
@param 	symmetry_string		symmetry code.
@param 	origin				transformation origin.
@param 	ref_view			reference view.
@param 	flags				1=find asu.
@return long					number of components (<0 means failure).

	Only the first model is processed.

**/
/*long		model_apply_point_group(Bmodel* model, Bsymmetry& sym,
					Vector3<double> origin, View2<double> ref_view, int flags)
{
	if ( ! model->comp ) return 0;
	
	Bsymmetry		sym(symmetry_string);
	
	long 			i, j, k;
	long 			ncomp(0), nlink(0), id(0);
	double			distance(1e30), d;
	
	int 			nunits(sym.order());
	if ( nunits < 2 ) return 0;

	if ( ref_view.vector_size() < 1e-10 ) ref_view = View2<double>(0, 0, 1, 0);

	Bcomponent*		comp;
	Bcomponent*		nu_comp = NULL;

	for ( ncomp=1, comp = model->comp; comp->next; comp = comp->next, ncomp++ ) {
		for ( nu_comp = comp->next; nu_comp; nu_comp = nu_comp->next ) {
			d = comp->location().distance(nu_comp->location());
			if ( d < distance ) distance = d;
		}
	}
	if ( ncomp < 2 ) distance = model->comp->location().length();
	distance /= 2;
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Applying symmetry " << sym.label() << ":" << endl;
		cout << "Number of asymmetric units:     " << nunits << endl;
		cout << "Origin for symmetrization:      " << origin << endl;
		cout << "Reference symmetry axis:        " << ref_view << endl;
		cout << "Reference rotation angle:       " << ref_view.angle()*180/M_PI << endl;
		cout << "Overlap distance cutoff:        " << distance << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Applying symmetry " << sym.label() << endl << endl;
	
	model->symmetry(sym.label());
	
	Matrix3			ref_mat = ref_view.matrix();
	Matrix3			mat(1), cmat;
	Vector3<double>	new_axis;
	
	double			clen;
	Vector3<double>	v;
	Blink*			link;
	Blink*			nu_link = NULL;
	Blink*			link_start = NULL;
	
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		if ( id < stol(comp->identifier()) ) id = stol(comp->identifier());
		if ( flags & 1 ) {
			clen = comp->location().length();
//			comp->view(View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0));
//			comp->view(find_asymmetric_unit_view(sym, comp->view()));
			View2<double>	tv(comp->location());
//			tv = find_asymmetric_unit_view(sym, tv);
			tv = sym.find_asymmetric_unit_view(tv);
			comp->view(View2<float>(tv[0], tv[1], tv[2], tv[3]));
			comp->location(comp->view().vector3());
			comp->scale(clen);
		}
	}
	
	for ( nlink = 0, link = model->link; link; link = link->next ) nlink++;
	
	for ( i=0; i<sym.operations(); i++ ) {
		new_axis = ref_mat * sym[i].axis();
		for ( j=1; j<sym[i].order(); j++ ) {
			mat = Matrix3(new_axis, j*TWOPI/sym[i].order());
			mat *= ref_mat;
			link_start = NULL;
			for ( k=0, link = model->link; k<nlink; link = link->next, k++ ) {
//				nu_link = link_add(&link, link->comp[0], link->comp[1], link->length, link->radius);
				nu_link = link->add(link);
				if ( !link_start ) link_start = nu_link;
//				nu_link->comp[0] = link->comp[0];
//				nu_link->comp[1] = link->comp[1];
//				nu_link->length(link->length());
//				nu_link->radius(link->radius());
//				nu_link->color(link->color());
			}
			for ( k=0, comp = model->comp; k<ncomp; comp = comp->next, k++ ) {
//				nu_comp = component_add(&comp, id);
//				component_copy(comp, nu_comp);
				nu_comp = comp->add(comp);
				nu_comp->identifier() = to_string(++id);
				v = comp->location() - origin;
				nu_comp->location((mat * v) + origin);
				cmat = comp->view().matrix();
				cmat = mat * cmat;
				nu_comp->view(View2<float>(cmat));
				for ( link = link_start; link; link = link->next ) {
					if ( link->comp[0] == comp ) link->comp[0] = nu_comp;
					if ( link->comp[1] == comp ) link->comp[1] = nu_comp;
				}
			}
		}
		ncomp *= sym[i].order();
		nlink *= sym[i].order();
	}
	
//	ncomp = model_average_overlapped_components(model, distance);

	return ncomp;
}
*/

/**
@brief 	Applying symmetry to model components.
@param 	*model				model parameters.
@param 	&sym				point group symmetry.
@param 	origin				transformation origin.
@param 	ref_view			reference view.
@param 	flags				1=find asu.
@return long					error code (<0 means failure).

	All models in the list are processed.

**/
long		models_apply_point_group(Bmodel* model, Bsymmetry& sym,
					Vector3<double> origin, View2<double> ref_view, int flags)
{
	long 			i, j, k, m, n, nmod, nmods, ncomp(model->component_count());
	Bmodel*			mp;
	Bmodel*			nu_mp;
	Bcomponent*		comp;

	if ( ref_view.vector_size() < 1e-10 ) ref_view = View2<double>(0, 0, 1, 0);

	// Set up the ids
	vector<string>	ids;
	for ( nmod = 0, mp = model; mp; mp = mp->next, ++nmod )
		ids.push_back(mp->identifier());

	if ( verbose & VERB_PROCESS ) {
		cout << "Applying symmetry " << sym.label() << ":" << endl;
		cout << "Number of models:               " << nmod << endl;
		cout << "Number of asymmetric units:     " << sym.order() << endl;
		cout << "Origin for symmetrization:      " << origin << endl;
		cout << "Reference symmetry axis:        " << ref_view << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Applying symmetry " << sym.label() << endl << endl;
	
	model->symmetry(sym.label());
	
	Matrix3			ref_mat = ref_view.matrix();
	Matrix3			mat(1), cmat;
	Vector3<double>	new_axis;
	Vector3<double>	v;
	string			id;

	for ( i=0, n=1, nmods = nmod; i<sym.operations(); ++i ) {
		if ( verbose )
			cout << "Symmetry operation " << i+1 << tab << sym[i].order() << endl;
		new_axis = ref_mat * sym[i].axis();
		for ( j=1; j<sym[i].order(); ++j ) {
			mat = Matrix3(new_axis, j*TWOPI/sym[i].order());
			for ( k=0, mp = model; mp && k<nmods; mp = mp->next, ++k ) {
				m = k%nmod;
				if ( m == 0 ) n++;
				id = ids[m] + to_string(n);
//				cout << j << tab << mp->identifier() << tab << id << endl;
				nu_mp = mp->copy(id);
				mp->add(nu_mp);
				for ( comp = nu_mp->comp; comp; comp = comp->next ) {
					v = comp->location() - origin;
					comp->location((mat * v) + origin);
					cmat = comp->view().matrix();
					cmat = mat * cmat;
					comp->view(View2<float>(cmat));
					ncomp++;
				}
			}
		}
		nmods = model->count();
	}
	
	if ( verbose )
		cout << "Number of models:               " << nmods << endl << endl;
		
	return ncomp;
}


/**
@brief 	Symmetrize a model.
@param 	*model				model parameters.
@param 	&sym				point group symmetry.
@return long					error code (<0 means failure).

	For each component, a new location is calculated from the average location
	of the closest symmetry-related components.

**/
/*int 		model_symmetrize2(Bmodel* model, Bstring& symmetry_string)
{
	Bsymmetry		sym(symmetry_string);
	
	int 			i, nunits(1);
	for ( i=0; i<sym.erations(); i++ ) nunits *= sym.[i].order;
	if ( nunits < 2 ) return 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Symmetrizing " << symmetry_string << ":" << endl;
		cout << "Number of asymmetric units:     " << nunits << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Symmetrizing " << symmetry_string << endl << endl;
	
	long 			ncomp(0);
	double			clen, d, dmin, R(0);
	Vector3<double>	loc, locsym;
	View*			views, *v, comp_view;

	Bcomponent*		comp;
	Bcomponent*		compsym;
	Bcomponent*		compsel = NULL;
	
	
	for ( comp = model->comp; comp; comp = comp->next ) {
		comp->force(0);
		comp->select(0);
	}
	for ( comp = model->comp; comp; comp = comp->next ) {
		clen = comp->location().length();
		comp_view = View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0);
		views = symmetry_get_all_views(sym, comp_view);
		for ( v = views; v; v = v->next ) {
			loc = Vector3<double>(v->x(), v->y(), v->z());
			loc *= clen;
			for ( dmin = 1e30, compsym = model->comp; compsym; compsym = compsym.next ) {
				d = loc.distance(compsym.loc);
				if ( dmin > d ) {
					dmin = d;
					compsel = compsym;
				}
			}
			compsel->vec += loc;
			compsel->sel++;
			R += dmin*dmin;
			ncomp++;
		}
	}
	for ( comp = model->comp; comp; comp = comp->next ) {
		comp->location(comp->vec/comp->select();
		comp->force(0);
	}
	
	R = sqrt(R/ncomp);

	if ( verbose & VERB_PROCESS )
		cout << "Symmetry deviation:             " << R << " A" << endl << endl;
		
	return ncomp;
}
*/
long		models_symmetrize(Bmodel* model, Bsymmetry& sym)
{
	int 			nunits(sym.order());
	if ( nunits < 2 ) return 0;
	
	if ( verbose & VERB_PROCESS ) {
		cout << endl << "Symmetrizing " << sym.label() << ":" << endl;
		cout << "Number of asymmetric units:     " << nunits << endl << endl;
	} else if ( verbose & VERB_LABEL )
		cout << endl << "Symmetrizing " << sym.label() << endl << endl;
	
	long 			ncomp(0), nasu;
	double			tol(1e-4), d, dmin, R(0);
	Vector3<double>	loc;
	View2<double>	view, view_asu;

	Bmodel*			mp;
	Bcomponent*		comp;
	Bcomponent*		compasu;
	Bcomponent*		compsel = NULL;
	
	// Select all components in the ASU
	for ( nasu=0, mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) {
			comp->select(0);
//			view = View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0);
//			view_asu = find_asymmetric_unit_view(sym, view);
			view = View2<double>(comp->location());
			view_asu = sym.find_asymmetric_unit_view(view);
			comp->force(Vector3<double>(view_asu[0], view_asu[1], view_asu[2]));
			if ( view.distance(view_asu) < tol ) {
				comp->location(comp->force() * comp->location().length());
				comp->select(1);
				nasu++;
			}
		}
	}
	
	for ( mp = model; mp; mp = mp->next ) {
		for ( comp = mp->comp; comp; comp = comp->next ) if ( !comp->select() ) {
			compsel = NULL;
			for ( dmin = 1e30, compasu = model->comp; compasu; compasu = compasu->next ) if ( compasu->select() ) {
				d = comp->force().distance(compasu->force());
				if ( dmin > d ) {
					dmin = d;
					compsel = compasu;
				}
			}
			if ( compsel ) {
				compsel->shift(comp->force() * comp->location().length());
				compsel->select_increment();
				dmin *= comp->location().length();
				R += dmin*dmin;
				ncomp++;
			}
		}
		if ( verbose & VERB_PROCESS ) {
			cout << "Model:                          " << mp->identifier() << endl;
			cout << "Asymmetric unit components:     " << nasu << endl;
			cout << "Comp\tCount" << endl;
		}
		
		for ( comp = mp->comp; comp; comp = comp->next ) if ( comp->select() ) {
			comp->scale(1.0/comp->select());
			comp->force(Vector3<float>(0,0,0));
			if ( verbose & VERB_PROCESS )
				cout << comp->identifier() << tab << comp->select() << endl;
		}
		if ( verbose & VERB_PROCESS )
			cout << endl;
	}
	
	
	R = sqrt(R/ncomp);

	if ( verbose & VERB_PROCESS )
		cout << "Symmetry deviation:             " << R << " A" << endl << endl;
		
	models_delete_non_selected(&model);
	
	Vector3<double>		origin;
	View2<double>		ref_view;
	
	models_apply_point_group(model, sym, origin, ref_view, 1);	

	return ncomp;
}

struct MassCOM {
	double	mass;
	Vector3<double>	com;
} ;

static int  QsortMassCOM(const void *x, const void *y)
{
	MassCOM*		c1 = (MassCOM *) x;
	MassCOM*		c2 = (MassCOM *) y;
	
	if ( fabs((c1->mass - c2->mass)/(c1->mass + c2->mass)) > 0.1 ) {
		if ( c1->mass < c2->mass ) return 1;
		else return -1;
	} else {
		if ( c1->com.length() < c2->com.length() ) return 1;
		else return -1;		// Sort from high to low
	}
}

/**
@brief	Searches for the standard view based on point group symmetry.
@param 	*model			linked list of models.
@param 	&sym			point group symmetry.
@param 	ref_view		reference view.
@return Transform			0.

	The set of models is first analyzed to identify the different
	chains and calculate their centers-of-mass and weights.
	The overall center-of-mass defines a point on at least the
	major symmetry axis (cyclic symmetries), or the likely intersection
	of symmetry axes.
	Note: This function does a reasonable job of orienting the model set,
	but it may be off by up to an angstrom!!!

**/
Transform 	model_find_standard_view(Bmodel* model, Bsymmetry& sym, View2<double> ref_view)
{
	Transform		t;

	Vector3<double>	coc = models_center_of_coordinates(model);
	models_shift(model, -coc);
	
	int				g, i, j, k, nmol, nunits(sym.order()), nmassif;
	Bmodel*			mp;
	
	nmol = model->count();
	
	if ( nmol < nunits ) {
		cerr << "Error: Not enough models (" << nmol << ") to fit into " << nunits << " asymmetric units!" << endl << endl;
		return t;
	}

	// Get all chains' centers of mass and weight
	double			themass, themass2, maxmass(0);
	MassCOM*		c = new MassCOM[nmol];
	for ( i=0; i<nmol; i++ ) c[i].com = c[i].mass = 0;
	
	nmassif = 0;
	for ( i=0, mp = model; mp; mp = mp->next, i++ ) {
		c[i].com = mp->center_of_coordinates();
		c[i].mass = mp->component_count();
		if ( maxmass < c[i].mass ) maxmass = c[i].mass;
	}
	
	// Sort the chains based on weight and distance from COM
	qsort((void *) c, nmol, sizeof(MassCOM), (int (*)(const void *, const void *)) QsortMassCOM);
	
	if ( verbose & VERB_FULL ) {
		cout << endl << "Ordered list of centers-of-mass:" << endl;
		cout << "    x\t    y\t    z\tMass (Da)" << endl;
		for ( i=0; i<nmol; i++ )
			cout << c[i].com[0] << tab << c[i].com[1] << tab << c[i].com[2] << tab << c[i].mass << endl;
		cout << endl << endl;
	}
	
	// Determine the number of groups of chains, each group likely symmetry-related
	// A group is defined as a set of chains (number = symmetry order),
	// that differ little in weight (stdev < 5%)
	int				ngroup(0);
	Vector3<double>	com, translate;
	
	if ( verbose &  VERB_FULL )
		cout << endl << "Group\tMass\tStDev" << endl;
	for ( i=0; i<nmol-nunits+1; i+=nunits ) {   // Loop over all groups of size defined by symmetry
		themass = themass2 = 0;
		com = 0;
		for ( j=i; j<i+nunits; j++ ) {
			themass += c[j].mass;
			themass2 += c[j].mass*c[j].mass;
			com += c[j].com;	// Center of group on symmetry axis
		}
		themass /= nunits;
		themass2 = themass2/nunits - themass*themass;
		if ( themass2 > 0 ) themass2 = sqrt(themass2);
		else themass2 = 0;
		if ( themass > 100 && themass2/themass < 0.05 ) {
			ngroup++;
			com *= -1.0/nunits;	// Translation to symmetry axis
			translate += com;
		}
		if ( verbose &  VERB_FULL )
			cout << ngroup << tab << themass << tab << themass2 << endl;
	}
	translate /= ngroup;
	
	// Subtract the global center-of-mass from the chain COM's
	com = 0;
	for ( i=0; i<nmassif; i++ ) com += c[i].com;
	com /= nmassif;
	for ( i=0; i<nmassif; i++ ) c[i].com -= com;
	
	// Generate pairs of vectors and calculate normals
	int				n(0), nvec = (ngroup*nunits*(nunits-1)*(nunits-2))/6;
	double			a;
	Vector3<double>	axis, v1, v2, d;
	Vector3<double>*	v = new Vector3<double>[nvec];
	
	if ( verbose &  VERB_FULL )
		cout << endl << "Group\tv1\tv2\tv3\tx\ty\tz" << endl;
	for ( g=0; g<ngroup; g++ ) {
		// Calculate normal vectors
		for ( i=2+nunits*g; i<nunits*(g+1); i++ ) {
			for ( j=1+nunits*g; j<i; j++ ) {
				v1 = c[i].com - c[j].com;
				for ( k=nunits*g; k<j; k++ ) {
					v2 = c[i].com - c[k].com;
//					v[n] = v1.cross(v2) * c[i].mass;
					v[n] = v1.cross(v2);
					v[n].normalize();
					if ( verbose &  VERB_FULL )
						cout << g << tab << k << tab << j << tab << i << tab
							<< v[n][0] << tab << v[n][1] << tab << v[n][2] << endl;
					n++;
				}
			}
		}
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG molgroup_find_standard_view: nvec=" << nvec << " n=" << n << endl;
		
	// Find a common normal when comparing pairs of vectors and assign it as the major axis
	int				vbest = -1, m = 100;
	int				nx = 2*m, ny = 2*m, nz = 2;
	int				sv = nx * ny * nz;
	vector<int>		nv(sv,0);
	
	for ( i=1; i<n; i++ ) {
		for ( j=0; j<i; j++ ) {
			v1 = v[i];
			v1.normalize();
			a = fabs(v[i].angle(v[j]));
			if ( a > M_PI_2 ) {
				a = M_PI - a;
				v1 = -v1;
			}
			v1 *= m;
			if ( a < 0.05 ) {
				k = (int) (nx*(v1[1]+m)) + (int) (v1[0]+m);
				if ( k > nx*ny ) k = 0;
				if ( v1[2] > 0 ) k += nx*ny;
				nv[k]++;
				if ( vbest < nv[k] ) {
					vbest = nv[k];
					axis = v1;
				}
			}
		}
	}
	
	for ( i=0; i<n; i++ ) {
		a = fabs(axis.angle(v[i]));
		if ( a > M_PI_2 ) {
			a = M_PI - a;
			v[i] = -v[i];
		}
		if ( a < 0.05 ) axis += v[i];
	}
	
	axis.normalize();
	
	// Find the 2-fold on the x-axis as the sum vector of two COM's
	// that has a difference vector that is most parallel with the major axis
	double		ang(0), ang_min = 1000;
	if ( sym.point() > 200 ) {
		for ( i=1; i<nunits; i++ ) {
			for ( j=0; j<i; j++ ) {
				d = c[i].com - c[j].com;
				ang = axis.angle(d);
				if ( fabs(ang_min) > fabs(ang) ) {
					ang_min = ang;
					v2 = c[i].com + c[j].com;
					v2.normalize();
				}
			}
		}
	}
	
	delete[] c;
	delete[] v;
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Angle with major axis: " << ang_min*180/M_PI << endl;
		cout << "Dihedral axis:         " << v2 << endl;
		cout << "Angle between axes:    " << (axis.angle(v2))*180/M_PI << endl;
	}
	
	if ( verbose & VERB_LABEL ) {
		cout << "Translate:                      " << translate << endl;
		cout << "Major symmetry axis:            " << axis << endl;
		if ( sym.point() > 200 )
			cout << "Minor symmetry axis:            " << v2 << endl;
		cout << endl;
	}
	
	// The first matrix rotates it to the major symmetry axis
	// The second does the in-plane rotation
	View2<double>	view(-axis[0], -axis[1], axis[2], 0);
	Matrix3			mat = view.matrix();
	if ( sym.point() > 200 ) {
		v1 = 0;
		v1[0] = 1;
		v2 = mat * v2;
		ang = fmod(v1.angle(v2), M_PI*2.0/sym[0].order());
		view = View2<double>(0, 0, 1, ang);
		mat = view.matrix() * mat;
	}
	
	if ( verbose & VERB_DEBUG ) {
		cout << "Dihedral axis:                      " << v2 << endl;
		cout << "Angle:                              " << ang*180.0/M_PI << endl;
	}
	
	t = Transform(mat);
	t.trans = translate;
	t.origin = coc;
	
	return t;
}

/**
@brief	Searches for the standard view based on point group symmetry and transform.
@param 	*model 		set of linked models..
@param 	&sym		point group symmetry.
@param 	ref_view	reference view (default should be 0,0,1,0).
@return int 			0.

	Each pair of chains in the molecule groupis tested for sequence
	identity to find symmetry-related molecules. For each pair of matched
	molecules, the transformation to superimpose the one onto the other
	is determined and the symmetry axis and translation calculated.
	The collection of symmetry axes are clustered with a k-means
	algorithm and the predominant class assigned to the major
	symmetry axis. For dihedral point groups, a minor axis is also
	assigned (randomly at this time). The molecule group is then
	transformed to orient it with the major axis on {0,0,1} and
	the minor axis on {1,0,0}, and the symmetry center at {0,0,0}.
	Note: This function has not been extensively tested with all
	point groups!!!

**/
int 		model_orient_to_standard_view(Bmodel* model, Bsymmetry& sym,
				View2<double> ref_view)
{
	random_seed();
	
	if ( sym.point() < 102 ) return 0;
	
	long			i, j, k, m, nmol, nres_cut, cut_max(30);
	long			nunit(sym.order()), ngroup, nid(0), n2(0), nn(0), no;
	Vector3<double>	a2[1000], an[1000];
	Bmodel*			mp, *mp2, *symmp;
	
	// Find the highest symmetry order
	for ( i=m=0; i<sym.operations(); i++ )
		if ( m < sym[i].order() ) m = sym[i].order();
	
	double			a, refang(TWOPI/m);
	
	nmol = model->count();
	ngroup = nmol/nunit;
	
	if ( ngroup < 1 ) {
		cerr << "Error: Not enough models (" << nmol << ") to fit into " << nunit << " asymmetric units!" << endl << endl;
		return -1;
	}

	Vector3<double>	coc = models_center_of_coordinates(model);
	models_shift(model, -coc);

	if ( verbose & VERB_PROCESS ) {
		cout << "Orienting to the standard view for symmetry " << sym.label() << endl;
		cout << "Number of models:               " << nmol << endl;
		cout << "Number of groups:               " << ngroup << endl;
	}
	
	vector<int>		g(nmol,0);	// Group indices
	Transform		t;
	bool			notfound;
	
	// Calculate a matrix of rotations between pairs of models
	for ( k=0, i=1; k<nmol; k++, i++, i=(i>ngroup)?1:i ) g[k] = i;
	for ( i=0, mp = model; mp->next; mp = mp->next, i++ ) {
		nres_cut = mp->component_count()/2;
		if ( nres_cut > cut_max ) for ( j=i+1, mp2 = mp->next; mp2; mp2 = mp2->next, j++ ) {
			if ( nid > cut_max || nid >= nres_cut ) {
				g[j] = g[i];
				t = model_find_transform(mp, mp2);
				symmp = mp->copy();
				model_rotate(symmp, t);
				t.fom = model_compare_corresponding(mp2, symmp);
				delete symmp;
				angle_set_negPI_to_PI(t.angle);
				a = t.angle/refang;
				no = (long) (a + 0.5);
				// Check for 2-fold axis
				if ( fabs(t.angle - M_PI) < 0.1 ) {
					for ( k=0, notfound=1; k<n2 && notfound; k++ ) {
						if ( a2[k].angle(t.axis) < 0.1 ) {
							a2[k] += t.axis;
							notfound = 0;
						}
					}
					if ( notfound ) {
						a2[n2] = t.axis;
						n2++;
					}
				}
				// Check for n-fold axis
				if ( fabs(a - no) < 0.03 ) {
					for ( k=0, notfound=1; k<nn && notfound; k++ ) {
						if ( an[k].angle(t.axis) < 0.1 ) {
							an[k] += t.axis;
							notfound = 0;
						}
					}
					if ( notfound ) {
						an[nn] = t.axis;
						nn++;
					}
				}
				m++;
			}
		}
	}

	if ( verbose & VERB_PROCESS ) {
		cout << "Num\tMol\tGroup\tLength" << endl;
		for ( i=0, mp = model; mp; mp = mp->next, i++ )
			cout << i+1 << tab << mp->identifier() << tab << g[i] << tab << mp->component_count() << endl;
		cout << endl;
	}
	
	double			d, dm;
	Vector3<double>	axis1, axis2;
	
	if ( sym.point() == 202 || sym.point() == 320 || sym.point() == 532 ) {
		// Pick the longest 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			d = a2[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = a2[k];
			}
		}
		// Pick the longest perpendicular 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			a = axis1.angle(a2[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = a2[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = a2[k];
				}
			}
		}
	} else if ( sym.point() == 432 ) {
		// Pick the longest 4-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			d = an[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = an[k];
			}
		}
		// Pick the longest perpendicular 4-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			a = axis1.angle(an[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = an[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = an[k];
				}
			}
		}
	} else {
		// Pick the longest n-fold axis
		for ( k=0, dm=0; k<nn; k++ ) {
			d = an[k].length();
			if ( dm < d ) {
				dm = d;
				axis1 = an[k];
			}
		}
		// Pick the longest perpendicular 2-fold axis
		for ( k=0, dm=0; k<n2; k++ ) {
			a = axis1.angle(a2[k]);
			if ( fabs(a - M_PI_2) < 0.1 ) {
				d = a2[k].length();
				if ( dm < d ) {
					dm = d;
					axis2 = a2[k];
				}
			}
		}
	}
	
	if ( axis1.length() < 0.9 ) {
		cerr << "Error: Major axis not found!" << endl;
	}
	
	if ( sym.point() > 200 && axis2.length() < 0.9 ) {
		cerr << "Error: Minor axis not found!" << endl;
	}
	
	axis1.normalize();
	axis2.normalize();
		
	// Translation required is the negative of the COM
	//	Vector3<double>	translate = -(molgroup_center_of_mass(molgroup));
	
	if ( verbose & VERB_LABEL ) {
		if ( sym.point() < 200 ) {
			cout << "Symmetry axis:                  " << axis1 << endl;
		} else {
			cout << "Major symmetry axis:            " << axis1 << endl;
			cout << "Minor symmetry axis:            " << axis2 << endl;
		}
	}
	
	// The first matrix rotates it to the major symmetry axis
	// The second does the in-plane rotation
	Vector3<double>	v1(0,0,1), v2(1,0,0);
	Matrix3         mat = Matrix3(axis1, v1);
	if ( sym.point() > 200 ) {
		v1 = mat * axis2;
		mat = Matrix3(v1, v2) * mat;
	}
	
	models_rotate(model, mat);
	
	return 0;
}

/*
vector<Vector3<double>>	model_symmetry_axes(Bmodel* model)
{
	long					i, j, nmod(model->count());
	double					angle;
	vector<string>			ids;
	vector<Vector3<double>>	axes;
	vector<double>			angles;
	vector<Vector3<double>>	com(nmod);
	Vector3<double>			axis, comall;
	Bmodel*					mp;
	
	if ( verbose )
		cout << "Calculating axes from model centers:" << endl;

	for ( i=0, mp = model; mp; mp = mp->next, ++i ) {
		com[i] = mp->center_of_coordinates();
		comall += com[i];
		ids.push_back(mp->identifier());
	}
	comall /= nmod;
	for ( i=0; i<nmod; ++i ) com[i] -= comall;
	
	if ( verbose )
		cout << "#\t#\tAxis\tAngle" << endl;
	for ( i=1; i<nmod; ++i ) {
		for ( j=0; j<i; ++j ) {
			angle = com[i].angle(com[j]);
			axis = com[i].cross(com[j]);
			axis.normalize();
			if ( axis[2] < 0 ) {
				axis = -axis;
				angle = -angle;
			}
			axes.push_back(axis);
			angles.push_back(angle);
			if ( verbose )
				cout << i << tab << j << tab << axes.back() << tab << angles.back()*180.0/M_PI << endl;
		}
	}

	return axes;
}
*/
vector<Transform>	model_symmetry_axes(Bmodel* model)
{
	Transform				t;
	vector<Transform>		tv;
	Bmodel*					mp, *mp2, *mpc;
	
	if ( verbose )
		cout << "Calculating axes from model fits:" << endl;

	if ( verbose )
		cout << "#\t#\tAxis\tAngle\tRMSD" << endl;
	for ( mp = model; mp->next; mp = mp->next ) {
		for ( mp2 = mp->next; mp2; mp2 = mp2->next ) {
			if ( model_select_corresponding_CA(mp, mp2) ) {
				t = model_find_transform(mp2, mp);
				tv.push_back(t);
				mpc = mp2->copy();
				model_rotate(mpc, t);
				t.fom = model_compare_corresponding(mp, mpc);
				if ( verbose )
					cout << mp->identifier() << tab << mp2->identifier() << tab
						<< t.axis << tab << t.angle*180.0/M_PI << tab << t.fom << endl;
			}
		}
	}

	return tv;
}

/**
@brief	Calculates the RMSD for symmetry in the standard orientation.
@param 	*model			linked list of models.
@param 	*sym			point group symmetry.
@param 	ref_view		reference view.
@return double			RMSD.

	The centers-of-mass for all the molecules are calculated as a reduced
	representation of the molecule group. All symmetry operations are
	imposed on the centers, with the RMSD defined as minimum distance
	between an original cneter and a transformed center.

**/
double		model_symmetry_RMSD(Bmodel* model, Bsymmetry& sym, View2<double> ref_view)
{
	int				i, j, k, l, lm, n;
	double			d, mind, R(0);
	Vector3<double>	coms, comall;
	vector<string>	ids;
	Matrix3			mat(1);
	Matrix3			ref_mat = ref_view.matrix();
	
	Bmodel*			mp;
	
	long			nmol = model->count();
	
	vector<Vector3<double>>	com(nmol);
	
	for ( i=0, mp=model; mp; mp=mp->next, i++ ) {
		com[i] = mp->center_of_coordinates();
		comall += com[i];
		ids.push_back(mp->identifier());
	}
	comall /= nmol;
	for ( i=0; i<nmol; ++i ) com[i] -= comall;

	if ( verbose ) {
		cout << "Calculating RMSD for symmetry " << sym.label() << endl;
		cout << "Reference view:                 " << ref_view << endl;
		cout << "Center of coordinates:          " << comall << endl;
	}
	
	if ( verbose & VERB_PROCESS )
		cout << "Order\tModel\tModel\tDistance" << endl;
	for ( i=n=0; i<sym.operations(); i++ ) {
		for ( j=1; j<sym[i].order(); j++ ) {
			mat = Matrix3(sym[i].axis(), j*TWOPI/sym[i].order());
			mat *= ref_mat;
			for ( k=0; k<nmol; k++ ) {
				coms = mat * com[k];
				for ( mind=1e30, l=0; l<nmol; l++ ) {
					d = coms.distance(com[l]);
					if ( mind > d ) {
						mind = d;
						lm = l;
					}
				}
				if ( verbose & VERB_PROCESS )
					cout << j << tab << ids[k] << tab << ids[lm] << tab << mind << endl;
				R += mind*mind;
				n++;
			}
		}
	}
	
	R = sqrt(R/n);
	
	if ( verbose )
		cout << "Symmetry RMSD:                  " << R << " A (" << n << ")" << endl << endl;
	
	return R;
}


/**
@brief	Calculates the B factors from symmetry-related models.
@param 	*model			linked list of models.
@param 	*sym			point group symmetry.
@return double			RMSD.

	The centers-of-mass for all the models are calculated as a reduced
	representation of the set. All symmetry operations are
	imposed on the centers, with the RMSD defined as minimum distance
	between an original cneter and a transformed center.

**/
double		model_symmetry_B(Bmodel* model, Bsymmetry& sym, View2<double> ref_view)
{
	long			i, j, nmol(0), nsym(0);
	double			d, mind, R(0);
	Vector3<double>	coms, coor;
	Matrix3			m(1);
	Matrix3			ref_mat = ref_view.matrix();

	Bmodel*			mp, *mpr;
	Bcomponent*		comp, *compr;
	
	nmol = model->count();
	
	if ( verbose )
		cout << "Calculating B-factors for symmetry " << sym.label() << endl;
	
	vector<Vector3<double>>	com(nmol);
	vector<Matrix3>			mat = sym.matrices();
	nsym = mat.size();

	for ( auto& m: mat ) m *= ref_mat;

	for ( i=0, mp=model; mp; mp=mp->next, i++ ) {
		com[i] = mp->center_of_coordinates();
		for ( mind=1e30, j=0; j<nsym; j++ ) {
			coms = mat[j] * com[i];
			d = com[0].distance(coms);
			if ( mind > d ) {
				mind = d;
				mp->select(j);
			}
		}
		cout << mp->identifier() << endl << mat[mp->select()] << endl;
	}
	
	// Initialize the vectors for the first molecule
	mpr = model;
	for( comp = model->comp; comp; comp = comp->next ) {
		comp->FOM(0);
		comp->velocity() = Vector3<double>(0,0,0);
	}
	
	// Accumulate the sums and square sums in the first molecule atoms
	for ( i=0, mp = model; mp; mp = mp->next, i++ ) {
		for( comp = mp->comp, compr = mpr->comp; comp; comp = comp->next, compr = compr->next ) {
			coor = mat[mp->select()] * comp->location();
			compr->velocity() += coor;
			compr->FOM(compr->FOM() + coor.length2());
		}
	}
	
	// Calculate the B factors in the first molecule
	for( comp = model->comp; comp; comp = comp->next )
		comp->FOM((comp->FOM() - comp->velocity().length2()/nmol)/nmol);
	
	// Distribute the B factors to all molecules
	for ( mp = model; mp; mp = mp->next )
		for( comp = mp->comp, compr = mpr->comp; comp; comp = comp->next, compr = compr->next )
			comp->FOM(compr->FOM());
	
	return R;
}

/**
@brief 	Generates unit cells from a set of coordinates.
@param 	*model		molecule group.
@param 	uc			unit cell dimensions.
@param 	lattice		number of unit cells in each lattice direction.
@return int 			0, <0 if error.

	The input model is replicated to generate the requested number
	of copies in each lattice direction.

**/
int 		models_generate_lattice(Bmodel* model, UnitCell uc, Vector3<long> lattice)
{
	if ( lattice.volume() < 2 ) return 0;

	int				i, x, y, z, nmodel(0);
	Vector3<double>	d;
	Matrix3			mat = uc.skew_matrix_inverse();
	Bmodel*			mp = model;
	Bmodel*			nu_mp = mp;
	Bcomponent		*comp;
	
	for ( nmodel=0, mp = model; mp; mp = mp->next ) {
		nmodel++;
		nu_mp = mp;
	}
	
	if ( verbose ) {
		cout << "Generating new unit cells for " << nmodel << " models:" << endl;
		cout << "Lattice:                        " << lattice << " = " << (long)lattice.volume() << endl;
		cout << mat << endl;
	}
	
	for ( z=0; z<lattice[2]; z++ ) {
		for ( y=0; y<lattice[1]; y++ ) {
			for ( x=0; x<lattice[0]; x++ ) {
				d = Vector3<double>(x,y,z);
				d = mat * d;
				if ( verbose & VERB_FULL )
					cout << "Generating unit cell:           " << x << " " << y << " " << z << tab << d << endl;
				if ( x+y+z > 0 ) for ( i=0, mp = model; i<nmodel; mp = mp->next, i++ ) {
					nu_mp->next = mp->copy();
					nu_mp = nu_mp->next;
					for ( comp = nu_mp->comp; comp; comp = comp->next )
						comp->location() = comp->location() + d;
				}
			}
		}
	}
	
	return 0;
}

