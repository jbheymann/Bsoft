/**
@file	model_symmetry.cpp
@brief	Library routines used for model symmetry operations
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20221115
**/

#include "model_util.h"
#include "model_transform.h"
#include "model_select.h"
#include "model_compare.h"
#include "mol_transform.h"
#include "mol_compare.h"
#include "mol_util.h"
#include "symmetry.h"
#include "Matrix3.h"
#include "random_numbers.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Set model component locations within the asymmetric unit.
@param 	*model				model parameters.
@param 	symmetry_string		symmetry code.
@return long					number of components (<0 means failure).

	Only the first model is processed.

**/
long		model_find_asymmetric_unit(Bmodel* model, string symmetry_string)
{
	Bsymmetry		sym(symmetry_string);
	
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
long		model_apply_point_group(Bmodel* model, string symmetry_string,
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
//				nu_link = (Blink *) add_item((char **) &link, sizeof(Blink));
				nu_link = link->add(link);
				if ( !link_start ) link_start = nu_link;
//				nu_link->comp[0] = link->comp[0];
//				nu_link->comp[1] = link->comp[1];
//				nu_link->length(link->length());
//				nu_link->radius(link->radius());
//				nu_link->color(link->color());
			}
			for ( k=0, comp = model->comp; k<ncomp; comp = comp->next, k++ ) {
//				nu_comp = (Bcomponent *) add_item((char **) &comp, sizeof(Bcomponent));
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

/**
@brief 	Applying symmetry to model components.
@param 	*model				model parameters.
@param 	symmetry_string		symmetry code.
@param 	origin				transformation origin.
@param 	ref_view			reference view.
@param 	flags				1=find asu.
@return long					error code (<0 means failure).

	All models in the list are processed.

**/
long		models_apply_point_group(Bmodel* model, string symmetry_string,
					Vector3<double> origin, View2<double> ref_view, int flags)
{
	long 			ncomp(0);
	Bmodel*			mp;
	
	for ( mp = model; mp; mp = mp->next )
		ncomp += model_apply_point_group(mp, symmetry_string, origin, ref_view, flags);
		
	return ncomp;
}

/**
@brief 	Symmetrize a model.
@param 	*model				model parameters.
@param 	symmetry_string		symmetry code.
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
long		model_symmetrize(Bmodel* model, string symmetry_string)
{
	Bsymmetry		sym(symmetry_string);
	
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

	Bcomponent*		comp;
	Bcomponent*		compasu;
	Bcomponent*		compsel = NULL;
	
	// Select all components in the ASU
	for ( nasu=0, comp = model->comp; comp; comp = comp->next ) {
		comp->select(0);
//		view = View2<float>(comp->location()[0], comp->location()[1], comp->location()[2], 0);
//		view_asu = find_asymmetric_unit_view(sym, view);
		view = View2<double>(comp->location());
		view_asu = sym.find_asymmetric_unit_view(view);
		comp->force(Vector3<double>(view_asu[0], view_asu[1], view_asu[2]));
		if ( view.distance(view_asu) < tol ) {
			comp->location(comp->force() * comp->location().length());
			comp->select(1);
			nasu++;
		}
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( !comp->select() ) {
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
		cout << "Model:                          " << model->identifier() << endl;
		cout << "Asymmetric unit components:     " << nasu << endl;
		cout << "Comp\tCount" << endl;
	}
	
	for ( comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		comp->scale(1.0/comp->select());
		comp->force(Vector3<float>(0,0,0));
		if ( verbose & VERB_PROCESS )
			cout << comp->identifier() << tab << comp->select() << endl;
	}
	if ( verbose & VERB_PROCESS )
		cout << endl;
	
	R = sqrt(R/ncomp);

	if ( verbose & VERB_PROCESS )
		cout << "Symmetry deviation:             " << R << " A" << endl << endl;
		
	model_delete_non_selected(&model);
	
	Vector3<double>		origin;
	View2<double>		ref_view;
	
	model_apply_point_group(model, symmetry_string, origin, ref_view, 1);	

	return ncomp;
}

/**
@brief 	Generates a list of symmetry-related models.
@param 	*model				model parameters.
@param 	symmetry_string	symmetry code.
@return long				total number of components.

	For each component, a new location is calculated from the average location
	of the closest symmetry-related components.

**/
/*long		model_symmetry_related(Bmodel* model, string symmetry_string)
{
	Bsymmetry		sym(symmetry_string);
	
	View2<float>	v2;
	View*			v;
	View			ref(0,0,1,0);
	View*			views = symmetry_get_all_views(sym, ref);
	
	long			i(1), ncomp(model->component_count());
	string			id;
	Bmodel*			mp = model;

	for ( v = views->next; v; v = v->next ) {
//		mp->next = model_copy(model);
		mp->next = model->copy();
		mp = mp->next;
		id = model->identifier() + "_" + to_string(++i);
		mp->identifier() = id;
		v2 = View2<float>((*v)[0],(*v)[1],(*v)[2],(*v)[3]);
		model_rotate(mp, v2);
	}
	
	ncomp *= i;
	
	return ncomp;
}*/
long		model_symmetry_related(Bmodel* model, string symmetry_string)
{
	Bsymmetry		sym(symmetry_string);
	
//	View2<double>	v2;
	View2<double>	ref(0,0,1,0);
	vector<View2<double>> views = sym.get_all_views(ref);
	
	long			i(1), ncomp(model->component_count());
	string			id;
	Bmodel*			mp = model;

	for ( auto v: views ) {
		mp->next = model->copy();
		mp = mp->next;
		id = model->identifier() + "_" + to_string(++i);
		mp->identifier() = id;
//		v2 = View2<float>((*v)[0],(*v)[1],(*v)[2],(*v)[3]);
		model_rotate(mp, v);
	}
	
	ncomp *= i;
	
	return ncomp;
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
int 		model_generate_lattice(Bmodel* model, UnitCell uc, Vector3<long> lattice)
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

