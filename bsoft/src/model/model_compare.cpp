/**
@file	model_compare.cpp
@brief	Functions to compare models and components
@author 	Bernard Heymann
@date	Created: 20060908
@date	Modified: 20250615
**/

#include "model_compare.h"
#include "model_select.h"
#include "model_transform.h"
#include "model_util.h"
//#include "Bimage.h"
#include "Matrix3.h"
//#include "rwresprop.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief	Compares the numbers of components in two models.
@param 	*model1			first model.
@param 	*model2			second model.
@return long			difference in the number of components.

	Only the first models in the lists are compared.

**/
long		model_component_number_difference(Bmodel* model1, Bmodel* model2)
{
	long			n1, n2;
	Bcomponent*		comp1, *comp2;
	
	for ( n1 = 0, comp1 = model1->comp; comp1; comp1 = comp1->next ) n1++;

	for ( n2 = 0, comp2 = model2->comp; comp2; comp2 = comp2->next ) n2++;
	
	if ( verbose & VERB_PROCESS ) if ( n1 - n2 )
		cout << "Component numbers are different:   " << n1 << " != " << n2 << endl;

	return n1 - n2;
}

/**
@brief	Calculates the largest number of components in a model.
@param 	*model			model.
@return long			largest number of components in a model.
**/
long		model_maxnum_components(Bmodel* model)
{
	long	ncomp, n;
	Bmodel*			mp;
	Bcomponent*		comp;
	
	for ( ncomp = 0, mp = model; mp; mp = mp->next ) {
		for ( n=0, comp = mp->comp; comp; comp = comp->next ) n++;
		if ( ncomp < n ) ncomp = n;
	}

	return ncomp;
}

/**
@brief	Compares two models by corresponding selected compoents.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@return double			RMSD.

	Only the first models in the linked lists are compared.

**/
double		model_compare(Bmodel* model1, Bmodel* model2)
{
	long			n(0);
	double			d, R(0);
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	if ( verbose & VERB_PROCESS )
		cout << "Comparing " << model1->identifier() << " with " << model2->identifier() << ":" << endl;

	// The first selected component on both models sets the alignment mapping
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) break; 	
	for ( comp2 = model2->comp; comp1 && comp2; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( ; comp2; comp2 = comp2->next ) if ( comp2->select() ) break;
		if ( !comp2 ) {
			cerr << "Error: No second component selected!" << endl;
			bexit(-1);
		}	
		if ( verbose & VERB_FULL )
			cout << "Mapping component " << comp1->identifier() << " to " << comp2->identifier() << endl;
		d = comp1->location().distance(comp2->location());
		R += d*d;
		n++;
		comp2 = comp2->next;
		if ( verbose & VERB_FULL )
			cout << comp1->description()[4] << tab << comp1->description()[2] << tab << 
				comp2->description()[4] << tab << comp2->description()[2] << tab << d << endl;
	}
	
	R = sqrt(R/n);

	if ( verbose & VERB_PROCESS )
		cout << "RMSD:                           " << R << " (" << n << ")" << endl;
	
	return R;
}

/**
@brief	Compares two models by closest distance between selected components.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@return double			RMSD.

	Only the first models in the linked lists are compared.

**/
double		model_compare_by_distance(Bmodel* model1, Bmodel* model2)
{
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	Bcomponent*		compsel;
	
	long			n(0), n1(0), n2(0);
	double			d, dmin, R(0);

	if ( verbose & VERB_PROCESS )
		cout << "Comparing " << model1->identifier() << " with " << model2->identifier() << ":" << endl;
	
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		dmin = 1e30;
		compsel = NULL;
		for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			if ( dmin > d ) {
				dmin = d;
				compsel = comp2;
			}
		}
		if ( compsel ) {
			comp1->select(stol(compsel->identifier()));
			comp1->FOM(dmin);
		}
		n1++;
	}

	for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
		dmin = 1e30;
		compsel = NULL;
		for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
			d = comp1->location().distance(comp2->location());
			if ( dmin > d ) {
				dmin = d;
				compsel = comp1;
			}
		}
		if ( compsel ) {
			comp2->select(stol(compsel->identifier()));
			comp2->FOM(dmin);
		}
		n2++;
	}

	if ( verbose & VERB_FULL )
		cout << "Comp1\tComp2\tDmin" << endl;
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = model2->comp; comp2 && comp1->select() != stoi(comp2->identifier()); comp2 = comp2->next ) ;
		if ( comp2 && comp2->select() == stoi(comp1->identifier()) ) {
			R += comp1->FOM()*comp1->FOM();
			n++;
			if ( verbose & VERB_FULL )
				cout << comp1->identifier() << tab << comp2->identifier() << tab << comp1->FOM() << endl;
		}
	}

	R = sqrt(R/n);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Components compared:            " << n << endl;
		cout << "Model 1 components:             " << n1 << " (" << n*100.0/n1 << "%)" << endl;
		cout << "Model 2 components:             " << n2 << " (" << n*100.0/n2 << "%)" << endl;
		cout << "RMSD:                           " << R << endl << endl;
	}
	
	return R;
}

/**
@brief	Compares two models by corresponding selected components.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@return double			RMSD.

	Only the first models in the linked lists are compared.

**/
double		model_compare_corresponding(Bmodel* model1, Bmodel* model2)
{
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	long			n(0);
	double			d, R(0);

	if ( verbose & VERB_PROCESS )
		cout << "Comparing " << model1->identifier() << " with " << model2->identifier() << ":" << endl;

	if ( model1->component_count_selected() != model2->component_count_selected() ) {
		cerr << "Different selected component counts!" << endl;
		return 0;
	}
		
	for ( comp1 = model1->comp, comp2 = model2->comp; comp1 && comp2; comp1 = comp1->next, comp2 = comp2->next ) 
		if ( comp1->select() && comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			R += d*d;
			n++;
	}

	if ( n ) R = sqrt(R/n);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Components compared:            " << n << endl;
		cout << "RMSD:                           " << R << endl << endl;
	}
	
	return R;
}

/**
@brief	Calculates the RMSD between two corresponding sets of models.
@param 	*model1			first model set.
@param 	*model2			second model set.
@return vector<double>	array of RMSDs.

	The number of selected models and their order must correspond.

**/
vector<double>	models_compare_corresponding(Bmodel* model1, Bmodel* model2)
{
	long			n(0);
	vector<double>	rmsd;

	if ( !model1 || !model2 ) return rmsd;
	
	long			nmod1 = model1->count_selected();
	long			nmod2 = model2->count_selected();
	if ( nmod1 != nmod2 ) {
		cerr << "The two sets of models must ahve the same number selected!" << endl;
		return rmsd;
	}
	
	rmsd.resize(nmod1);
	
	if ( verbose ) cout << "Comparing corresponding models:" << endl;
	
	Bmodel*			mp1;
	Bmodel*			mp2;
	
	if ( verbose )
		cout << "#\tModel1\tModel2\tRMSD" << endl;
	for ( mp1 = model1, mp2 = model2; mp1 && mp2; mp1 = mp1->next, mp2 = mp2->next ) 
		if ( mp1->select() && mp2->select() ) {
			rmsd[n] = model_compare_corresponding(mp1, mp2);
			if ( verbose )
				cout << n << tab << mp1->identifier() << tab << mp2->identifier() << tab << rmsd[n] << endl;
			n++;
	}

	return rmsd;
}


/**
@brief	Counts the number of components overlapping between two sets of models.
@param 	*model1			first model set.
@param 	*model2			second model set.
@param 	dcut			distance between components to define an interface.
@return Matrix			matrixs of interface counts.

	Only the selected models in the linked lists are compared.

**/
Matrix		models_interfaces(Bmodel* model1, Bmodel* model2, double dcut)
{
	long			n(0), n1(0), n2(0);
	double			d;
	Matrix			mat;

	if ( !model2 ) model2 = model1;
	
	if ( !model1 ) return mat;
	
	long			nmod1 = model1->count_selected();
	long			nmod2 = model2->count_selected();
	if ( nmod1<1 || nmod2<1 ) return mat;
	
	mat = Matrix(nmod1,nmod2);
	
	Bmodel*			mp1;
	Bmodel*			mp2;
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	if ( verbose & VERB_PROCESS )
		cout << "#1\tModel1\t#2\tModel2\tOverlap" << endl;
	for ( mp1 = model1; mp1; mp1 = mp1->next ) if ( mp1->select() ) {
		n2 = 0;
		for ( mp2 = model2; mp2; mp2 = mp2->next ) if ( mp2->select() ) {
			if ( mp1 != mp2 ) {
				n = 0;
				for ( comp1 = mp1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
					for ( comp2 = mp2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
						d = comp1->location().distance(comp2->location());
						if ( d <= dcut ) n++;
					}
				}
				if ( verbose & VERB_PROCESS )
					cout << n1 << tab << mp1->identifier() << tab << n2 << tab << mp2->identifier() << tab << n << endl;
				mat[n1][n2] = n;
			}
			n2++;
		}
		n1++;
	}

	return mat;
}

/**
@brief	Finds the interfaces between two models.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@param 	dcut			distance between components to define the interface.
@return long				number of components in interfaces.

	Only the first models in the linked lists are compared.

**/
long		model_interface(Bmodel* model1, Bmodel* model2, double dcut)
{
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	long			n(0), n1(0), n2(0);
	double			d;

	if ( verbose )
		cout << "Interface between " << model1->identifier() << " and " << model2->identifier() << ":" << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Comp1\tType1\tComp2\tType2\tDist" << endl;
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			d = comp1->location().distance(comp2->location());
			if ( d <= dcut ) {
				if ( verbose & VERB_PROCESS )
					cout << comp1->description()[4] << tab << comp1->description()[2] << tab 
						<< comp2->description()[4] << tab << comp2->description()[2] << tab << d << endl;
				comp1->select(2);
				comp2->select(2);
				n++;
			}
		}
	}

	if ( verbose )
		cout << "Comp1\tType1" << endl;
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() == 2 ) {
		cout << comp1->description()[4] << tab << comp1->description()[2] << endl; 
		n1++;
	}

	if ( verbose )
		cout << "Comp2\tType2" << endl;
	for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() == 2 ) {
		cout << comp2->description()[4] << tab << comp2->description()[2] << endl; 
		n2++;
	}

	if ( verbose ) {
		cout << "Components in interfaces:       " << n << endl;
		cout << "Components in model 1:          " << n1 << endl;
		cout << "Components in model 2:          " << n2 << endl;
	}
	
	return n;
}

/**
@brief	Finds the interfaces between two models.
@param 	*model1			first model structure.
@param 	*model2			second model structure.
@param 	res_prop		residue properties.
@return long				number of components in interfaces.

	Only the first models in the linked lists are compared.

**/
long		model_interface(Bmodel* model1, Bmodel* model2, map<string,Bresidue_type> res_prop)
{
	Bcomponent*		comp1;
	Bcomponent*		comp2;
	
	long			n(0), n1(0), n2(0);
	double			d, res1_ext, res2_ext;

	if ( verbose )
		cout << "Interface between " << model1->identifier() << " and " << model2->identifier() << ":" << endl;
	
	if ( verbose & VERB_PROCESS )
		cout << "Comp1\tType1\tComp2\tType2\tDist" << endl;
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() ) {
		res1_ext = res_prop[comp1->description()[2]].extension();
		for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			res2_ext = res_prop[comp2->description()[2]].extension();
			d = comp1->location().distance(comp2->location());
			if ( d <= res1_ext + res2_ext ) {
				if ( verbose & VERB_PROCESS )
					cout << comp1->description()[4] << tab << comp1->description()[2] << tab 
						<< comp2->description()[4] << tab << comp2->description()[2] << tab << d << endl;
				comp1->select(2);
				comp2->select(2);
				n++;
			}
		}
	}

	if ( verbose )
		cout << "Comp1\tType1" << endl;
	for ( comp1 = model1->comp; comp1; comp1 = comp1->next ) if ( comp1->select() == 2 ) {
		cout << comp1->description()[4] << tab << comp1->description()[2] << endl; 
		n1++;
	}

	if ( verbose )
		cout << "Comp2\tType2" << endl;
	for ( comp2 = model2->comp; comp2; comp2 = comp2->next ) if ( comp2->select() == 2 ) {
		cout << comp2->description()[4] << tab << comp2->description()[2] << endl; 
		n2++;
	}

	if ( verbose ) {
		cout << "Components in interfaces:       " << n << endl;
		cout << "Components in model 1:          " << n1 << endl;
		cout << "Components in model 2:          " << n2 << endl;
	}
	
	return n;
}

/**
@brief 	Constructs the adjacency matrix for a model.
@param 	*model			model structure.
@return Matrix			adjacency matrix.

	The matrix contains ones for adjacent components and zero elsewhere. 
	Only the first model in the linked list is used.
	The component selections are reset.

**/
Matrix		model_adjacency_matrix(Bmodel* model)
{
	int				i, j, k, ncomp;
	Bcomponent*		comp;
	
	for ( ncomp = 0, comp = model->comp; comp; comp = comp->next, ncomp++ )
		comp->select(ncomp);
	
	Matrix			mat(ncomp,ncomp);

	for ( comp = model->comp; comp; comp = comp->next ) {
		i = comp->select();
		for ( k=0; k<comp->link.size() && comp->link[k]; k++ ) {
			j = comp->link[k]->select();
			mat[i][j] = 1;
		}
	}

	for ( comp = model->comp; comp; comp = comp->next )
		comp->select(1);
	
	return mat;
}

/**
@brief 	Calculates the distance matrix for a model.
@param 	*model			model structure.
@param 	view_flag		flag to include views.
@return Matrix			distance matrix.

	The matrix is calculated from the pairwise distances between selected 
	components, with the option to include the views.
	When the views are included, the euclidean distances are rescaled to
	the maximum so that their ranges are similar to the view differences. 
	Only the first model in the linked list is used.

**/
Matrix		model_distance_matrix(Bmodel* model, int view_flag)
{
	long			n, i, j;
	double			d, dmax;
	Bcomponent*		comp;
	Bcomponent*		comp2;
	
	for ( n=0, dmax = 0, comp = model->comp; comp; comp = comp->next ) if ( comp->select() ) {
		d = comp->location().length();
		if ( dmax < d ) dmax = d;
		n++;
	}
	
	if ( verbose )
		cout << "Calculating a " << n << " x " << n << " distance matrix:" << endl << endl;

	Matrix			dmat(n, n);
	
	for ( i=0, comp = model->comp; comp->next; comp = comp->next ) if ( comp->select() ) {
		for ( j=i+1, comp2 = comp->next; comp2; comp2 = comp2->next ) if ( comp2->select() ) {
			if ( view_flag )
				dmat[i][j] = dmat[j][i] =
					comp->location().distance(comp2->location())/dmax + comp->view().residual(comp2->view());
			else
				dmat[i][j] = dmat[j][i] = comp->location().distance(comp2->location());
			j++;
		}
		i++;
	}

	if ( verbose & VERB_FULL )
		cout << dmat << endl;
	
	return dmat;
}

/**
@brief 	Calculates the distance matrix between two models.
@param 	*m1				first model structure.
@param 	*m2				second model structure.
@return Matrix			distance matrix.

	The matrix is calculated from the pairwise distances between selected 
	components. 
	Only the first model in each linked list is used.

**/
Matrix		model_distance_matrix(Bmodel* m1, Bmodel* m2)
{
	long			i, j, n1(0), n2(0);
	Bcomponent		*c1, *c2;
	
	n1 = m1->component_count();
	n2 = m2->component_count();
	
	Matrix			mat(n1, n2);

	for ( i=0, c1 = m1->comp; c1; c1 = c1->next, ++i )
		for ( j=0, c2 = m2->comp; c2; c2 = c2->next, ++j )
			mat[i][j] = c1->location().distance(c2->location());
	
	return mat;
}

/**
@brief 	Consolidates close components within a model.
@param 	*model			model structure.
@param 	distance		cutoff distance to consider components to be the same.
@return long			number of components retained.

	A matrix is calculated from the pairwise distances between selected components.
	Components closer to each other than the given distance are consolidated. 
	Only the first model in the linked list is used.

**/
long		model_consolidate(Bmodel* model, double distance)
{
	long			i, j, k, n(0);
	Bcomponent		*c1, *c2;
	Bcomponent*		c = NULL;
	Bcomponent*		c_list = NULL;
	string			cid;
	Matrix			mat;

	n = model->component_count();

	mat = model_distance_matrix(model, 0);
	
	if ( verbose ) {
		cout << model->identifier() << tab << n << endl;
		mat.show_below_cutoff(distance);
	}
	
	for ( i=0; i<mat.rows(); ++i ) mat[i][i] = 2*distance;

	for ( i=n=0, c1 = model->comp; c1; c1 = c1->next, ++i ) if ( mat[i][i] ) {
		cid = to_string(++n);
//		c = component_add(&c_list, cid);
//		component_copy(c1, c);
		if ( c_list ) c = c_list->add(c1);
		else c = c_list = new Bcomponent(c1);
		c->identifier() = cid;
		for ( j=i+1, k=1, c2 = c1->next; c2; c2 = c2->next, ++j ) {
			if ( mat[i][j] <= distance ) {
				c->shift(c2->location());
				mat[j][j] = 0;
				k++;
			}
		}
		if ( k > 1 ) c->scale(1.0L/k);
		c->select(k);
	}
	
//	component_list_kill(model->comp);
	model->clear_components();
	model->comp = c_list;

	if ( verbose )
		cout << "Number of components:     " << n << endl << endl;
	
	return n;
}

/**
@brief 	Determines a consensus between models in a list.
@param 	*model			model list.
@param 	distance		cutoff distance to consider components to be the same.
@return long			number of components retained.

	A matrix is calculated from the pairwise distances between selected components
	from each apir of models.
	Components closer to each other than the given distance in different
	models are consolidated.
	The component selection field contains the number of contributing components.
	A new model containing the result is returned.

**/
Bmodel*		models_consensus(Bmodel* model, double distance)
{
	long			nmod(0), i, j, k, n;
	Bmodel			*m1, *m2;
	Bcomponent		*c1, *c2;
	
	for ( m1 = model; m1; m1 = m1->next ) nmod++;
	
	for ( m1 = model; m1; m1 = m1->next )
		model_consolidate(m1, distance);
	
	long			nmat(factorial(nmod)/(2*factorial(nmod-2)));
	
	Matrix*			mat = new Matrix[nmat];
				
	for ( i=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( m2 = m1->next; m2; m2 = m2->next, ++i ) {
			mat[i] = model_distance_matrix(m1, m2);
			if ( verbose & VERB_FULL ) {
				cout << m1->identifier() << " vs " << m2->identifier() << endl;
				mat[i].show_below_cutoff(distance);
			}
		}
	}
	
	models_select_all(model);

	for ( i=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( m2 = m1->next; m2; m2 = m2->next, ++i ) {
			for ( j=0, c1 = m1->comp; c1; c1 = c1->next, ++j ) if ( c1->select() ) {
//				c1->select(1);
				for ( k=0, c2 = m2->comp; c2; c2 = c2->next, ++k ) if ( c2->select() ) {
//					c2->select(1);
					if ( mat[i][j][k] <= distance ) {
						c1->shift(c2->location());
						c1->select_increment();
						c2->select(0);
					}
				}
//				c1->loc /= c1->select();
			}
		}
	}
	
	delete[] mat;
	
	Bmodel*			numod = new Bmodel(model->identifier());
	string			cid;
	RGBA<float>		color;
	int*			ns = new int[nmod+1];
	for ( i=0; i<=nmod; i++ ) ns[i] = 0;
	
	numod->mapfile() = model->mapfile();
	numod->image_number(model->image_number());
	numod->comment(model->comment());
	Bcomptype*		ct = numod->add_type(model->type->identifier());

	for ( n=0, m1 = model; m1->next; m1 = m1->next ) {
		for ( c1 = m1->comp; c1; c1 = c1->next ) if ( c1->select() ) {
			cid = to_string(++n);
			c2 = numod->add_component(c1);
			c2->type(ct);
			c2->identifier(cid);
			c2->location(c2->location()/c2->select());
			c2->FOM(c2->select()*1.0L/nmod);
			ns[c2->select()]++;
			color.spectrum(c2->select(), nmod - 0.5, 0.5);
			c2->color(color);
			if ( verbose & VERB_FULL )
				cout << n << tab << c2->location() << tab << c2->select() << endl;
		}
	}
	
	if (verbose ) {
		cout << "Number\tCount\t%" << endl;
		for ( i=0; i<=nmod; i++ )
			cout << i << tab << ns[i] << tab << ns[i]*100.0/n << endl;
		cout << endl;
	}
	
	delete[] ns;
	
	return numod;
}

/**
@brief 	Fits a model to a reference model.
@param 	model			model structure.
@param 	refmod			reference model.
@param 	id				model id to fit.
@return int				0.

**/
int			model_fit(Bmodel* model, Bmodel* refmod, string id)
{
	Bmodel*			mp = model;
	
	if ( id.length() )
		while ( mp && mp->identifier() != id ) mp = mp->next;
	
	if ( !mp ) {
		cerr << "Error: Model with id " << id << " not found!" << endl;
		bexit(-1);
	}
	
	Transform		t = model_find_transform(mp, refmod);

	models_rotate(model, t);
	
	return 0;
}


