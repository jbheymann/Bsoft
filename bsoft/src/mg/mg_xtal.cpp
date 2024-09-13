/**
@file	mg_xtal.cpp
@brief	Functions to process crystallographic data
@author	Bernard Heymann
@date	Created: 20061110
@date	Modified: 20240306
**/

#include "mg_processing.h"
#include "mg_xtal.h"
#include "matrix_linear.h"
#include "Matrix.h"
#include "linked_list.h"
#include "Complex.h"
#include "utilities.h"
#include "timer.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates the unit cell vectors for a 2D crystal.
@param 	*mg			micrograph.
@return int			0, <0 on error.

	Finds the unit cell vectors u and v by solving the equation:
		x = uh + vk
	where x is the location of the reflection or structure factor,
	and h and k are the associated Miller indices.

**/
int			mg_unitcell_vectors(Bmicrograph* mg)
{
	long			i;
	vector<double>	bx(2,0), by(2,0);
	Matrix			a(2,2);
	Bstrucfac*		sf;
	
	for ( i = 0, sf = mg->sf; sf; sf = sf->next, i++ ) {
		a[0][0] += sf->index[0]*sf->index[0];
		a[0][1] += sf->index[0]*sf->index[1];
		a[1][1] += sf->index[1]*sf->index[1];
		bx[0] += sf->index[0]*sf->loc[0];
		bx[1] += sf->index[1]*sf->loc[0];
		by[0] += sf->index[0]*sf->loc[1];
		by[1] += sf->index[1]*sf->loc[1];
	}
	a[1][0] = a[0][1];
	
	if ( i < 2 ) {
		cerr << "Error: Too few structure factors to calculate unit cell vectors!" << endl;
		return -1;
	}
	
	a.LU_decomposition();
	a.multiply_in_place(bx);
	a.multiply_in_place(by);
	
	mg->hvec[0] = bx[0];
	mg->hvec[1] = by[0];
	mg->kvec[0] = bx[1];
	mg->kvec[1] = by[1];
		
	return 0;
}

/**
@brief 	Generates reflections given the unit cell vectors.
@param 	*mg				micrograph.
@param	real_size		physical image size.
@param 	resolution		resolution limit.
@return long			number of reflections generated, <0 on error.

	The structure factor location is given by:
		x = uh + vk
	where u and v are the unit cell vectors,
	and h and k are the associated Miller indices.

**/
long		mg_generate_reflections(Bmicrograph* mg, Vector3<double> real_size, double resolution)
{
	if ( resolution < 2*mg->pixel_size[0] )
		resolution = 2*mg->pixel_size[0];
	
	Vector3<double>		lim2 = real_size/resolution + 0.001;
	double				k_lim = lim2[0];
	if ( k_lim < lim2[1] ) k_lim = lim2[1];
	lim2 *= lim2;
//	double				hkl_lim = real_size[0]/resolution;
//	if ( real_size[1] > real_size[0] ) hkl_lim = real_size[1]/resolution;
//	double				hkl_lim2 = hkl_lim*hkl_lim + 0.001;
	
	if ( verbose )
		cout << "Generating reflections to " << resolution << " Å" << endl;
		
	long				n(0), h, k, l, hmin, kmin, lmin, hmax, kmax, lmax;
	Vector3<double>		loc, loc2;
	Bstrucfac*			sf = NULL;
	
	hmax = (long) (k_lim/mg->hvec.length()+1);
	hmin = -hmax;

	kmax = (long) (k_lim/mg->kvec.length()+1);
	kmin = -kmax;
	
	lmin = lmax = 0;
	
	kill_list((char *) mg->sf, sizeof(Bstrucfac));
	mg->sf = NULL;
	
	if ( verbose ) {
		cout << "real_size=" << real_size << endl;
		cout << "hvec=" << mg->hvec << endl;
		cout << "kvec=" << mg->kvec << endl;
		cout << "hmax=" << hmax << " kmax=" << kmax << endl;
		cout << "mg->origin=" << mg->origin << endl;
	}
	
	for ( l=lmin; l<=lmax; l++ ) {
		for ( k=kmin; k<=kmax; k++ ) {
			for ( h=hmin; h<=hmax; h++ ) {
				loc = mg->hvec * h + mg->kvec * k + mg->lvec * l;
				loc2 = (loc*loc)/lim2;
				if ( loc2.length2() < 1 ) {
					sf = (Bstrucfac *) add_item((char **) &sf, sizeof(Bstrucfac));
					if ( !mg->sf ) mg->sf = sf;
					sf->index[0] = h;
					sf->index[1] = k;
					sf->index[2] = l;
					sf->loc = loc;
					sf->fom = sf->sel = 1;
					n++;
				}
			}
		}
	}
	
	return n;
}

/**
@brief 	Masks the image using the list of reflections.
@param 	*p			complex image.
@param 	*sflist		reflection list.
@param 	radius		radius around reflection to mask.
@return int			error code.
**/
int			img_mask_reflections(Bimage* p, Bstrucfac* sflist, double radius)
{
	long			i, datasize = (long) p->size().volume();
//	Vector3<long>	mid(p->size()/2);
	Vector3<double>	loc;
	Bstrucfac*		sf = NULL;

	Bimage*			pmask = new Bimage(Float, p->compound_type(), p->size(), 1);
	Complex<double>	cv;

	if ( p->compound_type() == TSimple ) {	// Power spectrum
		for ( sf = sflist; sf; sf = sf->next )
			pmask->sphere(sf->loc, radius, 0, FILL_USER, 1);
		for ( i=0; i<datasize; i++ )
			if ( (*pmask)[i] < 1 ) p->set(i, 0);
	} else if ( p->compound_type() == TComplex ) {	// Standard Fourier transform
		for ( sf = sflist; sf; sf = sf->next ) {
			loc = vector3_set_PBC(sf->loc, p->size());
			pmask->sphere(loc, radius, 0, FILL_USER, 1);
		}
		for ( i=0; i<datasize; i++ )
			if ( (*pmask)[i] < 1 ) p->set(i, cv);;
	}
	
	delete pmask;
	
    return 0;
}

Vector3<double>	sinc_fit(vector<Vector3<double>> coor, vector<double> value, double& R)
{
	R = 0;
	
	long			i;
	double			k(2), d, R1;
	Vector3<double>	ori, bori, off;
	
	for ( ori[1]=-k; ori[1]<=k; ori[1]+=0.1 ) {
		for ( ori[0]=-k; ori[0]<=k; ori[0]+=0.1 ) {
			for ( i=0, R1=0; i<coor.size(); ++i ) {
				off = coor[i] - ori;
				d = value[i]*sinc(off[0])*sinc(off[1]);
				R1 += d;
			}
			R1 = sqrt(R1/coor.size());
//			cout << ori << tab << R << endl;
			if ( R < R1 ) {
				R = R1;
				bori = ori;
			}
		}
	}
	
	if ( verbose & VERB_FULL ) {
		cout << "Best origin:            " << bori << tab << R << endl;
	}
	
	return bori;
}

Vector3<double>	img_sinc_fit(Bimage* p, Vector3<double> loc, double bkg, double peak, double& R)
{
	long					i, xx, yy, k(2);
//	double					sum(0);
	vector<Vector3<double>>	coor;
	vector<double>			value;
	
	for ( yy=loc[1]-k; yy<=loc[1]+k; ++yy ) {
		for ( xx=loc[0]-k; xx<=loc[0]+k; ++xx ) {
			coor.push_back(Vector3<double>(xx,yy,0)-loc);
			i = p->index(xx, yy);
			value.push_back(((*p)[i]-bkg)/peak);
//			sum += value.back();
		}
	}
	
	loc += sinc_fit(coor, value, R);
	
	return loc;
}

/**
@brief 	Given a spatial frequency, find reflections in a power spectrum.
@param 	*p						power spectrum.
@param 	ref_res					corresponding resolution at reflection.
@param 	threshold				detection threshold.
@param 	kedge					radius around reflection to include in intensity summation.
@return vector<Vector3<double>>	list of reflections.
**/
vector<Vector3<double>>	img_find_reflections(Bimage* p, double ref_res, double threshold, long kedge)
{
	if ( threshold < 0.1 ) threshold = 20;
	if ( kedge < 1 ) kedge = 1;
	
	long				i, j, na;
	long				ks(kedge/2);				// Summation kernel half edge and volume
	double				a;
	double				s(1/ref_res);				// Spatial frequency
	Vector3<double>		k(p->real_size()*s);		// Frequency space pixel distance
	k[2] = 0;
	double				rd(1+5/k.length());			// Relative displacement for reference kernel
	double				angsam(ks);					// Angular sampling in pixels
	double				da(angsam/k.length());		// Angle increment
	double				pmax, ravg, bkg(0), tR;
	Vector3<double>		tloc, rloc, sloc;
	vector<Vector3<double>> loc;
	vector<double>		R, Rfit;
	
	if ( verbose ) {
		cout << "Finding reflections:" << endl;
		cout << "Power spectrum size:            " << p->size() << endl;
		cout << "Spatial frequency:              " << s << " 1/A" << endl;
		cout << "Threshold:                      " << threshold << " e/px" << endl;
		cout << "Kernel edge size:               " << 2*ks+1 << endl;
		cout << "Frequency space pixels:         " << k << endl;
		cout << "Frequency space pixel size:     " << 1/p->real_size()[0] << " 1/Å" << endl;
		cout << "Angular increment:              " << da*180.0/M_PI << " degrees" << endl << endl;
	}

	for ( a=0, j=na=0; a<TWOPI; a+=da, ++na ) {
//		tloc = Vector3<double>(p->image->origin()[0] + k[0]*cos(a), p->image->origin()[1] + k[1]*sin(a), 0);
		tloc = Vector3<double>(k[0]*cos(a), k[1]*sin(a), 0);
		rloc = tloc * rd;
		tloc += p->image->origin();
		rloc += p->image->origin();
//		cout << na << tab << tloc << tab << rloc << endl;
		i = p->kernel_max(p->index(tloc,0), ks);
		pmax = (*p)[i];
		ravg = p->kernel_average(p->index(rloc,0), ks, 0, 1e30);
		bkg += ravg;
		if ( pmax > threshold ) {
			tloc = p->coordinates(i);
			tloc = img_sinc_fit(p, tloc, ravg, pmax - ravg, tR);
			if ( loc.size() && tloc.distance(loc.back()) < 2 ) {
				if ( pmax > R.back() ) {
					loc.back() = tloc;
					R.back() = pmax;
					Rfit.back() = tR;
				}
			} else {
				loc.push_back(tloc);
				R.push_back(pmax);
				Rfit.push_back(tR);
			}
		}
	}
	
	bkg /= na;
	
	if ( loc.size() < 1 ) {
		cerr << "Error: No reflections found, check the threshold." << endl;
		return loc;
	}
	
	// Handle wrap-around reflections
	if ( loc[0].distance(loc.back()) < 2 ) {
		if ( R[0] < R.back() ) {
			loc[0] = loc.back();
			R[0] = R.back();
			Rfit[0] = Rfit.back();
		}
		loc.pop_back();
		R.pop_back();
		Rfit.pop_back();
	}
	
	if ( verbose )
		cout << "#\tx\ty\tz\ta\tR\tRfit" << endl;
	for ( i=0; i<loc.size(); ++i ) {
		sloc = loc[i] - p->image->origin();
		if ( verbose )
			cout << setprecision(2) << i+1 << tab << loc[i] << tab <<
				atan2(sloc[1], sloc[0])*180.0/M_PI << tab << setprecision(3) << R[i] << tab << Rfit[i] << endl;
		R[i] = p->kernel_sum(p->index(loc[i],0), 1) - 9*bkg;
		loc[i] = sloc/p->real_size();
	}
	
	if ( verbose ) {
		cout << "Number of reflections found:    " << loc.size() << endl;
		cout << "Background intensity:           " << bkg << " e/px" << endl << endl;
	}
	
	return loc;
}

/**
@brief 	Find symmetry-related reflections and group them.
@param 	loc						input set of reflections.
@param 	sym						cyclic symmetry to assess.
@param 	symdist					threshold distance to consider reflections as symmetry-related.
@return vector<Vector3<double>>	list of re-ordered reflections.
**/
vector<Vector3<double>>	symmetrize_reflections(vector<Vector3<double>> loc, long sym, double symdist)
{
	bool				found;
	long				i, j, n, ns(0);
	double				sf;
	Vector3<double>		rloc, sloc;
	Matrix3				mat(Vector3<double>(0,0,1), TWOPI/sym);
	vector<int>			sel(loc.size(), 0);
	vector<Vector3<double>> locsym;
	
	if ( verbose ) {
		cout << "Symmetrizing reflections:" << endl;
		cout << "x\ty\tz\tAngle\t#\tResolution" << endl;
	}
	
	for ( i=0; i<loc.size(); ++i ) if ( !sel[i] ) {
		ns++;
		sel[i] = ns;
		rloc = sloc = loc[i];
		sf = sloc.length();
		locsym.push_back(sloc);
		cout << setprecision(2) << loc[i] << tab << atan2(sloc[1], sloc[0])*180.0/M_PI << tab
			<< sel[i] << tab << setprecision(4) << 1/sf << endl;
		for ( n=1; n<sym; ++n ) {
			rloc = mat*rloc;
			found = 0;
			for ( j=0; j<loc.size() && !found; ++j ) if ( !sel[j] ) {
				sloc = loc[j];
				if ( sloc.distance(rloc) < symdist ) {
					found = 1;
					sel[j] = sel[i];
					locsym.push_back(sloc);
					sf = sloc.length();
					if ( verbose )
						cout << setprecision(2) << loc[j] << tab << atan2(sloc[1], sloc[0])*180.0/M_PI << tab
							<< sel[j] << tab << setprecision(4) << 1/sf << endl;
				}
			}
			if ( !found )
				locsym.push_back(rloc);
		}
		if ( verbose )
			cout << endl;
	}

	return locsym;
}

/**
@brief 	Calculate the real space scale from the difference between the reflections and the nominal spatial frequency.
@param 	loc						input set of reflections.
@param 	s						nominal spatial frequency.
@return Vector3<double>			real space scale.
**/
Vector3<double>	analyze_reflections(vector<Vector3<double>>& loc, double s)
{
//	long				i;
	double				u2(0), v2(0), iso(0);
	double				sf, sf_avg(0), sf_std(0);
	Matrix				amat(2,2);
	vector<double>		f2(2,0);

	if ( verbose )
		cout << "Analyzing reflections:" << endl;

	for ( auto& l1: loc ) {
		sf = l1.length();
		sf_avg += sf;
		sf_std += sf*sf;
		u2 = l1[0]*l1[0];
		v2 = l1[1]*l1[1];
		f2[0] += u2;
		f2[1] += v2;
		amat[0][0] += u2*u2;
		amat[1][1] += v2*v2;
		amat[0][1] += u2*v2;
		amat[1][0] += u2*v2;
	}
	
	f2[0] *= s*s;
	f2[1] *= s*s;
	amat.singular_value_decomposition(f2);
	iso = sqrt((f2[0] + f2[1])/2);
	f2[0] = sqrt(f2[0]);
	f2[1] = sqrt(f2[1]);
	
	sf_avg /= loc.size();
	sf_std = sqrt(sf_std/loc.size() - sf_avg*sf_avg);
	
	if ( verbose ) {
		cout << "Spatial frequency average:      " << sf_avg << " 1/A" << tab << "(" << sf_std << ")" << endl;
		cout << "Resolution average:             " << 1/sf_avg << " A" << endl;
		cout << "Scale in u and v:               " << f2 << endl;
		cout << "Isotropic scale:                " << iso << tab << "(" << loc.size() << ")" << endl << endl;
	}
	
	Vector3<double>		scale(f2[0],f2[1],1);
	
	return scale;
}

