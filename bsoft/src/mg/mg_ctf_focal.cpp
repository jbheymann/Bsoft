/**
@file	mg_ctf_focal.cpp
@brief	Processing focal series
@author 	Bernard Heymann
@date	Created: 20220808
@date	Modified: 20240527
**/

#include "rwimg.h"
#include "mg_ctf.h"
#include "mg_ctf_fit.h"
#include "mg_ctf_focal.h"
#include "simplex.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

/**
@brief 	Calculates an array of CTF parameters.
@param 	cp			source CTF parameters.
@param 	dstart		Starting defocus.
@param 	dend		Ending defocus.
@param 	dinc		Defocus increment (can be negative).
@return vector<CTFparam>	array of CTF parameters.

**/
vector<CTFparam>	ctf_series(CTFparam& cp, double& dstart, double& dend, double& dinc)
{
	if ( dstart > dend && dinc > 0 ) dinc = -dinc;
	if ( dstart < dend && dinc < 0 ) dinc = -dinc;
	
	if ( verbose )
		cout << "Focus series:                     " << dstart << " - " << dend << " @ " << dinc << endl;

	double				def_avg;
	vector<CTFparam>	cpa;

	long			n, np((dend-dstart)/dinc+1);
	if ( np < 1 ) np = 1;
	if ( verbose )
		cout << "#\tDefocus" << endl;
	for ( n=1, def_avg = dstart; n <= np; ++n, def_avg += dinc ) {
		cpa.push_back(cp);
		cpa.back().identifier(n);
		cpa.back().defocus_average(def_avg);
		if ( verbose )
			cout << n << tab << def_avg << endl;
	}

	return cpa;
}

vector<double> 		focus_series(double def_start, double def_inc, long nfoc)
{
	vector<double> 		dfocus;
	for ( long i=0; i<nfoc; ++i )
		dfocus.push_back(def_start+def_inc*i);
	return dfocus;
}

vector<double> 		focus_series(double def_start, double def_end, double def_inc)
{
	if ( def_start > def_end && def_inc > 0 ) def_inc = -def_inc;
	if ( def_start < def_end && def_inc < 0 ) def_inc = -def_inc;

	vector<double> 		dfocus;
	for ( double def = def_start; def <= def_end; def += def_inc )
		dfocus.push_back(def);
	return dfocus;
}

/**
@brief 	Calculates a focal series of aberration images.
@param 	cp				CTF & aberration parameters.
@param 	def_start		Starting defocus.
@param 	def_end			Ending defocus.
@param 	def_inc			Defocus increment (can be negative).
@param 	size			new image size.
@param 	sam				new image pixel size.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param	flag			0=complex, 1=sine, 2=envelope, 4=combined.
@return Bimage*			new complex CTF function image.

	If the third size dimension is one, multiple 2D images will be calculated,
	otherwise mulyiple slices in a single image will be calculated.
	Note: Defocus is positive for underfocus and negative for overfocus.

**/
Bimage*		img_ctf_focal_series(CTFparam& cp, double def_start, double def_end, double def_inc,
				Vector3<long> size, Vector3<double> sam, double hires, double lores, int flag)
{
	vector<double> 		dfocus = focus_series(def_start, def_end, def_inc);
	
	if ( dfocus.size() < size[2] )
		cerr << "Error in img_ctf_focal_series: The series size " << dfocus.size() << " is smaller than the image z size " << size[2] << endl;

	return img_ctf_focal_series(cp, dfocus, size, sam, lores, hires, flag);
}

Bimage*		img_ctf_focal_series(CTFparam& cp, long nimg, Vector3<long> size, Vector3<double> sam, double hires, double lores, int flag)
{
	vector<double> 		dfocus = focus_series(cp.focus_step()*(-nimg/2), cp.focus_step(), nimg);

	return img_ctf_focal_series(cp, dfocus, size, sam, lores, hires, flag);
}

Bimage*		img_ctf_focal_series(CTFparam& cp, vector<double>& dfocus,
				Vector3<long> size, Vector3<double> sam, double hires, double lores, int flag)
{
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;
	
	bool			env(flag&2);
	bool			sine(flag&1);
	bool			comb(flag&4);
	
	long			nimg(dfocus.size());
	if ( size[2] > 1 ) {
		nimg = 1;
		size[2] = dfocus.size();
	}
	
	double			shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	double			shi2(shi*shi), slo2(slo*slo);
	
	Bimage*			p = new Bimage(Float, TComplex, size, nimg);
	if ( sam.volume() > 0 ) p->sampling(sam);
	p->origin(0,0,size[2]/2);
	p->fourier_type(Standard);
	
	long 			i, x, y, z, n;
	double			sx, sy, s, s2, a, dphi, w;
	Complex<double>	cv(1,0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	
	double			def(cp.defocus_average());
	double			fac(0.5*cp.lambda()*fabs(dfocus.back()-dfocus[0]));

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF focal series:" << endl;
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Defocus start, end, increment:  " << def+dfocus[0] << " - " << def+dfocus.back() << " ∆ " << dfocus[1] - dfocus[0] << endl;
		cout << "First sinc node:                " << sqrt(fabs(dfocus.back()-dfocus[0])*cp.lambda()/2) << " A" << endl;
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		cout << "CTF flag:                       " << flag << endl;
		cp.show_envelope();
		cout << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[1]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[1]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}

	for ( i=n=0; n<nimg; ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			cp.defocus_average(def+dfocus[n+z]);
//			if ( verbose & VERB_FULL )
				cout << n+z << tab << cp.defocus_average() << tab << dfocus[n+z] << endl;
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); ++x, ++i ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy;
					if ( s2 >= slo2 && s2 <= shi2 ) {
						s = sqrt(s2);
						a = atan2(sy,sx);
						if ( sine ) {
							cv = cp.aberration_odd_complex(s2, a);
							dphi = cp.calculate_aberration_even(s2, a);
							cv = cv.conj() * sinl(dphi);
						} else if ( comb ) {
							w = 1 - fac*s2;
/*							if ( w < 0 ) {
								cv = cp.calculate_complex(s, a);	// Note: the phase is shifted by pi/2 to make the even terms real
							} else {
								dphi = -M_PI_2 + cp.calculate_aberration(s, a);
								cv = Complex<double>(cosl(dphi),sinl(dphi));
								dphi *= -2;
								cv += Complex<double>(cosl(dphi),sinl(dphi)) * w;
							}*/
							cv = cp.calculate_complex(s, a);	// Note: the phase is shifted by pi/2 to make the even terms real
							if ( w > 0 ) cv += cv.conj() * w;
						} else {
							cv = cp.calculate_complex(s, a);	// Note: the phase is shifted by pi/2 to make the even terms real
						}
						if ( env ) cv *= cp.coherence_envelope(s);
						p->set(i, cv);
					}
				}
			}
		}
	}
	
	cp.defocus_average(def);

	return p;
}
/*
Bimage*		img_ctf_focal_series(CTFparam& cp, Vector3<long> size, Vector3<double> sam, double hires, double lores)
{
	if ( lores < 0 ) lores = 0;
	if ( hires <= 0 ) hires = sam[0];
	if ( lores > 0 && lores < hires ) swap(lores, hires);
	if ( size[2] == 1 ) sam[2] = 1;
	
	double			shi(1/hires);
	double			slo = (lores > 0)? 1/lores: 0;
	double			shi2(shi*shi), slo2(slo*slo);
	
	Bimage*			p = new Bimage(Float, TComplex, size, 1);
	if ( sam.volume() > 0 ) p->sampling(sam);
	p->origin(0,0,size[2]/2);
	p->fourier_type(Standard);
	
	long 			i, x, y, z, n;
	double			sx, sy, s, s2, a, dphi;
	Complex<double>	cv(1,0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	
	double			def(cp.defocus_average());

	if ( verbose & ( VERB_LABEL | VERB_PROCESS ) ) {
		cout << "Calculating a CTF focal series:" << endl;
	}
	if ( verbose & VERB_PROCESS ) {
		cp.show();
		cout << "Defocus average and increment:  " << def << " ∆ " << cp.focus_step() << endl;
		cout << "First sinc node:                " << sqrt(fabs(p->sizeZ()*cp.focus_step())*cp.lambda()/2) << " A" << endl;
		cout << "Resolution range:               " << hires << " - ";
		if ( lores > 0 ) cout << lores << " A" << endl;
		else cout << "inf A" << endl;
		cout << "Frequency range:                " << slo << " - " << shi << " 1/A" << endl;
		cp.show_envelope();
		cout << endl;
		double		rel_size = cp.lambda()*cp.defocus_average()/(sam[0]*sam[1]);
		if ( rel_size > 500 ) {
			cerr << "Warning: The oscillations are too high and create artifacts!" << endl;
			cerr << tab << "Either decrease the defocus below " << 1e-4*500*sam[0]*sam[1]/cp.lambda() << " um" << endl;
			cerr << tab << "or increase the pixel size above " << sqrt(cp.lambda()*cp.defocus_average()/500) << " Å" << endl << endl;
		}
	}

	for ( i=n=0; n<p->images(); ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			cp.defocus_average(def+cp.focus_step()*(z-p->sizeZ()/2));
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); ++x, ++i ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy;
					if ( s2 >= slo2 && s2 <= shi2 ) {
						s = sqrt(s2);
						a = atan2(sy,sx);
						dphi = cp.calculate_aberration_even(s2, a);
						cv = cp.aberration_odd_complex(s2, a);
						cv = cv.conj() * (sinl(dphi) * cp.calc_envelope(s));
						p->add(i, cv);
					}
				}
			}
		}
	}
	
//	p->multiply(1.0L/n);

	cp.defocus_average(def);

	return p;
}
*/

/**
@brief 	Extracts a transverse section from focal series power spectra.
@param 	*p				Focal series.
@param 	which			0=x, 1=y.
@return Bimage*			transverse section.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_extract_section(Bimage* p, int which)
{
	long			i, j, xx, yy, zn, depth;
	Vector3<long>	size(1,p->sizeZ(),1);
	if ( which ) size[0] = p->sizeY();
	else size[0] = p->sizeX();
	if ( p->sizeZ() > 1 ) depth = p->sizeZ();
	else depth = p->images();
	
	Bimage*			ps = new Bimage(Float, TSimple, size, 1);
	if ( which ) ps->sampling(p->sampling(0)[1], p->sampling(0)[2], 1);
	else ps->sampling(p->sampling(0)[0], p->sampling(0)[2], 1);
	ps->origin(ps->size()/2);

	if ( which ) {
		for ( j=zn=0; zn<depth; ++zn ) {
			xx = p->sizeX()/2;
			i = zn*p->sizeY()*p->sizeX() + xx;
			for ( yy=0; yy<p->sizeY(); ++yy, i+=p->sizeX(), ++j )
				ps->set(j, (*p)[i]);
		}
	} else {
		for ( j=zn=0; zn<depth; ++zn ) {
			yy = p->sizeY()/2;
			i = (2*zn+1)*yy*p->sizeX();
			for ( xx=0; xx<p->sizeX(); ++xx, ++i, ++j )
				ps->set(j, (*p)[i]);
		}
	}
	
	return ps;
}

/**
@brief 	Extracts a wedge average from one or focal series power spectra.
@param 	*p				Powere spectra.
@param 	dir				wedge direction (radians).
@param 	width			wedge width (radians).
@return Bimage*			radial average of wedge.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_extract_wedge(Bimage* p, double dir, double width)
{
	long			i, j, xx, yy, zz, nz;
	long			hx(p->sizeX()/2);
	double			dx, dy, dx2, dy2, d, dphi, hw(width/2);
	Vector3<long>	size(p->sizeX(),p->sizeZ(),1);
	
	Bimage*			ps = new Bimage(Float, TSimple, size, 1);
	ps->sampling(p->sampling(0)[1], p->sampling(0)[2], 1);
	ps->origin(ps->size()/2);
	
	vector<long>	num(ps->image_size(), 0);
	
//	cerr << p->image->origin() << endl;
//	cerr << dir << tab << hw << tab << hx << endl;

	for ( i=zz=0; zz<p->sizeZ(); ++zz ) {
		nz = zz*ps->sizeX();
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			dy = double(yy) - p->image->origin()[1];
			dy2 = dy*dy;
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				dx = double(xx) - p->image->origin()[0];
				dx2 = dx*dx;
				d = sqrt(dx2 + dy2);
				if ( d < hx ) {
					dphi = fabs(angle_set_negPI_to_PI(atan2(dy,dx) - dir));
					if ( dphi <= hw || dphi >= M_PI - hw ) {
						if ( dphi > hw ) d = -d;
						j = nz + d + hx;
						if ( j < 0 || j >= ps->image_size() ) {
							cerr << xx << tab << yy << tab << zz << tab << d << tab << dphi << tab << j << endl;
							bexit(-1);
						}
						ps->add(j, (*p)[i]);
						num[j]++;
					}
				}
			}
		}
	}

	for ( i=0; i<ps->image_size(); ++i ) if ( num[i] ) ps->set(i, (*ps)[i]/num[i]);

	return ps;
}

/**
@brief 	Calculates a fit for a defocus value to a section from focal series modified power spectra.
@param 	*p				Section.
@param 	cp				CTF parameters.
@param 	def				Defocus to test for.
@param	&df				Vector of focus offsets.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_ctf_section_fit(Bimage* p, CTFparam cp, double def, vector<double>& df, double hires, double lores)
{
//	double			df = p->sampling(0)[1];
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), v, c, cs(0), rs(0), CC(0);
//	long			h(p->sizeX()/2);
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
//		cp.defocus_average(def - df*(yy - p->sizeY()/2));
		cp.defocus_average(def + df[yy]);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				c = 2*c*c - 1;
				v = (*p)[i];
				CC += c*v;
				cs += c*c;
				rs += v*v;
			}
		}
	}

	CC /= sqrt(cs*rs);

	return CC;
}

/**
@brief 	Calculates a fit for a defocus value to a section from focal series modified power spectra.
@param 	*p				Section.
@param 	cp				CTF parameters.
@param 	def				Defocus to test for.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_ctf_section_fit(Bimage* p, CTFparam cp, double def, double hires, double lores)
{
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), v, c, cs(0), rs(0), CC(0);
//	long			h(p->sizeX()/2);
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		cp.defocus_average(def + cp.focus_step()*(yy - p->sizeY()/2));
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				c = 2*c*c - 1;
				v = (*p)[i];
				CC += c*v;
				cs += c*c;
				rs += v*v;
			}
		}
	}

	CC /= sqrt(cs*rs);

	return CC;
}

/**
@brief 	Determines average defocus from a transverse section of focal series power spectra.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param	&df				Vector of focus offsets.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			Correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_find_section_defocus(Bimage* p, CTFparam& cp, vector<double>& df, double hires, double lores)
{
	double			def, def_best(0), CC(0), CC_best(0);
	double			def_start(cp.defocus_average() - 1000);
	double			def_end(cp.defocus_average() + 1000);
	double			dinc(10);
	if ( df.size() > 1 ) dinc = fabs(df[1]-df[0]);
	
	for ( def=def_start; def<=def_end; def+=dinc ) {
		CC = img_ctf_section_fit(p, cp, def, df, hires, lores);
		if ( verbose & VERB_FULL )
			cout << def << tab << CC << endl;
		if ( CC_best < CC ) {
			CC_best = CC;
			def_best = def;
		}
	}
	
	if ( verbose ) {
		cout << "Best defocus:                   " << def_best << endl;
		cout << "Correlation coefficient:        " << CC_best << endl << endl;
	}
	
	cp.defocus_average(def_best);

	return CC_best;
}

/**
@brief 	Determines average defocus from a transverse section of focal series power spectra.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@return double			Correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
double		img_find_section_defocus(Bimage* p, CTFparam& cp, double hires, double lores)
{
	double			def, def_best(0), CC(0), CC_best(0);
	double			def_start(cp.defocus_average() - 1000);
	double			def_end(cp.defocus_average() + 1000);
	double			dinc(10);
	
	for ( def=def_start; def<=def_end; def+=dinc ) {
		CC = img_ctf_section_fit(p, cp, def, hires, lores);
		if ( verbose & VERB_FULL )
			cout << def << tab << CC << endl;
		if ( CC_best < CC ) {
			CC_best = CC;
			def_best = def;
		}
	}
	
	if ( verbose ) {
		cout << "Best defocus:                   " << def_best << endl;
		cout << "Correlation coefficient:        " << CC_best << endl << endl;
	}
	
	cp.defocus_average(def_best);

	return CC_best;
}


/**
@brief 	Calculates a transverse section of focal series power spectra from CTF parameters.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param	&df				Vector of focus offsets.
@param 	res				High resolution limit.
@return Bimage*			Transverse section.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_ctf_section_calc(Bimage* p, CTFparam& cp, vector<double>& df, double res)
{
	double			def(cp.defocus_average());
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), c;
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smax = (res)? 1/res: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose )
		cout << "Calculating the section CTF for defocus " << def << " A" << endl << endl;

	Bimage*			pc = p->copy();
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		cp.defocus_average(def + df[yy]);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				pc->set(i, c*c);
			} else {
				pc->set(i, 0);
			}
		}
	}
	
	cp.defocus_average(def);

	return pc;
}

/**
@brief 	Calculates a transverse section of focal series power spectra from CTF parameters.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param 	res				High resolution limit.
@return Bimage*			Transverse section.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.

**/
Bimage*		img_ctf_section_calc(Bimage* p, CTFparam& cp, double res)
{
	double			def(cp.defocus_average());
	long			i, xx, yy;
	double			u(0), s2, dphi, a(0), c;
	double			fspace_scale(1.0/p->real_size()[0]);
	double			smax = (res)? 1/res: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose )
		cout << "Calculating the section CTF for defocus " << def << " A" << endl << endl;

	Bimage*			pc = p->copy();
	
	for ( i=yy=0; yy<p->sizeY(); ++yy ) {
		cp.defocus_average(def + cp.focus_step()*(yy-p->sizeY()/2));
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
//			u = (xx < h)? xx: xx-p->sizeX();
			u = xx - p->image->origin()[0];
			u *= fspace_scale;
			s2 = u*u;
			if ( s2 <= smax2 ) {
				dphi = cp.calculate_aberration_even(s2, a);
				c = sinl(dphi);
				pc->set(i, c*c);
			} else {
				pc->set(i, 0);
			}
		}
	}
	
	cp.defocus_average(def);

	return pc;
}

double		focus_cs_amp_section_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			fac = -simp.constant(0)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i, ++j ) {
		v = simp.parameter(0);				// Amp
		v += simp.parameter(1) * x[j];		// Defocus
		v += simp.parameter(2)*x[j]*x[j];	// Cs
		w = exp(fac*x[j]);					// Envelope
		v += x[++j];						// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focus_cs_amp_section_fit_R2(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			fac = -simp.constant(0)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i, ++j ) {
		v = simp.parameter(0);				// Amp
		v += simp.parameter(1) * x[j];		// Defocus
		v += simp.parameter(3)*x[j]*x[j];	// Cs
		w = exp(fac*x[j]);					// Envelope
		v += simp.parameter(2)*x[++j];		// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focus_cs_section_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			amp = simp.constant(0);		// Amplitude contrast
	double			fac = -simp.constant(1)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i, ++j ) {
		v = amp;				// Amp
		v += simp.parameter(0) * x[j];		// Defocus
		v += simp.parameter(1)*x[j]*x[j];	// Cs
		w = exp(fac*x[j]);					// Envelope
		v += x[++j];						// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focus_step_cs_section_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			amp = simp.constant(0);		// Amplitude contrast
	double			fac = -simp.constant(1)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i, ++j ) {
		v = amp;				// Amp
		v += simp.parameter(0) * x[j];		// Defocus
		v += simp.parameter(2)*x[j]*x[j];	// Cs
		w = exp(fac*x[j]);					// Envelope
		v += simp.parameter(1)*x[++j];		// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focus_fit_R(Bsimplex& simp)
{
	long			i, j;
	double			v, v2(0), f2(0), w(1), CC(0);
	double			fac = -simp.constant(2)/4;	// B-factor
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		v = simp.parameter(0) * x[j++];	// Average defocus
		w = exp(fac*x[j]);				// Envelope
		v += simp.constant(0) + simp.constant(1)*x[j]*x[j];	// Amp & Cs
		v += simp.parameter(1) * x[j++];	// Astigmatism
		v += simp.parameter(2) * x[j++];	// Astigmatism
		v += x[j++];					// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focal_aberration_fit_R(Bsimplex& simp)
{
	long			i, j, k;
	double			v, v2(0), f2(0), w(1), CC(0);
	long			iw = simp.constant(0);		// Index for s2
	double			fac = -simp.constant(1)/4;	// B-factor
	double			amp = simp.constant(2);		// Amplitude contrast
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		w = exp(fac*x[j+iw]);			// Envelope
		if ( simp.parameters() < 5 ) v = amp;
		else v = 0;
		for ( k=0; k<simp.parameters(); ++k )
			v += simp.parameter(k) * x[j++];
		v += x[j++];					// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

double		focal_aberration_fit_R2(Bsimplex& simp)
{
	long			i, j, k;
	double			v, v2(0), f2(0), w(1), CC(0);
	long			iw = simp.constant(0);		// Index for s2
	double			fac = -simp.constant(1)/4;	// B-factor
	double			amp = simp.constant(2);		// Amplitude contrast
	vector<double>&	f = simp.dependent_values();
	vector<double>&	x = simp.independent_values();
	
	for ( i=j=0; i<simp.points(); ++i ) {
		w = exp(fac*x[j+iw]);			// Envelope
		if ( simp.parameters() < 6 ) v = amp;
		else v = 0;
		for ( k=0; k<simp.parameters(); ++k )
			v += simp.parameter(k) * x[j++];
//		v += x[j++];					// Relative focus
		v = sinl(v);
		v = 2*v*v - 1;
		v *= w;
		CC += f[i]*v;
		v2 += v*v;
		f2 += f[i]*f[i];
	}
	
	CC /= sqrt(v2*f2);
		
	return 1-CC;
}

/**
@brief 	Fits 3 CTF parameters to a transverse section of focal series power spectra.
@param 	*p				Transverse section.
@param 	&cp				CTF parameters.
@param	&df				Vector of focus offsets.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.

	The x axis is the original x or y axis from the power spectra.
	The y axis is the original series with focus change specified in sampling.
	
	The 3 parameters are the isotropic aberrations:
		constant phase shift (amplitude contrast)
		defocus
		spherical aberration
		

**/
double		img_ctf_fit_section_3p(Bimage* p, CTFparam& cp, vector<double>& df, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting section:" << endl;

//	cout << "origin = " << p->image->origin() << endl;
	
	long			i, xx, yy;
	double			wl(cp.lambda());
	double			dd, s2;
	vector<double>	t, terms, val;
		
	for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
		dd = -M_PI*wl*df[yy];
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s2 = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s2 >= 0.5 ) s2 -= 1;
			s2 /= p->image->sampling()[0];
			s2 *= s2;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				terms.push_back(s2);
				terms.push_back(dd*s2);
				val.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 3, 1, val.size(), terms, val);
	simp.constant(0, Bfactor);
	simp.parameter(0, cp.aberration_weight(0,0));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, cp.aberration_weight(4,0));
	simp.limits(0, cp.aberration_weight(0,0)-0.2, cp.aberration_weight(0,0)+0.2);
	simp.limits(1, cp.aberration_weight(2,0)-50, cp.aberration_weight(2,0)+50);
	simp.limits(2, cp.aberration_weight(4,0)-50, cp.aberration_weight(4,0)+50);

	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_cs_amp_section_fit_R);
	CC = 1 - iR;
	if ( verbose )
		cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_cs_amp_section_fit_R, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(0,0,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(4,0,simp.parameter(2));

	if ( verbose ) {
		cout << "Amplitude phase:                " << cp.amp_shift()*180.0/M_PI << " degrees" << endl;
		cout << "Defocus average:                " << cp.defocus_average()*1e-4 << " um" << endl;
		cout << "Cs:                             " << cp.Cs()*1e-7 << " mm" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

double		img_ctf_fit_section_3p_df(Bimage* p, CTFparam& cp, vector<double>& df, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting section:" << endl;

//	cout << "origin = " << p->image->origin() << endl;
	
	long			i, xx, yy;
	double			wl(cp.lambda());
	double			dd, s2;
	vector<double>	t, terms, val;
		
	for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
//		dd = -M_PI*wl*df[yy];
		dd = -M_PI*wl*(yy-p->sizeY()/2);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s2 = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s2 >= 0.5 ) s2 -= 1;
			s2 /= p->image->sampling()[0];
			s2 *= s2;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				terms.push_back(s2);
				terms.push_back(dd*s2);
				val.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 4, 1, val.size(), terms, val);
	simp.constant(0, Bfactor);
	simp.parameter(0, cp.aberration_weight(0,0));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, df[1] - df[0]);
	simp.parameter(3, cp.aberration_weight(4,0));
	simp.limits(0, cp.aberration_weight(0,0)-0.2, cp.aberration_weight(0,0)+0.2);
	simp.limits(1, cp.aberration_weight(2,0)-50, cp.aberration_weight(2,0)+50);
	simp.limits(2, 0.8*(df[1] - df[0]), 1.2*(df[1] - df[0]));
	simp.limits(3, cp.aberration_weight(4,0)-50, cp.aberration_weight(4,0)+50);

	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_cs_amp_section_fit_R2);
	CC = 1 - iR;
	if ( verbose )
		cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_cs_amp_section_fit_R2, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(0,0,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(4,0,simp.parameter(3));

	if ( verbose ) {
		cout << "Amplitude phase:                " << cp.amp_shift()*180.0/M_PI << " degrees" << endl;
		cout << "Defocus average:                " << cp.defocus_average()*1e-4 << " um" << endl;
		cout << "Focus step size:                " << simp.parameter(2) << " A" << endl;
		cout << "Cs:                             " << cp.Cs()*1e-7 << " mm" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

double		img_ctf_fit_section_df(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting section:" << endl;

//	cout << "origin = " << p->image->origin() << endl;
	
	long			i, xx, yy;
	double			wl(cp.lambda());
	double			dd, s2;
	vector<double>	t, terms, val;
		
	for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
		dd = -M_PI*wl*(yy-p->sizeY()/2);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s2 = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s2 >= 0.5 ) s2 -= 1;
			s2 /= p->image->sampling()[0];
			s2 *= s2;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				terms.push_back(s2);
				terms.push_back(dd*s2);
				val.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 3, 2, val.size(), terms, val);
	simp.constant(0, cp.aberration_weight(0,0));
	simp.constant(1, Bfactor);
	simp.parameter(0, cp.aberration_weight(2,0));
	simp.parameter(1, cp.focus_step());
	simp.parameter(2, cp.aberration_weight(4,0));
	simp.limits(0, cp.aberration_weight(2,0)-50, cp.aberration_weight(2,0)+50);
	simp.limits(1, 0.7*cp.focus_step(), 1.3*cp.focus_step());
	simp.limits(2, cp.aberration_weight(4,0)-50, cp.aberration_weight(4,0)+50);

	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_step_cs_section_fit_R);
	CC = 1 - iR;
	if ( verbose )
		cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_step_cs_section_fit_R, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(2,0,simp.parameter(0));
	cp.aberration_weight(4,0,simp.parameter(2));
	cp.focus_step(simp.parameter(1));
		
	if ( verbose ) {
		cout << "Amplitude phase:                " << cp.amp_shift()*180.0/M_PI << " degrees" << endl;
		cout << "Defocus average:                " << cp.defocus_average()*1e-4 << " um" << endl;
		cout << "Focus step size:                " << cp.focus_step() << " A" << endl;
		cout << "Cs:                             " << cp.Cs()*1e-7 << " mm" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

double		img_ctf_fit_section(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting section:" << endl;

//	cout << "origin = " << p->image->origin() << endl;
	
	long			i, xx, yy;
	double			wl(cp.lambda());
	double			dd, s2;
	vector<double>	t, terms, val;
		
	for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
		dd = -M_PI*wl*cp.focus_step()*(yy-p->sizeY()/2);
		for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
			s2 = (double(xx) - p->image->origin()[0])/p->sizeX();
			if ( s2 >= 0.5 ) s2 -= 1;
			s2 /= p->image->sampling()[0];
			s2 *= s2;
			if ( s2 >= smin2 && s2 <= smax2 ) {
				terms.push_back(s2);
				terms.push_back(dd*s2);
				val.push_back((*p)[i]);
			}
		}
	}

	Bsimplex			simp(2, 3, 2, val.size(), terms, val);
	simp.constant(0, cp.aberration_weight(0,0));
	simp.constant(1, Bfactor);
	simp.parameter(0, cp.aberration_weight(2,0));
	simp.parameter(1, cp.aberration_weight(4,0));
	simp.limits(0, cp.aberration_weight(2,0)-50, cp.aberration_weight(2,0)+50);
	simp.limits(1, cp.aberration_weight(4,0)-50, cp.aberration_weight(4,0)+50);

	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_cs_section_fit_R);
	CC = 1 - iR;
	if ( verbose )
		cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_cs_section_fit_R, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(2,0,simp.parameter(0));
	cp.aberration_weight(4,0,simp.parameter(1));
	
	if ( verbose ) {
		cout << "Amplitude phase:                " << cp.amp_shift()*180.0/M_PI << " degrees" << endl;
		cout << "Defocus average:                " << cp.defocus_average()*1e-4 << " um" << endl;
		cout << "Cs:                             " << cp.Cs()*1e-7 << " mm" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits defocus and astigmatism to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param	&df				Vector of focus offsets.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.

	The constant phase shift and spherical aberration are fixed at initial values.
	
**/
double		img_ctf_fit_astigmatism(Bimage* p, CTFparam& cp, vector<double>& df, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting astigmatism:" << endl;

	map<pair<long,long>,double>		wa;
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);

//	for ( auto w: wa )
//		cout << w.first.first << tab << w.first.second << tab << w.second << endl;
//	cout << "origin = " << p->image->origin() << endl;
//	cout << "focal step = " << p->image->sampling()[2] << endl;
	
	long			nt(wa.size()), nc(3);
	
	long			i, j, xx, yy, zn, depth;
	double			s2, dd;
	double			wl(cp.lambda());
	Vector3<double>	s;
	vector<double>	t, terms, val;
	
	if ( p->sizeZ() > 1 ) depth = p->sizeZ();
	else depth = p->images();
		
	for ( i=0, zn=0; zn<depth; ++zn ) {
		dd = -M_PI*wl*df[zn];
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
			if ( s[1] >= 0.5 ) s[1] -= 1;
			s[1] /= p->image->sampling()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
				if ( s[0] >= 0.5 ) s[0] -= 1;
				s[0] /= p->image->sampling()[0];
				s2 = s.length2();
				if ( s2 >= smin2 && s2 <= smax2 ) {
					t = aberration_terms(wa, s[0], s[1]);
					for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
					terms.push_back(dd*s2);
					val.push_back((*p)[i]);
				}
			}
		}
	}

//	cout << "independent variable size = " << terms.size() << endl;
//	cout << "dependent variable size =   " << val.size() << endl;
/*
	for ( i=j=0; i<val.size(); ++i ) {
		cout << i;
		for ( long k=0; k<4; ++k ) cout << tab << terms[j++];
		cout << tab << val[i] << endl;
	}
*/
	double				dlim(100);
	Bsimplex			simp(nt+1, nt, nc, val.size(), terms, val);
	simp.constant(0, cp.aberration_weight(0,0));
	simp.constant(1, cp.aberration_weight(4,0));
	simp.constant(2, Bfactor);
	simp.parameter(0, cp.aberration_weight(2,-2));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, cp.aberration_weight(2,2));
	simp.limits(0, -dlim, dlim);
	simp.limits(1, cp.aberration_weight(2,0)-10*dlim, cp.aberration_weight(2,0)+10*dlim);
	simp.limits(2, -dlim, dlim);

	if ( verbose ) {
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
		j = 0;
		for ( auto w: wa ) {
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
			j++;
		}
		cout << endl;
	}
	
	simp.show();
	
	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_fit_R);
	CC = 1 - iR;
	cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_fit_R, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(2,-2,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(2,2,simp.parameter(2));

	if ( verbose ) {
		cout << "Defocus average:                " << cp.defocus_average() << " A" << endl;
		cout << "Defocus deviation:              " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180.0/M_PI << " °" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

double		img_ctf_fit_astigmatism(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting astigmatism:" << endl;

	map<pair<long,long>,double>		wa;
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);

//	for ( auto w: wa )
//		cout << w.first.first << tab << w.first.second << tab << w.second << endl;
//	cout << "origin = " << p->image->origin() << endl;
//	cout << "focal step = " << p->image->sampling()[2] << endl;
	
	long			nt(wa.size()), nc(3);
	
	long			i, j, xx, yy, zn, depth;
	double			s2, dd;
	double			wl(cp.lambda());
	Vector3<double>	s;
	vector<double>	t, terms, val;
	
	if ( p->sizeZ() > 1 ) depth = p->sizeZ();
	else depth = p->images();
		
	for ( i=0, zn=0; zn<depth; ++zn ) {
//		dd = -M_PI*wl*df[zn];
		dd = -M_PI*wl*cp.focus_step()*(zn-depth/2);
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
			if ( s[1] >= 0.5 ) s[1] -= 1;
			s[1] /= p->image->sampling()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
				if ( s[0] >= 0.5 ) s[0] -= 1;
				s[0] /= p->image->sampling()[0];
				s2 = s.length2();
				if ( s2 >= smin2 && s2 <= smax2 ) {
					t = aberration_terms(wa, s[0], s[1]);
					for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
					terms.push_back(dd*s2);
					val.push_back((*p)[i]);
				}
			}
		}
	}

//	cout << "independent variable size = " << terms.size() << endl;
//	cout << "dependent variable size =   " << val.size() << endl;
/*
	for ( i=j=0; i<val.size(); ++i ) {
		cout << i;
		for ( long k=0; k<4; ++k ) cout << tab << terms[j++];
		cout << tab << val[i] << endl;
	}
*/
	double				dlim(100);
	Bsimplex			simp(nt+1, nt, nc, val.size(), terms, val);
	simp.constant(0, cp.aberration_weight(0,0));
	simp.constant(1, cp.aberration_weight(4,0));
	simp.constant(2, Bfactor);
	simp.parameter(0, cp.aberration_weight(2,-2));
	simp.parameter(1, cp.aberration_weight(2,0));
	simp.parameter(2, cp.aberration_weight(2,2));
	simp.limits(0, -dlim, dlim);
	simp.limits(1, cp.aberration_weight(2,0)-10*dlim, cp.aberration_weight(2,0)+10*dlim);
	simp.limits(2, -dlim, dlim);

	if ( verbose ) {
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
		j = 0;
		for ( auto w: wa ) {
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
			j++;
		}
		cout << endl;
	}
	
	simp.show();
	
	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focus_fit_R);
	CC = 1 - iR;
	cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focus_fit_R, maxiter/10);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
	}

	cp.aberration_weight(2,-2,simp.parameter(0));
	cp.aberration_weight(2,0,simp.parameter(1));
	cp.aberration_weight(2,2,simp.parameter(2));

	if ( verbose ) {
		cout << "Defocus average:                " << cp.defocus_average() << " A" << endl;
		cout << "Defocus deviation:              " << cp.defocus_deviation() << " A" << endl;
		cout << "Astigmatism angle:              " << cp.astigmatism_angle()*180.0/M_PI << " °" << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param	&focus_step		Focus step size.
@param 	wl				Electron wavelength.
@param 	&wa				Aberration weights (replaced).
@param 	wd				Aberration weight limits for fitting.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	amp				amplitude contrast.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_aberration(Bimage* p, double& focus_step, double wl, map<pair<long,long>,double>& wa,
				map<pair<long,long>,double> wd, double hires, double lores,
				double Bfactor, double amp, long maxiter)
{
	double			smin = (lores)? 1/lores: 1e-3;
	double			smin2(smin*smin);
	double			smax = (hires)? 1/hires: 0.5/p->sampling(0)[0];
	double			smax2(smax*smax);
	
	if ( verbose & VERB_FULL )
		cout << "Fitting aberrations:" << endl;

	long			nt(wa.size());
	
	long			i, j, xx, yy, zn, depth;
	double			s2, dd;
	Vector3<double>	s;
	vector<double>	t, terms, val;
		
	if ( p->sizeZ() > 1 ) depth = p->sizeZ();
	else depth = p->images();
		
	for ( i=0, zn=0; zn<depth; ++zn ) {
		dd = -M_PI*wl*(zn-depth/2);
		for ( yy=0; yy<p->sizeY(); ++yy ) {
			s[1] = (double(yy) - p->image->origin()[1])/p->sizeY();
			if ( s[1] >= 0.5 ) s[1] -= 1;
			s[1] /= p->image->sampling()[1];
			for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
				s[0] = (double(xx) - p->image->origin()[0])/p->sizeX();
				if ( s[0] >= 0.5 ) s[0] -= 1;
				s[0] /= p->image->sampling()[0];
				s2 = s.length2();
				if ( s2 >= smin2 && s2 <= smax2 ) {
					t = aberration_terms(wa, s[0], s[1]);
					for ( j=0; j<nt; ++j ) terms.push_back(t[j]);
					terms.push_back(dd*s2);
					val.push_back((*p)[i]);
				}
			}
		}
	}

	Bsimplex			simp(nt+1, nt+1, 3, val.size(), terms, val);
	simp.constant(0, 1);							// Index of s2 for weighting function
	if ( wa.size() > 4 ) simp.constant(0, 2);		// Index of s2 for weighting function
	simp.constant(1, Bfactor);						// B-factor of weighting function
	simp.constant(2, amp);	// Constant amplitude contrast
	
	i = 0;
	for ( auto w: wa )
		simp.parameter(i++, w.second);
	simp.parameter(i, focus_step);

	i = 0;
	for ( auto w: wd ) {
		simp.limits(i, simp.parameter(i)-w.second, simp.parameter(i)+w.second);
		i++;
	}
	simp.limits(i, focus_step*0.7, focus_step*1.3);

	if ( verbose ) {
		cout << "Aberration parameters:\nn\tm\tw\twmin\twmax" << endl;
		j = 0;
		for ( auto w: wa ) {
			cout << w.first.first << tab << w.first.second << tab <<
				simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
			j++;
		}
		cout << "∆f" << tab << "-" << tab <<
			simp.parameter(j) << tab << simp.limit_low(j) << tab << simp.limit_high(j) << endl;
		cout << endl;
	}
	
	simp.show();
	
	double			R(0), iR(0), CC(0);
	
	iR = simp.R(focal_aberration_fit_R2);
	CC = 1 - iR;
	if ( verbose )
		cout << "Correlation coefficient (start): " << CC << endl << endl;

	R = simp.run(maxiter, 0.01, focal_aberration_fit_R2, maxiter/10);

	i = 0;
	for ( auto& w: wa )
		w.second = simp.parameter(i++);
	focus_step = simp.parameter(i);

	CC = 1 - R;

	string		ws = "[" + concatenate(simp.parameter_vector()) + "]";

	if ( verbose & VERB_FULL ) {
		cout << ws << endl;
		cout << "iR: " << iR << endl;
		cout << "R: " << R << endl;
		cout << "Correlation coefficient:        " << CC << endl << endl;
	}

	return CC;
}

/**
@brief 	Fits 4 even aberration parameters to focal series power spectra, keeping the amplitude contrast constant.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_even4(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double							focus_step(cp.focus_step());
	map<pair<long,long>,double>		wa, wd;
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);
	wa[{4,0}] = cp.aberration_weight(4,0);
	wd[{2,-2}] = 100;
	wd[{2,0}] = 1000;
	wd[{2,2}] = 100;
	wd[{4,0}] = 100;

	double		CC = img_ctf_fit_aberration(p, focus_step, cp.lambda(), wa, wd, hires, lores, Bfactor, cp.aberration_weight(0,0), maxiter);
	
	cp.aberration_weights(wa);
	cp.focus_step(focus_step);

	return CC;
}

/**
@brief 	Fits 5 even aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_even5(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double							focus_step(cp.focus_step());
	map<pair<long,long>,double>		wa, wd;
	wa[{0,0}] = cp.aberration_weight(0,0);
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);
	wa[{4,0}] = cp.aberration_weight(4,0);
	wd[{0,0}] = 0.1;
	wd[{2,-2}] = 100;
	wd[{2,0}] = 1000;
	wd[{2,2}] = 100;
	wd[{4,0}] = 100;

	double		CC = img_ctf_fit_aberration(p, focus_step, cp.lambda(), wa, wd, hires, lores, Bfactor, 0, maxiter);
	
	cp.aberration_weights(wa);

	return CC;
}

/**
@brief 	Fits 9 even aberration parameters to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
**/
double		img_ctf_fit_even9(Bimage* p, CTFparam& cp, double hires, double lores, double Bfactor, long maxiter)
{
	double							focus_step(cp.focus_step());
	map<pair<long,long>,double>		wa, wd;
	wa[{0,0}] = cp.aberration_weight(0,0);
	wa[{2,-2}] = cp.aberration_weight(2,-2);
	wa[{2,0}] = cp.aberration_weight(2,0);
	wa[{2,2}] = cp.aberration_weight(2,2);
	wa[{4,-4}] = cp.aberration_weight(4,-4);
	wa[{4,-2}] = cp.aberration_weight(4,-2);
	wa[{4,0}] = cp.aberration_weight(4,0);
	wa[{4,2}] = cp.aberration_weight(4,2);
	wa[{4,4}] = cp.aberration_weight(4,4);
	wd[{0,0}] = 0.1;
	wd[{2,-2}] = 100;
	wd[{2,0}] = 1000;
	wd[{2,2}] = 100;
	wd[{4,-4}] = 10;
	wd[{4,-2}] = 10;
	wd[{4,0}] = 100;
	wd[{4,2}] = 10;
	wd[{4,4}] = 10;

	double		CC = img_ctf_fit_aberration(p, focus_step, cp.lambda(), wa, wd, hires, lores, Bfactor, 0, maxiter);
	
	cp.aberration_weights(wa);

	return CC;
}

/*
@brief 	Fits the CTF to focal series power spectra.
@param 	*p				Focal series.
@param 	&cp				CTF parameters.
@param	&dfocus			Vector of focus offsets.
@param 	hires			High resolution limit.
@param 	lores			Low resolution limit.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			Maximum number of iterations.
@return double			Correlation coefficient.
	
	The image should be 3D with the third sampling interval the
	change in focus per 2D image.
*/
/*Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, vector<double>& dfocus, double hires, double lores,
				double tmax, double Bfactor, long maxiter)
{
	p->set_hi_lo_resolution(hires, lores);
	
	Bimage*			pf = img_ctf_fit_prepare(p, tmax, Bfactor);
	pf->sampling(p->sampling(0));

	if ( verbose ) {
		cout << "Fitting focal series power spectra:" << endl;
		cout << "Focus change:                   " << p->sampling(0)[2] << " A" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "Truncation maximum:             " << tmax << endl;
		cout << "B-factor:                       " << Bfactor << " A2" << endl;
		cp.show();
	}
	
	Bimage*			ps = NULL;

//	if ( pf->sizeZ() > 1 || pf->images() > 1 )
//		ps = img_extract_section(pf, 0);
//	else
		ps = img_extract_wedge(pf, 0, M_PI/20);

	double			Rx = img_find_section_defocus(ps, cp, dfocus, hires, lores);

	Rx = img_ctf_fit_section3(ps, cp, dfocus, hires, lores, Bfactor, maxiter);

	double			dx = cp.defocus_average();
	double			ddfx = dfocus[1] - dfocus[0];

	Bimage*			psc = img_ctf_section_calc(ps, cp, dfocus, hires);

	ps->replace_half(psc);
	
	write_img("px.grd", ps, 0);
	
	delete psc;
	delete ps;
	
//	if ( pf->sizeZ() > 1 )
//		ps = img_extract_section(pf, 1);
//	else
		ps = img_extract_wedge(pf, M_PI_2, M_PI/20);

	double			Ry = img_find_section_defocus(ps, cp, dfocus, hires, lores);
	
	Ry = img_ctf_fit_section3(ps, cp, dfocus, hires, lores, Bfactor, maxiter);

	double			dy = cp.defocus_average();
	double			ddfy = dfocus[1] - dfocus[0];

	psc = img_ctf_section_calc(ps, cp, dfocus, hires);

	ps->replace_half(psc);
	
	write_img("py.grd", ps, 0);
	
	delete psc;
	delete ps;
	
	if ( verbose ) {
		cout << "Best defocus in x:              " << dx << " A (" << Rx << ")" << endl;
		cout << "Best defocus in y:              " << dy << " A (" << Ry << ")" << endl;
		cout << "Best focus step in x:           " << ddfx << " A" << endl;
		cout << "Best focus step in y:           " << ddfy << " A" << endl;
		cout << "Average defocus:                " << (dx+dy)/2 << " A" << endl;
		cout << "Average focus step:             " << (ddfx+ddfy)/2 << " A" << endl << endl;
	}
	
	long			i, nfoc(dfocus.size());
	cp.defocus_average(0.5*(dx+dy));
	for ( i=0; i<nfoc; ++i )
		dfocus[i] = 0.5*(ddfx + ddfy)*(i-nfoc/2);
	
//	cout << "amp = " << cp.aberration_weight(0,0) << endl;
//	bexit(0);

	img_ctf_fit_astigmatism(pf, cp, dfocus, hires, lores, Bfactor, maxiter);
	
//	img_ctf_fit_even4(pf, cp, dfocus, hires, lores, Bfactor, maxiter);
	cp.focus_step(dfocus[1]-dfocus[0]);
	img_ctf_fit_even4(pf, cp, hires, lores, Bfactor, maxiter);
	for ( i=0; i<nfoc; ++i )
		dfocus[i] = cp.focus_step()*(i-nfoc/2);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	img_ctf_fit_even5(pf, cp, dfocus, hires, lores, Bfactor, maxiter);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	img_ctf_fit_even9(pf, cp, dfocus, hires, lores, Bfactor, maxiter);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	ps = img_ctf_focal_series(cp, dfocus, p->size(), p->image->sampling(), 0, hires/2);

	ps->complex_to_intensities();
	
	ps->shift_wrap(Vector3<double>(p->sizeX()/2, p->sizeY()/2, 0));
	ps->multiply(2);
	ps->add(-1);
	
	pf->replace_half(ps);
	
	delete ps;
	
	return pf;
}
*/

/**
@brief 	Fits the CTF to focal series power spectra.
@param	*p				focal series.
@param 	&cp				CTF parameters.
@param 	hires			high resolution limit.
@param 	lores			low resolution limit.
@param	tmax			truncation maximum.
@param 	Bfactor			B-factor for weighting.
@param 	maxiter			maximum number of iterations.
@return double			correlation coefficient.
	
	The image should be 3D with 2D power spectra packed into the slices and
	with the third sampling interval the change in focus per 2D image.
**/
Bimage*		img_ctf_focal_fit(Bimage* p, CTFparam& cp, double hires, double lores,
				double tmax, double Bfactor, long maxiter)
{
	p->set_hi_lo_resolution(hires, lores);
	
	Bimage*			pf = img_ctf_fit_prepare(p, tmax, Bfactor);
	pf->sampling(p->sampling(0));

	if ( verbose ) {
		cout << "Fitting focal series power spectra:" << endl;
		cout << "Focus step:                     " << cp.focus_step() << " A" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "Truncation maximum:             " << tmax << endl;
		cout << "B-factor:                       " << Bfactor << " A2" << endl;
		cp.show();
	}
	
	Bimage*			ps = NULL;

//	if ( pf->sizeZ() > 1 || pf->images() > 1 )
//		ps = img_extract_section(pf, 0);
//	else
		ps = img_extract_wedge(pf, 0, M_PI/20);

	double			Rx = img_find_section_defocus(ps, cp, hires, lores);

	Rx = img_ctf_fit_section_df(ps, cp, hires, lores, Bfactor, maxiter);

	double			dx = cp.defocus_average();
	double			dfx = cp.focus_step();

	Bimage*			psc = img_ctf_section_calc(ps, cp, hires);

	ps->replace_half(psc);
	
	write_img("px.grd", ps, 0);
	
	delete psc;
	delete ps;
	
//	if ( pf->sizeZ() > 1 )
//		ps = img_extract_section(pf, 1);
//	else
		ps = img_extract_wedge(pf, M_PI_2, M_PI/20);

	double			Ry = img_find_section_defocus(ps, cp, hires, lores);
	
	Ry = img_ctf_fit_section_df(ps, cp, hires, lores, Bfactor, maxiter);

	double			dy = cp.defocus_average();
	double			dfy = cp.focus_step();

	psc = img_ctf_section_calc(ps, cp, hires);

	ps->replace_half(psc);
	
	write_img("py.grd", ps, 0);
	
	delete psc;
	delete ps;
	
	if ( verbose ) {
		cout << "Best defocus in x:              " << dx << " A (" << Rx << ")" << endl;
		cout << "Best defocus in y:              " << dy << " A (" << Ry << ")" << endl;
		cout << "Best focus step in x:           " << dfx << " A" << endl;
		cout << "Best focus step in y:           " << dfy << " A" << endl;
		cout << "Average defocus:                " << (dx+dy)/2 << " A" << endl;
		cout << "Average focus step:             " << (dfx+dfy)/2 << " A" << endl;
	}
	
	cp.defocus_average(0.5*(dx+dy));
	cp.focus_step(0.5*(dfx+dfy));
	
//	cout << "amp = " << cp.aberration_weight(0,0) << endl;
//	bexit(0);

//	img_ctf_fit_astigmatism(pf, cp, hires, lores, Bfactor, maxiter);
	
	img_ctf_fit_even4(pf, cp, hires, lores, Bfactor, maxiter);

	if ( verbose ) {
		cout << "CTF parameters:" << endl;
		cp.show();
	}

	ps = img_ctf_focal_series(cp, p->images()*p->sizeZ(), p->size(), p->image->sampling(), 0, hires/2, 1);
	
	ps->shift_wrap(Vector3<double>(p->sizeX()/2, p->sizeY()/2, 0));

	ps->complex_to_real();
	ps->square();
	ps->multiply(2);
	ps->add(-1);
	
	pf->replace_half(ps);
	
	delete ps;
	
	return pf;
}

/**
@brief 	Weighs a 3D transform with a focal coherent sphere corresponding to a given acceleration voltage.
@param 	*p				3D Fourier transform.
@param 	volt			Acceleration voltage (V).
@return Bimage*			Image with data from the sphere.
	
**/
//Bimage*		img_fspace_weigh_sphere(Bimage* p, double volt)
int			img_fspace_weigh_sphere(Bimage* p, double volt)
{
	if ( p->compound_type() != TComplex )
		p->fft();

	long			i, xx, yy, zz, nn, xh(p->sizeX()/2), yh(p->sizeY()/2), zh(p->sizeZ()/2);
	double			u, v, w, ws, s2, f, t(p->real_size()[2]);
	double			wl = electron_wavelength(volt);
	
	if ( verbose )
		cout << "Weighing with a focal coherent sphere for a wavelength of " << wl << " A" << endl << endl;
	
/*	Bimage*			ps = new Bimage(Float, TSimple, p->size(), p->images());
	ps->sampling(p->image->sampling());
	ps->origin(p->image->origin());
	ps->fourier_type(Standard);
*/
	for ( nn=i=0; nn<p->images(); ++nn ) {
		for ( zz=0; zz<p->sizeZ(); ++zz ) {
			w = ( zz < zh )? zz/p->real_size()[2]: (zz-p->sizeZ())/p->real_size()[2];
			for ( yy=0; yy<p->sizeY(); ++yy ) {
				v = ( yy < yh )? yy/p->real_size()[1]: (yy-p->sizeY())/p->real_size()[1];
				for ( xx=0; xx<p->sizeX(); ++xx, ++i ) {
					u = ( xx < xh )? xx/p->real_size()[0]: (xx-p->sizeX())/p->real_size()[0];
					s2 = u*u + v*v;
					ws = (wl/2)*s2;
					f = 0.5*(sinc((w-ws)*t) + sinc((w+ws)*t));
					p->set(i, p->complex(i) * f);
				}
			}
		}
	}
	
//	return ps;
	return 0;
}

/**
@brief 	Extracts a sphere corresponding to a given acceleration voltage.
@param 	*p				Focal series 2D Fourier transforms.
@param 	cp				CTF parameters.
@param	&dfocus			Vector of focus offsets.
@return Bimage*			Image with data from the sphere.

	The transforms are expected in standard format with the origin at {0,0}.
	
**/
Bimage*		img_fspace_extract_sphere(Bimage* p, CTFparam cp, vector<double>& dfocus)
{
	if ( dfocus.size() != p->images()*p->sizeZ() ) {
		cerr << "Error: The number of focal values (" << dfocus.size() <<
			") must equal the number of images (" << p->images() <<
			") or the number of slices (" << p->sizeZ() << ")" << endl;
		return NULL;
	}

	long			i, j, xx, yy, zz, nn;
//	long			xy(p->sizeX()*p->sizeY());
	double			sx, sy, s2, a, dphi;
	double			wl(cp.lambda());
	double			def(cp.defocus_average());
	Complex<double>	cv;
	
	if ( verbose ) {
		cout << "Extracting a focal coherent wave:" << endl;
		cout << "Wavelength:                       " << wl << " A" << endl;
		cout << "Focus step size:                  " << dfocus[1] - dfocus[0] << " A" << endl << endl;
	}
	
//	Bimage*			ps = new Bimage(Float, TComplex, p->sizeX(), p->sizeY(), p->images()*p->sizeZ(), 1);
	Bimage*			ps = new Bimage(Float, TComplex, p->sizeX(), p->sizeY(), 1, 1);
	ps->sampling(p->image->sampling());
	ps->origin(p->image->origin());
	ps->fourier_type(Standard);
	
	for ( j=nn=0; nn<p->images(); ++nn ) {
		for ( zz=0; zz<p->sizeZ(); ++zz ) {
//			dd = M_PI*wl*dfocus[zz+nn];
			cp.defocus_average(def+dfocus[zz+nn]);
			if ( verbose )
				cout << zz+nn << tab << cp.defocus_average() << endl;
//			for ( i=(nn+zz)*xy, yy=0; yy<p->sizeY(); ++yy ) {
			for ( i=0, yy=0; yy<p->sizeY(); ++yy ) {
				sy = (yy < p->sizeY()/2 )? yy: yy - p->sizeY();
				sy /= p->real_size()[1];
				for ( xx=0; xx<p->sizeX(); ++xx, ++i, ++j ) {
					sx = (xx < p->sizeX()/2 )? xx: xx - p->sizeX();
					sx /= p->real_size()[0];
					s2 = sx*sx + sy*sy;
//					s = sqrt(s2);
					a = atan2(sy,sx);
					dphi = cp.calculate_aberration_even(s2, a);
					cv = cp.aberration_odd_complex(s2, a);
//					cv = cv.conj() * (sinl(dphi) * cp.calc_envelope(s));
					cv = cv.conj() * sinl(dphi);
					cv = p->complex(j) * cv;
					ps->add(i, cv);
				}
			}
		}
	}
	
	return ps;
}

/**
@brief 	Applies CTF to a focal series .
@param 	*p				Focal series (images or transforms).
@param 	cp				CTF parameters with average defocus.
@param	&dfocus			Vector of focus offsets.
@return Bimage*			Resultant complex image.
	
**/
Bimage*		img_ctf_apply_to_focal_series(Bimage* p, CTFparam cp, vector<double>& dfocus)
{
	if ( p->compound_type() != TComplex )
		p->fftxy();

	double			hires(p->sampling(0)[0]);
	
	Bimage*			pr = img_ctf_focal_series(cp, dfocus,
						p->size(), p->image->sampling(), 0, hires, 1);

	pr->complex_product(p);
//	pr->complex_conjugate_product(p);

	return pr;
}

double		img_focal_wave_intensity(Bimage* p, double wl, double focus_step, double hires, double lores)
{
	if ( p->compound_type() != TComplex )
		p->fftxy();

	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	if ( lores < hires ) lores = p->real_size()[0];
	
	long			i, j, x, y, z, n, depth;
	double			slo2(1/(lores*lores)), shi2(1/(hires*hires));
	double			sx, sy, sy2, s2;
	double			dd, pwr, fcw(0);
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	Complex<double>	cv;
	if ( p->sizeZ() > 1 ) depth = p->sizeZ();
	else depth = p->images();
	
	Bimage*			pr = new Bimage(Float, TComplex, p->sizeX(), p->sizeY(), 1, 1);
	pr->fourier_type(Standard);
	pr->sampling(p->sampling(0));

	for ( i=n=0; n<p->images(); ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			dd = -M_PI*wl*focus_step*(z+n-depth/2);
			for ( j=y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				sy2 = sy*sy;
				for ( x=0; x<p->sizeX(); ++x, ++i, ++j ) {
					if ( sy2 <= shi2 ) {
						sx = x;
						if ( x > h[0] ) sx -= p->sizeX();
						sx *= freq_scale[0];
						s2 = sx*sx + sy*sy;
						if ( s2 >= slo2 && s2 <= shi2 ) {
							cv = p->complex(i);
							cv.shift_phi(-dd*s2);
							pr->add(j, cv);
						}
					}
				}
			}
		}
	}

	for ( i=n=0; i<pr->sizeX()*pr->sizeY(); ++i ) {
		pwr = pr->complex(i).power();
		if ( pwr ) {
			fcw += pwr;
			n++;
		}
	}
	
	if ( n ) fcw /= n;
	
	delete pr;

	return fcw;
}

/**
@brief 	Determines the focal step size from the focal coherent wave.
@param 	*p				focal series (images or transforms).
@param 	wl				wavelength (angstrom).
@param	&focus_step		focus step size to be refined (angstrom).
@param	focus_inc		initial focus step size (angstrom).
@param	hires			high resolution limit.
@param	lores			low resolution limit.
@param	linear			linear search extent.
@return double			focal coherent wave intensity.
	
**/
double		img_ctf_refine_focus_step(Bimage* p, double wl, double& focus_step, double focus_inc, double hires, double lores, double linear)
{
	if ( p->compound_type() != TComplex )
		p->fftxy();

	if ( lores < 0.1 ) lores = p->real_size()[0];
	if ( hires > lores ) swap(hires, lores);
	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	
	double			f(focus_step), fcw(0), fcw_best(0);
	double			flo(focus_step-linear*focus_inc), fhi(focus_step+linear*focus_inc);
//	if ( fabs(flo) < fhi && flo < 0 ) flo = 0;
//	if ( fasb(fhi) < fabs(flo) && fhi < 0 )
	
	if ( verbose ) {
		cout << "Determining the focus step size:" << endl;
		cout << "Resolution limits:              " << hires << " - " << lores << " A" << endl;
		cout << "Initial focus step size:        " << focus_step << " A" << endl;
		cout << "Focus increment:                " << focus_inc << " A" << endl;
		if ( linear )
			cout << "Linear search:                  " << flo << " - " << fhi << endl;
		else
			cout << "Convergent search" << endl;
	}
	
	if ( verbose )
		cout << "Step\tIntensity" << endl;
		
	if ( linear ) {
		for ( f=flo; f<=fhi; f+=focus_inc ) {
			fcw = img_focal_wave_intensity(p, wl, f, hires, lores);
			if ( fcw_best < fcw ) {
				fcw_best = fcw;
				focus_step = f;
			}
			if ( verbose )
				cout << f << tab << fcw << endl;
		}
	} else {
		fcw_best = img_focal_wave_intensity(p, wl, f, hires, lores);
		if ( verbose )
			cout << f << tab << fcw_best << endl;

		while ( fabs(focus_inc) > 0.01 ) {
			f += focus_inc;
			fcw = img_focal_wave_intensity(p, wl, f, hires, lores);
			if ( fcw_best < fcw ) {
				fcw_best = fcw;
				focus_step = f;
				if ( verbose )
					cout << f << tab << fcw << endl;
			} else {
				f -= focus_inc;
				focus_inc = -0.9*focus_inc;
			}
		}
	}

	if ( verbose ) {
		cout << "Final focus step size:          " << focus_step << " A" << endl;
		cout << "Focal coherent wave intensity:  " << fcw_best << endl << endl;
	}

	return fcw_best;
}


/**
@brief 	Reconstructs an image from a focal series.
@param 	*p				focal series (images or transforms).
@param 	cp				CTF parameters with average defocus.
@param	&dfocus			vector of focus offsets.
@param	hires			high resolution limit.
@param	flag			0=no correction, 1=phase flip, 2=flatten sphere, 3=both, 4=phase flip after, 8=average.
@return Bimage*			resultant complex image.
	
**/
Bimage*		img_ctf_focal_reconstruct(Bimage* p, CTFparam cp, vector<double>& dfocus, double hires, int flag)
{
	if ( p->compound_type() != TComplex )
		p->fftxy();

	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	
	long			i, j, x, y, z, n;
	double			shi2(1/(hires*hires));
	double			sx, sy, s, s2, a;
	double			def(cp.defocus_average()), pil(M_PI*cp.lambda());
	Vector3<double>	freq_scale(1.0L/p->real_size());
	Vector3<double>	h((p->size() - 1)/2);
	Complex<double>	cv;
	Complex<double>	ctf;
	
	if ( verbose ) {
		cout << "Reconstructing a focal series:" << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Reconstruction flag:            " << flag << endl;
		cp.show();
	}
	
	Bimage*			pr = new Bimage(Float, TComplex, p->size(), p->images());
	pr->fourier_type(Standard);
	pr->sampling(p->sampling(0));
	
	if ( verbose )
		cout << "Image\tDefocus\tDifference" << endl;
	for ( i=j=n=0; n<p->images(); ++n ) {
		for ( z=0; z<p->sizeZ(); ++z ) {
			if ( flag & 1 ) cp.defocus_average(def+dfocus[n+z]);
			if ( verbose )
				cout << n+z << tab << cp.defocus_average() << tab << dfocus[n+z] << endl;
			for ( y=0; y<p->sizeY(); ++y ) {
				sy = y;
				if ( y > h[1] ) sy -= p->sizeY();
				sy *= freq_scale[1];
				for ( x=0; x<p->sizeX(); ++x, ++i, ++j ) {
					sx = x;
					if ( x > h[0] ) sx -= p->sizeX();
					sx *= freq_scale[0];
					s2 = sx*sx + sy*sy;
					if ( s2 <= shi2 ) {
						s = sqrt(s2);
						cv = p->complex(i);
						if ( flag & 5 ) {
							a = atan2(sy,sx);
//							dphi = cp.calculate_aberration_even(s2, a);
//							if ( dphi > 0 ) cv = -cv;
//							ctf = cp.calculate_complex(s2, a);
							ctf = cp.calculate_complex(s, a);
							cv *= ctf.conj();
						}
						if ( flag & 2 )
							cv.shift_phi(-pil*s2*dfocus[n+z]);
						pr->add(j, cv);
					}
				}
			}
		}
	}

	cp.defocus_average(def);
	
	if ( flag & 8 ) {
		pr->slices_to_images();
		pr->average_images();
	}

	return pr;
}

/**
@brief 	Reconstructs an image from a focal series.
@param 	*p				focal series (images or transforms).
@param 	cp				CTF parameters with average defocus.
@param	hires			high resolution limit.
@param	lores			low resolution limit.
@param	flag			0=no correction, 1=phase flip, 2=flatten sphere, 3=both, 4=phase flip after, 8=average.
@return Bimage*			resultant complex image.
	
**/
Bimage*		img_ctf_focal_reconstruct(Bimage* p, CTFparam cp, double hires, double lores, int flag)
{
	if ( p->compound_type() != TComplex )
		p->fftxy();
		
//		p->information();

	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	
	Bimage*			pfs = img_ctf_focal_series(cp, p->images()*p->sizeZ(), p->size(), p->image->sampling(), lores, hires, flag);

	if ( verbose )
		cout << "Reconstructing a focal series" << endl << endl;
	
	p->complex_conjugate_product(pfs);
	
	delete pfs;

//	p->average_images();

	return p->fspace_sum(0);
//	return p->copy();
}

/**
@brief 	Reconstructs the exit wave from a focal series.
@param 	*p				focal series (images or transforms).
@param 	cp				CTF parameters with average defocus.
@param	hires			high resolution limit.
@param	max_iter		maximum number of iterations.
@param	tol				tolerance: minimum change to accept to set a break point.
@param	ctf_flag		CTF form: 0=complex, 1=sine, 2=envelope, 4=combined.
@param	exit_flag		reference images: 0=intensities, 1=sqrt(int), 2=real.
@return Bimage*			resultant complex image.

	The input must be a real space stack of images.
	Reference: Allen et al. Ultramicroscopy 100 (2004) 91–104.
	
**/
Bimage*		img_ctf_focal_exit_wave_reconstruct(Bimage* p, CTFparam cp,
				double hires, long max_iter, double tol, int ctf_flag, int exit_flag)
{
	if ( hires < p->sampling(0)[0] ) hires = p->sampling(0)[0];
	if ( tol < 1e-30 ) tol = 1e-6;
	
	double			ehs(1.0/p->sampling(0)[0]);		// Envelope half maximum
	vector<double>	ce = cp.coherence_envelope(1000, 1.0/p->real_size()[0]);
	for ( long i=0; i<ce.size(); ++i ) if ( ce[i] < 0.5 ) {
		ehs = i*1.0/p->real_size()[0];
		break;
	}
	
	p->slices_to_images();
	
	if ( exit_flag&1 ) p->square_root();

	Bimage*			pfs = img_ctf_focal_series(cp, p->images()*p->sizeZ(),
						p->size(), p->image->sampling(), 0, hires, ctf_flag);

	if ( verbose ) {
		cout << "Reconstructing the exit wave from a focal series:" << setprecision(3) << endl;
		cout << "Focus step size:                " << cp.focus_step() << " A" << endl;
		cout << "Resolution limit:               " << hires << " A" << endl;
		cout << "Tolerance:                      " << setprecision(10) << tol << endl;
		if ( exit_flag&1 ) cout << "Fitting to the square root of the input images" << endl;
		cout << "Coherence envelope half max:    " << setprecision(3) << ehs << " /A" << endl << endl;;
	}

	Bimage*			pc = p->copy();

	fft_plan		planf = pc->fft_setup(FFTW_FORWARD, 0);
	fft_plan		planb = pc->fft_setup(FFTW_BACKWARD, 0);

	long			iter, nn;
	double			R(10), dR(1);
	
	Bimage*			pew = new Bimage(Float, TComplex, p->size(), 1);
	pew->fourier_type(Standard);
	pew->sampling(p->sampling(0));
		
	if ( verbose )
		cout << setprecision(7) << "Iter\tR\tdR" << endl;
	for ( iter=0; iter < max_iter && fabs(dR) > tol; ++iter ) {
		pc->fftp(planf);
		pc->complex_conjugate_product(pfs);
		pew->clear();
		for ( nn=0; nn<pc->images(); ++nn ) pew->add(pc, nn);
		pew->multiply(1.0/nn);
		for ( nn=0; nn<pc->images(); ++nn ) pc->replace(nn, pew);
		pc->complex_product(pfs);
		pc->fftp(planb);
		dR = R;
		R = pc->merge_amplitudes_and_phases(p);
		dR -= R;
		if ( verbose )
			cout << iter+1 << tab << R << tab << dR << endl;
	}
	
	delete pfs;
	delete pc;

	return pew;
}

