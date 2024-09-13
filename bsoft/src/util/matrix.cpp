/**
@file	matrix.cpp
@brief	Matrix manipulation functions
@author Bernard Heymann 
@date	Created: 20000501
@date	Modified: 20240312
**/
 
#include "Matrix.h" 
#include "matrix_linear.h"
#include "math_util.h"
#include "random_numbers.h"
#include "utilities.h" 

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

ostream& operator<<(ostream& output, Matrix& mat) {
	output.setf(ios::scientific, ios::floatfield);
//	output.precision(4);
	long		i, j;
	for ( i=0; i<mat.rows(); i++ ) {
		for ( j=0; j<mat.columns(); j++ )
			output << "\t" << mat[i][j];
		output << endl;
	}
	return output;
}

/**
@brief 	Matrix inversion by LU decomposition.
@return double	 		determinant.

	This inverts matrix A by LU decomposition.
	The matrix A must be square and is converted to and replaced by its inverse.
	Note: The matrix is modified.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double 		Matrix::LU_decomposition()
{
	if ( rows() != columns() ) {
		cerr << "Error: rows:" << rows() << " != columns:" << columns() << endl;
		return 0;
	}
	
	long 			i, j, k, ibig, m(rows()), n(columns());
	double			big, det(1);
	vector<long>	ind(m,0);
	vector<double>	x(m,0);
	Matrix			a(m,m);

	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Matrix::LU_decomposition: rows = " << m << endl;
	
	// Find the largest element in each row for scaling and pivoting
	// Additionally test for singularity
	for ( i=0; i<n; i++ ) {
		ind[i] = i; 		// Initialize the indices for tracking pivoting
		big = 0;
		for ( j=0; j<n; j++ )
			if ( big < fabs(d[i][j]) ) big = fabs(d[i][j]);
		if ( big < 1e-37 ) {
			cerr << "Error: Singular matrix in Matrix::LU_decomposition" << endl;
//			if ( verbose & VERB_DEBUG ) {
				cerr << *this << endl;
//				bexit(-1);
//			}
			return -1e30;
		} else
			x[i] = 1.0/big;
	}
	
	// The LU decomposition
	for ( j=0; j<n; j++ ) {
		for ( i=0; i<j; i++ ) {			// Calculating the upper triangle
			for ( k=0; k<i; k++ )		//  not the diagonal
				d[i][j] -= d[i][k]*d[k][j];
		}
		big = 0;			// Largest element for pivoting
		ibig = -1;
		for ( i=j; i<n; i++ ) { 		// Calculating the lower triangle
			for ( k=0; k<j; k++ )		//  and diagonal (i==j)
				d[i][j] -= d[i][k]*d[k][j];
			if ( x[i]*fabs(d[i][j]) >= big ) { 	// Get the biggest element
				big = x[i]*fabs(d[i][j]);		//   for pivoting
				ibig = i;
			}
		}
		if ( j != ibig ) {						// Interchange rows if necessary
			for ( k=0; k<n; k++ ) swap(d[ibig][k], d[j][k]);
			swap(x[ibig], x[j]); 		// Switch scales
			swap(ind[ibig], ind[j]);	// Switch indices
			det = -det; 				// Change the sign of the determinant
		}
		if ( d[j][j] == 0 ) d[j][j] = 1e-37;
		for ( i=j+1; i<n; i++ )	d[i][j] /= d[j][j];	// Divide by pivot element
	}
	
	// Calculate the determinant
	for ( i=0; i<n; i++ ) det *= d[i][i];
	
	// Invert matrix
	for ( k=0; k<n; k++ ) {
		for ( i=0; i<n; i++ ) x[i] = 0;
		x[k] = 1;
		for ( i=0; i<n; i++ )			// Forward substitution
			for ( j=0; j<i; j++ ) x[i] -= d[i][j]*x[j];
		for ( i=m-1; i>=0; i-- ) {		// Backward substitution
			for ( j=i+1; j<n; j++ ) x[i] -= d[i][j]*x[j];
			x[i] /= d[i][i];
			a[i][k] = x[i];
		}
	}
	
	// Reorder the matrix back to the original row order
	for ( i=0; i<n; i++ ) {
		j = ind[i]; 					// Find the old row at this position
		for ( k=0; k<n; k++ ) 			// Pack it back in order into matrix A
			d[k][j] = a[k][i];
	}
	
	if ( verbose & VERB_DEBUG )
		cout << "DEBUG Matrix::LU_decomposition: done" << endl << *this << endl;
	
	return det;
}

/**
@brief 	Singular value decomposition of a matrix .
@return double*			0.

The matrix A is replaced by the matrix U.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
double		Matrix::singular_value_decomposition()
{
	long			flag, i, its, j, jj, k, l, nm, m(rows()), n(columns());
	double			anorm, c, f, g, h, s, scale, x, y, z, t;

	vector<double>	w(m,0);
	vector<double>	tmp(m,0);
	Matrix			v(m,m);

	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+1;
		tmp[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<m;k++) scale += fabs(d[k][i]);
			if (scale) {
				for (k=i;k<m;k++) {
					d[k][i] /= scale;
					s += d[k][i]*d[k][i];
				}
				f=d[i][i];
				g = -copysign(sqrt(s),f);
				h=f*g-s;
				d[i][i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += d[k][i]*d[k][j];
					f=s/h;
					for (k=i;k<m;k++) d[k][j] += f*d[k][i];
				}
				for (k=i;k<m;k++) d[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if ( i < m && i != n ) {
			for (k=l;k<n;k++) scale += fabs(d[i][k]);
			if (scale) {
				for (k=l;k<n;k++) {
					d[i][k] /= scale;
					s += d[i][k]*d[i][k];
				}
				f=d[i][l];
				g = -copysign(sqrt(s),f);
				h=f*g-s;
				d[i][l]=f-g;
				for (k=l;k<n;k++) tmp[k]=d[i][k]/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += d[j][k]*d[i][k];
					for (k=l;k<n;k++) d[j][k] += s*tmp[k];
				}
				for (k=l;k<n;k++) d[i][k] *= scale;
			}
		}
		t = fabs(w[i])+fabs(tmp[i]);
		if ( anorm < t ) anorm = t;
	}
	
	for (i=m-1;i>=0;i--) {
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++)
					v[j][i]=(d[i][j]/d[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += d[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=tmp[i];
		l=i;
	}
	
	for (i=( m < n )? m-1: n-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) d[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += d[k][i]*d[k][j];
				f=(s/d[i][i])*g;
				for (k=i;k<m;k++) d[k][j] += f*d[k][i];
			}
			for (j=i;j<m;j++) d[j][i] *= g;
		} else for (j=i;j<m;j++) d[j][i]=0.0;
		++d[i][i];
	}

	for (k=m-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(fabs(tmp[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ( l==0 ) break;	// JBH
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
//			printf("k=%d  its=%d  nm=%d  flag=%d\n", k, its, nm, flag);
			if (flag) {	// Examples never get into this block
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*tmp[i];
					tmp[i]=c*tmp[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=sqrt(f*f + g*g);	// JBH
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=d[j][nm];
						z=d[j][i];
						d[j][nm]=y*c+z*s;
						d[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) cerr << "Error: no convergence in 30 svdcmp iterations" << endl;
			x=w[l];
			nm=k-1;
			y = g = 0;			// JBH
			if ( nm >= 0 ) {	// JBH
				y=w[nm];		// JBH
				g=tmp[nm];		// JBH
			}
			h=tmp[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=sqrt(f*f + 1.0);	// JBH
			f=((x-z)*(x+z)+h*((y/(f+copysign(g,f)))-h))/x;
			c=s=1.0;
//			cout << "k=%d  its=%d  l=%d  nm=%d\n", k, its, l, nm);
			for (j=l;j<=nm;j++) {	// if nm<0 => no loop
				i=j+1;
				g=tmp[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=sqrt(f*f + h*h);	// JBH
				tmp[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=sqrt(f*f + h*h);	// JBH
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=d[jj][j];
					z=d[jj][i];
					d[jj][j]=y*c+z*s;
					d[jj][i]=z*c-y*s;
				}
			}
			tmp[l]=0.0;
			tmp[k]=f;
			w[k]=x;
		}
	}
	
	// Setting singular values to zero
	double		wmax, wmin;
	for ( wmax=w[0], l=0; l<n; l++ )
		if ( wmax < w[l] ) wmax = w[l];
	wmin = 1e-6*wmax;
	for ( l=0; l<n; l++ )
		if ( w[l] < wmin ) w[l] = 0;
		else w[l] = 1/w[l];

	// Invert
	Matrix		u(*this);
	
	for ( j=0; j<n; j++ )
		for ( i=0; i<n; i++)
			v[j][i] *= w[i];
		
	for ( j=0; j<n; j++ )
		for ( i=0; i<m; i++ ) {
			d[j][i] = 0;
			for ( k=0; k<n; k++ )
				d[j][i] += v[j][k]*u[i][k];
		}

	return 0;
}

int			Matrix::jrotate(double s, double tau, long i, long j, long k, long l)
{
	double			g, h;
	
	g = d[i][j];
	h = d[k][l];
	d[i][j] = g - s*(h + g*tau);
	d[k][l] = h + s*(g - h*tau);
	
	return 0;
}

/**
@brief 	Computes all eigenvalues and eigenvectors of a real symmetric matrix.
@return	double* val			eigenvalues.

The eigenvectors are returned in the columns of the input matrix.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
vector<double>	Matrix::jacobi_rotation()
{
	long			j, iq, ip, i, m(rows()), n(columns());
	double			tresh, theta, tau, t, sm, s, h, g, c;
	vector<double>	val(n);
	Matrix			vec(m,n);
	
	for ( ip=i=j=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++, i++ ) {
			j += ( fabs(d[ip][iq]) < 1e-12 );	// Check if small to prevent too many iterations bailout
		}
		vec[ip][ip] = 1.0;
		val[ip] = d[ip][ip];
	}

	if ( j == n*n ) return val;			// The matrix is zero

	vector<double>	b(n);
	vector<double>	z(n);

	for ( ip=0; ip<n; ip++ ) {
		b[ip] = val[ip];
		z[ip] = 0;
	}

	for ( i=1; i<=50; i++ ) {
		sm=0.0;
		for ( ip=0; ip<m-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ )
				sm += fabs(d[ip][iq]);
		}
		if ( sm == 0.0 ) {
			for ( ip=0; ip<m; ++ip )
				for ( iq=0; iq<n; ++iq )
					d[ip][iq] = vec.d[ip][iq];
			return val;						// Function exit
		}
		if ( i < 4 )
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for ( ip=0; ip<m-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ ) {
				g = 100.0*fabs(d[ip][iq]);
				if ( i > 4 && (double)(fabs(val[ip])+g) == (double)fabs(val[ip])
							&& (double)(fabs(val[iq])+g) == (double)fabs(val[iq]))
					d[ip][iq] = 0.0;
				else if ( fabs(d[ip][iq]) > tresh ) {
					h = val[iq] - val[ip];
					if ( (double)(fabs(h)+g) == (double)fabs(h))
						t = d[ip][iq]/h;
					else {
						theta = 0.5*h/d[ip][iq];
						t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*d[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					val[ip] -= h;
					val[iq] += h;
					d[ip][iq] = 0.0;
					for ( j=0; j<ip; j++ ) {
						jrotate(s, tau, j, ip, j, iq);
					}
					for ( j=ip+1; j<iq; j++ ) {
						jrotate(s, tau, ip, j, j, iq);
					}
					for ( j=iq+1; j<n; j++ ) {
						jrotate(s, tau, ip, j, iq, j);
					}
					for ( j=0; j<n; j++ ) {
						vec.jrotate(s, tau, j, ip, j, iq);
					}
//					nrot++;
				}
			}
		}
		for ( ip=0; ip<n; ip++ ) {
			b[ip] += z[ip];
			val[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	
	error_show("Error in Matrix::jacobi_rotation", __FILE__, __LINE__);
	cerr << "Too many iterations (" << i << ")" << endl;
	
	for ( ip=0; ip<m; ip++ ) {
		for ( iq=0; iq<n; iq++ ) d[ip][iq] = 0.0;
		d[ip][ip] = 1.0;
		val[ip] = 0;
	}

	return val;
}
/*
long	Matrix::jacobi_rotation_row(long ip, long n, long i, double* val, double* z, Matrix* vec)
{
	long		j, iq, nrot(0);
	double		tresh, theta, tau, t, sm, s, h, g, c;

	for ( iq=ip+1; iq<n; iq++ ) {
		g = 100.0*fabs(d[ip][iq]);
		if ( i > 4 && (double)(fabs(val[ip])+g) == (double)fabs(val[ip])
					&& (double)(fabs(val[iq])+g) == (double)fabs(val[iq]))
			d[ip][iq] = 0.0;
		else if ( fabs(d[ip][iq]) > tresh ) {
			h = val[iq] - val[ip];
			if ( (double)(fabs(h)+g) == (double)fabs(h))
				t = d[ip][iq]/h;
			else {
				theta = 0.5*h/d[ip][iq];
				t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
				if (theta < 0.0) t = -t;
			}
			c = 1.0/sqrt(1+t*t);
			s = t*c;
			tau = s/(1.0+c);
			h = t*d[ip][iq];
			z[ip] -= h;
			z[iq] += h;
			val[ip] -= h;
			val[iq] += h;
			d[ip][iq] = 0.0;
			for ( j=0; j<ip; j++ ) {
				jrotate(s, tau, j, ip, j, iq);
			}
			for ( j=ip+1; j<iq; j++ ) {
				jrotate(s, tau, ip, j, j, iq);
			}
			for ( j=iq+1; j<n; j++ ) {
				jrotate(s, tau, ip, j, iq, j);
			}
			for ( j=0; j<n; j++ ) {
				vec->jrotate(s, tau, j, ip, j, iq);
			}
			nrot++;
		}
	}
	
	return nrot;
}
double*		Matrix::jacobi_rotation_parallel()
{
	long		j, iq, ip, i, nrot(0);
	double		tresh, theta, tau, t, sm, s, h, g, c;
	double*		val = new double[n];
	Matrix*		vec = new Matrix(m,n);
	
	for ( ip=i=j=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++, i++ ) {
			j += ( fabs(d[ip][iq]) < 1e-12 );	// Check if small to prevent too many iterations bailout
		}
		(*vec)[ip][ip] = 1.0;
		val[ip] = d[ip][ip];
	}

	if ( j == n*n ) return val;			// The matrix is zero

	double*		b = new double[n];
	double*		z = new double[n];

	for ( ip=0; ip<n; ip++ ) {
		b[ip] = val[ip];
		z[ip] = 0;
	}

	for ( i=1; i<=50; i++ ) {
		sm=0.0;
		for ( ip=0; ip<n-1; ip++ ) {
			for ( iq=ip+1; iq<n; iq++ )
				sm += fabs(d[ip][iq]);
		}
		if ( sm == 0.0 ) {
			delete[] b;
			delete[] z;
			for ( ip=0; ip<len; ip++ ) d[ip] = vec->d[ip];
			return val;						// Function exit
		}
		if ( i < 4 )
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		
		for ( ip=0; ip<n-1; ip++ ) {
#ifdef HAVE_GCD
			dispatch_apply(n - 1, dispatch_get_global_queue(0, 0), ^(size_t ip){
				jacobi_rotation_row(ip, n, i, val, z, vec);
			});
#else
#pragma omp parallel for
			for ( ip=0; ip<n-1; ip++ ) {
				jacobi_rotation_row(ip, n, i, val, z, vec);
			}
#endif
		}

		for ( ip=0; ip<n; ip++ ) {
			b[ip] += z[ip];
			val[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	
	error_show("Error in Matrix::jacobi_rotation", __FILE__, __LINE__);
	cerr << "Too many iterations (" << i << ")" << endl;
	
	for ( ip=0; ip<n; ip++ ) {
		for ( iq=0; iq<n; iq++ ) d[ip][iq] = 0.0;
		d[ip][ip] = 1.0;
		val[ip] = 0;
	}

	return val;
}
*/

/**
@brief 	Sorts eigenvalues into descending order and rearranges matrix columns accordingly.
@param 	val					eigenvalues.

	The eigenvectors are in the columns.
	This method uses straight insertion.
Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
void 	Matrix::eigen_sort(vector<double>& val)
{
	long			k,j,i;
	double			t;

	for ( i=0; i<rows()-1; i++ ) {
		t = val[k=i];
		for ( j=i+1; j<columns(); j++ )
			if ( val[j] >= t ) t = val[k=j];
		if ( k != i ) {
			val[k] = val[i];
			val[i] = t;
			for ( j=0; j<columns(); j++ )
				swap(d[j][i], d[j][k]);
		}
	}
}

/**
@brief 	Calculates the principal axes of 3D coordinates.
@param 	avg				average of vectors.
@param 	avg2			average of squared vectors.
@param 	avgx			average of cross products {xy, xz, yz}.
@param 	*eigenvec		3 return eigen vectors (can be NULL).
@return Vector3<double>	principal axes.

Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
Vector3<double> 	principal_axes(Vector3<double> avg, Vector3<double> avg2, Vector3<double> avgx, Vector3<double>* eigenvec)
{
	Matrix			a(3,3);
	
	a[0][0] = avg2[0] - avg[0]*avg[0];
	a[1][1] = avg2[1] - avg[1]*avg[1];
	a[2][2] = avg2[2] - avg[2]*avg[2];
	a[0][1] = a[1][0] = avgx[0] - avg[0]*avg[1];
	a[0][2] = a[2][0] = avgx[1] - avg[0]*avg[2];
	a[1][2] = a[2][1] = avgx[2] - avg[1]*avg[2];
	
	vector<double>	d = a.jacobi_rotation();
	a.eigen_sort(d);
	
	Vector3<double>	eigenval(d[0], d[1], d[2]);
	Vector3<double>	pax = (eigenval * 3.0).square_root();
	
	long			i, j;
	if ( eigenvec ) {
		for ( i=0; i<3; i++ )
			for ( j=0; j<3; j++ )
				eigenvec[i][j] = a[j][i];
	}
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Center of mass:                 " << avg << " A" << endl;
		cout << "Principal axis lengths:         " << pax << " A" << endl;
		cout << "Volume:                         " << (4.0*M_PI/3.0)*pax.volume() << " A3" << endl;
		cout << "Eigenvalues:                    " << eigenval << endl; 
		cout << "Eigenvectors:" << endl;
		for ( i=0; i<3; i++ )
			cout << tab << a[0][i] << tab << a[1][i] << tab << a[2][i] << endl;
		cout << endl;
	}
	
	return pax;
}

/**
@brief 	Calculates the principal axes of 3D coordinates.
@param 	coor			vector of coordinates.
@param 	&a				matrix with eigen vectors.
@return Vector3<double>	principal axes.

Reference: 	Press W.H. et al (1992) Numerical Recipes in C.

**/
Vector3<double> principal_axes(vector<Vector3<double>>& coor, Matrix& a)
{
	long			n(0);
	Vector3<double>	s, s2, sx, pax;
	
	for ( auto loc: coor ) {
		s += loc;					// Sums
		s2 += loc*loc;				// Square sums
		sx[0] += loc[0]*loc[1];		// Cross-term sums
		sx[1] += loc[0]*loc[2];
		sx[2] += loc[1]*loc[2];
		n++;
	}
	
	if ( n < 1 ) return pax;
	
	s /= n;
	s2 /= n;
	sx /= n;

	a[0][0] = s2[0] - s[0]*s[0];
	a[1][1] = s2[1] - s[1]*s[1];
	a[2][2] = s2[2] - s[2]*s[2];
	a[0][1] = a[1][0] = sx[0] - s[0]*s[1];
	a[0][2] = a[2][0] = sx[1] - s[0]*s[2];
	a[1][2] = a[2][1] = sx[2] - s[1]*s[2];
	
	vector<double>	d = a.jacobi_rotation();
	a.eigen_sort(d);
	
	a = a.transpose();
	
	Vector3<double>	eigenval(d[0], d[1], d[2]);
	
	pax = Vector3<double>(sqrt(3.0*d[0]), sqrt(3.0*d[1]), sqrt(3.0*d[2]));
	
	double		sph(0), pln(0);
	sph = pax[2]/pax[0];
	pln = (1 - sph)*((pax[1] - pax[2])/pax[0]);
	
	if ( verbose & VERB_PROCESS ) {
		cout << "Center of mass:                 " << s << " A" << endl;
		cout << "Principal axis lengths:         " << pax << " A" << endl;
		cout << "Sphericity:                     " << sph << endl;
		cout << "Planarity:                      " << pln << endl;
		cout << "Volume:                         " << (4.0*M_PI/3.0)*pax.volume() << " A3" << endl;
		cout << "Eigenvalues:                    " << eigenval << endl;
		cout << "Eigenvectors:" << endl;
		cout << a << endl;
		cout << endl;
	}
	
	return pax;
}

