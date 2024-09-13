/**
@file	Matrix.h 
@brief	Generalized matrix class
@author Bernard Heymann
@date	Created: 20000501
@date	Modified: 20240313
**/

#include "string_util.h"
#include "random_numbers.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen

#ifndef _Matrix_
#define _Matrix_

class Matrix
{
private:
	vector<vector<double>>	d;
	void		resize(long nr, long nc, double fill) {
		d.resize(nr);
		for ( auto& r: d ) r.resize(nc, fill);
	}
	bool		check_same_size(Matrix& mat) {
		if ( rows() == mat.rows() && columns() == mat.columns() )
			return 1;
		cerr << "Matrices not the same size!" << endl;
		return 0;
	}
public:
	Matrix() { }
	Matrix(Matrix& mat) {
		d.resize(mat.rows());
		for ( long r=0; r<d.size(); ++r )
			d[r] = mat[r];
	}
	Matrix(long rows, long cols) {
		resize(rows, cols, 0);
	}
	Matrix(Bstring& filename) {
		long			m(0), n(0);
		ifstream		fmat(filename.c_str());
		if ( fmat.fail() ) return;
		string			s;
		long		 	i(0), j;
		vector<string>	token;
		while ( getline(fmat, s) && !fmat.eof() && s.find("matrix_") == string::npos ) ;
		if ( verbose & VERB_DEBUG )
			cout << s << endl;
		while ( getline(fmat, s) ) {
			token = split(s);
			if ( token.size() ) {
				if ( token[0] == "_rows" ) {
					m = to_integer(token[1]);
				} else if ( token[0] == "_columns" ) {
					n = to_integer(token[1]);
					if ( m && n ) { resize(m, n, 0); }
//					cout << "Columns: " << n << endl;
				} else if ( i < m ) {
					if ( verbose & VERB_DEBUG )
						cout << "Row " << i+1 << ": " << token.size() << endl;
					if ( token.size() >= n ) {
						for ( j=0; j<n; ++j ) d[i][j] = to_real(token[j]);
						i++;
					}
				}
			}
		}
		fmat.close();
		if ( m < 1 || n < 1 ) {
			cerr << "Error: The matrix size is zero!" << endl;
			bexit(-1);
		}
		if ( verbose & VERB_DEBUG )
			cout << filename << " read" << endl;
	}
	
	void	write(Bstring& filename) {
		ofstream        fmat(filename.c_str());
		if ( fmat.fail() ) return;
		long	 		j;
		fmat << "# Matrix written by Bsoft\n\ndata_\n\nmatrix_" << endl;
		fmat << "_rows" << tab << d.size() << endl;
		fmat << "_columns" << tab << d[0].size() << endl;
		for ( auto& r: d ) {
			fmat << r[0];
			for ( j=1; j<columns(); ++j ) fmat << " " << r[j];
			fmat << endl;
		}
		fmat << endl;
		fmat.close();
	}
	
	long	rows() { return d.size(); }
	long	columns() { return d[0].size(); }
	
	vector<double>&	operator[](long i) { return d[i]; }

/*	Matrix	operator=(Matrix& mat) {
		return Matrix(mat);
	}
	
	Matrix	operator=(Matrix& mat) {
		d.resize(mat.rows());
		for ( long r=0; r<d.size(); ++r )
			d[r] = mat[r];
		return *this;
	}
	Matrix	operator=(const Matrix& mat) {
		d.resize(mat.rows());
		for ( long r=0; r<d.size(); ++r )
			d[r] = mat[r];
		return *this;
	}*/
	Matrix	operator+=(Matrix& mat) {
		if ( !check_same_size(mat) ) return *this;
		long		i, j;
		for ( i=0; i<rows(); i++ )
			for ( j=0; j<columns(); j++ )
				d[i][j] += mat[i][j];
		return *this;
	}
	Matrix	operator+(Matrix& mat) {
		Matrix		numat(*this);
		numat += mat;
		return numat;
	}
	Matrix	operator-=(Matrix& mat) {
		if ( !check_same_size(mat) ) return *this;
		long		i, j;
		for ( i=0; i<rows(); i++ )
			for ( j=0; j<columns(); j++ )
				d[i][j] -= mat[i][j];
		return *this;
	}
	Matrix	operator-(Matrix& mat) {
		Matrix		numat(*this);
		numat -= mat;
		return numat;
	}
	Matrix	operator-() {
		Matrix		mat(*this);
		for ( auto& r: mat.d )
			for ( auto& v: r )
				v = -v;
		return mat;
	}
/*	Matrix	operator*=(Matrix mat) {
		return *this * mat;
	}*/
	Matrix	operator*=(double d) {
		Matrix		mat(*this);
		for ( auto& r: mat.d )
			for ( auto& v: r )
				v *= d;
		return mat;
	}
	Matrix	operator*(Matrix& mat) {
		long		i, j, k;
		Matrix		numat(mat.rows(), columns());
		if ( rows() != mat.columns() ) {
			cerr << "Matrix rows not equal to second vector columns!" << endl;
			return numat;
		}
		for ( i=0; i<columns(); ++i )
			for ( j=0; j<mat.rows(); ++j )
				for ( k=0; k<rows(); ++k ) numat[i][j] += (*this)[i][k]*mat[k][j];
		return numat;
	}
	vector<double>	operator*(vector<double>& vec) {
		long			i, j;
		vector<double>	nuvec(columns(),0);
		if ( columns() != vec.size() ) {
			cerr << "Matrix columns not equal to vector size!" << endl;
			return nuvec;
		}
		for ( i=0; i<rows(); i++ )
			for ( j=0; j<columns(); j++ )
				nuvec[i] += (*this)[i][j]*vec[j];
		return nuvec;
	}

	void	show_below_cutoff(double d) {
		for ( long i=0; i<rows(); ++i ) {
			cout << i;
			for ( long j=0; j<columns(); ++j )
				if ( (*this)[i][j] <= d )
					cout << tab << j << "," << (*this)[i][j];
			cout << endl;
		}
	}

	void	swap_rows_columns(long rc1, long rc2) {
		long			i;
		for ( i=0; i<rows(); ++i ) swap((*this)[i][rc1], (*this)[i][rc2]);
		for ( i=0; i<columns(); ++i ) swap((*this)[rc1][i], (*this)[rc2][i]);
	}
	
	Matrix	delete_row_column(long rc) {
		long			i, j, k, l;
		Matrix			mat(rows()-1,columns()-1);
		for ( i=k=0; i<rows(); ++i ) {
			for ( j=l=0; j<columns(); ++j ) {
				if ( i != rc && j != rc ) mat[k][l] = d[i][j];
				if ( j != rc ) l++;
			}
			if ( i != rc ) k++;
		}
		return mat;
	}

	Matrix	transpose() {
		Matrix			mat(columns(),rows());
		for ( long i=0; i<columns(); i++ )
			for ( long j=0; j<rows(); j++ )
				mat[i][j] = d[j][i];
		return mat;
	}

	void	fill(double f) {
		for ( auto& r: d )
			for ( auto& v: r )
				v = f;
	}
	int		check_for_singularity() {
		int				c(0);
		long			i, j;
		double			max;
		for ( i=0; i<rows(); ++i ) {
			for ( j=0, max=0; j<columns(); ++j ) if ( max < fabs(d[i][j]) ) max = fabs(d[i][j]);
			if ( max < 1e-37 ) c |= 1 << i;
		}
		return c;
	}
	/**
		The rows and columns are alternatively iteratively normalized until
		the error is small enough.
	**/
	void			normalize()
	{
		long			i, r, c;
		double			err(1);
		vector<double>	rw(rows());
		vector<double>	cw(columns());
	
		if ( verbose & VERB_FULL )
			cout << "Cycle\tError" << endl;
		for ( i=0; i<100 && err > 1e-20; i++ ) {
			for ( r=0; r<rows(); r++ ) rw[r] = 0;	// Row scaling
			for ( r=0; r<rows(); r++ )
				for ( c=0; c<columns(); c++ )
					rw[r] += d[r][c];
			for ( r=0; r<rows(); r++ ) rw[r] = 1/rw[r];
			for ( r=0; r<rows(); r++ )
				for ( c=0; c<columns(); c++ )
					d[r][c] *= rw[r];

			for ( c=0; c<columns(); c++ ) cw[c] = 0;	// Column scaling
			for ( r=0; r<rows(); r++ )
				for ( c=0; c<columns(); c++ )
					cw[c] += d[r][c];
			for ( c=0; c<columns(); c++ ) cw[c] = 1/cw[c];
			for ( r=0; r<rows(); r++ )
				for ( c=0; c<columns(); c++ )
					d[r][c] *= cw[c];

			for ( err=0, r=0; r<rows(); r++ ) err += rw[r];
			for ( c=0; c<columns(); c++ ) err += cw[c];
			err = fabs(err/(rows() + columns()) - 1);
			if ( verbose & VERB_FULL )
				cout << i+1 << "\t" << err << endl;
		}
	}
	void		randomize() {
		random_seed();
		double		rm = INT_MAX/4, irm = 10.0L/INT_MAX;
		for ( auto& r: d )
			for ( auto& v: r )
				v = irm*(random() - rm);
	}
	int			multiply_in_place(vector<double>& vec) {
		if ( vec.size() != columns() ) return -1;
		long			i, j;
		vector<double>	t(vec.size(),0);
		for ( i=0; i<vec.size(); i++ )
			for ( j=0; j<vec.size(); j++ )
				t[i] += d[i][j]*vec[j];
		vec = t;
		return 0;
	}
	double		determinant() { return LU_decomposition(); }
	double 		LU_decomposition();
	double 		LU_decomposition(vector<double>& b) {
		double	det = LU_decomposition();
		multiply_in_place(b);
		return det;
	}
	double 		singular_value_decomposition();
	double 		singular_value_decomposition(vector<double>& b) {
		singular_value_decomposition();
		multiply_in_place(b);
		return 0;
	}
	int			jrotate(double s, double tau, long i, long j, long k, long l);
	vector<double>	jacobi_rotation();
	void 		eigen_sort(vector<double>& val);
};

ostream& operator<<(ostream& output, Matrix& mat);

#endif

// Function prototypes 
Vector3<double> principal_axes(Vector3<double> avg, Vector3<double> avg2, Vector3<double> avgx, Vector3<double>* eigenvec);
Vector3<double> principal_axes(vector<Vector3<double>>& coor, Matrix& eigenvec);
