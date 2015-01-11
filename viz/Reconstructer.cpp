//
// Original authors: Andrew Dowsey <andrew.dowsey@manchester.ac.uk>
//
// Copyright 2014 CADET Bioinformatics Laboratory
//                University of Manchester, UK
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "Reconstructer.hpp"

#include <omp.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <math.h>

#include <Eigen/SparseCore>


#define DEBUG 0


namespace bspline {

// quartic b-spline basis function
inline double
qbs0(const double& b4, const double& b3, const double& b2, const double& b)
{
	return  1.0/24.0*b4;
}

inline double
qbs1(const double& b4, const double& b3, const double& b2, const double& b)
{
	return -4.0/24.0*b4 +  4.0/24.0*b3 + 6.0/24.0*b2 +  4.0/24.0*b +  1.0/24.0;
}

inline double
qbs2(const double& b4, const double& b3, const double& b2, const double& b)
{
	return  6.0/24.0*b4 - 12.0/24.0*b3 - 6.0/24.0*b2 + 12.0/24.0*b + 11.0/24.0;
}

inline double
qbs3(const double& b4, const double& b3, const double& b2, const double& b)
{
	return -4.0/24.0*b4 + 12.0/24.0*b3 - 6.0/24.0*b2 - 12.0/24.0*b + 11.0/24.0;
}

inline double
qbs4(const double& b4, const double& b3, const double& b2, const double& b)
{
	return  1.0/24.0*b4 -  4.0/24.0*b3 + 6.0/24.0*b2 -  4.0/24.0*b +  1.0/24.0;
}

// i-spline basis function (culminative cubic b-spline basis function) [we can do faster]
double
iqbs(const double& x)
{
	if (x <= 0.0) return 0.0;
	if (x >= 4.0) return 1.0;

	int ix = (int) x;
	double b = x - ix;
	double b2 = b*b;
	double b3 = b*b2;
	double b4 = b*b3;

	double v =        qbs0(b4,b3,b2,b);
	if (ix >= 1) v += qbs1(b4,b3,b2,b);
	if (ix >= 2) v += qbs2(b4,b3,b2,b);
	if (ix >= 3) v += qbs3(b4,b3,b2,b);
	//if (ix >= 4) v += qbs3(b4,b3,b2,b); // will always be 0
	return v;
}

int
factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// uniform knot-insertion filter for cubic b-spline
double cbsf[5] = {0.0625, 0.25, 0.375, 0.25, 0.0625};

}


bool 
Coef::
operator()(Coef& a, Coef& b)
{
	if (a.ly < b.ly)
	{
		return true;
	}
	else if (a.ly == b.ly)
	{
		if (a.lx < b.lx)
		{
			return true;
		}
		else if (a.lx == b.lx)
		{
			if (a.y < b.y)
			{
				return true;
			}
			else if (a.y == b.y)
			{
				if (a.x < b.x)
				{
					return true;
				}
			}
		}
	}
	return false;
}


Reconstructer::
Reconstructer(const string& out_dir, int width, int height) :
	out_path(out_dir),
	img(height, width),
	stream_index(0),
	read_time(0.0),
	viz_time(0.0)
{
	if (filesystem::exists(out_path))
	{
		for (filesystem::directory_iterator end_dir_it, it(out_path); it!=end_dir_it; ++it)
		{
			filesystem::remove_all(it->path());
		}
	}
	else
	{
		filesystem::create_directories(out_path);
	}

	//setNbThreads(1);
	//cout << "Using " << nbThreads() << " threads" << endl;
}

void
Reconstructer::
next_stream(double _mz_min, double _mz_max, double _rt_min, double _rt_max, double _max_counts)
{
	mz_min = _mz_min;
	mz_max = _mz_max;
	rt_min = _rt_min;
	rt_max = _rt_max;

	//un-normalise to display resolution
	double x_res = (mz_max - mz_min) / img.cols() * 60.0/1.0033548378;
	double y_res = (rt_max - rt_min) /img.rows();
	max_counts = _max_counts / (x_res * y_res); 

	ostringstream oss; oss << stream_index;
	stream_path = out_path / oss.str();
	filesystem::create_directory(stream_path);

	filesystem::path file_path = stream_path / "stream.csv";	
	stream_ofs.open(file_path.string().c_str());
    stream_ofs << "0,0.0,0.0" << endl;

	img.setZero();
	stream_index++;
	chunk_index = 0;
	chunk_count = 0;

	read_start = omp_get_wtime();

	nx = 0; ny = 0; nxy = 0;
}

void
Reconstructer::
next_chunk(const vector<Coef>& cs)
{
	// time chunk
	chunk_count += cs.size();
	read_time += omp_get_wtime() - read_start;
	stream_ofs << chunk_count << "," << read_time << ",";

	// write chunk to csv file
	ostringstream txtfile_oss; txtfile_oss << setfill('0') << setw(8) << chunk_index << ".txt";
	filesystem::path txtfile_path = stream_path / txtfile_oss.str();
	ofstream txtfile_ofs(txtfile_path.string().c_str());
	for (size_t i = 0; i < cs.size(); i++)
	{
		txtfile_ofs << cs[i].v << ":[" << cs[i].lx << "," << cs[i].ly << "]:[" << cs[i].x << "," << cs[i].y << "]" << endl;
	}

	// time display
	viz_start = omp_get_wtime();

	// CONSTRUCT VISUALISATION //////////////////////////////////////////////////////////

	// determine lowest and highest lx and ly
	int lx0 = numeric_limits<int>::max(); 
	int ly0 = numeric_limits<int>::max(); 
	int lx1 = numeric_limits<int>::min(); 
	int ly1 = numeric_limits<int>::min(); 
	for (size_t i = 0; i < cs.size(); i++)
	{
		lx0 = lx0 < cs[i].lx ? lx0 : cs[i].lx; 
		ly0 = ly0 < cs[i].ly ? ly0 : cs[i].ly; 
		lx1 = lx1 > cs[i].lx ? lx1 : cs[i].lx; 
		ly1 = ly1 > cs[i].ly ? ly1 : cs[i].ly; 
	}

	// what is the finest lx and ly appropriate for our output resolution?
	int cutoff = 1; // 2 is also ok
	int mzr = (int) (-log((mz_max - mz_min) / img.cols() * 60.0/1.0033548378) / log(2.0)) - cutoff;
	int rtr = (int) (-log((rt_max - rt_min) / img.rows()) / log(2.0)) - cutoff;

	lx1 = lx1 < mzr ? lx1 : mzr;
	ly1 = ly1 < rtr ? ly1 : rtr;
	int nlx = lx1 - lx0 + 1;
	int nly = ly1 - ly0 + 1;

	// for basis functions smaller than our cutoff, add them to individually
	int processed = 0;
	for (int i = 0; i < cs.size(); i++)
	if (cs[i].lx > lx1 || cs[i].ly > ly1)
	{
		// bounding box for this basis function in mz/rt space
		double mz = pow(2.0, -cs[i].lx) * 1.0033548378/60.0;
		double mz0 = (cs[i].x - 3) * mz;
		double mz1 = (cs[i].x + 1) * mz;
		double rt = pow(2.0, -cs[i].ly);
		double rt0 = (cs[i].y - 3) * rt;
		double rt1 = (cs[i].y + 1) * rt;

		// project into viewport x/y space
		double x0f = img.cols() * (mz0 - mz_min) / (mz_max - mz_min);
		double x1f = img.cols() * (mz1 - mz_min) / (mz_max - mz_min);
		double y0f = img.rows() * (rt0 - rt_min) / (rt_max - rt_min);
		double y1f = img.rows() * (rt1 - rt_min) / (rt_max - rt_min);

		// bounding box intersection
		int x0 = (int) floor(x0f); x0 = x0 > 0 ? x0 : 0;
		int x1 = (int) ceil(x1f); x1 = x1 < img.cols() ? x1 : img.cols();
		int y0 = (int) floor(y0f); y0 = y0 > 0 ? y0 : 0;
		int y1 = (int) ceil(y1f); y1 = y1 < img.rows() ? y1 : img.rows();

		vector<double> bx(x1-x0+1);
		for (size_t x = 0; x < bx.size(); ++x)
		{
			bx[x] = bspline::iqbs(4.0 * (x + x0 - x0f) / (x1f - x0f));
		}

		vector<double> by(y1-y0+1);
		for (int y = 0; y < by.size(); ++y)
		{
			by[y] = bspline::iqbs(4.0 * (y + y0 - y0f) / (y1f - y0f));
		}

		// add separable b-spline integral to image
		for (int y = y0; y < y1; ++y)
		for (int x = x0; x < x1; ++x)
		{
			img(y, x) += cs[i].v * (bx[x-x0+1] - bx[x-x0]) * (by[y-y0+1] - by[y-y0]);
		}

		processed++;
	}
	//cout << processed << "/" << cs.size() << endl;

	// if any left, use the b-spline recursive subdivision algorithm to efficiently display the remaining big basis functions
	if (processed < cs.size())
	{
		// determine x bounds for each lx
		vector<int> x0(nlx);
		vector<int> x1(nlx);
		vector<int> nx(nlx);
		for (int lx = 0; lx < nlx; lx++)
		{
			x0[lx] = (int) (pow(2.0, lx0+lx) * mz_min * 60.0/1.0033548378);
			x1[lx] = (int) (pow(2.0, lx0+lx) * mz_max * 60.0/1.0033548378) + 3;
			nx[lx] = x1[lx] - x0[lx] + 1;
		}

		// determine y bounds for each ly
		vector<int> y0(nly);
		vector<int> y1(nly);
		vector<int> ny(nly);
		for (int ly = 0; ly < nly; ly++)
		{
			y0[ly] = (int) (pow(2.0, ly0+ly) * rt_min);
			y1[ly] = (int) (pow(2.0, ly0+ly) * rt_max) + 3;
			ny[ly] = y1[ly] - y0[ly] + 1;
		}

		// index of mz levels (lx) in first dimension synthesis (ly)
		vector<int> ix(nlx+1);
		ix[0] = 0;
		for (int lx = 0; lx < nlx; lx++)
		{
			ix[lx+1] = ix[lx] + nx[lx]; 
		}

		// build sparse matrix triplet list of coefs Cy at each ly
		vector< SparseMatrix<double> > matCy(nly);
		{
			vector< std::vector< Triplet<double> > > listCy(nly);
			for (int i = 0; i < cs.size(); i++)
			if (cs[i].lx <= lx1 && cs[i].ly <= ly1)
			{
				// these are big basis functions compared to our viewport resolution, so we do the recursive subdivision
				int _i = cs[i].y-y0[cs[i].ly-ly0];
				int _j = ix[cs[i].lx-lx0] + cs[i].x-x0[cs[i].lx-lx0];
				double _v = cs[i].v;
				// some basis functions are not actually in the viewport - boundary condition problem
				if (_i >= 0 && _j >= 0 && _i < ny[cs[i].ly-ly0] && _j < ix.back()) // temp hack
				{
					listCy[cs[i].ly-ly0].push_back(Triplet<double>(_i,_j,_v));
				}
			}
	
			for (int ly = 0; ly < matCy.size(); ly++)
			{
				matCy[ly].resize(ny[ly], ix.back());
				matCy[ly].setFromTriplets(listCy[ly].begin(), listCy[ly].end());

				//SparseMatrix<double> tmp(ny[ly], ix.back());
				//tmp.setFromTriplets(listCy[ly].begin(), listCy[ly].end());
				//matCy[ly] = tmp;


				if (DEBUG)
				{
					ostringstream pngfile_oss; pngfile_oss << "Cy_" << ly0+ly << ".png";	
					filesystem::path pngfile_path = stream_path / pngfile_oss.str();
					writer.write(pngfile_path, matCy[ly], true);
				}
			}
		}

		// perform first dimension synthesis
		for (int ly = 0; ly < matCy.size()-1; ly++)
		{
			SparseMatrix<double> matAy(matCy[ly+1].rows(), matCy[ly].rows());
			vector< Triplet<double> > listAy;
			int offset = 3 + ((y0[ly] + 1) % 2);
			for (int j = 0; j < matAy.cols(); j++)
			for (int i = 0; i < 5; i++)
			{
				int _i = 2 * j + i - offset;
				if (_i >= 0 && _i < matAy.rows())
				{
					int _j = j;
					double _v = bspline::cbsf[i];
					listAy.push_back(Triplet<double>(_i,_j,_v));
				}
			}
			matAy.setFromTriplets(listAy.begin(), listAy.end());	

			matCy[ly+1] += matAy * matCy[ly];	

			if (DEBUG)
			{
				cout << "Cy[" << ly0+ly+1 << "](" << matCy[ly+1].rows() << "," << matCy[ly+1].cols() << ") += ";
				cout << "Ay[" << ly0+ly << "](" << matAy.rows() << "," << matAy.cols() << ") * ";
				cout << "Cy[" << ly0+ly << "](" << matCy[ly].rows() << "," << matCy[ly].cols() << ")" << endl;

				ostringstream pngfile_oss; pngfile_oss << "Ay_" << ly0+ly << "_" << ly0+ly+1 << ".png";
				filesystem::path pngfile_path = stream_path / pngfile_oss.str();
				writer.write(pngfile_path, matAy, true);

				ostringstream pngfile_oss2; pngfile_oss2 << "Fy_" << ly0+ly+1 << ".png";
				filesystem::path pngfile_path2 = stream_path / pngfile_oss2.str();
				writer.write(pngfile_path2, matCy[ly+1], true);
			}
		}
		if (DEBUG) cout << endl;

		// split Cy.back() into Cx for each lx, for second dimension synthesis
		vector< SparseMatrix<double> > matCx(nlx);
		for (int lx = 0; lx < nlx; lx++)
		{
			matCx[lx] = matCy.back().middleCols(ix[lx], ix[lx+1]-ix[lx]).transpose();

			if (DEBUG)
			{
				cout << "Cx[" << lx0+lx << "](" << matCx[lx].rows() << "," << matCx[lx].cols() << ") = ";
				cout << "t(Cy[" << ly0+(int)matCy.size()-1 << "](" << matCy.back().rows() << "," << ix[lx] << ":" <<  ix[lx+1]-1 << "))" << endl;

				ostringstream pngfile_oss; pngfile_oss << "Cx_" << lx0+lx << ".png";
				filesystem::path pngfile_path = stream_path / pngfile_oss.str();
				writer.write(pngfile_path, matCx[lx], true);
			}
		}
		if (DEBUG) cout << endl;

		// create sparse matrices Ax for knot insertion from lx[i] -> lx[i+1]
		for (int lx = 0; lx < matCx.size()-1; lx++)
		{
			SparseMatrix<double> matAx(matCx[lx+1].rows(), matCx[lx].rows());
			vector< Triplet<double> > listAx;
			int offset = 3 + ((x0[lx] + 1) % 2);
			for (int j = 0; j < matAx.cols(); j++)
			for (int i = 0; i < 5; i++)
			{
				int _i = 2 * j + i - offset;
				if (_i >= 0 && _i < matAx.rows())
				{
					int _j = j;
					double _v = bspline::cbsf[i];
					listAx.push_back(Triplet<double>(_i,_j,_v));
				}
			}
			matAx.setFromTriplets(listAx.begin(), listAx.end());

			matCx[lx+1] += matAx * matCx[lx];	

			if (DEBUG)
			{
				cout << "Cx[" << lx0+lx+1 << "](" << matCx[lx+1].rows() << "," << matCx[lx+1].cols() << ") += ";
				cout << "Ax[" << lx0+lx << "](" << matAx.rows() << "," << matAx.cols() << ") * ";
				cout << "Cx[" << lx0+lx << "](" << matCx[lx].rows() << "," << matCx[lx].cols() << ")" << endl;

				ostringstream pngfile_oss; pngfile_oss << "Ax_" << lx0+lx << "_" << lx0+lx+1 << ".png";
				filesystem::path pngfile_path = stream_path / pngfile_oss.str();
				writer.write(pngfile_path, matAx, true);

				ostringstream pngfile_oss2; pngfile_oss2 << "Fx_" << lx0+lx+1 << ".png";
				filesystem::path pngfile_path2 = stream_path / pngfile_oss2.str();
				writer.write(pngfile_path2, matCx[lx+1], true);
			}
		}
		if (DEBUG) cout << endl;

		// create sparse matrix for b-spline rebin Cx to display Dx dimension
		SparseMatrix<double> matDx(img.cols(), matCx.back().cols());
		{
			SparseMatrix<double> matBx(img.cols(), matCx.back().rows());
			vector< Triplet<double> > listBx;
			for (int j = 0; j < matBx.cols(); j++)
			{
				// bounding box for this basis function in mz/rt space
				double mz = pow(2.0, -lx1) * 1.0033548378/60.0;
				double mz0 = (x0.back()+j - 3) * mz;
				double mz1 = (x0.back()+j + 1) * mz;

				// project into viewport x/y space
				double x0f = img.cols() * (mz0 - mz_min) / (mz_max - mz_min);
				double x1f = img.cols() * (mz1 - mz_min) / (mz_max - mz_min);

				// bounding box intersection
				int x0 = (int) floor(x0f); x0 = x0 > 0 ? x0 : 0;
				int x1 = (int) ceil(x1f); x1 = x1 < img.cols() ? x1 : img.cols();			

				vector<double> bx(x1-x0+1);
				for (size_t x = 0; x < bx.size(); ++x)
				{
					bx[x] = bspline::iqbs(4.0 * (x + x0 - x0f) / (x1f - x0f));
				}

				for (int i = 0; i < x1-x0; i++)
				{
					int _i = x0+i;
					int _j = j;
					double _v = bx[i+1] - bx[i];
					listBx.push_back(Triplet<double>(_i,_j,_v));
				}
				//cout << mz0 << ":" << mz1 << " " << x0f << ":" << x1f << " " << x0 << ":" << x1 << endl;
			}		
			matBx.setFromTriplets(listBx.begin(), listBx.end());

			//SparseMatrix<double> tmp(img.cols(), matCx.back().cols());
			//tmp = matBx * matCx.back();
			//matDx = tmp;
			matDx = matBx * matCx.back();	

			if (DEBUG)
			{
				cout << "Dx(" << matDx.rows() << "," << matDx.cols() << ") = ";
				cout << "Bx(" << matBx.rows() << "," << matBx.cols() << ") * ";
				cout << "Cx[" << lx0+(int)matCx.size()-1 << "](" << matCx.back().rows() << "," << matCx.back().cols() << ")" << endl;

				ostringstream pngfile_oss; pngfile_oss << "Bx.png";
				filesystem::path pngfile_path = stream_path / pngfile_oss.str();
				writer.write(pngfile_path, matBx, true);

				ostringstream pngfile_oss2; pngfile_oss2 << "Dx.png";
				filesystem::path pngfile_path2 = stream_path / pngfile_oss2.str();
				writer.write(pngfile_path2, matDx, true);
			}
		}
		if (DEBUG) cout << endl;

		// create sparse matrix for b-spline rebin Cy to display Dy dimension
		{
			SparseMatrix<double> matBy(img.rows(), matDx.cols());
			vector< Triplet<double> > listBy;
			for (int j = 0; j < matBy.cols(); j++)
			{
				// bounding box for this basis function in mz/rt space
				double rt = pow(2.0, -ly1);
				double rt0 = (y0.back()+j - 3) * rt;
				double rt1 = (y0.back()+j + 1) * rt;

				// project into viewport x/y space
				double y0f = img.rows() * (rt0 - rt_min) / (rt_max - rt_min);
				double y1f = img.rows() * (rt1 - rt_min) / (rt_max - rt_min);

				// bounding box intersection
				int y0 = (int) floor(y0f); y0 = y0 > 0 ? y0 : 0;
				int y1 = (int) ceil(y1f); y1 = y1 < img.rows() ? y1 : img.rows();

				vector<double> by(y1-y0+1);
				for (int y = 0; y < by.size(); ++y)
				{
					by[y] = bspline::iqbs(4.0 * (y + y0 - y0f) / (y1f - y0f));
				}

				for (int i = 0; i < y1-y0; i++)
				{
					int _i = y0+i;
					int _j = j;
					double _v = by[i+1] - by[i];
					listBy.push_back(Triplet<double>(_i,_j,_v));
				}
				//cout << rt0 << ":" << rt1 << " " << y0f << ":" << y1f << " " << y0 << ":" << y1 << endl;
			}				
			matBy.setFromTriplets(listBy.begin(), listBy.end());

			SparseMatrix<double> matDy(img.rows(), img.cols());
			matDy = matBy * matDx.transpose();
			img += matDy;	
			//img += matBy * matDx.transpose();

			if (DEBUG)
			{
				cout << "Dy(" << img.rows() << "," << img.cols() << ") = ";
				cout << "By(" << matBy.rows() << "," << matBy.cols() << ") * ";
				cout << "t(Dx(" << matDx.rows() << "," << matDx.cols() << "))" << endl;

				ostringstream pngfile_oss; pngfile_oss << "By.png";
				filesystem::path pngfile_path = stream_path / pngfile_oss.str();
				writer.write(pngfile_path, matBy, true);

				ostringstream pngfile_oss2; pngfile_oss2 << "Dy.png";
				filesystem::path pngfile_path2 = stream_path / pngfile_oss2.str();
				writer.write(pngfile_path2, matDy, true);
			}
		}
		if (DEBUG) cout << endl;
	}

	// END CONSTRUCT VISUALISATION ////////////////////////////////////////////

	// write display timing
	viz_time += omp_get_wtime() - viz_start;
	stream_ofs << viz_time << endl;

	// write png
	ostringstream pngfile_oss; pngfile_oss << setfill('0') << setw(8) << chunk_index << ".png";
	filesystem::path pngfile_path = stream_path / pngfile_oss.str();
	writer.write(pngfile_path, img, max_counts);

	// increment chunk index
	chunk_index++;
	read_start = omp_get_wtime();

	//cout << nx << "," << ny << "," << nxy << endl;
	if (DEBUG) exit(0);
}
