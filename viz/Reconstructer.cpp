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

#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <math.h>


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

}


Reconstructer::
Reconstructer(const string& out_dir, int width, int height) :
	out_path(out_dir),
	img(width * height),
	w(width),
	h(height),
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
	double x_res = (mz_max - mz_min)/w/1.0033548378/60.0;
	double y_res = (rt_max - rt_min)/h;
	//cout << x_res << ":" << y_res << endl;
	max_counts = _max_counts / (x_res * y_res); 

	ostringstream oss; oss << stream_index;
	stream_path = out_path / oss.str();
	filesystem::create_directory(stream_path);

	filesystem::path file_path = stream_path / "stream.csv";	
	stream_ofs.open(file_path.string().c_str());

	img.assign(w*h, 0.0);
	stream_index++;
	chunk_index = 0;
	chunk_count = 0;

	read_start = omp_get_wtime();
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

	// add new basis functions to image
	for (size_t i = 0; i < cs.size(); i++)
	{
		// bounding box for this basis function in mz/rt space
		double mz = pow(2.0, -cs[i].lx) * 1.0033548378/60.0;
		double mz0 = (cs[i].x - 3) * mz;
		double mz1 = (cs[i].x + 1) * mz;
		double rt = pow(2.0, -cs[i].ly);
		double rt0 = (cs[i].y - 3) * rt;
		double rt1 = (cs[i].y + 1) * rt;

		// project into viewpoint x/y space
		double x0f = w * (mz0 - mz_min) / (mz_max - mz_min);
		double x1f = w * (mz1 - mz_min) / (mz_max - mz_min);
		double y0f = h * (rt0 - rt_min) / (rt_max - rt_min);
		double y1f = h * (rt1 - rt_min) / (rt_max - rt_min);

		// i-spline in x-axis 
		int x0 = (int) floor(x0f); x0 = x0 > 0 ? x0 : 0;
		int x1 = (int) ceil(x1f); x1 = x1 < w ? x1 : w;
		
		vector<double> bx(x1-x0+1);
		for (size_t x = 0; x < bx.size(); ++x)
		{
			bx[x] = bspline::iqbs(4.0 * (x + x0 - x0f) / (x1f - x0f));
		}

		// i-spline in y-axis
		int y0 = (int) floor(y0f); y0 = y0 > 0 ? y0 : 0;
		int y1 = (int) ceil(y1f); y1 = y1 < h ? y1 : h;
		
		vector<double> by(y1-y0+1);
		for (int y = 0; y < by.size(); ++y)
		{
			by[y] = bspline::iqbs(4.0 * (y + y0 - y0f) / (y1f - y0f));
		}

		// add separable b-spline integral to image
		for (int y = y0; y < y1; ++y)
		for (int x = x0; x < x1; ++x)
		{
			img[x+y*w] += cs[i].v * (bx[x-x0+1] - bx[x-x0]) * (by[y-y0+1] - by[y-y0]);
		}
	}

	// write display timing
	viz_time += omp_get_wtime() - viz_start;
	stream_ofs << viz_time << endl;

	// write png
	ostringstream pngfile_oss; pngfile_oss << setfill('0') << setw(8) << chunk_index << ".png";
	filesystem::path pngfile_path = stream_path / pngfile_oss.str();
	writer.write(pngfile_path, img, w, h, max_counts);

	// increment chunk index
	chunk_index++;
	read_start = omp_get_wtime();
}
