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

#pragma once

#include <iostream>
#include <sstream>
#include <math.h>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>


using namespace std;
using namespace boost;
using namespace Eigen;


// fp is the desired floating point precision
typedef float fp; 
//typedef double fp; 


template <class fp>
class Triplets
{
protected:
	vector< Triplet<fp> > trips;

public:
	Triplets() {}

	void push_back(int i, int j, fp v)
	{
		trips.push_back(Triplet<fp>(i,j,v));
	}

	void init(SparseMatrix<fp>& mat)
	{
		mat.setFromTriplets(trips.begin(), trips.end());
	}
};


template <class fp> 
class PNGWriter
{
public:
	PNGWriter() {}

	void write(const filesystem::path& path, Matrix<fp,Dynamic,Dynamic>& mat, fp max_counts, bool show_sparsity = false);
	void write(const filesystem::path& path, Matrix<fp,Dynamic,Dynamic>& mat, bool show_sparsity = false);
	
	void write(const filesystem::path& path, SparseMatrix<fp>& mat, fp max_counts, bool show_sparsity = false);
	void write(const filesystem::path& path, SparseMatrix<fp>& mat, bool show_sparsity = false);
};


///////////////////////////////////////////////////////////////////////////////


typedef struct {
    fp r,g,b;
} COLOUR;

COLOUR GetColour(fp v,fp vmin,fp vmax);


template <class fp> 
void
PNGWriter<fp>::
write(const filesystem::path& filename, Matrix<fp,Dynamic,Dynamic>& mat, fp max_counts, bool show_sparsity)
{
	fp min_intensity = log(2.0f*sqrt(3.0f/8.0f));
	fp max_intensity = log(2.0f*sqrt(max_counts+3.0f/8.0f));

	gil::rgb8_image_t png_img(mat.cols(), mat.rows());
	gil::rgb8_view_t png_view = gil::view(png_img);
	gil::rgb8_pixel_t bg(0, 0, show_sparsity ? 0 : 255);
    gil::fill_pixels(png_view, bg);
	
	for (int k = 0; k < mat.outerSize(); ++k)
	for (Matrix<fp,Dynamic,Dynamic>::InnerIterator it(mat,k); it; ++it)
	{
		COLOUR c = GetColour(log(2.0f*sqrt(it.value()+3.0f/8.0f)), min_intensity, max_intensity);
		png_view(it.col(), it.row()) = gil::rgb8_pixel_t((unsigned char) (c.r*255.9999f), (unsigned char) (c.g*255.9999f), (unsigned char) (c.b*255.9999f));
	}
    gil::png_write_view(filename.string(), gil::const_view(png_img));
}


template <class fp> 
void
PNGWriter<fp>::
write(const filesystem::path& filename, Matrix<fp,Dynamic,Dynamic>& mat, bool show_sparsity)
{
	fp max_counts = 0.0;
	for (int k = 0; k < mat.outerSize(); ++k)
	for (Matrix<fp,Dynamic,Dynamic>::InnerIterator it(mat,k); it; ++it)
	{
		max_counts = max_counts > it.value() ? max_counts : it.value();
	}
	write(filename, mat, max_counts, show_sparsity);
}


template <class fp> 
void
PNGWriter<fp>::
write(const filesystem::path& filename, SparseMatrix<fp>& mat, fp max_counts, bool show_sparsity)
{
	fp max_intensity = log(2.0f*sqrt(max_counts+3.0f/8.0f));

	gil::rgb8_image_t png_img(mat.cols(), mat.rows());
	gil::rgb8_view_t png_view = gil::view(png_img);
	gil::rgb8_pixel_t bg(0, 0, show_sparsity ? 0 : 255);
    gil::fill_pixels(png_view, bg);
	
	for (int k = 0; k < mat.outerSize(); ++k)
	for (SparseMatrix<fp>::InnerIterator it(mat,k); it; ++it)
	{
		COLOUR c = GetColour(log(2.0f*sqrt(it.value()+3.0f/8.0f)), 0.0, max_intensity);
		png_view(it.col(), it.row()) = gil::rgb8_pixel_t((unsigned char) (c.r*255.9999f), (unsigned char) (c.g*255.9999f), (unsigned char) (c.b*255.9999f));
	}
    gil::png_write_view(filename.string(), gil::const_view(png_img));
}


template <class fp> 
void
PNGWriter<fp>::
write(const filesystem::path& filename, SparseMatrix<fp>& mat, bool show_sparsity)
{
	fp max_counts = 0.0;
	for (int k = 0; k < mat.outerSize(); ++k)
	for (SparseMatrix<fp>::InnerIterator it(mat,k); it; ++it)
	{
		max_counts = max_counts > it.value() ? max_counts : it.value();
	}
	write(filename, mat, max_counts, show_sparsity);
}