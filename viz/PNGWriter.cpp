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

#include "PNGWriter.hpp"

#include <iostream>
#include <sstream>
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <math.h>

// http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale

typedef struct {
    double r,g,b;
} COLOUR;

COLOUR GetColour(double v,double vmin,double vmax)
{
   COLOUR c = {1.0,1.0,1.0}; // white
   double dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      c.r = 0;
      c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
      c.r = 0;
      c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
      c.r = 4 * (v - vmin - 0.5 * dv) / dv;
      c.b = 0;
   } else {
      c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      c.b = 0;
   }

   return(c);
}

// end http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale


void
PNGWriter::
write(const filesystem::path& filename, Matrix<double,Dynamic,Dynamic>& mat, double max_counts, bool show_sparsity)
{
	double min_intensity = log(2.0*sqrt(3.0/8.0));
	double max_intensity = log(2.0*sqrt(max_counts+3.0/8.0));

	gil::rgb8_image_t png_img(mat.cols(), mat.rows());
	gil::rgb8_view_t png_view = gil::view(png_img);
	gil::rgb8_pixel_t bg(0, 0, show_sparsity ? 0 : 255);
    gil::fill_pixels(png_view, bg);
	
	for (int k = 0; k < mat.outerSize(); ++k)
	for (Matrix<double,Dynamic,Dynamic>::InnerIterator it(mat,k); it; ++it)
	{
		COLOUR c = GetColour(log(2.0*sqrt(it.value()+3.0/8.0)), min_intensity, max_intensity);
		png_view(it.col(), it.row()) = gil::rgb8_pixel_t(c.r*255.9999, c.g*255.9999, c.b*255.9999);
	}
    gil::png_write_view(filename.string(), gil::const_view(png_img));
}


void
PNGWriter::
write(const filesystem::path& filename, Matrix<double,Dynamic,Dynamic>& mat, bool show_sparsity)
{
	double max_counts = 0.0;
	for (int k = 0; k < mat.outerSize(); ++k)
	for (Matrix<double,Dynamic,Dynamic>::InnerIterator it(mat,k); it; ++it)
	{
		max_counts = max_counts > it.value() ? max_counts : it.value();
	}
	write(filename, mat, max_counts, show_sparsity);
}


void
PNGWriter::
write(const filesystem::path& filename, SparseMatrix<double>& mat, double max_counts, bool show_sparsity)
{
	double max_intensity = log(2.0*sqrt(max_counts+3.0/8.0));

	gil::rgb8_image_t png_img(mat.cols(), mat.rows());
	gil::rgb8_view_t png_view = gil::view(png_img);
	gil::rgb8_pixel_t bg(0, 0, show_sparsity ? 0 : 255);
    gil::fill_pixels(png_view, bg);
	
	for (int k = 0; k < mat.outerSize(); ++k)
	for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
	{
		COLOUR c = GetColour(log(2.0*sqrt(it.value()+3.0/8.0)), 0.0, max_intensity);
		png_view(it.col(), it.row()) = gil::rgb8_pixel_t(c.r*255.9999, c.g*255.9999, c.b*255.9999);
	}
    gil::png_write_view(filename.string(), gil::const_view(png_img));
}


void
PNGWriter::
write(const filesystem::path& filename, SparseMatrix<double>& mat, bool show_sparsity)
{
	double max_counts = 0.0;
	for (int k = 0; k < mat.outerSize(); ++k)
	for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
	{
		max_counts = max_counts > it.value() ? max_counts : it.value();
	}
	write(filename, mat, max_counts, show_sparsity);
}