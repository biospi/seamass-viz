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

#include <sstream>
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <math.h>

void
PNGWriter::
write(const filesystem::path& filename, const vector<double>& img, int w, int h, double max_counts)
{
	double max_intensity = log(2.0*sqrt(max_counts+3.0/8.0));

	gil::rgb8_image_t png_img(w, h);
	for (int y = 0; y < h; ++y)
	{
		gil::rgb8_view_t::x_iterator it = gil::view(png_img).row_begin(y);
		for (int x = 0; x < w; ++x)
		{
			double s = log(2.0*sqrt(img[x+y*w]+3.0/8.0)) / max_intensity;
			//s = s > 0.0 ? s : 0.0;
			//s = s < 1.0 ? s : 1.0;
			s = 255 - floor(255 * s);

			it[x][0] = s;
			it[x][1] = s;
			it[x][2] = s;
		}
	}
    gil::png_write_view(filename.string(), gil::const_view(png_img));
}
