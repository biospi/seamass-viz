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

#include "aux_Eigen.hpp"

#include <iostream>
#include <sstream>
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>
#include <math.h>

// http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale

COLOUR GetColour(fp v,fp vmin,fp vmax)
{
   COLOUR c = {1.0f,1.0f,1.0f}; // white
   fp dv;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25f * dv)) {
      c.r = 0;
      c.g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5f * dv)) {
      c.r = 0;
      c.b = 1 + 4 * (vmin + 0.25f * dv - v) / dv;
   } else if (v < (vmin + 0.75f * dv)) {
      c.r = 4 * (v - vmin - 0.5f * dv) / dv;
      c.b = 0;
   } else {
      c.g = 1 + 4 * (vmin + 0.75f * dv - v) / dv;
      c.b = 0;
   }

   return(c);
}

// end http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale
