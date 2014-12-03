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

#include <string>
#include <vector>
#include <SpatialIndex.h>

#include "Reconstructer.hpp"

using namespace std;
using namespace SpatialIndex;

class SMVStreamer
{
protected:
	string basename;
	IStorageManager* diskfile;
	StorageManager::IBuffer* file;
	ISpatialIndex* tree;

	double mz_min, mz_max, rt_min, rt_max, counts_max;

public:
	SMVStreamer(const string& filename);
	~SMVStreamer();
    
    void stream_cs(double mz_min, double mz_max, double rt_min, double rt_max,
		           size_t chunk_size, Reconstructer* recon);
};
