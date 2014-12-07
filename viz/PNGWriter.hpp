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

#include <vector>
#include <boost/filesystem.hpp>
#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace std;
using namespace boost;
using namespace Eigen;

class PNGWriter
{
public:
	PNGWriter() {}

	void write(const filesystem::path& path, const Matrix<double,Dynamic,Dynamic>& mat, double max_counts, bool show_sparsity = false);
	void write(const filesystem::path& path, const Matrix<double,Dynamic,Dynamic>& mat, bool show_sparsity = false);
	
	void write(const filesystem::path& path, const SparseMatrix<double>& mat, double max_counts, bool show_sparsity = false);
	void write(const filesystem::path& path, const SparseMatrix<double>& mat, bool show_sparsity = false);
};