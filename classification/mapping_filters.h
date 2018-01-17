//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <cstdint>

#include "graphs/graph_mapping.h"

bool CheckIfStartAndEndMapWell(GraphMapping mapping,
                               int32_t num_bases_to_rescore,
                               const int32_t kHighMatchCount = 8,
                               const double kHighMatchProp = 0.8,
                               const double kLowMatchProp = 0.6);

bool CheckIfStartMapsWell(GraphMapping mapping, int32_t num_bases_to_rescore,
                          const int32_t kHighMatchCount,
                          const double kHighMatchProp,
                          const double kLowMatchProp);

bool CheckIfEndMapsWell(GraphMapping mapping, int32_t num_bases_to_rescore,
                        const int32_t kHighMatchCount,
                        const double kHighMatchProp,
                        const double kLowMatchProp);