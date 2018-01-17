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

#include "classification/mapping_filters.h"

#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"
#include "graphs/path.h"

using std::string;

TEST(CheckingIfStartAndEndMapWell, MappingWithIffyEnd_FilterFailed) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TAAT", "CCG", "CCTT");
  //                   FFFRRRFFF
  const string read = "AATCCGCCT";
  const GraphMapping spanning_mapping =
      DecodeFromString(1, "0[3M]1[3M]2[3M]", read, graph_ptr);

  ASSERT_FALSE(CheckIfStartAndEndMapWell(spanning_mapping, 4, 3));
}

TEST(CheckingIfStartAndEndMapWell, MappingWithIffyStart_FilterFailed) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TCTG", "CCG", "CCTT");
  //                   FFFRRRFFFF
  const string read = "CTGCCGCCTT";
  const GraphMapping spanning_mapping =
      DecodeFromString(1, "0[3M]1[3M]2[3M]", read, graph_ptr);

  ASSERT_FALSE(CheckIfStartAndEndMapWell(spanning_mapping, 4, 3));
}

TEST(CheckingIfStartAndEndMapWell, MappingWithConfidentEnds_FilterPass) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TAAT", "CCG", "TAAT");
  //                   FFFRRRFFF
  const string read = "AATCCGTAA";
  const GraphMapping spanning_mapping =
      DecodeFromString(1, "0[3M]1[3M]2[3M]", read, graph_ptr);

  ASSERT_TRUE(CheckIfStartAndEndMapWell(spanning_mapping, 4, 3));
}
