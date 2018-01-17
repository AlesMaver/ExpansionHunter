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

#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "graphs/gapless_aligner.h"
#include "graphs/path.h"
#include "graphs/path_operations.h"

using std::list;
using std::string;
using std::vector;

bool CheckIfStartAndEndMapWell(GraphMapping mapping,
                               int32_t num_bases_to_rescore,
                               const int32_t kHighMatchCount,
                               const double kHighMatchProp,
                               const double kLowMatchProp) {
  const bool is_start_maps_well =
      CheckIfStartMapsWell(mapping, num_bases_to_rescore, kHighMatchCount,
                           kHighMatchProp, kLowMatchProp);
  const bool is_end_maps_well =
      CheckIfEndMapsWell(mapping, num_bases_to_rescore, kHighMatchCount,
                         kHighMatchProp, kLowMatchProp);
  return is_start_maps_well && is_end_maps_well;
}

static double CalcProportionOfMatches(int32_t num_matches, size_t total_len) {
  return static_cast<double>(num_matches) / total_len;
}

static bool IsOneScoreHighAndRestLow(vector<int32_t>& num_matches_at_ends,
                                     int32_t max_num_matches,
                                     const int32_t kHighMatchCount,
                                     const double kHighMatchProp,
                                     const double kLowMatchProp) {
  std::sort(num_matches_at_ends.begin(), num_matches_at_ends.end());

  const int32_t max_matches = num_matches_at_ends.back();

  int32_t second_most_matches = 0;
  if (num_matches_at_ends.size() > 1) {
    second_most_matches = *(num_matches_at_ends.end() - 2);
  }

  const bool best_mapping_has_many_matches = max_matches >= kHighMatchCount;
  const bool best_mapping_has_high_prop_matches =
      CalcProportionOfMatches(max_matches, max_num_matches) >= kHighMatchProp;

  const bool second_best_mapping_has_low_prop_matches =
      CalcProportionOfMatches(second_most_matches, max_num_matches) <
      kLowMatchProp;

  return best_mapping_has_many_matches && best_mapping_has_high_prop_matches &&
         second_best_mapping_has_low_prop_matches;
}

vector<int32_t> AlignToEachPath(const list<GraphPath>& paths,
                                const string& bases) {
  vector<int32_t> num_matches_at_ends;
  for (const GraphPath& path : paths) {
    GraphMapping mapping = AlignWithoutGaps(path, bases);
    num_matches_at_ends.push_back(mapping.NumMatches());
  }
  return num_matches_at_ends;
}

bool CheckIfStartMapsWell(GraphMapping mapping, int32_t num_bases_to_rescore,
                          const int32_t kHighMatchCount,
                          const double kHighMatchProp,
                          const double kLowMatchProp) {
  list<GraphPath> left_endings =
      ComputeLeftEndings(mapping.Path(), num_bases_to_rescore - 1);
  if (left_endings.empty()) {
    return false;
  }

  const string& query = mapping.Query();
  const string query_piece = query.substr(0, num_bases_to_rescore);

  vector<int32_t> num_matches_at_ends =
      AlignToEachPath(left_endings, query_piece);

  return IsOneScoreHighAndRestLow(num_matches_at_ends, query_piece.length(),
                                  kHighMatchCount, kHighMatchProp,
                                  kLowMatchProp);
}

bool CheckIfEndMapsWell(GraphMapping mapping, int32_t num_bases_to_rescore,
                        const int32_t kHighMatchCount,
                        const double kHighMatchProp,
                        const double kLowMatchProp) {
  list<GraphPath> right_endings =
      ComputeRightEndings(mapping.Path(), num_bases_to_rescore - 1);
  if (right_endings.empty()) {
    return false;
  }

  const string& query = mapping.Query();
  const string query_piece =
      query.substr(query.length() - num_bases_to_rescore, num_bases_to_rescore);

  vector<int32_t> num_matches_at_ends =
      AlignToEachPath(right_endings, query_piece);

  return IsOneScoreHighAndRestLow(num_matches_at_ends, query_piece.length(),
                                  kHighMatchCount, kHighMatchProp,
                                  kLowMatchProp);
}
