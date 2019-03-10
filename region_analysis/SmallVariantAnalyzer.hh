//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#pragma once

#include <boost/optional.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "classification/ClassifierOfAlignmentsToVariant.hh"
#include "genotyping/AlleleChecker.hh"
#include "genotyping/SmallVariantGenotyper.hh"
#include "region_analysis/VariantAnalyzer.hh"
#include "region_spec/VariantSpecification.hh"

namespace ehunter
{

class SmallVariantAnalyzer : public VariantAnalyzer
{
public:
    SmallVariantAnalyzer(
        std::string variantId, VariantSubtype variantSubtype, AlleleCount expectedAlleleCount,
        const graphtools::Graph& graph, std::vector<graphtools::NodeId> nodeIds,
        boost::optional<graphtools::NodeId> optionalRefNode, const GenotyperParameters& params)
        : VariantAnalyzer(std::move(variantId), expectedAlleleCount, graph, std::move(nodeIds))
        , variantSubtype_(variantSubtype)
        , optionalRefNode_(optionalRefNode)
        , alignmentClassifier_(nodeIds_)
        , allelePresenceChecker_(params.errorRate, params.likelihoodRatioThreshold)
        , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))
    {
        // Only indels are allowed
        assert(nodeIds_.size() <= 2);
    }

    ~SmallVariantAnalyzer() = default;

    std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) const override;

    void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

protected:
    int countReadsSupportingNode(graphtools::NodeId nodeId) const;

    VariantSubtype variantSubtype_;
    boost::optional<graphtools::NodeId> optionalRefNode_;

    ClassifierOfAlignmentsToVariant alignmentClassifier_;
    AlleleChecker allelePresenceChecker_;

    std::shared_ptr<spdlog::logger> console_;
};

}
