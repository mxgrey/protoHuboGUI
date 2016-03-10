/*
 * Copyright (c) 2016, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Michael X. Grey <mxgrey@gatech.edu>
 *
 * Georgia Tech Graphics Lab and AMBER Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Aaron Ames
 * <karenliu@cc.gatech.edu> <ames@gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HUBO_DRCMODEL_HPP
#define HUBO_DRCMODEL_HPP

#include <deque>

#include <dart/dynamics/Skeleton.h>
#include <dart/simulation/World.h>

#include "hubo/config.hpp"

namespace hubo {
namespace DrcModel {

//==============================================================================
dart::dynamics::SkeletonPtr create(
    const std::string& modelFile = "drchubo.urdf",
    const std::string& modelPath = DRCHUBO_MODEL_PATH,
    bool lockedFeet = false);

//==============================================================================
void goodPosture(const dart::dynamics::SkeletonPtr& hubo);

//==============================================================================
void laxPosture(const dart::dynamics::SkeletonPtr& hubo);

//==============================================================================
void disableAdjacentPairs(const dart::dynamics::SkeletonPtr& hubo,
                          const dart::simulation::WorldPtr& world);

//==============================================================================
void removeFingers(const dart::dynamics::SkeletonPtr& hubo);

//==============================================================================
void removeHead(const dart::dynamics::SkeletonPtr& hubo);

//==============================================================================
std::vector<size_t> getUsefulDofIndices(
    const dart::dynamics::ConstSkeletonPtr& hubo);

//==============================================================================
std::deque<Eigen::VectorXd> getUsefulSeeds(const dart::dynamics::SkeletonPtr& robot);

} // namespace DrcModel
} // namespace hubo

#endif // HUBO_DRCMODEL_HPP
