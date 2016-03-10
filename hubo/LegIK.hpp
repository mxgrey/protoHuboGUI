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

#ifndef HUBO_LEGIK_HPP
#define HUBO_LEGIK_HPP

#include <dart/dynamics/InverseKinematics.h>

namespace hubo {

using InverseKinematics = dart::dynamics::InverseKinematics;

class LegIK : public InverseKinematics::Analytical
{
public:

  /// baseLink should be Body_LHY or Body_RHY
  LegIK(InverseKinematics* _ik, const std::string& baseLinkName,
            const Analytical::Properties& properties = Analytical::Properties());

  std::unique_ptr<GradientMethod> clone(InverseKinematics* _newIK) const override;

  const std::vector<Solution>& computeSolutions(
      const Eigen::Isometry3d& _desiredBodyTf) override;

  const std::vector<size_t>& getDofs() const override;

  const double zeroSize = 1e-8;

protected:

  void configure() const;

  mutable double L4, L5, L6;
  mutable Eigen::Isometry3d waist;
  mutable Eigen::Isometry3d hipRotation;
  mutable Eigen::Isometry3d footTfInv;
  mutable Eigen::Matrix<int, 8, 3> alternatives;

  mutable std::vector<size_t> mDofs;

  mutable bool configured;

  std::string mBaseLinkName;

  mutable dart::dynamics::WeakBodyNodePtr mBaseLink;

};

} // namespace hubo

#endif // HUBO_LEGIK_HPP
