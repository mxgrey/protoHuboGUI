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

#ifndef HUBO_ARMIK_HPP
#define HUBO_ARMIK_HPP

#include <dart/dynamics/InverseKinematics.h>

namespace hubo {

using InverseKinematics = dart::dynamics::InverseKinematics;

class ArmIK : public InverseKinematics::Analytical
{
public:

  ArmIK(InverseKinematics* _ik, const std::string& baseLinkName,
            const Analytical::Properties& properties = Analytical::Properties());

  std::unique_ptr<GradientMethod> clone(InverseKinematics* _newIK) const override;

  const std::vector<Solution>& computeSolutions(
      const Eigen::Isometry3d& _desiredBodyTf);

  const std::vector<size_t>& getDofs() const override;

  const double zeroSize = 1e-8;

protected:

  void configure() const;

  mutable bool configured;

  mutable Eigen::Isometry3d shoulderTf;
  mutable Eigen::Isometry3d wristTfInv;
  mutable Eigen::Isometry3d mNodeOffsetTfInv;
  mutable double L3, L4, L5;

  mutable Eigen::Matrix<int, 8, 3> alterantives;

  mutable std::vector<size_t> mDofs;

  std::string mBaseLinkName;
  mutable dart::dynamics::WeakBodyNodePtr mBaseLink;

  mutable dart::dynamics::JacobianNode* mWristEnd;
};

}

#endif // HUBO_ARMIK_HPP
