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

#include <dart/dynamics/BodyNode.h>
#include <dart/dynamics/DegreeOfFreedom.h>

#include "hubo/ArmIK.hpp"
#include "hubo/utils.hpp"

using namespace dart;
using namespace dynamics;

namespace hubo {

//==============================================================================
ArmIK::ArmIK(InverseKinematics* _ik, const std::string& baseLinkName,
             const Analytical::Properties& properties)
  : Analytical(_ik, "HuboArmIK_"+baseLinkName, properties),
    configured(false),
    mBaseLinkName(baseLinkName)
{
  // Do nothing
}

//==============================================================================
std::unique_ptr<InverseKinematics::GradientMethod> ArmIK::clone(
    InverseKinematics* _newIK) const
{
  return std::unique_ptr<GradientMethod>(
        new ArmIK(_newIK, mBaseLinkName, getAnalyticalProperties()));
}

//==============================================================================
const std::vector<ArmIK::Solution>&
ArmIK::computeSolutions(const Eigen::Isometry3d& _desiredBodyTf)
{
  mSolutions.clear();
  mSolutions.reserve(8);

  if(!configured)
  {
    configure();

    if(!configured)
    {
      dtwarn << "[HuboArmIK::computeSolutions] This analytical IK was not able "
             << "to configure properly, so it will not be able to compute "
             << "solutions\n";
      return mSolutions;
    }
  }

  const BodyNodePtr& base = mBaseLink.lock();
  if(nullptr == base)
  {
    dterr << "[HuboArmIK::computeSolutions] Attempting to perform an IK on a "
          << "limb that no longer exists [" << getMethodName() << "]!\n";
    assert(false);
    return mSolutions;
  }

  if(nullptr == mWristEnd)
  {
    dterr << "[HuboArmIK::computeSolutions] Attempting to perform IK without "
          << "a wrist!\n";
    assert(false);
    return mSolutions;
  }

  const size_t SP = 0;
  const size_t SR = 1;
  const size_t SY = 2;
  const size_t EP = 3;
  const size_t WY = 4;
  const size_t WP = 5;

  const SkeletonPtr& skel = base->getSkeleton();

  Eigen::Isometry3d B =
      base->getParentBodyNode()->getWorldTransform().inverse()
      * _desiredBodyTf * mWristEnd->getTransform(mIK->getNode());

  Eigen::Isometry3d shoulder_from_wrist = shoulderTf.inverse() * B;
  Eigen::Vector3d p = shoulder_from_wrist.inverse().translation();

  const double a2 = L5*L5 + L4*L4;
  const double b2 = L3*L3 + L4*L4;
  const double a = sqrt(a2);
  const double b = sqrt(b2);

  const double alpha = atan2(L5, L4);
  const double beta = atan2(L3, L4);

  bool startValid = checkDist(p, a, b);

  double c2 = p.dot(p);
  double x = p.x();
  double y = p.y();
  double z = p.z();

  for(size_t i = 0; i < 8; ++i)
  {
    const int flipEP = alterantives(i,0);
    const int incWY = alterantives(i,1);
    const int flipShoulder = alterantives(i,2);

    Eigen::Vector6d testQ;
    bool isValid = startValid;

    double cosGamma = (a2 + b2 - c2) / (2*a*b);
    clamp_sincos(cosGamma, isValid);

    double gamma = flipEP * acos( cosGamma );
    double theta3 = alpha + beta + gamma - 2*M_PI;

    testQ(EP) = theta3;

    double c3 = cos(theta3);
    double s3 = sin(theta3);

    double numer = -y;
    double denom = (-L4*c3 - L3*s3 + L4);

    double s2, theta2;

    if(std::abs(denom) < zeroSize)
    {
      isValid = false;
      const double& prevWY = skel->getPosition(mDofs[WY]);
      theta2 = incWY ? prevWY : M_PI - prevWY;
      s2 = sin(theta2);
    }
    else
    {
      s2 = numer / denom;
      clamp_sincos(s2, isValid);
      theta2 = incWY ? M_PI - asin(s2) : asin(s2);
    }

    testQ(WY) = theta2;

    double c2 = cos(theta2);

    double r =  L4*c2 - L4*c2*c3 - L3*s3*c2;
    double q = -L4*s3 + L3*c3 + L5;

    double det = -(q*q + r*r);

    if(std::abs(det) < zeroSize)
      isValid = false;

    double k = det < 0 ? -1 : 1;

    double ks1 = k*( q*x - r*z );
    double kc1 = k*(-r*x - q*z );

    double theta1 = atan2(ks1, kc1);
    testQ(WP) = theta1;

    Eigen::Quaterniond Rlower =
      Eigen::Quaterniond(Eigen::AngleAxisd(testQ(EP), Eigen::Vector3d::UnitY())) *
      Eigen::Quaterniond(Eigen::AngleAxisd(testQ(WY), Eigen::Vector3d::UnitZ())) *
      Eigen::Quaterniond(Eigen::AngleAxisd(testQ(WP), Eigen::Vector3d::UnitY()));

    Eigen::Matrix3d Rupper = B.rotation() * Rlower.inverse().matrix();

    Eigen::Vector3d euler = Rupper.eulerAngles(1, 0, 2);

    if(flipShoulder)
      euler = flipEuler3Axis(euler);

    testQ(SP) = euler[0];
    testQ(SR) = euler[1];
    testQ(SY) = euler[2];

    for(size_t j=0; j < 6; ++j)
    {
      testQ[j] = dart::math::wrapToPi(testQ[j]);
      if(std::abs(testQ[j]) < zeroSize)
        testQ[j] = 0.0;
    }

    int validity = isValid? VALID : OUT_OF_REACH;
    mSolutions.push_back(Solution(testQ, validity));
  }

  checkSolutionJointLimits();

  return mSolutions;
}

//==============================================================================
const std::vector<size_t>& ArmIK::getDofs() const
{
  if(!configured)
    configure();

  return mDofs;
}

//==============================================================================
void ArmIK::configure() const
{
  configured = false;

  mBaseLink = mIK->getNode()->getSkeleton()->getBodyNode(mBaseLinkName);

  BodyNode* base = mBaseLink.lock();
  if(nullptr == base)
  {
    dterr << "[HuboArmIK::configure] base link is a nullptr\n";
    assert(false);
    return;
  }

  const SkeletonPtr& skel = base->getSkeleton();
  const BodyNodePtr& pelvis = skel->getBodyNode("Body_TSY");
  if(nullptr == pelvis)
  {
    dterr << "[HuboArmIK::configure] Could not find Hubo's pelvis "
          << "(Body_TSY)\n";
    assert(false);
    return;
  }

  Eigen::Vector6d saved_q;

  DegreeOfFreedom* dofs[6];
  BodyNode* bn = base;
  for(size_t i=0; i < 6; ++i)
  {
    Joint* joint = bn->getParentJoint();
    if(joint->getNumDofs() != 1)
    {
      dterr << "[HuboArmIK::configure] Invalid number of DOFs ("
            << joint->getNumDofs() << ") in the Joint [" << joint->getName()
            << "]\n";
      assert(false);
      return;
    }

    dofs[i] = joint->getDof(0);
    saved_q[i] = dofs[i]->getPosition();
    dofs[i]->setPosition(0.0);
    bn = bn->getChildBodyNode(0);
  }

  BodyNode* elbow = dofs[3]->getChildBodyNode();
  L3 = std::abs(elbow->getTransform(dofs[2]->getParentBodyNode()).translation()[2]);
  L4 = std::abs(elbow->getTransform(dofs[3]->getParentBodyNode()).translation()[0]);

  BodyNode* wrist = dofs[5]->getChildBodyNode();
  Eigen::Isometry3d wrist_tf = wrist->getTransform(elbow);
  L5 = std::abs(wrist_tf.translation()[2]);

  shoulderTf = Eigen::Isometry3d::Identity();
  shoulderTf.translate(dofs[3]->getParentBodyNode()->getTransform(pelvis)
      .translation()[0] * Eigen::Vector3d::UnitX());
  shoulderTf.translate(dofs[2]->getParentBodyNode()->getTransform(pelvis)
      .translation()[1] * Eigen::Vector3d::UnitY());
  shoulderTf.translate(dofs[2]->getParentBodyNode()->getTransform(pelvis)
      .translation()[2] * Eigen::Vector3d::UnitZ());

  mWristEnd = dofs[5]->getChildBodyNode();

  alterantives <<
       1,  1,  1,
       1,  1,  0,
       1,  0,  1,
       1,  0,  0,
      -1,  1,  1,
      -1,  1,  0,
      -1,  0,  1,
      -1,  0,  0;

  for(size_t i=0; i < 6; ++i)
  {
    dofs[i]->setPosition(saved_q[i]);
    mDofs.push_back(dofs[i]->getIndexInSkeleton());
  }

  configured = true;
}

} // namespace hubo
