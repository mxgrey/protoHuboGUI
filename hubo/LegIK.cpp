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

#include "hubo/LegIK.hpp"
#include "hubo/utils.hpp"

using namespace dart;
using namespace dynamics;

namespace hubo {

//==============================================================================
LegIK::LegIK(InverseKinematics* _ik, const std::string& baseLinkName,
          const Analytical::Properties& properties)
  : Analytical(_ik, "HuboLegIK_"+baseLinkName, properties),
    configured(false),
    mBaseLinkName(baseLinkName)
{
  // Do nothing
}

//==============================================================================
std::unique_ptr<InverseKinematics::GradientMethod> LegIK::clone(
    InverseKinematics* _newIK) const
{
  return std::unique_ptr<GradientMethod>(
        new LegIK(_newIK, mBaseLinkName, getAnalyticalProperties()));
}

//==============================================================================
const std::vector<LegIK::Solution>& LegIK::computeSolutions(
    const Eigen::Isometry3d& _desiredBodyTf)
{
  mSolutions.clear();
  mSolutions.reserve(8);

  if(!configured)
  {
    configure();

    if(!configured)
    {
      dtwarn << "[HuboLegIK::computeSolutions] This analytical IK was not able "
            << "to configure properly, so it will not be able to compute "
            << "solutions\n";
      return mSolutions;
    }
  }

  const BodyNodePtr& base = mBaseLink.lock();
  if(nullptr == base)
  {
    dterr << "[HuboLegIK::computeSolutions] Attempting to perform IK on a "
          << "limb that no longer exists!\n";
    assert(false);
    return mSolutions;
  }

  double nx, ny, sx, sy, ax, ay, az, px, py, pz;
  double q1, q2, q3, q4, q5, q6;
  double S2, S4, S6;
  double C2, C4, C5, C6;
  double C45, psi, q345;
  std::complex<double> radical;
  std::complex<double> sqrt_radical;
  Eigen::Isometry3d B, Binv;

  Eigen::Vector6d testQ;

  B = (base->getParentBodyNode()->getWorldTransform() * waist).inverse()
      * _desiredBodyTf * footTfInv;
  Binv = B.inverse();

  nx = Binv(0,0); sx = Binv(0,1); ax = Binv(0,2); px = Binv(0,3);
  ny = Binv(1,0); sy = Binv(1,1); ay = Binv(1,2); py = Binv(1,3);
                                  az = Binv(2,2); pz = Binv(2,3);

  for(size_t i=0; i < 8; ++i)
  {
    bool isValid = true;

    C4 = ((px+L6)*(px+L6) - L4*L4 - L5*L5 + py*py + pz*pz)/(2*L4*L5);
    radical = 1-C4*C4;
    sqrt_radical = std::sqrt(radical);
    if(sqrt_radical.imag() != 0)
        isValid = false;
    q4 = atan2(alternatives(i,0)*sqrt_radical.real(), C4);

    S4 = sin(q4);
    psi = atan2(S4*L4, C4*L4+L5);
    radical = (px+L6)*(px+L6) + py*py;
    sqrt_radical = std::sqrt(radical);
    if(sqrt_radical.imag() != 0)
        isValid = false;

    q5 = dart::math::wrapToPi(atan2(-pz, alternatives(i,1)*sqrt_radical.real())-psi);

    q6 = atan2(py, -(px+L6));
    C45 = cos(q4+q5);
    C5 = cos(q5);
    if( C45*L4 + C5*L5 < 0 )
        q6 = dart::math::wrapToPi(q6+M_PI);

    S6 = sin(q6);
    C6 = cos(q6);

    S2 = C6*ay + S6*ax;
    radical = 1-S2*S2;
    sqrt_radical = std::sqrt(radical);
    if(sqrt_radical.imag() != 0)
        isValid = false;
    q2 = atan2(S2, alternatives(i,2)*sqrt_radical.real());

    q1 = atan2(C6*sy + S6*sx, C6*ny + S6*nx);
    C2 = cos(q2);
    if( C2 < 0 )
        q1 = dart::math::wrapToPi(q1+M_PI);

    q345 = atan2(-az/C2, -(C6*ax - S6*ay)/C2);
    q3 = dart::math::wrapToPi(q345 - q4 - q5);

    testQ[0]=q1; testQ[1]=q2; testQ[2]=q3; testQ[3]=q4; testQ[4]=q5; testQ[5]=q6;

    for(int k=0; k<testQ.size(); ++k)
        if( fabs(testQ[k]) < zeroSize )
            testQ[k] = 0;

    int validity = isValid? VALID : OUT_OF_REACH;
    mSolutions.push_back(Solution(testQ, validity));
  }

  checkSolutionJointLimits();

  return mSolutions;
}

//==============================================================================
const std::vector<size_t>& LegIK::getDofs() const
{
  if(!configured)
    configure();

  return mDofs;
}

//==============================================================================
void LegIK::configure() const
{
  configured = false;

  mBaseLink = mIK->getNode()->getSkeleton()->getBodyNode(mBaseLinkName);

  BodyNode* base = mBaseLink.lock();
  if(nullptr == base)
  {
    dterr << "[HuboLegIK::configure] base link is a nullptr\n";
    assert(false);
    return;
  }

  const SkeletonPtr& skel = mIK->getNode()->getSkeleton();
  BodyNode* pelvis = skel->getBodyNode("Body_TSY");
  if(nullptr == pelvis)
  {
    dterr << "[HuboLegIK::configure] Could not find Hubo's pelvis "
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
      dterr << "[HuboLegIK::configure] Invalid number of DOFs ("
            << joint->getNumDofs() << ") in the Joint [" << joint->getName()
            << "]\n";
      assert(false);
      return;
    }

    dofs[i] = joint->getDof(0);
    saved_q[i] = dofs[i]->getPosition();
    dofs[i]->setPosition(0.0);

    if(bn->getNumChildBodyNodes() > 0)
      bn = bn->getChildBodyNode(0);
  }

  L4 = std::abs(dofs[3]->getChildBodyNode()->
      getRelativeTransform().translation()[2]);

  L5 = std::abs(dofs[4]->getChildBodyNode()->
      getRelativeTransform().translation()[2]);

  // This offset will be taken care of with footTfInv
  L6 = 0.0;

  hipRotation = Eigen::Isometry3d::Identity();
  hipRotation.rotate(Eigen::AngleAxisd(90*M_PI/180.0,
                                       Eigen::Vector3d::UnitZ()));

  waist = dofs[2]->getChildBodyNode()->getTransform(
            dofs[0]->getParentBodyNode()) * hipRotation;

  footTfInv = Eigen::Isometry3d::Identity();
  footTfInv.rotate(Eigen::AngleAxisd(-90*M_PI/180.0,
                                     Eigen::Vector3d::UnitY()));
  footTfInv = footTfInv * mIK->getNode()->getTransform(dofs[5]->getChildBodyNode());
  footTfInv = footTfInv.inverse();

  alternatives <<
       1,  1,  1,
       1,  1, -1,
       1, -1,  1,
       1, -1, -1,
      -1,  1,  1,
      -1,  1, -1,
      -1, -1,  1,
      -1, -1, -1;

  for(size_t i=0; i < 6; ++i)
  {
    dofs[i]->setPosition(saved_q[i]);
    mDofs.push_back(dofs[i]->getIndexInSkeleton());
  }

  configured = true;
}



} // namespace hubo
