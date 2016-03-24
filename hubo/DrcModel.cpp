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

#include <dart/dynamics/DegreeOfFreedom.h>
#include <dart/constraint/BalanceConstraint.h>
#include <dart/utils/urdf/DartLoader.h>
#include <dart/collision/CollisionDetector.h>
#include <dart/constraint/ConstraintSolver.h>

#include "hubo/DrcModel.hpp"
#include "hubo/ArmIK.hpp"
#include "hubo/LegIK.hpp"
#include "hubo/RelaxedPosture.hpp"

#include "hubo/utils.hpp"

using namespace dart;
using namespace dynamics;

namespace hubo {
namespace DrcModel {

//==============================================================================
static void setStartupConfiguration(const SkeletonPtr& hubo)
{
  hubo->getDof("LHP")->setPosition(-45.0*M_PI/180.0);
  hubo->getDof("LKP")->setPosition( 90.0*M_PI/180.0);
  hubo->getDof("LAP")->setPosition(-45.0*M_PI/180.0);

  hubo->getDof("RHP")->setPosition(-45.0*M_PI/180.0);
  hubo->getDof("RKP")->setPosition( 90.0*M_PI/180.0);
  hubo->getDof("RAP")->setPosition(-45.0*M_PI/180.0);

  hubo->getDof("LSP")->setPosition(  30.0*M_PI/180.0);
//  hubo->getDof("LSP")->setPosition(  90.0*M_PI/180.0);
  hubo->getDof("LEP")->setPosition(-120.0*M_PI/180.0);
//  hubo->getDof("LEP")->setPosition( -60*M_PI/180.0);

  hubo->getDof("RSP")->setPosition(  30.0*M_PI/180.0);
//  hubo->getDof("RSP")->setPosition(  90.0*M_PI/180.0);
  hubo->getDof("REP")->setPosition(-120.0*M_PI/180.0);
//  hubo->getDof("REP")->setPosition( -60*M_PI/180.0);


  hubo->getDof("LSY")->setPositionLowerLimit(-90.0*M_PI/180.0);
  hubo->getDof("LSY")->setPositionUpperLimit( 90.0*M_PI/180.0);
  hubo->getDof("LWY")->setPositionLowerLimit(-90.0*M_PI/180.0);
  hubo->getDof("LWY")->setPositionUpperLimit( 90.0*M_PI/180.0);

  hubo->getDof("RSY")->setPositionLowerLimit(-90.0*M_PI/180.0);
  hubo->getDof("RSY")->setPositionUpperLimit( 90.0*M_PI/180.0);
  hubo->getDof("RWY")->setPositionLowerLimit(-90.0*M_PI/180.0);
  hubo->getDof("RWY")->setPositionUpperLimit( 90.0*M_PI/180.0);
}

//==============================================================================
static void setupEndEffectors(const SkeletonPtr& hubo, bool lockedFeet)
{
  Eigen::VectorXd rootjoint_weights = Eigen::VectorXd::Ones(7);
  rootjoint_weights = 0.01*rootjoint_weights;

  double extra_error_clamp = 0.1;

  Eigen::Vector3d linearBounds =
      Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());

  Eigen::Vector3d angularBounds =
      Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());

  Eigen::Isometry3d tf_hand(Eigen::Isometry3d::Identity());
  tf_hand.translate(Eigen::Vector3d(0.0, 0.0, -0.09));

  EndEffector* l_hand = hubo->getBodyNode("Body_LWR")->
      createEndEffector("l_hand");
  l_hand->setDefaultRelativeTransform(tf_hand, true);

  l_hand->getIK(true)->useWholeBody();

  l_hand->getIK()->setGradientMethod<ArmIK>("Body_LSP");

  l_hand->getIK()->getAnalytical()->setExtraDofUtilization(
        IK::Analytical::PRE_ANALYTICAL);

  l_hand->getIK()->getAnalytical()->setExtraErrorLengthClamp(extra_error_clamp);

  l_hand->getIK()->getGradientMethod().setComponentWeights(rootjoint_weights);

  l_hand->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds, linearBounds);

  l_hand->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);


  EndEffector* r_hand = hubo->getBodyNode("Body_RWR")->
      createEndEffector("r_hand");
  r_hand->setDefaultRelativeTransform(tf_hand, true);

  r_hand->getIK(true)->useWholeBody();

  r_hand->getIK()->setGradientMethod<ArmIK>("Body_RSP");

  r_hand->getIK()->getAnalytical()->setExtraDofUtilization(
        IK::Analytical::PRE_ANALYTICAL);

  r_hand->getIK()->getAnalytical()->setExtraErrorLengthClamp(extra_error_clamp);

  r_hand->getIK()->getGradientMethod().setComponentWeights(rootjoint_weights);

  r_hand->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds, linearBounds);

  r_hand->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);


  double y_range = 0.01;
//  double y_range = 0.005;
//  double left_y_offset = -0.02;
  double left_y_offset = -0.01;
  dart::math::SupportGeometry left_foot_support;
  left_foot_support.push_back(Eigen::Vector3d(-0.11,  y_range+left_y_offset, 0.0));
  left_foot_support.push_back(Eigen::Vector3d(-0.15,  y_range+left_y_offset, 0.0));
  left_foot_support.push_back(Eigen::Vector3d(-0.15, -y_range+left_y_offset, 0.0));
  left_foot_support.push_back(Eigen::Vector3d(-0.11, -y_range+left_y_offset, 0.0));

//  double right_y_offset = -left_y_offset + 0.005;
  double right_y_offset = -left_y_offset + 0.01;
  dart::math::SupportGeometry right_foot_support = left_foot_support;
  for(size_t i=0; i < right_foot_support.size(); ++i)
    right_foot_support[i][1] += right_y_offset - left_y_offset;


  Eigen::Isometry3d tf_foot(Eigen::Isometry3d::Identity());
  double ground_dist = 0.00;
  tf_foot.translation() = Eigen::Vector3d(0.14, 0.0, -0.136+ground_dist);

  if(lockedFeet)
  {
    linearBounds = 1e-8*Eigen::Vector3d::Ones();
    angularBounds = 1e-8*Eigen::Vector3d::Ones();
  }
  else
  {
    linearBounds[2] = 1e-8;

    angularBounds[0] = 1e-8;
    angularBounds[1] = 1e-8;
  }

  Eigen::Vector3d ground_offset = ground_dist * Eigen::Vector3d::UnitZ();

  EndEffector* l_foot = hubo->getBodyNode("Body_LAR")->
      createEndEffector("l_foot");
  l_foot->setDefaultRelativeTransform(tf_foot, true);

  l_foot->createIK();

  l_foot->getIK()->setHierarchyLevel(1);

  l_foot->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds + ground_offset, linearBounds + ground_offset);
  l_foot->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);

  l_foot->getIK()->setGradientMethod<LegIK>("Body_LHY");

  l_foot->getSupport(true)->setGeometry(left_foot_support);
  l_foot->getSupport()->setActive();


  EndEffector* r_foot = hubo->getBodyNode("Body_RAR")->
      createEndEffector("r_foot");
  r_foot->setDefaultRelativeTransform(tf_foot, true);

  r_foot->createIK();

  r_foot->getIK()->setHierarchyLevel(1);

  r_foot->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds + ground_offset, linearBounds + ground_offset);
  r_foot->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);

  r_foot->getIK()->setGradientMethod<LegIK>("Body_RHY");

  r_foot->getSupport(true)->setGeometry(right_foot_support);
  r_foot->getSupport()->setActive();


  dart::math::SupportGeometry peg_support;
  peg_support.push_back(Eigen::Vector3d::Zero());

  linearBounds = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
  angularBounds = linearBounds;

  Eigen::Isometry3d tf_peg(Eigen::Isometry3d::Identity());
  tf_peg.translation() = Eigen::Vector3d(0.0, 0.0, 0.09);

  EndEffector* l_peg = hubo->getBodyNode("Body_LWP")->createEndEffector("l_peg");
  l_peg->setDefaultRelativeTransform(tf_peg, true);

  l_peg->getIK(true)->setGradientMethod<ArmIK>("Body_LSP");

  l_peg->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds, linearBounds);
  l_peg->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);

  l_peg->getSupport(true)->setGeometry(peg_support);

  EndEffector* r_peg = hubo->getBodyNode("Body_RWP")->createEndEffector("r_peg");
  r_peg->setDefaultRelativeTransform(tf_peg, true);

  r_peg->getIK(true)->setGradientMethod<ArmIK>("Body_RSP");

  r_peg->getIK()->getErrorMethod().setLinearBounds(
        -linearBounds, linearBounds);
  r_peg->getIK()->getErrorMethod().setAngularBounds(
        -angularBounds, angularBounds);

  r_peg->getSupport(true)->setGeometry(peg_support);

  double heightChange = -r_foot->getWorldTransform().translation()[2]+ground_dist;
  hubo->getDof("rootJoint_pos_z")->setPosition(heightChange);

  l_foot->getIK()->getTarget()->setTransform(l_foot->getTransform());
  r_foot->getIK()->getTarget()->setTransform(r_foot->getTransform());
}

//==============================================================================
#define SET_IDEAL_POSTURE_LIMITS( X )\
{ size_t index = hubo->getDof( #X )->getIndexInSkeleton();\
  lower_ideal[index] = ideal[index] - tol;\
  upper_ideal[index] = ideal[index] + tol;\
}

//==============================================================================
static void setupWholeBodySolver(const SkeletonPtr& hubo)
{
  std::shared_ptr<dart::optimizer::GradientDescentSolver> solver =
      std::dynamic_pointer_cast<dart::optimizer::GradientDescentSolver>(
        hubo->getIK(true)->getSolver());

  size_t nDofs = hubo->getNumDofs();

  const double tol = 1e-2;
  const double tolFactor = 0.05;

  Eigen::VectorXd ideal = hubo->getPositions();

  double default_weight = 0.01;
  Eigen::VectorXd weights = default_weight * Eigen::VectorXd::Ones(nDofs);
  weights[2] = 0.0;
  weights[3] = 0.0;
  weights[4] = 0.0;

  Eigen::VectorXd lower_posture = Eigen::VectorXd::Constant(nDofs,
        -std::numeric_limits<double>::infinity());
  lower_posture[0] = -0.35;
  lower_posture[1] = -0.35;
  lower_posture[5] =  0.55;
  Eigen::VectorXd lower_ideal = tolFactor*(lower_posture - ideal) + ideal;

  Eigen::VectorXd upper_posture = Eigen::VectorXd::Constant(nDofs,
         std::numeric_limits<double>::infinity());
  upper_posture[0] =  0.35;
  upper_posture[1] =  0.50;
  upper_posture[5] =  0.95;
  Eigen::VectorXd upper_ideal = tolFactor*(upper_posture - ideal) + ideal;

  SET_IDEAL_POSTURE_LIMITS( TSY );

  SET_IDEAL_POSTURE_LIMITS( RSP );
  SET_IDEAL_POSTURE_LIMITS( RSR );
  SET_IDEAL_POSTURE_LIMITS( RSY );
  SET_IDEAL_POSTURE_LIMITS( REP );
  SET_IDEAL_POSTURE_LIMITS( RWY );
  SET_IDEAL_POSTURE_LIMITS( RWP );
  SET_IDEAL_POSTURE_LIMITS( RWR );

  SET_IDEAL_POSTURE_LIMITS( LSP );
  SET_IDEAL_POSTURE_LIMITS( LSR );
  SET_IDEAL_POSTURE_LIMITS( LSY );
  SET_IDEAL_POSTURE_LIMITS( LEP );
  SET_IDEAL_POSTURE_LIMITS( LWY );
  SET_IDEAL_POSTURE_LIMITS( LWP );
  SET_IDEAL_POSTURE_LIMITS( LWR );

  std::shared_ptr<RelaxedPosture> objective = std::make_shared<RelaxedPosture>(
        ideal, lower_ideal, upper_ideal, lower_posture, upper_posture, weights);

  hubo->getIK()->setObjective(objective);

  std::shared_ptr<dart::constraint::BalanceConstraint> balance =
      std::make_shared<dart::constraint::BalanceConstraint>(hubo->getIK());
  hubo->getIK()->getProblem()->addEqConstraint(balance);

  solver->setNumMaxIterations(1000);

  objective->enforceIdealPosture = true;
  balance->setBalanceMethod(dart::constraint::BalanceConstraint::SHIFT_SUPPORT);
  balance->setErrorMethod(dart::constraint::BalanceConstraint::OPTIMIZE_BALANCE);

  hubo->getIK()->solve();
  solver->setNumMaxIterations(200);

  objective->enforceIdealPosture = false;
  balance->setErrorMethod(dart::constraint::BalanceConstraint::FROM_CENTROID);
  balance->setBalanceMethod(dart::constraint::BalanceConstraint::SHIFT_SUPPORT);
//  balance->setBalanceMethod(dart::constraint::BalanceConstraint::SHIFT_COM);

//  hubo->getEndEffector("l_foot")->getIK()->setHierarchyLevel(0);
//  hubo->getEndEffector("r_foot")->getIK()->setHierarchyLevel(0);

  hubo->enableSelfCollision();
}

//==============================================================================
void removeFingers(const SkeletonPtr& hubo)
{
  for(size_t i=0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    if(bn->getName().substr(0, 7) == "Body_LF"
       || bn->getName().substr(0, 7) == "Body_RF")
    {
//      bn->removeAllCollisionShapes();
      bn->remove();
      --i;
    }
  }
}

//==============================================================================
void removeHead(const SkeletonPtr& hubo)
{
  for(size_t i=0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    if(bn->getName().substr(0, 7) == "Body_NK")
    {
      bn->remove();
      --i;
    }
  }
}

//==============================================================================
SkeletonPtr create(const std::string& modelFile, const std::string& modelPath,
                   bool lockedFeet)
{
  dart::utils::DartLoader loader;
  loader.addPackageDirectory("drchubo", modelPath);
  SkeletonPtr hubo =
      loader.parseSkeleton(modelPath + modelFile);

  removeFingers(hubo);
  removeHead(hubo);

  setStartupConfiguration(hubo);
  setupEndEffectors(hubo, lockedFeet);
  setupWholeBodySolver(hubo);

  for(size_t i=0; i<2; ++i)
    hubo->getDof(i)->setPositionLimits(-M_PI/4.0, M_PI/4.0);
  hubo->getDof(2)->setPositionLimits(-M_PI, M_PI);

  const double z = hubo->getDof(5)->getPosition();
  for(size_t i=3; i < 5; ++i)
    hubo->getDof(i)->setPositionLimits(-2.5, 2.5);
  hubo->getDof(5)->setPositionLimits(z-0.5, z+0.5);
  hubo->getDof(5)->setRestPosition(1.0);

  return hubo;
}

//==============================================================================
void goodPosture(const SkeletonPtr& hubo)
{
  std::shared_ptr<RelaxedPosture> posture =
      std::dynamic_pointer_cast<RelaxedPosture>(
        hubo->getIK(true)->getObjective());

  if(posture)
    posture->enforceIdealPosture = true;

  const std::shared_ptr<dart::optimizer::Problem>& problem =
      hubo->getIK(true)->getProblem();
  std::shared_ptr<dart::constraint::BalanceConstraint> balance;
  for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
  {
    balance = std::dynamic_pointer_cast<dart::constraint::BalanceConstraint>(
          problem->getEqConstraint(i));
    if(balance)
      break;
  }

  if(balance)
    balance->setErrorMethod(dart::constraint::BalanceConstraint::OPTIMIZE_BALANCE);
}

//==============================================================================
void laxPosture(const SkeletonPtr& hubo)
{
  std::shared_ptr<RelaxedPosture> posture =
      std::dynamic_pointer_cast<RelaxedPosture>(
        hubo->getIK(true)->getObjective());

  if(posture)
    posture->enforceIdealPosture = false;

  const std::shared_ptr<dart::optimizer::Problem>& problem =
      hubo->getIK(true)->getProblem();
  std::shared_ptr<dart::constraint::BalanceConstraint> balance;
  for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
  {
    balance = std::dynamic_pointer_cast<dart::constraint::BalanceConstraint>(
          problem->getEqConstraint(i));
    if(balance)
      break;
  }

  if(balance)
    balance->setErrorMethod(dart::constraint::BalanceConstraint::FROM_CENTROID);
}

//==============================================================================
#define DISABLE_PAIR( X1, X2 ) cd->disablePair(hubo->getBodyNode( #X1 ), hubo->getBodyNode( #X2 ))

//==============================================================================
void disableAdjacentPairs(const dart::dynamics::SkeletonPtr& hubo,
                          const dart::simulation::WorldPtr& world)
{
  dart::collision::CollisionDetector* cd =
      world->getConstraintSolver()->getCollisionDetector();

  DISABLE_PAIR(Body_LHY,Body_LHP);
  DISABLE_PAIR(Body_LKP,Body_LAR);
  DISABLE_PAIR(Body_RHY,Body_RHP);
  DISABLE_PAIR(Body_RKP,Body_RAR);
  DISABLE_PAIR(Body_LSP,Body_LSY);
  DISABLE_PAIR(Body_RSP,Body_RSY);
}

//==============================================================================
#define ADD_INDEX( X ) indices.push_back(hubo->getDof( #X )->getIndexInSkeleton());

//==============================================================================
std::vector<size_t> getUsefulDofIndices(const ConstSkeletonPtr& hubo)
{
  std::vector<size_t> indices;

  // Add the root
  for(size_t i=0; i<6; ++i)
    indices.push_back(i);

  // Add the important joints
  ADD_INDEX( TSY );

  ADD_INDEX( RHY );
  ADD_INDEX( RHR );
  ADD_INDEX( RHP );
  ADD_INDEX( RKP );
  ADD_INDEX( RAP );
  ADD_INDEX( RAR );

  ADD_INDEX( LHY );
  ADD_INDEX( LHR );
  ADD_INDEX( LHP );
  ADD_INDEX( LKP );
  ADD_INDEX( LAP );
  ADD_INDEX( LAR );

  ADD_INDEX( RSP );
  ADD_INDEX( RSR );
  ADD_INDEX( RSY );
  ADD_INDEX( REP );
  ADD_INDEX( RWY );
  ADD_INDEX( RWP );
  ADD_INDEX( RWR );

  ADD_INDEX( LSP );
  ADD_INDEX( LSR );
  ADD_INDEX( LSY );
  ADD_INDEX( LEP );
  ADD_INDEX( LWY );
  ADD_INDEX( LWP );
  ADD_INDEX( LWR );

  return indices;
}

//==============================================================================
static void lockHand(EndEffector* hand)
{
  const Eigen::Vector6d bounds(Eigen::Vector6d::Constant(1e-6));
  hand->getIK()->getErrorMethod().setBounds(-bounds, bounds);
}

//==============================================================================
static void unlockHand(EndEffector* hand)
{
  const Eigen::Vector6d bounds =
      Eigen::Vector6d::Constant(std::numeric_limits<double>::infinity());
  hand->getIK()->getErrorMethod().setBounds(-bounds, bounds);
}

//==============================================================================
Eigen::Isometry3d projectToXY(const Eigen::Isometry3d& tf)
{
  Eigen::Isometry3d result(Eigen::Isometry3d::Identity());
  Eigen::Vector3d v = tf.linear().block<3,1>(0,2);
  Eigen::Vector3d z(Eigen::Vector3d::UnitZ());
  Eigen::Vector3d vp = v - v.dot(z)*z;
  assert(vp.norm() > 1e-4);
  vp.normalize();

  Eigen::AngleAxisd R = findRotation(v, vp);
  result.pretranslate(tf.translation());
  result.rotate(R);
  result.rotate(tf.rotation());
  return result;
}

//==============================================================================
static Eigen::VectorXd flattenHand(
    const Eigen::VectorXd& q, dart::dynamics::EndEffector* hand)
{
  const dart::dynamics::SkeletonPtr& robot = hand->getSkeleton();
  robot->setPositions(q);
  Eigen::Isometry3d tf = hand->getWorldTransform();
  hand->getIK(true)->getTarget()->setTransform(projectToXY(tf));
  hand->getIK(true)->solve();
  return robot->getPositions();
}

//==============================================================================
// TODO: Replace this with a config file rather than hardcoding it
std::deque<Eigen::VectorXd> getUsefulSeeds(
    const dart::dynamics::SkeletonPtr& robot)
{
  std::deque<Eigen::VectorXd> seeds;
  Eigen::VectorXd seed(33);

  dart::dynamics::EndEffector* hand = robot->getEndEffector("l_hand");
  lockHand(hand);

  seed << 0.0506239, 0.220517, 0.0674325, -0.023945, 0.0148171, 0.760608,
      -0.0210546, -0.0489457, -1.08156, 1.549, -0.687142, -0.00523577,
      -0.0246972, -0.0569363, -1.10521, 1.59255, -0.707247, 0.00335567,
      -0.379456, -0.371874, 0.219561, 0.336045, -1.69848, 0.145026, 0.35897,
      -0.0541734, 1.39029, -0.129258, -0.00157321, -2.08574, -0.00111973,
      0.000392012, -3.40398e-05;
  seeds.push_back(flattenHand(seed, hand));


  seed << -0.00400777, 0.5, 0.0368712, -0.010865, -0.00126129, 0.886177,
      -0.0316179, 0.0179162, -1.06899, 0.830162, -0.261053, -0.00575012,
      -0.0317096, 0.0177484, -1.06561, 0.844042, -0.278313, -0.0055589, 0,
      -0.501077, 0.0810193, -0.0311712, -1.57298, -0.157108, 0.0222733,
      0.209454, 1.08349, -0.150689, -0.186461, -2.14418, -0.0953182, -0.493902,
      -0.150918,
  seeds.push_back(flattenHand(seed, hand));


  seed << -0.00441859, 0.49995, 0.0365283, -0.0101918, -0.00109831, 0.886141,
      -0.0312797, 0.0180589, -1.06851, 0.831575, -0.262887, -0.00556233,
      -0.0313703, 0.0178929, -1.0648, 0.844716, -0.279737, -0.00537319,
      -0.000342878, -1.41277, 0.875731, -1.56433, -1.71245, 0.102654, -0.545887,
      0.884793, 1.08349, -0.150689, -0.186461, -2.14418, -0.0953182, -0.493902,
      -0.150918;
  seeds.push_back(flattenHand(seed, hand));


  seed << -0.000425537, -0.0174812, -0.000117557, -0.0153074,
      9.48786e-05, 0.763267, 0.000245928, 0.00025312, -0.849937, 1.5792,
      -0.711781, 0.000175712, -0.000263399, -0.000987922, -0.849768, 1.57888,
      -0.711628, 0.00140765, -1.4335e-05, 0.526692, -9.07258e-05, 0.000259322,
      -2.09439, -3.3246e-05, -1.41967e-05, -9.63282e-07, 0.526692, 0.000120258,
      -0.000236058, -2.09439, 3.32741e-05, -1.41967e-05, 9.93573e-07;
  seeds.push_back(flattenHand(seed, hand));


  seed << 0.00148122, -0.00144614, 0.00132295, 0.000555056, 0.00115728,
      0.767155, -1.18855e-06, 0.00040551, -0.805603, 1.61041, -0.803356,
      -0.00188577, -1.18978e-06, 0.000406358, -0.806154, 1.61151, -0.803906,
      -0.00188662, 0.00132295, 0.0209502, 0.144105, 0.857714, -1.36985,
      0.282422, -0.128234, -0.197028, 1.663, -0.210131, -0.114765, -2.19723,
      -0.118105, -0.424089, -0.27985;
  seeds.push_back(flattenHand(seed, hand));


  seed << 0, 0, 0, 0, 0, 0.76662, 0, 0, -0.806609, 1.61321, -0.806602, 0, 0, 0,
      -0.806609, 1.61321, -0.806602, 0, 0.874904, 0.475794, 0.389718, 0.862845,
      -1.3429, 0.210355, -0.189995, 0.431563, 1.22701, -0.11455, -0.439105,
      -2.29864, -0.26014, -0.336643, -0.277946;
  seeds.push_back(flattenHand(seed, hand));

//  seed << 0.0928673, 0.135991, 0.0862818, -0.0153082, 9.49365e-05, 0.763268,
//      -0.0116168, -0.0295227, -0.877085, 1.53289, -0.788849, -0.0675422,
//      -0.0130428, -0.0304285, -0.914682, 1.60774, -0.826224, -0.0664553,
//      0.307651, 1.37671, -0.012133, 0.00428279, -2.09283, 0.000175332,
//      5.73918e-05, 1.01263e-05, 0.00716002, -0.273675, -0.498219, -2.03159,
//      -0.113259, 0.463986, -0.114782;
//  seeds.push_back(seed);

  unlockHand(hand);

  return seeds;
}

} // namespace DrcModel
} // namespace hubo
