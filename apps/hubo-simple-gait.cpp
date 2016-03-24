/*
 * Copyright (c) 2016, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Michael X. Grey <mxgrey@gatech.edu>
 *
 * Georgia Tech Graphics Lab and AMBER Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Aaron D. Ames
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

#include <vector>
#include <Eigen/Core>
#include <dart/dynamics/Skeleton.h>
#include <dart/dynamics/EndEffector.h>
#include <dart/dynamics/DegreeOfFreedom.h>
#include <dart/simulation/World.h>
#include <dart/constraint/BalanceConstraint.h>
#include <osgDart/osgDart.h>
#include <hubo/DrcModel.hpp>
#include <hubo/RelaxedPosture.hpp>

#include <HuboPath/Operator.hpp>

using Trajectory = std::vector<Eigen::VectorXd>;
using namespace dart::dynamics;
using namespace dart::simulation;

const double DefaultStepLength = 0.13;

const double DefaultSwayOffset_X = -0.14;
const double DefaultSwayOffset_Y = 0.0;

const size_t DefaultNumSteps = 4;


const double frequency = 200.0;
//const double frequency = 10.0;

//==============================================================================
std::vector<Eigen::VectorXd> shiftOver(
    const dart::dynamics::SkeletonPtr& robot,
    const Eigen::Vector2d& swayOffset,
    bool leftStep,
    bool finish)
{
  std::vector<Eigen::VectorXd> walk;
  walk.push_back(robot->getPositions());

  // Shift to the stance side
  dart::dynamics::EndEffector* stance = leftStep?
        robot->getEndEffector("r_foot") : robot->getEndEffector("l_foot");
  dart::dynamics::EndEffector* swing = leftStep?
        robot->getEndEffector("l_foot") : robot->getEndEffector("r_foot");

  Eigen::VectorXd start_q = robot->getPositions();
  Eigen::VectorXd end_q = start_q;

  if(finish)
  {
    end_q.block<2,1>(3,0) =
        (stance->getWorldTransform().translation().block<2,1>(0,0)
        + swing->getWorldTransform().translation().block<2,1>(0,0)) / 2.0;
  }
  else
  {
    end_q.block<2,1>(3,0) =
        stance->getWorldTransform().translation().block<2,1>(0,0);
    end_q[3] += swayOffset[0];
    end_q[4] += leftStep? -swayOffset[1] : swayOffset[1];
  }

  stance->getIK(true)->getTarget()->setTransform(
        stance->getWorldTransform());

  swing->getIK(true)->getTarget()->setTransform(
        swing->getWorldTransform());

  const Eigen::Vector6d x0 = start_q.block<6,1>(0,0);
  const Eigen::Vector6d xf = end_q.block<6,1>(0,0);
  const Eigen::Vector6d dx = xf - x0;

  stance->getSupport()->setActive(true);
  swing->getSupport()->setActive(true);

  const size_t res = 10;
  for(size_t i=0; i < res+1; ++i)
  {
    const double t = (double)(i)/(double)(res);

    const Eigen::Vector6d x = x0 + t*dx;
    std::cout << "x,y init:  " << x.block<2,1>(3,0).transpose() << std::endl;
    robot->getJoint(0)->setPositions(x);


    stance->getIK(true)->solve();
    swing->getIK(true)->solve();

    std::cout << "x,y final: " << robot->getJoint(0)->getPositions().block<2,1>(3,0).transpose()
        << "\n" << std::endl;

//    if(!robot->getIK(true)->solve())
//    {
//      walk.push_back(robot->getPositions());
//      std::cout <<" ===== COULD NOT SOLVE FOOT STEP AT INDEX " << walk.size()-1 << std::endl;
//      return walk;
//    }

    walk.push_back(robot->getPositions());
  }

  return walk;
}

//==============================================================================
std::vector<Eigen::VectorXd> stepForward(
    const dart::dynamics::SkeletonPtr& robot,
    double stepDist,
    bool leftStep,
    bool finish)
{
  std::vector<Eigen::VectorXd> walk;
  walk.push_back(robot->getPositions());

  std::shared_ptr<hubo::RelaxedPosture> posture =
      std::dynamic_pointer_cast<hubo::RelaxedPosture>(
        robot->getIK(true)->getObjective());

  dart::dynamics::EndEffector* stance = leftStep?
        robot->getEndEffector("r_foot") : robot->getEndEffector("l_foot");
  dart::dynamics::EndEffector* swing = leftStep?
        robot->getEndEffector("l_foot") : robot->getEndEffector("r_foot");

  stance->getSupport()->setActive(true);
  swing->getSupport()->setActive(false);

  const double tol = 1e-4;
  dart::math::SupportGeometry geom = stance->getSupport()->getGeometry();
  Eigen::Vector3d com = robot->getCOM(stance);
  com[2] = 0.0;
  geom[0] = com + Eigen::Vector3d( tol,  tol, 0);
  geom[1] = com + Eigen::Vector3d(-tol,  tol, 0);
  geom[2] = com + Eigen::Vector3d(-tol, -tol, 0);
  geom[3] = com + Eigen::Vector3d( tol, -tol, 0);

  stance->getSupport()->setGeometry(geom);
  stance->getIK(true)->getTarget()->setTransform(
        stance->getWorldTransform());
  const Eigen::Isometry3d stanceTf = stance->getWorldTransform();

  const Eigen::Isometry3d swingStart = swing->getWorldTransform();
  Eigen::Isometry3d swingGoal = swingStart;
  swingGoal.translation()[0] = stanceTf.translation()[0];
  if(!finish)
    swingGoal.translation()[0] += stepDist;

  const Eigen::Vector2d v =
      (swingGoal.translation()-swingStart.translation()).block<2,1>(0,0)/(2*M_PI);
  const double D = v.norm();

  Eigen::AngleAxisd aa(swingStart.linear().transpose()*swingGoal.linear());
  const double R = aa.angle()/(2*M_PI);
  const Eigen::Vector3d axis = aa.axis();

  const size_t res = 10;
  for(size_t i=0; i < res+1; ++i)
  {
    const double t = (double)(i)/(double)(res)*2*M_PI;
    const Eigen::Vector2d x = v*(t - sin(t));
    const double z = D*(1 - cos(t));
    const double theta = t*R;

    Eigen::Isometry3d swingTarget = swingStart;
    swingTarget.translation().block<2,1>(0,0) += x;
    swingTarget.translation()[2] += z;
    swingTarget.rotate(Eigen::AngleAxisd(theta, axis));

    swing->getIK(true)->getTarget()->setTransform(swingTarget);

    if(posture)
      posture->enforceIdealPosture = true;

    robot->getIK(true)->solve();

    if(posture)
      posture->enforceIdealPosture = false;

    if(!robot->getIK(true)->solve())
    {
      walk.push_back(robot->getPositions());
      std::cout <<" ===== COULD NOT SOLVE FOOT STEP AT INDEX " << walk.size()-1 << std::endl;
      return walk;
    }

    walk.push_back(robot->getPositions());
  }

  return walk;
}

//==============================================================================
std::vector<std::vector<Eigen::VectorXd>> interpolateSteps(
    const dart::dynamics::SkeletonPtr& robot,
    double stepDist, const Eigen::Vector2d& swayOffset,
    size_t numSteps)
{
  std::shared_ptr<dart::constraint::BalanceConstraint> balance;
  const std::shared_ptr<dart::optimizer::Problem>& problem =
      robot->getIK(true)->getProblem();
  for(size_t i=0; i < problem->getNumEqConstraints(); ++i)
  {
    balance = std::dynamic_pointer_cast<dart::constraint::BalanceConstraint>(
          problem->getEqConstraint(i));
    if(balance)
      break;
  }

  std::vector<std::vector<Eigen::VectorXd>> walk;
  for(size_t i=0; i < numSteps; ++i)
  {
    std::cout << "Generating Step #" << i << " ("
              << (double)(i)/(double)(numSteps) << "%)" << std::endl;
    problem->removeEqConstraint(balance);
    std::cout << "Shift" << std::endl;
    walk.push_back(shiftOver(robot, swayOffset, i%2==0, false));

    problem->addEqConstraint(balance);
    std::cout << "Forward" << std::endl;
    walk.push_back(stepForward(robot, stepDist, i%2==0, i==numSteps-1));
  }

  problem->removeEqConstraint(balance);
  walk.push_back(shiftOver(robot, swayOffset, numSteps%2==0, true));

  return walk;
}

//==============================================================================
class TrajectoryDisplayWorld : public osgDart::WorldNode
{
public:

  TrajectoryDisplayWorld(dart::simulation::WorldPtr world,
                         const std::vector<Eigen::VectorXd>& traj)
    : osgDart::WorldNode(world), mTrajectory(traj)
  {
    hubo = world->getSkeleton(0);
    LSR = hubo->getDof("LSR")->getIndexInSkeleton();
    RSR = hubo->getDof("RSR")->getIndexInSkeleton();

    firstLoop = true;
  }

  void customPreRefresh() override
  {
    if(mTrajectory.size() == 0)
    {
      if(firstLoop)
      {
        std::cerr << "No trajectory was generated!" << std::endl;
        firstLoop = false;
      }
      return;
    }

    if(firstLoop)
    {
      mTimer.setStartTick();
      firstLoop = false;
    }

    double time = mTimer.time_s();
    size_t count = (size_t)(floor(frequency*time)) % mTrajectory.size();

    double dt = 1.0/frequency;

    size_t last_1 = count <= 0? 0 : count - 1;
    size_t last_2 = count <= 1? 0 : count - 2;
    size_t last_3 = count <= 2? 0 : count - 3;
    size_t last_4 = count <= 3? 0 : count - 4;

    size_t next_1 = count < mTrajectory.size() - 1? count + 1 : count;
    size_t next_2 = count < mTrajectory.size() - 2? count + 2 : count;
    size_t next_3 = count < mTrajectory.size() - 3? count + 3 : count;
    size_t next_4 = count < mTrajectory.size() - 4? count + 4 : count;

    Eigen::VectorXd velocities =
        (0.5*mTrajectory[next_1] - 0.5*mTrajectory[last_1])/dt;

    Eigen::Matrix<double, 5, 1> c;
    c << -205.0/72.0, 8.0/5.0, -1.0/5.0, 8.0/315.0, -1.0/560.0;

    Eigen::VectorXd accelerations =
        (  c[1]*mTrajectory[next_1] + c[2]*mTrajectory[next_2] + c[3]*mTrajectory[next_3] + c[4]*mTrajectory[next_4]
         + c[0]*mTrajectory[count]
         + c[1]*mTrajectory[last_1] + c[2]*mTrajectory[last_2] + c[3]*mTrajectory[last_3] + c[4]*mTrajectory[last_4])/pow(dt,2);


//    accelerations.head<6>().setZero();

    hubo->setPositions(mTrajectory[count]);
    hubo->setVelocities(velocities);
    hubo->setAccelerations(accelerations);

//    hubo->setPositions(mMapping, positions);
//    clone->setPositions(mMapping, mRaw[count]);
  }

protected:

  SkeletonPtr hubo;
  std::vector<Eigen::VectorXd> mTrajectory;
  size_t LSR;
  size_t RSR;

  bool firstLoop;
  osg::Timer mTimer;
};

//==============================================================================
int main()
{
  SkeletonPtr robot = hubo::DrcModel::create();
  robot->getIK(true)->getSolver()->setNumMaxIterations(500);


  Eigen::Vector6d bounds(Eigen::Vector6d::Constant(
                          dart::dynamics::DefaultIKTolerance));
  robot->getEndEffector("l_foot")->getIK(true)->getErrorMethod().setBounds(-bounds, bounds);
  robot->getEndEffector("r_foot")->getIK(true)->getErrorMethod().setBounds(-bounds, bounds);

  std::vector<std::vector<Eigen::VectorXd>> walk = interpolateSteps(
        robot, DefaultStepLength,
        Eigen::Vector2d(DefaultSwayOffset_X, DefaultSwayOffset_Y),
        DefaultNumSteps);

  std::cout << "waypoints: " << walk.size() << std::endl;

  HuboPath::Operator mOperator;
  std::vector<size_t> opIndices;
  std::vector<std::string> indexNames;
  for(size_t i=6; i < robot->getNumDofs(); ++i)
  {
    DegreeOfFreedom* dof = robot->getDof(i);
    opIndices.push_back(i);
    indexNames.push_back(dof->getName());
  }
  mOperator.setJointIndices(indexNames);

  HuboPath::Trajectory trajectory = mOperator.getCurrentTrajectory();
  trajectory.elements.clear();

  for(size_t s=0; s < walk.size(); ++s)
  {
    std::cout << "interpolating walk component #" << s << std::endl;
    const std::vector<Eigen::VectorXd>& step = walk[s];
    mOperator.clearWaypoints();
    for(size_t i=0; i < step.size(); ++i)
    {
      robot->setPositions(step[i]);
      mOperator.addWaypoint(robot->getPositions(opIndices));
    }

    mOperator.setInterpolationMode(HUBO_PATH_OPTIMAL);
    HuboPath::Trajectory step_trajectory = mOperator.getCurrentTrajectory();
    step_trajectory.elements.erase(step_trajectory.elements.begin());

    HuboCan::HuboDescription& desc = step_trajectory.desc;
    for(size_t i=0; i < desc.getJointCount(); ++i)
    {
      hubo_joint_info_t& info = desc.joints[i]->info;
      hubo_joint_limits_t& limits = info.limits;

      limits.nominal_speed /= 5.0;
      limits.nominal_accel /= 5.0;
    }
    step_trajectory.interpolate();

    for(size_t j=0; j < step_trajectory.size(); ++j)
      trajectory.elements.push_back(step_trajectory.elements[j]);
  }







//  for(size_t s=0; s < walk.size(); ++s)
//  {
//    const std::vector<Eigen::VectorXd>& step = walk[s];
//    for(size_t i=0; i < step.size(); ++i)
//    {
//      robot->setPositions(step[i]);
//      mOperator.addWaypoint(robot->getPositions(opIndices));
//    }
//  }
//  HuboPath::Trajectory trajectory = mOperator.getCurrentTrajectory();
//  trajectory.elements.erase(trajectory.elements.begin());







  bool operate = true;
  operate = false;

  if(operate)
  {
    mOperator.sendNewTrajectory(trajectory);
  }
  else
  {

    std::vector<Eigen::VectorXd> fullwalk;

    std::cout << "Converting back to trajectory (" << trajectory.size() << ")..."
              << std::endl;
    for(size_t i=0; i < trajectory.size(); ++i)
    {
      Eigen::VectorXd q(robot->getNumDofs());
      for(size_t j=0; j < 6; ++j)
        q[j] = 0.0;

      const IndexArray& indexMap = mOperator.getIndexMap();
      for(size_t j=6; j < robot->getNumDofs(); ++j)
        q[j] = trajectory[i].references[indexMap[j-6]];

      fullwalk.push_back(q);
    }
    std::cout << "... converted!" << std::endl;

    auto world = std::make_shared<World>();
    world->addSkeleton(robot);

    osg::ref_ptr<TrajectoryDisplayWorld> display =
        new TrajectoryDisplayWorld(world, fullwalk);

    osgDart::Viewer viewer;
    viewer.addWorldNode(display);
    viewer.allowSimulation(false);

    robot->getEndEffector("l_foot")->getSupport()->setActive();
    robot->getEndEffector("r_foot")->getSupport()->setActive();
    viewer.addAttachment(new osgDart::SupportPolygonVisual(robot, -0.97+0.02));

    viewer.setUpViewInWindow(0, 0, 1280, 960);

    // Set up the default viewing position
    viewer.getCameraManipulator()->setHomePosition(osg::Vec3( 5.34,  3.00, 1.00),
                                                   osg::Vec3( 0.00,  0.00, 0.00),
                                                   osg::Vec3(-0.20, -0.08, 0.98));

    // Reset the camera manipulator so that it starts in the new viewing position
    viewer.setCameraManipulator(viewer.getCameraManipulator());
    std::cout << "Launching viewer" << std::endl;
    viewer.run();
  }
}























