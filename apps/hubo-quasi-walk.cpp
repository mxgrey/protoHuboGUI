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
#include <osgDart/osgDart.h>
#include <hubo/DrcModel.hpp>
#include <hubo/RelaxedPosture.hpp>

#include <HuboPath/Operator.hpp>

using Trajectory = std::vector<Eigen::VectorXd>;
using namespace dart::dynamics;
using namespace dart::simulation;

//const double DefaultStepLength = 0.13;
const double DefaultStepLength = 0.155;
const double frequency = 200.0;

//==============================================================================
struct StepContext
{
  StepContext(EndEffector* stance, EndEffector* swing,
              size_t index, double t0, double dist)
    : stance_foot(stance),
      swing_foot(swing),
      route_index(index),
      start_t(t0),
      distance(dist),
      attempts(0)
  {
    // Do nothing
  }


  EndEffector* stance_foot;
  EndEffector* swing_foot;
  size_t route_index;
  double start_t;

  Eigen::VectorXd start_q;
  Eigen::VectorXd end_q;

  double distance;
  size_t attempts;
};

//==============================================================================
std::tuple<size_t, double> getNextRouteLocation(
    const std::vector<Eigen::VectorXd>& route,
    double stepDist,
    size_t start_index,
    double start_t)
{
  assert(start_t < 1.0);

  double remaining_D = stepDist;
  while(remaining_D > 0)
  {
    if(start_index == route.size()-1)
    {
      return std::tuple<size_t, double>(start_index, 0);
    }

    Eigen::VectorXd q0 = route[start_index];
    Eigen::VectorXd q1 = route[start_index+1];

    Eigen::Vector2d p0 = q0.block<2,1>(3,0);
    Eigen::Vector2d p1 = q1.block<2,1>(3,0);

    double D = (p1-p0).norm();
    double t = start_t + remaining_D/D;
    if(t >= 1)
    {
      remaining_D = remaining_D - D*(1.0-start_t);
      ++start_index;
      start_t = 0.0;
    }
    else
    {
      start_t = t;
      remaining_D = 0.0;
    }
  }

  return std::tuple<size_t, double>(start_index, start_t);
}

//==============================================================================
Eigen::VectorXd getStepConfig(const std::vector<Eigen::VectorXd>& route,
                              const StepContext& step)
{
  Eigen::VectorXd q = route[step.route_index];
  if(step.route_index < route.size()-1 && step.start_t > 0.0)
    q += step.start_t*(route[step.route_index+1]-route[step.route_index]);

  return q;
}

//==============================================================================
bool confirmStepReachability(const dart::dynamics::SkeletonPtr& robot,
                             const std::vector<Eigen::VectorXd>& route,
                             const StepContext& step0,
                             StepContext& step1)
{
//  std::cout << step0.route_index << ":" << step0.start_t
//            << " | " << step1.route_index << ":" << step1.start_t << std::endl;

  std::shared_ptr<hubo::RelaxedPosture> posture =
      std::dynamic_pointer_cast<hubo::RelaxedPosture>(robot->getIK(true)->getObjective());

  Eigen::VectorXd q0 = getStepConfig(route, step0);
  robot->setPositions(q0);
  const Eigen::Isometry3d tf_foot0 = step0.stance_foot->getWorldTransform();
  step0.stance_foot->getIK(true)->getTarget()->setTransform(tf_foot0);
  step0.stance_foot->getIK(true)->setActive(true);

//  std::cout << "\n------------------\n"
//            << step0.stance_foot->getName() <<"\n"
//            << tf_foot0.matrix() << std::endl;

  Eigen::VectorXd q1 = getStepConfig(route, step1);
  robot->setPositions(q1);
  const Eigen::Isometry3d tf_foot1 = step1.stance_foot->getWorldTransform();
  step1.stance_foot->getIK(true)->getTarget()->setTransform(tf_foot1);
  step1.stance_foot->getIK(true)->setActive(true);

//  std::cout << "\n" << step1.stance_foot->getName() << "\n"
//            << tf_foot1.matrix() << std::endl;

  assert(step0.stance_foot != step1.stance_foot);

  // Solve while leaning on first foot
  step0.stance_foot->getSupport()->setActive(true);
  step1.stance_foot->getSupport()->setActive(false);
  robot->setPositions(q0);

  if(posture)
    posture->enforceIdealPosture = true;

  robot->getIK(true)->solve();
//  if(posture)
//    posture->enforceIdealPosture = false;
//  if(!robot->getIK(true)->solve())
//    return false;

  step1.start_q = robot->getPositions();

  // Solve while leaning on second foot
  step1.stance_foot->getSupport()->setActive(true);
  step0.stance_foot->getSupport()->setActive(false);
  robot->setPositions(q1);

  if(posture)
    posture->enforceIdealPosture = true;

  robot->getIK(true)->solve();
//  if(posture)
//    posture->enforceIdealPosture = false;

//  if(!robot->getIK(true)->solve())
//    return false;

  step1.end_q = robot->getPositions();

  return true;
}

//==============================================================================
bool testNextStepContext(const dart::dynamics::SkeletonPtr& robot,
                         const std::vector<Eigen::VectorXd>& route,
                         const StepContext& lastStep,
                         StepContext& nextStep)
{
  size_t start_index = lastStep.route_index;

  if(start_index >= route.size())
    return true;

  std::tie(nextStep.route_index, nextStep.start_t) =
      getNextRouteLocation(route, nextStep.distance,
                           start_index, lastStep.start_t);

  return confirmStepReachability(robot, route, lastStep, nextStep);
}

//==============================================================================
bool getNextStepContext(const dart::dynamics::SkeletonPtr& robot,
                        const std::vector<Eigen::VectorXd>& route,
                        const StepContext& lastStep,
                        StepContext& nextStep)
{
  const size_t maxAttempts = 20;
  size_t attempts = 0;
  while(!testNextStepContext(robot, route, lastStep, nextStep))
  {
    nextStep.distance /= 2.0;

    ++attempts;
//    std::cout << "Attempt #" << attempts << " out of " << maxAttempts << std::endl;
    if(attempts > maxAttempts)
    {
//      std::cout << "[getNextStepContext] ERROR: Tried " << attempts << " times "
//                << "to find the next step, but it keeps failing! Reverting!"
//                << std::endl;

      return false;
    }
  }

  return true;
}

//==============================================================================
std::vector<StepContext> generateSteps(
    const dart::dynamics::SkeletonPtr& robot,
    const std::vector<Eigen::VectorXd>& route,
    double stepDistance = DefaultStepLength)
{
  robot->setPositions(route.front());

  StepContext lastStep(robot->getEndEffector("l_foot"),
                       robot->getEndEffector("r_foot"),
                       0, 0, stepDistance);
  lastStep.start_q = robot->getPositions();
  lastStep.end_q = lastStep.start_q;

  std::vector<StepContext> steps;
  steps.push_back(lastStep);

  std::swap(lastStep.stance_foot, lastStep.swing_foot);
  steps.push_back(lastStep);

  const size_t maxAttempts = 3;
  while(steps.back().route_index < route.size()-1
        && steps.size() > 1)
  {
    std::cout << "Steps: " << steps.size() << " (" <<
              100*steps.back().route_index/route.size() << "%)" << std::endl;
    StepContext& previousStep = steps[steps.size()-2];
    StepContext& currentStep = steps.back();
    if(!getNextStepContext(robot, route, previousStep, currentStep))
    {
      previousStep.distance /= 2.0;
      ++previousStep.attempts;
      steps.pop_back();
//      std::cout << " === stepping failed " << std::endl;
//      return steps;

      if(previousStep.attempts >= maxAttempts)
      {
//        steps.pop_back();
//        steps.back().distance /= 2.0;
//        ++steps.back().attempts;

        std::cout << " ==== Failing to generate steps!" << std::endl;

//        osgDart::Viewer viewer;
//        viewer.addWorldNode(new osgDart::WorldNode(world));
//        viewer.setUpViewInWindow(0, 0, 640, 480);
//        viewer.run();

        return steps;
      }
      continue;
    }

    StepContext nextStep(currentStep);
    std::swap(nextStep.stance_foot, nextStep.swing_foot);
    nextStep.distance = stepDistance;

    steps.push_back(nextStep);
  }

  if(steps.size() < 2)
  {
//    std::cout << " ======= Failed to generate steps! ===== " << std::endl;
    return steps;
  }

  // The last step is extraneous
  steps.pop_back();

//  std::cout << " --- Generating finish --- " << std::endl;
//  std::cout << "Last " << steps.back().route_index << " : " << steps.back().start_t << std::endl;
  // One last finishing step
  StepContext nextStep(steps.back());
  std::swap(nextStep.stance_foot, nextStep.swing_foot);
  nextStep.distance = stepDistance;
  getNextStepContext(robot, route, steps.back(), nextStep);
  nextStep.end_q = route.back();
  steps.push_back(nextStep);

  return steps;
}

//==============================================================================
void interpolateSteps(const dart::dynamics::SkeletonPtr& robot,
                      std::vector<Eigen::VectorXd>& walk,
                      std::vector<size_t>& stepRanges,
                      const StepContext& stepFrom,
                      const StepContext& stepTo,
                      bool excludeStep = false)
{
  std::shared_ptr<hubo::RelaxedPosture> posture =
      std::dynamic_pointer_cast<hubo::RelaxedPosture>(robot->getIK(true)->getObjective());

  dart::dynamics::EndEffector* stance = stepFrom.stance_foot;
  dart::dynamics::EndEffector* swing = stepFrom.swing_foot;

  if(!excludeStep)
  {
    stance->getSupport()->setActive(true);
    swing->getSupport()->setActive(false);

    robot->setPositions(stepFrom.end_q);
    stance->getIK(true)->getTarget()->setTransform(
          stance->getWorldTransform());
    const Eigen::Isometry3d swingStart = swing->getWorldTransform();

    robot->setPositions(stepTo.start_q);
    const Eigen::Isometry3d swingGoal = swing->getWorldTransform();

    const Eigen::Vector2d v =
        (swingGoal.translation()-swingStart.translation()).block<2,1>(0,0)/(2*M_PI);
    const double D = v.norm();

    Eigen::AngleAxisd aa(swingStart.linear().transpose()*swingGoal.linear());
    const double R = aa.angle()/(2*M_PI);
    const Eigen::Vector3d axis = aa.axis();

    const size_t res = 10;
    stepRanges.push_back(walk.size());
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
        return;
      }

      walk.push_back(robot->getPositions());
    }
  }

  const Eigen::VectorXd& start_q = stepTo.start_q;
  const Eigen::VectorXd& end_q = stepTo.end_q;

  if(excludeStep)
  {
    robot->setPositions(start_q);
    stance->getIK(true)->getTarget()->setTransform(
          stance->getWorldTransform());

    swing->getIK(true)->getTarget()->setTransform(
          swing->getWorldTransform());
  }

  const Eigen::Vector6d x0 = start_q.block<6,1>(0,0);
  const Eigen::Vector6d xf = end_q.block<6,1>(0,0);
  const Eigen::Vector6d dx = xf - x0;

  // Consider disabling the IK from moving the hips during this stage
  stance->getSupport()->setActive(true);
  swing->getSupport()->setActive(true);

  const size_t res = 10;
  stepRanges.push_back(walk.size());
  for(size_t i=0; i < res+1; ++i)
  {
    const double t = (double)(i)/(double)(res);

    const Eigen::Vector6d x = x0 + t*dx;
    robot->getJoint(0)->setPositions(x);

    if(!robot->getIK(true)->solve())
    {
      walk.push_back(robot->getPositions());
      std::cout <<" ===== COULD NOT SOLVE FOOT STEP AT INDEX " << walk.size()-1 << std::endl;
      return;
    }

    walk.push_back(robot->getPositions());
  }
}

//==============================================================================
std::vector<Eigen::VectorXd> convertRouteToWalk(
    const dart::dynamics::SkeletonPtr& robot,
    const Trajectory& trajectory,
    std::vector<size_t>& stepRanges)
{
  const std::vector<Eigen::VectorXd>& route = trajectory;

  auto originalBoundsLeft =
      robot->getEndEffector("l_foot")->getIK(true)->getErrorMethod().getBounds();
  auto originalBoundsRight =
      robot->getEndEffector("r_foot")->getIK(true)->getErrorMethod().getBounds();

  Eigen::Vector6d bounds(Eigen::Vector6d::Constant(
                          dart::dynamics::DefaultIKTolerance));
  robot->getEndEffector("l_foot")->getIK(true)->getErrorMethod().setBounds(-bounds, bounds);
  robot->getEndEffector("r_foot")->getIK(true)->getErrorMethod().setBounds(-bounds, bounds);

  std::shared_ptr<hubo::RelaxedPosture> posture =
      std::dynamic_pointer_cast<hubo::RelaxedPosture>(
        robot->getIK(true)->getObjective());
  if(posture)
    posture->enforceIdealPosture = true;

//  std::cout << "Generate steps" << std::endl;
  std::vector<StepContext> steps = generateSteps(robot, route);

  std::vector<Eigen::VectorXd> walk;

  StepContext firstStep = steps[0];
  firstStep.start_q = route[0];
  firstStep.end_q = steps[1].start_q;
  interpolateSteps(robot, walk, stepRanges, firstStep, firstStep, true);
  stepRanges.clear();
  for(size_t i=0; i < steps.size()-1; ++i)
  {
    const StepContext& stepFrom = steps[i];
    const StepContext& stepTo = i < steps.size()-1? steps[i+1] : steps[i];
    interpolateSteps(robot, walk, stepRanges, stepFrom, stepTo);
  }

  robot->getEndEffector("l_foot")->getIK(true)->getErrorMethod().setBounds(originalBoundsLeft);
  robot->getEndEffector("r_foot")->getIK(true)->getErrorMethod().setBounds(originalBoundsRight);

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
  robot->getIK(true)->getSolver()->setNumMaxIterations(1000);
//  std::cout << robot->getIK(true)->getSolver()->getNumMaxIterations() << std::endl;
  std::vector<Eigen::VectorXd> locations;
  locations.push_back(robot->getPositions());
  locations.push_back(robot->getPositions());

//  locations[1][3] += 3*DefaultStepLength;
  locations[1][3] += 4*DefaultStepLength;
//  locations[1][3] += 2.0;

  std::vector<size_t> stepRanges;
  std::vector<Eigen::VectorXd> walk = convertRouteToWalk(robot, locations, stepRanges);
  std::cout << "waypoints: " << walk.size() << std::endl;
  std::cout << "Step Ranges: ";
  for(size_t i=0; i < stepRanges.size(); ++i)
    std::cout << stepRanges[i] << ", ";
  std::cout << std::endl;

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

  HuboPath::Trajectory trajectory;
  for(size_t r=0; r < stepRanges.size()+1; ++r)
  {
    size_t start_i = r==0? 0 : stepRanges[r-1]-1;
    size_t end_i = r==stepRanges.size()? walk.size() : stepRanges[r];
    std::cout << "step from " << start_i << " to " << end_i << std::endl;

    mOperator.clearWaypoints();
    for(size_t i=start_i; i < end_i; ++i)
    {
      robot->setPositions(walk[i]);
      mOperator.addWaypoint(robot->getPositions(opIndices));
    }

    HuboPath::Trajectory partial_trajectory = mOperator.getCurrentTrajectory();

    if(r > 0)
      partial_trajectory.elements.erase(partial_trajectory.elements.begin());

    HuboCan::HuboDescription& desc = partial_trajectory.desc;
    for(size_t i=0; i < desc.getJointCount(); ++i)
    {
      hubo_joint_info_t& info = desc.joints[i]->info;
      hubo_joint_limits_t& limits = info.limits;

      limits.nominal_speed /= 5.0;
      limits.nominal_accel /= 5.0;
    }

    std::cout << "Initial size of partial: " << partial_trajectory.elements.size() << std::endl;
    std::cout << "Interpolating #" << r << "..." << std::endl;
    if(!partial_trajectory.interpolate(HUBO_PATH_OPTIMAL))
    {
      break;
    }
    std::cout << "Size of partial: " << partial_trajectory.elements.size() << std::endl;

    if(r==0)
    {
      trajectory = partial_trajectory;
    }
    else
    {
      for(size_t k=0; k < partial_trajectory.size(); ++k)
        trajectory.elements.push_back(partial_trajectory.elements[k]);
    }
  }


  bool operate = true;
  operate = false;

  if(operate)
  {
    if(!trajectory.check_limits())
    {
      std::cout << "Not sending walk trajectory, because it violates limits!" << std::endl;
      return 1;
    }

    mOperator.sendNewTrajectory(trajectory);
  }
  else
  {
    if(!trajectory.check_limits())
    {
      std::cout << "Walk trajectory violates limits!" << std::endl;
    }

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























