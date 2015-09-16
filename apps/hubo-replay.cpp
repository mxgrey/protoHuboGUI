/*
 * Copyright (c) 2015, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Michael X. Grey <mxgrey@gatech.edu>
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

#include <memory>
#include <iostream>
#include <fstream>

#include <dart/dart.h>
#include <osgDart/osgDart.h>

#include <osg/Timer>

#include "config.h"

const double frequency = 1000.0;

using namespace dart::dynamics;
using namespace dart::simulation;
using namespace osgDart;
class InputHandler : public osgGA::GUIEventHandler
{
public:

  InputHandler(osgDart::Viewer* viewer)
    : mViewer(viewer),
      recording(false)
  {

  }

  virtual bool handle(const osgGA::GUIEventAdapter& ea,
                      osgGA::GUIActionAdapter&) override
  {

    if( osgGA::GUIEventAdapter::KEYDOWN == ea.getEventType() )
    {

      if( ea.getKey() == 'p')
      {
        if(recording)
        {
          mViewer->pauseRecording();
          recording = false;
        }
        else
        {
          mViewer->record(PROJECT_PATH"dump");
          recording = true;
        }

        return true;
      }

      if( ea.getKey() == 'o')
      {
        time_t now = time(0);
        std::string timestr = ctime(&now);
        timestr.erase(timestr.end()-1, timestr.end());
        mViewer->captureScreen(PROJECT_PATH"dump/screenshot - "+timestr+".png");
        return true;
      }

    }

    return false;
  }

protected:

  osgDart::Viewer* mViewer;

  bool recording;

};

class SimulationWorld : public osgDart::WorldNode
{
public:

  SimulationWorld(WorldPtr world, const std::vector<Eigen::VectorXd>& trajectory)
    : osgDart::WorldNode(world), mTrajectory(trajectory)
  {
    mHubo = world->getSkeleton("drchubo");

    double height = mHubo->getPosition(5);
    mHubo->setPositions(mTrajectory[0]);
    mHubo->setPosition(5, height);

    finishedDumping = false;
    q_dump.open(PROJECT_PATH"sim/trajectory.dat");
    vel_dump.open(PROJECT_PATH"sim/velocity.dat");
    com_dump.open(PROJECT_PATH"sim/com.dat");
    zmp_dump.open(PROJECT_PATH"sim/zmp.dat");
    time_dump.open(PROJECT_PATH"sim/time.dat");
    std::cout << "Dumping trajectory data... ";
    std::cout << std::flush;
  }

  void customPreStep() override
  {
    if(mWorld->getTime() < t0)
      return;

    double time = mWorld->getTime() - t0;
    size_t bot_index = floor(frequency*time);
    size_t top_index = ceil(frequency*time);
    double t = time - (double)(bot_index)/frequency;

    if(top_index >= mTrajectory.size()-1)
    {
//      return;
      top_index = mTrajectory.size()-2;
      bot_index = mTrajectory.size()-2;
    }

//    std::cout << "Index: " << top_index << std::endl;

    const Eigen::VectorXd qd = (mTrajectory[top_index]-mTrajectory[bot_index])*t + mTrajectory[bot_index];
    for(size_t j=6; j < mHubo->getNumDofs(); ++j)
    {
      DegreeOfFreedom* dof = mHubo->getDof(j);
      double velocity = qd[j] - dof->getPosition();
      velocity /= mHubo->getTimeStep();

      dof->setCommand(velocity);
    }

  }

  void customPostStep() override
  {
    if(mWorld->getTime() < t0)
      return;

    double traj_time = mWorld->getTime()-t0;

//    std::cout << traj_time << " : " << (double)(mTrajectory.size())/(double)(frequency) << std::endl;

    if(traj_time > (double)(mTrajectory.size())/(double)(frequency))
    {
      // If the current time exceeds the length of the trajectory, we should no
      // longer be dumping data
      if(!finishedDumping)
      {
        q_dump.close();
        vel_dump.close();
        com_dump.close();
        zmp_dump.close();
        time_dump.close();

        finishedDumping = true;
        std::cout << "finished dumping!" << std::endl;
      }

      return;
    }

//    mHubo->computeForwardDynamics();

    for(size_t j=0; j < mHubo->getNumDofs(); ++j)
    {
      const DegreeOfFreedom* dof = mHubo->getDof(j);
      q_dump << dof->getPosition() << "\t";
      vel_dump << dof->getVelocity() << "\t";
    }
    q_dump << "\n";
    vel_dump << "\n";

    Eigen::Vector3d com = mHubo->getCOM();
    for(size_t j=0; j < 3; ++j)
      com_dump << com[j] << "\t";
    com_dump << "\n";

    Eigen::Vector3d zmp = mHubo->getZMP();
    for(size_t j=0; j < 3; ++j)
      zmp_dump << zmp[j] << "\t";
    zmp_dump << "\n";

    time_dump << traj_time << "\n";
  }

protected:

  std::ofstream q_dump;
  std::ofstream vel_dump;
  std::ofstream com_dump;
  std::ofstream zmp_dump;
  std::ofstream time_dump;

  SkeletonPtr mHubo;

  std::vector<Eigen::VectorXd> mTrajectory;

  const double t0 = 1.0;

  bool finishedDumping;
};

SkeletonPtr createHubo()
{
  dart::utils::DartLoader loader;
  loader.addPackageDirectory("drchubo", DART_DATA_PATH"/urdf/drchubo");
  SkeletonPtr hubo =
      loader.parseSkeleton(DART_DATA_PATH"/urdf/drchubo/drchubo.urdf");

  for(size_t i = 0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    if(bn->getName().substr(0, 7) == "Body_LF"
       || bn->getName().substr(0, 7) == "Body_RF"
       || bn->getName().substr(0, 7) == "Body_NK")
    {
      bn->remove();
      --i;
    }
  }

  hubo->setPosition(5, 0.97);

  for(size_t i=1; i < hubo->getNumJoints(); ++i)
  {
    hubo->getJoint(i)->setActuatorType(Joint::VELOCITY);
  }

  for(size_t i=0; i < hubo->getNumBodyNodes(); ++i)
  {
    BodyNode* bn = hubo->getBodyNode(i);
    for(size_t j=0; j < bn->getNumVisualizationShapes(); ++j)
    {
      const ShapePtr& shape = bn->getVisualizationShape(j);
      shape->setColor(Eigen::Vector3d(0.2, 0.2, 0.2));

      if(MeshShapePtr mesh = std::dynamic_pointer_cast<MeshShape>(shape))
        mesh->setColorMode(MeshShape::SHAPE_COLOR);
    }
  }

  hubo->setName("drchubo");
  return hubo;
}

SkeletonPtr createGround()
{
  // Create a Skeleton to represent the ground
  SkeletonPtr ground = Skeleton::create("ground");
  Eigen::Isometry3d tf(Eigen::Isometry3d::Identity());
  double thickness = 0.01;
  tf.translation() = Eigen::Vector3d(0,0,-thickness/2.0);
  WeldJoint::Properties joint;
  joint.mT_ParentBodyToJoint = tf;
  ground->createJointAndBodyNodePair<WeldJoint>(nullptr, joint);
  ShapePtr groundShape =
      std::make_shared<BoxShape>(Eigen::Vector3d(10,10,thickness));
  groundShape->setColor(dart::Color::Blue(0.2));

  ground->getBodyNode(0)->addVisualizationShape(groundShape);
  ground->getBodyNode(0)->addCollisionShape(groundShape);

  return ground;
}

int main()
{
  WorldPtr world(new dart::simulation::World);
  world->setTimeStep(1.0/frequency);

  SkeletonPtr hubo = createHubo();
  world->addSkeleton(hubo);
  world->addSkeleton(createGround());

  std::string name = PROJECT_PATH"ideal/trajectory.dat";

  std::ifstream file;
  file.open(name);
  std::vector<Eigen::VectorXd> trajectory;
  if(file.is_open())
  {
    while(!file.eof())
    {
      trajectory.push_back(Eigen::VectorXd(hubo->getNumDofs()));
      Eigen::VectorXd& q = trajectory.back();
      for(size_t i=0; i < hubo->getNumDofs(); ++i)
        file >> q[i];
    }
  }
  else
  {
    std::cerr << "Could not open file: " << name << std::endl;
  }

  osg::ref_ptr<SimulationWorld> node = new SimulationWorld(world, trajectory);
  node->setNumStepsPerCycle(30);

  osgDart::Viewer viewer;
  viewer.addWorldNode(node);
//  viewer.simulate(true);

  viewer.addEventHandler(new InputHandler(&viewer));

  viewer.setUpViewInWindow(0, 0, 1280, 960);

  viewer.addAttachment(new osgDart::SupportPolygonVisual(hubo));

  // Set up the default viewing position
  viewer.getCameraManipulator()->setHomePosition(osg::Vec3( 5.34,  3.00, 1.91),
                                                 osg::Vec3( 0.00,  0.00, 0.50),
                                                 osg::Vec3(-0.20, -0.08, 0.98));

  viewer.setCameraManipulator(viewer.getCameraManipulator());
//  viewer.record(PROJECT_PATH"dump");

  viewer.run();
}
