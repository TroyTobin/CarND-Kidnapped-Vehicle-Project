/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */


/**********************************************
 *                 Includes
 **********************************************/
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <ctime>
#include <memory>
#include <cassert>

#include "particle_filter.h"
#include "map.h"

using namespace std;

/**********************************************
 *                  Constants
 **********************************************/
#define X_INDEX         (0)
#define Y_INDEX         (1)
#define THETA_INDEX     (2)

#define NUM_PARTICLES  (75)

/**
 * init Initializes particle filter by initializing particles to Gaussian
 *   distribution around first position and all the weights to 1.
 * @param x Initial x position [m] (simulated estimate from GPS)
 * @param y Initial y position [m]
 * @param theta Initial orientation [rad]
 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 */
ParticleFilter::ParticleFilter(double posX, 
			       double posY, 
			       double theta, 
			       double stdDev[3]) 
{
  vector<Particle>::iterator particleIter;
  
  // Seed the random number generator with the current time.
  pRandomNumGen = make_shared<default_random_engine>(time(NULL));
  
  // Create the particles
  particles = vector<Particle>(NUM_PARTICLES);
  
  // Set weights to 1
  weights = vector<double>(NUM_PARTICLES, 1);

  particleIter = particles.begin();
  
  // Set the particle values
  for (; particleIter != particles.end(); particleIter++)
  {
    // Set particle with position randomly chosen from the normal distributions
    particleIter->id = (particleIter - particles.begin());
    addNoise(*particleIter, posX, posY, theta, stdDev);
    particleIter->weight = weights[particleIter - particles.begin()];
  }
}

/**********************************************
 *          Function Implementations
 **********************************************/
void ParticleFilter::addNoise(Particle &particle,
                              double posX,
                              double posY, 
                              double theta, 
                              double stdDev[3])
{  
  // Expect random number generator to have been created already
  assert(pRandomNumGen);
  
  // Create normal distribution around the initial input positions (x, y, theta)
  // Particles will be randomly chosen from these distributions
  normal_distribution<double> NormalX(posX, stdDev[X_INDEX]);
  normal_distribution<double> NormalY(posY, stdDev[Y_INDEX]);
  normal_distribution<double> NormanTheta(theta, stdDev[THETA_INDEX]);
  
  // Choose a value from the normal distribution (i.e. has some noise defined by the distribution)
  particle.x     = NormalX(*pRandomNumGen);
  particle.y     = NormalY(*pRandomNumGen);
  particle.theta = NormanTheta(*pRandomNumGen);
}


/**
 * prediction Predicts the state for the next time step
 *   using the process model.
 * @param deltaT Time between time step t and t+1 in measurements [s]
 * @param stdDevPos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
 *   standard deviation of yaw [rad]]
 * @param velocity Velocity of car from t to t+1 [m/s]
 * @param yawRate Yaw rate of car from t to t+1 [rad/s]
 */
void ParticleFilter::prediction(double deltaT,
                                double stdDevPos[3], 
                                double velocity, 
                                double yawRate) 
{
  vector<Particle>::iterator particleIter = particles.begin();
  
  for (; particleIter != particles.end(); particleIter++)
  {
    double predictedX;
    double predictedY;
    double predictedTheta;

    // Need to handle prediction differently if yaw rate is 0
    if (fabs(yawRate) > numeric_limits<float>::epsilon())
    {
      double yaw = yawRate * deltaT;
      predictedX = particleIter->x + (velocity * (sin(particleIter->theta + yaw) - sin(particleIter->theta))) / yawRate;
      predictedY = particleIter->y + (velocity * (cos(particleIter->theta) - cos(particleIter->theta + yaw))) / yawRate;
    }
    else
    {
      predictedX = particleIter->x + (velocity * cos(particleIter->theta) * deltaT);
      predictedY = particleIter->y + (velocity * sin(particleIter->theta) * deltaT);
    }

    predictedTheta = particleIter->theta + yawRate * deltaT;

    addNoise(*particleIter, predictedX, predictedY, predictedTheta, stdDevPos);
  }
}

void ParticleFilter::dataAssociation(LandmarkObs& predicted,
				     LandmarkObs observation,
				     Map mapLandmarks)
{
  Map::MapIterator landmarkIter = mapLandmarks.begin();
  double closestDist = -1;
  
  for (; landmarkIter != mapLandmarks.end(); landmarkIter++)
  {
    double testDist = dist(observation.x, observation.y, landmarkIter->x_f, landmarkIter->y_f);
    if ((closestDist == -1) ||
	(testDist < closestDist))
    {
      predicted.id = landmarkIter->id_i;
      predicted.x  = landmarkIter->x_f;
      predicted.y  = landmarkIter->y_f; 
      closestDist = testDist;
    }
  }
}

LandmarkObs GetLandmarkByID(Map landmarks,
			    int id)
{
  Map::MapIterator landmarkIter = landmarks.begin();
  LandmarkObs obs;
  obs.id = -1;
 
  for (; landmarkIter != landmarks.end(); landmarkIter++)
  {
    if ((*landmarkIter).id_i == id)
    {
      obs.id = id;
      obs.x = landmarkIter->x_f;
      obs.y = landmarkIter->y_f;
      break;
    }
  }
  return obs;
}


void ParticleFilter::updateWeights(double sensorRange, 
                                   double stdLandmark[], 
		                   std::vector<LandmarkObs> observations, 
                                   Map mapLandmarks)
{
  vector<LandmarkObs>::iterator obsIter = observations.begin();
  for (; obsIter != observations.end(); obsIter++)
  {
    // Apply the observation to the particle
    vector<Particle>::iterator particleIter = particles.begin();
    for(; particleIter != particles.end(); particleIter++)
    {
      if (obsIter == observations.begin())
      {
	particleIter->weight = 1;
      }
            
      // Convert the observation to the particle coordinate system (global/Map coordinates)
      // This essentially creates a new observation in the global cooridnate system for this particle.
      LandmarkObs mapCoordObs;
      LandmarkObs closestLandmark;
      
      // Translation equation from http://planning.cs.uiuc.edu/node99.html
      mapCoordObs.x = (obsIter->x * cos(particleIter->theta)) - (obsIter->y * sin(particleIter->theta)) + particleIter->x;
      mapCoordObs.y = (obsIter->x * sin(particleIter->theta)) + (obsIter->y * cos(particleIter->theta)) + particleIter->y;
      
      dataAssociation(closestLandmark, mapCoordObs, mapLandmarks);

      // Update weights.
      // P(x, y) = (1/(2*pi*sig_x*sig_y)) exp (-((((x - mu_x)^2)/2*sig_x^2) + ((y - mu_y)^2)/2*sig_y^2))'
      double diffX  = mapCoordObs.x - closestLandmark.x;
      double diffX2 = diffX * diffX;

      double diffY  = mapCoordObs.y - closestLandmark.y;
      double diffY2 = diffY * diffY;

      double stdLandmark01  = stdLandmark[0] * stdLandmark[1];
      double stdLandmark0_2 = 2 * stdLandmark[0] * stdLandmark[0];
      double stdLandmark1_2 = 2 * stdLandmark[1] * stdLandmark[1];
      double div = 2 * M_PI * stdLandmark01;
     
      double Pxy = exp(-((diffX2 / stdLandmark0_2) + (diffY2 / stdLandmark1_2))) / div;
      
      particleIter->weight *= Pxy;
      weights[particleIter - particles.begin()] = particleIter->weight;
    }
  }
}


void ParticleFilter::resample()
{
  assert(pRandomNumGen);
  
  discrete_distribution<int> distribution(weights.begin(), weights.end());
  
  // Vector of particles to populate from re-sampled particles (based on current weights)
  vector<Particle> resampleParticles;
  
  vector<Particle>::iterator particleIter = particles.begin();
  
  for (; particleIter != particles.end(); particleIter++)
  {
    resampleParticles.push_back(particles[distribution(*pRandomNumGen)]);
  }
  
  particles = resampleParticles;
}


Particle ParticleFilter::setAssociations(Particle particle,
					 std::vector<int> associations,
					 std::vector<double> sense_x,
					 std::vector<double> sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  
  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();
  
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  
  return particle;
}

#define VECTOR_TO_STRING(In, Type, Name) do {				\
    stringstream strStream;						\
    copy(In.begin(), In.end(), ostream_iterator<Type>(strStream, ""));	\
    Name = strStream.str();							\
    Name = Name.substr(0, (Name.length() - 1));				\
  } while(0);

string ParticleFilter::getAssociations(Particle best)
{
  string strOutput;
  vector<int> v = best.associations;
  VECTOR_TO_STRING(v, int, strOutput);
  return strOutput;
}
 
string ParticleFilter::getSenseX(Particle best)
{
  string strOutput;
  vector<double> v = best.sense_x;
  VECTOR_TO_STRING(v, double, strOutput);
  return strOutput;
}
 
string ParticleFilter::getSenseY(Particle best)
{
  string strOutput;
  vector<double> v = best.sense_y;
  VECTOR_TO_STRING(v, double, strOutput);
  return strOutput;
}
