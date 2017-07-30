/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <random>

#include "helper_functions.h"

typedef struct Particle
{
  int id;
  double x;
  double y;
  double theta;
  double weight;
  std::vector<int> associations;
  std::vector<double> sense_x;
  std::vector<double> sense_y;
} tParticle;



class ParticleFilter
{  
  std::shared_ptr<std::default_random_engine> pRandomNumGen;
  
  // Vector of weights of all particles
  std::vector<double> weights;
  
 public:
  
  // Set of current particles
  std::vector<Particle> particles;
  
  // Constructor
  // @param M Number of particles
  ParticleFilter(double posX,
		 double posY,
		 double theta,
		 double stdDev[3]);


  // Destructor
  ~ParticleFilter() {}


  /**
   * @brief Add noise to the position (x, y, theta) based on guassian distribution
   *        and randomly chosen location on that distribution
   * @param particle The particle to populate
   * @param posX The X position to add noise to
   * @param posY The Y position to add noise to
   * @param theta The angle to add noise to
   * @param stdDev Array of standard deviations for X, Y, Theta
   */
  void addNoise(Particle &particle,
		double posX,
		double posY,
		double theta,
		double stdDev[3]);
  
  /**
   * prediction Predicts the state for the next time step
   *   using the process model.
   * @param deltaT Time between time step t and t+1 in measurements [s]
   * @param stdDevPos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
   *   standard deviation of yaw [rad]]
   * @param velocity Velocity of car from t to t+1 [m/s]
   * @param yawRate Yaw rate of car from t to t+1 [rad/s]
   */
  void prediction(double deltaT,
		  double stdDevPos[3],
		  double velocity,
		  double yawRate);
  
  /**
   * dataAssociation Finds which observations correspond to which landmarks (likely by using
   *   a nearest-neighbors data association).
   * @param predicted predicted landmark from observations
   * @param observation A landmark observations
   */
  void dataAssociation(LandmarkObs& predicted,
		       LandmarkObs observation,
		       Map mapLandMarks);
  
  /**
   * updateWeights Updates the weights for each particle based on the likelihood of the 
   *   observed measurements. 
   * @param sensor_range Range [m] of sensor
   * @param std_landmark[] Array of dimension 2 [standard deviation of range [m],
   *   standard deviation of bearing [rad]]
   * @param observations Vector of landmark observations
   * @param map Map class containing map landmarks
   */
  void updateWeights(double sensorRange,
		     double stdLandmark[],
		     std::vector<LandmarkObs> observations,
		     Map mapLandmarks);
  
  /**
   * resample Resamples from the updated set of particles to form
   *   the new set of particles.
   */
  void resample();
  
  /*
   * Set a particles list of associations, along with the associations calculated world x,y coordinates
   * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
   */
  Particle setAssociations(Particle particle,
			   std::vector<int> associations,
			   std::vector<double> sense_x,
			   std::vector<double> sense_y);
  
  std::string getAssociations(Particle best);
  std::string getSenseX(Particle best);
  std::string getSenseY(Particle best);
  
};

#endif /* PARTICLE_FILTER_H_ */
