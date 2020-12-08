/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>

#include "helper_functions.h"

using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles (40, 100)
  num_particles = 20;

  /*
   * Initialize all particles to first position (based on estimates of x, y, theta and their uncertainties
   * from GPS) and all weights to 1. Add random Gaussian noise to each particle.
   */

  // Create normal (Gaussian) distribution for x, y, theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i <= num_particles; ++i){
    particles.push_back({i, dist_x(random_engine), dist_y(random_engine), dist_theta(random_engine), 1, {}, {}, {}});
    weights.push_back(1.0);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std[],
                                double velocity, double yaw_rate) {
  /**
   * NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */


  for (Particle & particle : particles){
    /*
    update each particles' position estimates
    */
    // avoid division by zero
    if (fabs(yaw_rate) > 0.001) {
      particle.x += velocity/yaw_rate * ( sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta));
      particle.y += velocity/yaw_rate * ( cos(particle.theta) - cos(particle.theta + yaw_rate*delta_t));
    } else {
      particle.x += velocity * delta_t * cos(particle.theta);
      particle.y += velocity * delta_t * sin(particle.theta);
    }
    particle.theta = particle.theta + yaw_rate*delta_t;

    // Create normal (Gaussian) distribution for x, y, theta
    normal_distribution<double> dist_x(particle.x, std[0]);
    normal_distribution<double> dist_y(particle.y, std[1]);
    normal_distribution<double> dist_theta(particle.theta, std[2]);

    /*
    Add Gaussian noise by sampling from the Gaussian distribution with mean updated to
    the predicted particle position and standard deviation equal to the standard deviation of the measurements.
    */
    particle.x = dist_x(random_engine);
    particle.y = dist_y(random_engine);
    particle.theta = dist_theta(random_engine);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs>& predicted,
                                     vector<LandmarkObs>& measured_landmarks) {
  /*
    Perform nearest neighbor data association. Assign each sensor observation a predictred map landmark ID, associated with it.
  */

  for(auto& measured_landmark : measured_landmarks){
    double min_distance = 9999;
    for (auto& predicted_landmark : predicted){
      double distance = dist(predicted_landmark.x, predicted_landmark.y, measured_landmark.x, measured_landmark.y);
      if (distance <= min_distance){
        measured_landmark.id = predicted_landmark.id;
        min_distance = distance;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {

  for(auto &particle : particles){
    // Reset particle weight
    particle.weight = 1;

    // Predict measurements to all map landmarks within sensor range for each particle. Both measurements are in Map coordinate system.
    vector<LandmarkObs> landmarks_global_within_sensor_range = extractWithinRange(sensor_range, particle.x, particle.y, map_landmarks);
    vector<LandmarkObs> measured_landmarks_global = convert_to_map_system(observations, particle.x, particle.y, particle.theta);

    // Associate the sensor measurements to map landmarks.
    dataAssociation(landmarks_global_within_sensor_range, measured_landmarks_global);

    /*
    * Update the weights of each particle using a multi-variate Gaussian distribution. You can read more about this distribution here:
    * https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    */
    for(auto &measured_landmark : measured_landmarks_global){
      for (auto & map_landmark : landmarks_global_within_sensor_range){
        if(measured_landmark.id == map_landmark.id){
          particle.weight *= multivariate_probability(std_landmark, measured_landmark.x, measured_landmark.y, map_landmark.x, map_landmark.y);
        }
      }
    }
  }
  /**
  * Normalize the weights of all particles to the [0;1] range.
  */
  normalize_weights();
}

void ParticleFilter::resample() {
  /**
   * Resample particles with replacement with probability proportional to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::vector<Particle> particles_resampled;

  /*
  * Choose first index randomly for resampling wheel.
  */
  int index = rand() % static_cast<int>(num_particles + 1);
  double beta = 0.0;
  double max_weight = *std::max_element(weights.begin(), weights.end());
  for(int i = 0; i < num_particles; ++i){
    double random_var = ((double) rand() / (RAND_MAX)) +1;
    beta += random_var *2.0 *max_weight;
    while (beta > weights[index]){
      beta -= weights[index];
      index = ((index+1) % num_particles);
    }
    particles_resampled.push_back(particles[index]);
  }
  particles = particles_resampled;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

void ParticleFilter::normalize_weights(){
  weights.clear();

  double denominator = 0.0f;
  for(auto &particle : particles){
    denominator +=particle.weight;
  }

  for(auto &particle : particles){
    particle.weight = particle.weight/denominator;
    weights.push_back(particle.weight);
  }
}
