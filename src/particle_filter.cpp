/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;
static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 20;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	int i;
	for (i = 0; i < num_particles; i++) {
	  Particle particle;
	  particle.id = i;
	  particle.x = dist_x(gen);
	  particle.y = dist_y(gen);
	  particle.theta = dist_theta(gen);
	  particle.weight = 1.0;
	  
	  particles.push_back(particle);
	  weights.push_back(1.0);
	}
	is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	normal_distribution<double> N_x(0, std_pos[0]);
	normal_distribution<double> N_y(0, std_pos[1]);
	normal_distribution<double> N_theta(0, std_pos[2]);

	for (int i = 0; i < num_particles; i++) {

	
	if (fabs(yaw_rate) < 0.001) {  
	  particles[i].x = particles[i].x +  velocity * delta_t * cos(particles[i].theta) + N_x(gen);
	  particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta) + N_x(gen);
	} 
	else {
	  particles[i].x += (N_x(gen) + (velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta))));
	  particles[i].y += (N_x(gen) + (velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t))));
	  particles[i].theta +=  (N_x(gen) +  yaw_rate * delta_t);
	}
	}	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  for (unsigned int i = 0; i < observations.size(); i++) {
    
    double min_dist = numeric_limits<double>::max();

    int location = -1;
    
    for (unsigned int j = 0; j < predicted.size(); j++) {
      
      double cur_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        location = predicted[j].id;
      }
    }

    observations[i].id = location;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  int i, j;
  /*This variable is used for normalizing weights of all particles and bring them in the range
    of [0, 1]*/

  //filter landmarks wihthin range


	
    // transfor coordinate frames
  double normalizer = 0.0;

  for (i = 0; i < num_particles; i++) {

    // vector<LandmarkObs> within_range_landmarks;
    // for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
    //   float landmark_x = map_landmarks.landmark_list[j].x_f;
    //   float landmark_y = map_landmarks.landmark_list[j].y_f;
    //   int id = map_landmarks.landmark_list[j].id_i;
    //   double dx = x - landmark_x;
    //   double dy = y - landmark_y;
    //   if ( dx*dx + dy*dy <= sensor_range * sensor_range ) {
    //     within_range_landmarks.push_back(LandmarkObs{ id, landmark_x, landmark_y });
    //   }
    // }
    // filter within tange landmarks
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    vector<LandmarkObs> transformed_observations;
     vector<LandmarkObs> within_range_landmarks;
    for (j = 0; j < map_landmarks.landmark_list.size(); j++) {
      Map::single_landmark_s current_landmark = map_landmarks.landmark_list[j];
      float landmark_dx = fabs((x - current_landmark.x_f));
      float landmark_dy = fabs((y - current_landmark.y_f));
      if (( landmark_dx <= sensor_range) && (landmark_dy <= sensor_range)) {
        within_range_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
      }
    } 	


    // transform coordinate frame
    for (j = 0; j < observations.size(); j++) {
      LandmarkObs transformed_obs;
      transformed_obs.id = j;
      transformed_obs.x = x + (cos(theta) * observations[j].x) - (sin(theta) * observations[j].y);
      transformed_obs.y = y + (sin(theta) * observations[j].x) + (cos(theta) * observations[j].y);
      transformed_observations.push_back(transformed_obs);
    }
   //Associate observations
    dataAssociation(within_range_landmarks, transformed_observations);

    /*Step 4: Calculate the weight of each particle using Multivariate Gaussian distribution.*/
    //Reset the weight of particle to 1.0
    particles[i].weight = 1.0;

    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    double sigma_x_2 = pow(sigma_x, 2);
    double sigma_y_2 = pow(sigma_y, 2);
    int k, l;
    
    // Calculate weight of particle
    for (k = 0; k < transformed_observations.size(); k++) {
      double trans_obs_x = transformed_observations[k].x;
      double trans_obs_y = transformed_observations[k].y;
      double trans_obs_id = transformed_observations[k].id;
      double multi_prob = 1.0;

      for (l = 0; l < within_range_landmarks.size(); l++) {
        double pred_landmark_x = within_range_landmarks[l].x;
        double pred_landmark_y = within_range_landmarks[l].y;
        double pred_landmark_id = within_range_landmarks[l].id;

        if (trans_obs_id == pred_landmark_id) {
          multi_prob = (1.0/(2.0 * M_PI * sigma_x * sigma_y)) * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
          particles[i].weight *= multi_prob;
        }
      }
    }
    normalizer += particles[i].weight;
  }

  // Normalise
  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= normalizer;
    weights[i] = particles[i].weight;
}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::vector<Particle> new_particles;
  int index;

  // Initializes discrete distribution function
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int> weight_distribution(weights.begin(), weights.end());

  for (int i = 0; i < num_particles; i++) {
    index = weight_distribution(gen);
    new_particles.push_back(particles[index]);
  }
  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
