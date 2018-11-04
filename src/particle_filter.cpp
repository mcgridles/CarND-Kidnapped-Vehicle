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

	if (is_initialized) {
	    return;
	}

    num_particles = 100;

    // Create normal distributions
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i=0; i<num_particles; i++) {
        Particle p;

        p.id = i;

        // Sample from each distribution
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);

        p.weight = 1.0;

        particles.push_back(p);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for (Particle& p : particles) {
	    double x, y;
	    double theta = p.theta;

	    // Find new x, y, and theta (mean) values
	    if (fabs(yaw_rate) < 0.00001) {
            x = p.x + velocity * delta_t * cos(p.theta);
            y = p.y + velocity * delta_t * sin(p.theta);
	    } else{
            x = p.x + (velocity/yaw_rate) * (sin(p.theta + yaw_rate*delta_t) - sin(p.theta));
            y = p.y + (velocity/yaw_rate) * (cos(p.theta) - cos(p.theta + yaw_rate*delta_t));
            theta += yaw_rate * delta_t;
	    }

	    // Normal distributions centered around prediction value to account for noise
        normal_distribution<double> dist_x(x, std_pos[0]);
        normal_distribution<double> dist_y(y, std_pos[1]);
        normal_distribution<double> dist_theta(theta, std_pos[2]);

        // Sample distribution
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (LandmarkObs& obs : observations) {
        int map_id = -1;
        double min_dist = numeric_limits<double>::max();

	    for (LandmarkObs pred : predicted) {
	        double cur_dist = dist(obs.x, obs.y, pred.x, pred.y);

	        if (cur_dist < min_dist) {
	            min_dist = cur_dist;
                map_id = pred.id;
	        }
	    }

	    obs.id = map_id;
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

    double std_x = std_landmark[0];
    double std_y = std_landmark[1];

    for (Particle& p : particles) {
        double x = p.x;
        double y = p.y;
        double theta = p.theta;

        std::vector<LandmarkObs> visible_landmarks;
        for (Map::single_landmark_s landmark : map_landmarks.landmark_list) {
            double distance = dist(x, y, landmark.x_f, landmark.y_f);

            // Check if landmark is within sensor range of each particle
            if (fabs(distance) <= sensor_range) {
                visible_landmarks.push_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
            }
        }

        // Transform observations to map coordinates
        std::vector<LandmarkObs> transformed_observations;
        for (LandmarkObs obs : observations) {
            double t_x = x + cos(theta)*obs.x - sin(theta)*obs.y;
            double t_y = y + sin(theta)*obs.x + cos(theta)*obs.y;

            transformed_observations.push_back(LandmarkObs{obs.id, t_x, t_y});
	    }

	    dataAssociation(visible_landmarks, transformed_observations);

        // Calculate new weight
        p.weight = 1.0;
        for (LandmarkObs trans : transformed_observations) {
            double observation_x, observation_y, landmark_x, landmark_y;

            observation_x = trans.x;
            observation_y = trans.y;

            // Find associated landmark
            for (LandmarkObs landmark : visible_landmarks) {
                if (landmark.id == trans.id) {
                    landmark_x = landmark.x;
                    landmark_y = landmark.y;
                    break;
                }
            }

            // Calculate weight
            double normalize = 1 / (2 * M_PI * std_x * std_y);
            double exp_x = pow(observation_x - landmark_x, 2) / (2 * pow(std_x, 2));
            double exp_y = pow(observation_y - landmark_y, 2) / (2 * pow(std_y, 2));
            double weight = normalize * exp( -(exp_x + exp_y) );

            if (weight == 0) {
                p.weight *= 0.00001;
            } else {
                p.weight *= weight;
            }
        }
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::vector<Particle> new_particles;

    // Find max particle weight
    std::vector<double> weights;
    for (Particle p : particles) {
        weights.push_back(p.weight);
    }

    double max_weight = *max_element(weights.begin(), weights.end());

    std::uniform_int_distribution<int> dist_index(0,num_particles-1);
    std::uniform_real_distribution<double> dist_beta(0, 2*max_weight);
    double beta = 0.0;

    // Resample new particles
    for (int i=0; i<num_particles; i++) {
        int idx = dist_index(gen);
        beta += dist_beta(gen);

        while (weights[idx] < beta) {
            beta -= weights[idx];
            idx = (idx + 1) % num_particles;
        }

        new_particles.push_back(particles[idx]);
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
