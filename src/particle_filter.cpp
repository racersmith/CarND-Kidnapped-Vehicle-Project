/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */


#define _USE_MATH_DEFINES

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Number of particles
	num_particles = 100;

	// Gaussian noise generators
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std[0]);
	std::normal_distribution<double> dist_y(y, std[1]);
	std::normal_distribution<double> dist_theta(theta, std[2]);

	// Weight initialize value
	double init_weight = 1;

	// Generate particles around initial state estimate
	for (int i=0; i < num_particles; i++) {
		Particle temp_particle;
		temp_particle.id = i;
		temp_particle.weight = init_weight;
		temp_particle.x = dist_x(gen);
		temp_particle.y = dist_y(gen);

		// Sample and normalize heading angle
		double theta = dist_theta(gen);
		theta = NormalizeAngle(theta);
		temp_particle.theta = theta;

		// Push particle and weight into their vectors
		particles.push_back(temp_particle);
		weights.push_back(init_weight);
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// random sample generator
	std::default_random_engine gen;

	// apply motion model and noise to each particle
	for (int i = 0; i < num_particles; i++) {
		// extract particle state for readability
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		// define future state parameters
		double x_f = x;
		double y_f = y;
		double theta_f = theta_f;

		if (fabs(theta_f) > 0.0001) {
			// common calculations
			double vel_yaw_rate = velocity / yaw_rate;

			// apply motion model
			theta_f += yaw_rate*delta_t;
			x_f += vel_yaw_rate*( sin(theta_f) - sin(theta));
			y_f += vel_yaw_rate*(-cos(theta_f) + cos(theta));
		}
		else {
			// common calculation
			double v_dt = velocity * delta_t;

			// apply motion model
			x_f	+= v_dt * cos(theta);
			y_f += v_dt * sin(theta);
		}

		// Gaussian noise generators
		std::normal_distribution<double> dist_x(x_f, std_pos[0]);
		std::normal_distribution<double> dist_y(y_f, std_pos[1]);
		std::normal_distribution<double> dist_theta(theta_f, std_pos[2]);

		// Add noise to particle
		x_f += dist_x(gen);
		y_f += dist_y(gen);
		theta_f += dist_theta(gen);

		// Write predicted particle state
		particles[i].x = x_f;
		particles[i].y = y_f;
		particles[i].theta = theta_f;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Find nearest landmark to each observation
	for (int i = 0; i < observations.size(); i++) {
		// extract for readability
		double x_o = observations[i].x;
		double y_o = observations[i].y;
		
		// Search for minimum distance
		double range_min = -1;
		int nearest_id = -1;
		for (int j = 0; i < predicted.size(); j++) {
			// extract for readability
			double x_p = predicted[j].x;
			double y_p = predicted[j].y;

			// calculate range
			double dx = x_o - x_p;
			double dy = y_o - y_p;
			double range = std::sqrt(dx*dx + dy*dy);
			
			// Check for better match or first pass
			if (range < range_min || range_min < 0) {
				range_min = range;
				nearest_id = predicted[j].id;
			}
		}

		// set nearest landmark id
		observations[i].id = nearest_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}

double ParticleFilter::NormalizeAngle(double angle) {
	//std::cout << angle << " -> ";
	// Constrain to less than pi
	while (angle > M_PI) angle -= 2.0*M_PI;

	// Constrain to greater than -pi
	while (angle < -M_PI) angle += 2.0*M_PI;
	//std::cout << angle << std::endl;
	return angle;
}