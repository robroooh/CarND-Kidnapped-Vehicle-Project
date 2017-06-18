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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	std::cout << "Initializing ParticleFilter..." << std::endl;
	default_random_engine gen;

	num_particles = 1;
	double std_x = std[0];
	double std_y = std[1];
	double std_yaw = std[2];

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std_x);
	
	// TODO: Create normal distributions for y and psi
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_yaw(theta, std_yaw);

	for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_psi;

		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_psi = dist_yaw(gen);	 
		
		Particle p;
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_psi;
		p.weight = 1.0;
		weights.push_back(1.0);
		particles.push_back(p);
	}
	std::cout << "Initialized ParticleFilter!" << std::endl;
	std::cout << weights[0] << std::endl;
	std::cout << particles[0].weight << std::endl;
	is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
		default_random_engine gen;
		for (int i = 0; i < num_particles; ++i) {
			double std_x = std_pos[0];
			double std_y = std_pos[1];
			double std_yaw = std_pos[2];
			double x0 = particles[i].x;
			double y0 = particles[i].y;
			double theta0 = particles[i].theta;
			double vel_over_yaw = velocity/yaw_rate;
			double pred_x;
			double pred_y;
			double pred_theta;


			if (fabs(yaw_rate) > 0.00001) {
				pred_x = x0 + (vel_over_yaw)*(sin(theta0+yaw_rate*delta_t) - sin(theta0));
				pred_y = y0 + (vel_over_yaw)*(cos(theta0) - cos(theta0 + yaw_rate*delta_t));
				pred_theta = theta0 + (yaw_rate*delta_t);
			} else {
				pred_x = x0 + velocity*delta_t*cos(theta0);
				pred_y = y0 + velocity*delta_t*sin(theta0);
				pred_theta = yaw_rate;
			}
			normal_distribution<double> dist_x(pred_x, std_x);
			normal_distribution<double> dist_y(pred_y, std_y);
			normal_distribution<double> dist_yaw(pred_theta, std_yaw);

			particles[i].x = pred_x + dist_x(gen);
			particles[i].y = pred_y + dist_y(gen);
			particles[i].theta = pred_theta + dist_yaw(gen);
		}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	int min_id = -1;
	double min_x = 0;
	double min_y = 0;

	std::cout <<"data ASS ";
	for(int i = 0; i < observations.size(); i++){
		double min_dist = 1000000000;
		for(int j = 0; j < predicted.size(); j++){
			double x1 = predicted[j].x;
			double x2 = observations[i].x;
			double y1 = predicted[j].y;
			double y2 = observations[i].y;

			double tmp_dist = dist(x1,y1,x2,y2);
			if (tmp_dist <= min_dist) {
				min_id = predicted[j].id;
				min_x = x2;
				min_y = y2;
				min_dist = tmp_dist;
			}
		}
		observations[i].id = min_id;
		observations[i].x = min_x;
		observations[i].y = min_y;
		// std::cout << "Obs id,x,y: "<< observations[i].id
		// << " " << observations[i].x << " " << observations[i].y << std::endl;
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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	weights.clear();
	for(int i = 0; i < particles.size(); i++){

		std::vector<LandmarkObs> landmark_inrange;
		std::vector<LandmarkObs> observations_pspace;

		std::cout << "----------Landmark in range ---------------- "<< std::endl;

		// Get only landmark that is in range
		for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
			LandmarkObs lm_tmp;
			lm_tmp.id = map_landmarks.landmark_list[j].id_i;
			lm_tmp.x = map_landmarks.landmark_list[j].x_f;
			lm_tmp.y = map_landmarks.landmark_list[j].y_f;
			
			double dist_x = fabs(lm_tmp.x - particles[i].x);
			double dist_y = fabs(lm_tmp.y - particles[i].y);
			if (dist_x <= sensor_range && dist_y <= sensor_range) {
				std::cout << "Landmark id,x,y: "<< lm_tmp.id
				 << " " << lm_tmp.x << " " << lm_tmp.y << std::endl;
				landmark_inrange.push_back(lm_tmp);
			}
		}

		std::cout << "-------Obs tranform -------: "<< std::endl;
		// Transform Obs to particle space
		for(int j = 0; j < observations.size(); j++){
			LandmarkObs obs_tmp;
			obs_tmp.id = observations[j].id;
			obs_tmp.x = observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			obs_tmp.y = observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y;
			std::cout << "From " << observations[j].x << " " 
			<< observations[j].y << " to " << obs_tmp.x << " " << obs_tmp.y << std::endl;
			observations_pspace.push_back(obs_tmp);
		}

		// Parse to this function to see which observation belongs to which landmark
		dataAssociation(landmark_inrange, observations_pspace);
		std::cout << "---------after dataAss ----------" << std::endl;
		for(int j = 0; j < observations_pspace.size(); j++){
			LandmarkObs obs_tmp;
			obs_tmp.id = observations_pspace[j].id;
			obs_tmp.x = observations_pspace[j].x;
			obs_tmp.y = observations_pspace[j].y;
			std::cout << "Observations id,x,y: "<< obs_tmp.id << " " << obs_tmp.x << " " << obs_tmp.y << std::endl;
		}

		double std_x = std_landmark[0];
		double std_y = std_landmark[1];
		std::cout << "---------Particle's Weight update ----------" << std::endl;
		std::cout << "Particle's Weight Before update: "<< particles[i].weight << std::endl;
		
		// UpdateWeights(really)
		for(int j = 0; j < observations_pspace.size(); j++){
			LandmarkObs obs_tmp;
			obs_tmp.id = observations_pspace[j].id;
			obs_tmp.x = observations_pspace[j].x;
			obs_tmp.y = observations_pspace[j].y;
			// std::cout << "Observations id,x,y: "<< obs_tmp.id << " " << obs_tmp.x << " " << obs_tmp.y << std::endl;
			
			LandmarkObs pred_tmp;
			for(int k = 0; k < landmark_inrange.size(); k++){
				std::cout << landmark_inrange[k].id << " ";
				if (obs_tmp.id == landmark_inrange[k].id){
					pred_tmp.id = landmark_inrange[k].id;
					pred_tmp.x = landmark_inrange[k].x;
					pred_tmp.y = landmark_inrange[k].y;
				}
			}
			std::cout << std::endl;
			std::cout << "Match Landmark id,x,y: "<< pred_tmp.id << " " << pred_tmp.x << " " << pred_tmp.y << std::endl;

			// Update the weights using gaussian 2d
			double sqrt_2picov = 1/sqrt(2.0*3.14159*std_x*std_y);
			double x_diff = ((pred_tmp.x-obs_tmp.x)*(pred_tmp.x-obs_tmp.x))/(2*std_x*std_x);
			double y_diff = ((pred_tmp.y-obs_tmp.y)*(pred_tmp.y-obs_tmp.y))/(2*std_y*std_y);
			double obs_w = sqrt_2picov * exp(-(x_diff + y_diff));
			std::cout << "sqrt_2picov: "<< sqrt_2picov << " x_diff: "<< x_diff << " y_diff: " 
			<< y_diff << "  exp: " <<  exp(-(x_diff + y_diff)) << "obs_w: "<< obs_w << std::endl;
			std::cout << "Particle's Weight Updating : " << particles[i].weight << std::endl;
			particles[i].weight = particles[i].weight*obs_w;
		}
		weights.push_back(particles[i].weight);
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::discrete_distribution<> d(weights.begin(), weights.end());
	default_random_engine gen;
	std::vector<Particle> resampling_particles;

	for(int i = 0; i < num_particles; i++){
		resampling_particles.push_back(particles[d(gen)]);
	}
	particles = resampling_particles;
	weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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
