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
#include <limits>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    // set the number of particles
    num_particles = 100;
    cout << "Num particles set" << endl;
    
    // random number generator
    default_random_engine gen;
    
    // create normal distributions for x, y and theta based on the initial values and
    // standard deviations
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    // Initialize the particles using this distribution
    // and also set all the weights to 1
    particles.resize(num_particles);
    weights.resize(num_particles);
    for(unsigned int i = 0; i < num_particles; ++i) {
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1.0;
        weights[i] = 1.0;
    }
    cout << "Particles set" << endl;
    
    is_initialized = true;
    
    /// END OF INITIALIZATION
    

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    // random number generator
    default_random_engine gen;
    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);
    
    // for each particle predict the new position and heading
    // using the velocity and the yaw rate as inputs
    for(unsigned int i = 0; i < num_particles; ++i) {
        double x0 = particles[i].x;
        double y0 = particles[i].y;
        double theta0 = particles[i].theta;
        double xf, yf, thetaf;
        if(fabs(yaw_rate) > 0.001){
            xf = x0 + (velocity / yaw_rate)*(sin(theta0 + yaw_rate*delta_t) - sin(theta0));
            yf = y0 + (velocity / yaw_rate)*(cos(theta0) - cos(theta0 + yaw_rate*delta_t));
            thetaf = theta0 + yaw_rate*delta_t;
        }
        else
        {
            double distance = velocity * delta_t;
            xf = x0 + distance * cos(theta0);
            yf = y0 + distance * sin(theta0);
            thetaf = theta0;
        }
        
        // Apply some jitter (gaussian random noise) on the predicted
        // state for each particle using the noise characteristics
        // similar for the GPS
        // Update the particle position and heading
        particles[i].x = xf + dist_x(gen);
        particles[i].y = yf + dist_y(gen);
        particles[i].theta = thetaf + dist_theta(gen);
    }
    

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations, Particle& particle) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    
    // I believe predicted contains the landmarks in the map
    // And observations contains the measurements from the sensors of the vehicle
    
    // I also assume that these observations have already been transformed
    // to the map co-ordinate system before this function is called
 
    
    // So we loop through the observations vector of LandmarkObs
    // LandmarkObs structure is:
    // int id
    // double x (m)
    // double y (m)
    double curr_dist, best_dist;
    // Here we will also update the associations, sensex and sensey for each particle
    vector<int> associated_ids;
    vector<double> associated_x;
    vector<double> associated_y;
    // For each observation we find the closest associated landmark
    // in the predicted vector. And we assign the id based on that.
    // For each observation
    for(LandmarkObs currObs : observations) {
        // Initialize the best distance to some large value
        best_dist = std::numeric_limits<double>::max();
        int matched_id;
        double matched_x, matched_y;
        // For each prediction
        for(LandmarkObs pred : predicted) {
            // If the distance between the current observation and
            // prediction is less than the best, update the id and the distance
            curr_dist = dist(currObs.x, currObs.y, pred.x, pred.y);
            if(curr_dist < best_dist){
                best_dist = curr_dist;
                matched_id = pred.id;
                matched_x = pred.x;
                matched_y = pred.y;
            }
        }
        associated_ids.push_back(matched_id);
        associated_x.push_back(matched_x);
        associated_y.push_back(matched_y);
    }
    SetAssociations(particle, associated_ids, associated_x, associated_y);
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
    
    // This function will take the noisy observations in the vehicle
    // co-ordinate system and transform to the map co-ordinate
    // system for each particle.
    // This can be stored in vector<LandmarkObs> transformed_observations
    // For each particle we do the following steps:
    // 0. Find the landmarks in the map which are in range of the particle
    // 1. For each observation we transform from VEHICLE system
    // to MAP system using the particle's heading and position.
    // 2. We store the observation as a LandmarkObs vector called transformed_
    //observations
    // 3. We pass the vector to dataAssociation function, where for each
    //transformed observation we compare with the map and find the nearest one
    // 4. We then use this updated transformed_observations vector to update the weight for each particle.
    
    
    int particle_number = 0;
    for(Particle& particle : particles){
        
        // Loop through the landmarks
        // on the map and find which can be possible associations
        // based on the given sensor range
        vector<LandmarkObs> landmarks_in_range;
        for(Map::single_landmark_s single_landmark_temp : map_landmarks.landmark_list){
            
            if(dist(particle.x, particle.y, single_landmark_temp.x_f, single_landmark_temp.y_f) <= sensor_range){
                LandmarkObs landmark;
                landmark.x = single_landmark_temp.x_f;
                landmark.y = single_landmark_temp.y_f;
                landmark.id = single_landmark_temp.id_i;
                //cout << "Landmark " << landmark.id << " in range" << endl;
                landmarks_in_range.push_back(landmark);
            }
        }
        
        // Get the heading and position of the current particle
        double heading = particle.theta;
        double xp = particle.x;
        double yp = particle.y;
        
        // Now generate the transformed observations from the given observations
        // based on the particles position and heading
        vector<LandmarkObs> transformed_observations;
        for(LandmarkObs observation : observations){
            // Calculate the transformed observation
            // and store it in the vector
            LandmarkObs transformed_observation;
            transformed_observation.x = xp + observation.x*cos(heading) -
            observation.y*sin(heading);
            transformed_observation.y = yp + observation.x*sin(heading) +
            observation.y*cos(heading);
            //cout << "Adding transformed observation" << endl;
            transformed_observations.push_back(transformed_observation);
        }
        
        
        // Do the data association using nearest neighbor
        // This will update assign the associated ids, sensex and sensey
        // for the given particle based on the landmarks and observation.
        dataAssociation(landmarks_in_range, transformed_observations, particle);
        //cout << "Associated IDs and Positions Updated" << endl;
        
        
        // Loop through each transformed observation and update the weight
        double particle_weight = 1.0;
        double measurement_likelihood;
        double observed_x, observed_y;
        // For each observation
        for(unsigned int i = 0; i < transformed_observations.size(); ++i){
            // Find the associated landmark for this observation
            //cout << "Good till here" << endl;
            // Retrieve the observed measurement
            observed_x = transformed_observations[i].x;
            observed_y = transformed_observations[i].y;
            // Calculate the measurement likelihood
            measurement_likelihood = meas_prob(observed_x, observed_y, particle.sense_x[i], particle.sense_y[i], std_landmark[0], std_landmark[1]);
            // Update particle weight
            particle_weight *= measurement_likelihood;
        }
        //cout << "Weight updated" << endl;
        //cout << "weight = " << particle_weight << endl;
        // Assign the new weight to the particle
        particle.weight = particle_weight;
        //particles[particle_number].weight = particle_weight;
        //particle_number++;
        
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    // Update the list of weights
    for(unsigned int i = 0; i < num_particles; ++i){
        cout << "weight = " << particles[i].weight << endl;
        weights[i] = particles[i].weight;
    }
    // Create a vector to store the resampled particles
    vector<Particle> new_particles;
    random_device rd;
    std::mt19937 gen(rd());
    // Create a discrete distribution to resample the particles with
    // probability proportional to their weight
    discrete_distribution<> dist(weights.begin(), weights.end());
    for(int i = 0; i < num_particles; ++i){
        Particle new_particle = particles[dist(gen)];
        new_particles.push_back(new_particle);
    }
    particles = new_particles;
    
    
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
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
