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
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  std::default_random_engine e; //引擎
  std::normal_distribution<double> x_gauss(x, std[0]); //均值, 方差
  std::normal_distribution<double> y_gauss(y, std[1]); //均值, 方差
  std::normal_distribution<double> theta_gauss(theta, std[2]); //均值, 方差
  for(int i=0;i<num_particles;++i)
  {
    Particle particle_temp;
    particle_temp.id = i;
    particle_temp.x = x_gauss(e);
    particle_temp.y = y_gauss(e);
    particle_temp.theta = theta_gauss(e); 
    particle_temp.weight = 1.0;
    particles.push_back(particle_temp);
  }
  is_initialized = true;
  std::cout << "init finished "  << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine e; //引擎
  std::normal_distribution<double> x_gauss(0.0, std_pos[0]); //均值, 方差
  std::normal_distribution<double> y_gauss(0.0, std_pos[1]); //均值, 方差
  std::normal_distribution<double> theta_gauss(0.0, std_pos[2]); //均值, 方差
  for(int i=0;i<particles.size();++i)
  {
    if(fabs(yaw_rate)<0.00001)
    {
      particles[i].x += velocity*delta_t*cos(particles[i].theta)+x_gauss(e);
      particles[i].y += velocity*delta_t*sin(particles[i].theta)+y_gauss(e);
      particles[i].theta += yaw_rate*delta_t + theta_gauss(e);
      std::cout<<" yaw=====================0"<<std::endl;
    }
    else
    {
      particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;
      particles[i].x +=x_gauss(e);
      particles[i].y +=y_gauss(e);
      particles[i].theta +=theta_gauss(e);
    }
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  vector<LandmarkObs> observations_temp;
   for(int i=0;i<predicted.size();++i)
   {
     double error_dist=100000;
     LandmarkObs LandmarkObs_min;
     for(int j=0;j<observations.size();++j)
     {
       double dist_temp = dist(observations[j].x,observations[j].y,predicted[i].x,predicted[i].y);
       if (dist_temp<error_dist)
       {
          error_dist = dist_temp;
          LandmarkObs_min = observations[j];
       }
     }
     observations_temp.push_back(LandmarkObs_min);
   }
   observations = observations_temp;
  //  observations.clear();
  //  for(int i=0;i<observations_temp.size();++i)
  //  {
  //    observations.push_back(observations_temp[i]);
  //  }
}

double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) 
{
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double weight_sum = 0;
  for(int i=0;i<particles.size();++i)
  {
    vector<LandmarkObs> observations_map;
    for(int j=0;j<observations.size();++j)
    {
      LandmarkObs observations_temp;
      observations_temp.x = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
      observations_temp.y = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;
      observations_temp.id = observations[j].id;
      observations_map.push_back(observations_temp);    
    }
    vector<LandmarkObs> observations_in_map;
    for(int j=0;j<map_landmarks.landmark_list.size();++j)
    {
      double dist_temp = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f, particles[i].x, particles[i].y); 
      if(dist_temp < sensor_range)
      {
        LandmarkObs observations_temp;
        observations_temp.id = map_landmarks.landmark_list[j].id_i;
        observations_temp.x = map_landmarks.landmark_list[j].x_f;
        observations_temp.y = map_landmarks.landmark_list[j].y_f;
        observations_in_map.push_back(observations_temp);
      }
    }
    dataAssociation(observations_map,observations_in_map);
    for(int j=0;j<observations_map.size();++j)
    {
      double weight_temp = multiv_prob(std_landmark[0],std_landmark[1],observations_map[j].x,observations_map[j].y,observations_in_map[j].x,
                                       observations_in_map[j].y);
      particles[i].weight *= weight_temp;                               
    }
    weight_sum += particles[i].weight;
  }
  weights.clear();
  for(int i=0;i<particles.size();++i)
  {
    particles[i].weight = particles[i].weight/weight_sum;
    weights.push_back(particles[i].weight);
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    std::vector<Particle> particles_resample;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(std::begin(weights), std::end(weights));
    int map[1000];
    for(int n=0; n<particles.size(); ++n) 
    {
      map[n]=0;
    }
    for(int n=0; n<particles.size(); ++n) 
    {
        int index = d(gen);
        map[index]++;
        Particle sample = particles[index];
        particles_resample.push_back(sample);
    }
    particles = particles_resample;
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