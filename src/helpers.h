#include <cmath>
#include <numeric>
#include <vector>



/**
 * Calculate the average value of a vector
 */
float vecAvg(std::vector<float> &vec) {
    float numerator = std::accumulate(vec.begin(), vec.end(), 0.0);
    float denominator = (float) vec.size();
    return numerator/denominator;
}

/**
 * Calculate the difference between two angles in the range (-pi, pi).
 */
float angleDif(float angle1,float angle2) {
    return atan2(sin(angle1-angle2), cos(angle1-angle2));
}

/**
 * Calculate delta R for two points in eta phi space.
 */
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return sqrt( pow(eta2 - eta1, 2) + pow(angleDif(phi2, phi1), 2));
}
