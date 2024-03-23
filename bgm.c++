# author: yimin
# 
# update: as of 23rd March 2024

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

// Define the parameters of our BGM model
const double initialForwardRate = 0.05;
const double initialVolatility = 0.1;
const double kappa = 0.1; 
// Mean reversion speed
const double theta = 0.05; 
// Long-term mean of the forward rate
const double volatilitySigma = 0.2; 
// Volatility of the forward rate
const double maturity = 1.0; 
// Maturity of the interest rate derivative
const int numPaths = 10000; 
// Number of simulation paths
const int numSteps = 100; 
// Number of time steps

// Function to generate random numbers from a standard normal distribution
double generateRandomNumber() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::normal_distribution<> dist(0, 1);
    return dist(gen);
}

int main() {
    std::vector<double> forwardRates(numSteps + 1);
    std::vector<double> volatilities(numSteps + 1);

    // Initialize initial values
    forwardRates[0] = initialForwardRate;
    volatilities[0] = initialVolatility;

    // Simulate paths
    for (int i = 1; i <= numSteps; ++i) {
        double dt = maturity / numSteps;
        double dW1 = generateRandomNumber() * sqrt(dt);
        double dW2 = generateRandomNumber() * sqrt(dt);
        forwardRates[i] = forwardRates[i - 1] + kappa * (theta - forwardRates[i - 1]) * dt + volatilities[i - 1] * sqrt(dt) * dW1;
        volatilities[i] = volatilities[i - 1] + volatilitySigma * volatilities[i - 1] * sqrt(dt) * dW2;
    }

    // Perform derivative pricing or other analysis using simulated paths

    return 0;
}








