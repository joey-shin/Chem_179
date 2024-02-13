#include <iostream>
#include <cmath>

// Define the objective function (quadratic function)
double objectiveFunction(double x) {
    return x * x - 4 * x + 4;
}

// Calculate the derivative of the objective function
double derivative(double x) {
    return 2 * x - 4;
}

// Steepest Descent Algorithm
double steepestDescent(double initialGuess, double learningRate, int maxIterations, double tolerance) {
    double currentX = initialGuess;

    for (int i = 0; i < maxIterations; ++i) {
        double gradient = derivative(currentX);

        if (std::abs(gradient) < tolerance) {
            std::cout << "Converged after " << i + 1 << " iterations." << std::endl;
            break;
        }

        // Update the current value using steepest descent formula
        currentX = currentX - learningRate * gradient;

        std::cout << "Iteration " << i + 1 << ": x = " << currentX << ", f(x) = " << objectiveFunction(currentX) << std::endl;
    }

    return currentX;
}

int main() {
    // Set initial values
    double initialGuess = 0.0;
    double learningRate = 0.1;
    int maxIterations = 100;
    double tolerance = 1e-6;

    // Run the steepest descent algorithm
    double result = steepestDescent(initialGuess, learningRate, maxIterations, tolerance);

    std::cout << "Minimum found at x = " << result << ", f(x) = " << objectiveFunction(result) << std::endl;

    return 0;
}
