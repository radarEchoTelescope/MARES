#include <iostream>

#include "extern/c8_ray_tracer/include/c8_tracer/transcribed/RayTracer.hpp"
#include "extern/c8_ray_tracer/include/c8_tracer/environment.hpp"

using namespace c8_tracer;

int main() {
    // n(z) = n_deep - delta_n * exp((z - z0) / length_scale)
    CartesianSingleExponentialEnvironment env(
        1.78, 0.423, 77.0,      // n_deep, delta_n, length_scale
        Vec3(0, 0, 1),          // axis (points toward zenith)
        Vec3(0, 0, 0)           // reference point z0
    );

    RayTracer2D ray_tracer(Vec3(0, 0, 1), 1e-4, 1.0);  // axis, minStep, maxStep

    Point start(0, 0, -50);
    Point end(300, 0, -200);

    std::vector<SignalPath> paths = ray_tracer.GetSignalPathsBrent(start, end, env);

    for (size_t i = 0; i < paths.size(); ++i) {
        std::cout << "Solution " << i << ": " << paths[i].getNSegments()
                  << " segments, end=" << paths[i].getEnd().to_string() << "\n";
        for (auto const &pt : paths[i]) {
            // pt.x, pt.y, pt.z — the actual path points
        }
    }
}
