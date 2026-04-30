#include "util.hpp"
#include "detection.hpp"
#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;

#define RED   "\033[31m"
#define RESET "\033[0m"

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: ./hole_detection <input_mesh.ply>" << endl;
        return -1;
    }

    auto start_time = high_resolution_clock::now();

    pcl::PolygonMesh mesh = loadMesh(argv[1]);
    auto boundaryEdges = extractMeshBoundaries(mesh);
    cout << "Extracted " << boundaryEdges.size()
         << " boundary edges." << endl;

    auto outerBoundary = extractOuterBoundary(boundaryEdges);


    const int minClusterSize = 5;
    auto holeEdges = findHolesFromMesh(mesh, minClusterSize);

    if (holeEdges.empty()) {
        cout << RED << "No hole boundaries detected"
             << RESET << endl;
    } else {
        cout << "Detected " << holeEdges.size()
             << " hole boundary edges." << endl;

        saveHoleBoundariesToTXT(
            holeEdges,
            mesh,
            "final_hole_boundary.txt"
        );
    }

    auto end_time = high_resolution_clock::now();
    auto duration =
        duration_cast<milliseconds>(end_time - start_time);

    cout << "Total computation time: "
         << duration.count() << " ms" << endl;

    return 0;
}



