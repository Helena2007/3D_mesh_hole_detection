#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/point_types.h>

using PointT = pcl::PointXYZ;
using PointCloudT = pcl::PointCloud<PointT>;

inline pcl::PolygonMesh loadMesh(const std::string& filename) {
    pcl::PolygonMesh mesh;

    if (pcl::io::loadPolygonFilePLY(filename, mesh) == -1) {
        std::cerr << "Error: Could not load mesh file "
                  << filename << std::endl;
        return mesh;
    }

   // std::cout << "Loaded mesh successfully: " << filename << "\n"
     //         << "  - Faces    : " << mesh.polygons.size() << "\n"
       //       << "  - Vertices : "
         //     << mesh.cloud.width * mesh.cloud.height
           //   << std::endl;

    return mesh;
}

inline void saveHoleBoundariesToTXT(
    const std::vector<std::pair<int, int>>& holeEdges,
    const pcl::PolygonMesh& mesh,
    const std::string& filename) {

    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file "
                  << filename << std::endl;
        return;
    }

    // Convert mesh cloud ONCE (stack allocation, faster)
    PointCloudT cloud;
    pcl::fromPCLPointCloud2(mesh.cloud, cloud);

    const auto& points = cloud.points;

    for (const auto& edge : holeEdges) {
        const PointT& p1 = points[edge.first];
        const PointT& p2 = points[edge.second];

        file << p1.x << ' ' << p1.y << ' ' << p1.z << "  "
             << p2.x << ' ' << p2.y << ' ' << p2.z << '\n';
    }

    file.close();
    std::cout << "Hole boundary coordinates saved to "
              << filename << std::endl;
}

