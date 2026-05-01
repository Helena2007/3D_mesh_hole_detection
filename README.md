## 3D Mesh Hole Detection
This repository contains my research work on hole detection in 3D triangular mesh models using high curvature-based analysis.
The work is currently under submission.
The work describes a robust and efficient hole detection method that detects complex irregular hole boundaries by using high curvature analysis.

### Installation
git clone https://github.com/Helena2007/3D_mesh_hole_detection.git

cd 3D_mesh_hole_detection

sudo apt update

sudo apt install -y build-essential cmake libpcl-dev

## Dependencies
- C++ compiler with C++17 support (GCC ≥ 7 or Clang ≥ 5)
- CMake ≥ 3.10
- PCL (Point Cloud Library)

### Required PCL Modules:
- common
- io
- surface
- search
- kdtree

### Installed Automatically with PCL:
- Eigen (linear algebra)
- Boost (utility libraries)
- FLANN (nearest neighbor search)
- VTK (visualization backend)

## How to compile
   mkdir build
   cd build
   cmake ..
   make -j$(nproc)

## How to run the project
./HoleDetection input.ply
