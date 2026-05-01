# 3D Mesh Hole Detection
This repository contains my research work on 3D hole detection in triangular mesh models using high curvature-based analysis.

The work is currently under submission.

The work describes a robust and efficient hole detection method that detects complex irregular hole boundaries by using high curvature analysis.

**#Installation**
git clone https://github.com/Helena2007/3D_mesh_hole_detection.git

cd 3D_mesh_hole_detection

sudo apt update

sudo apt install -y build-essential cmake libpcl-dev

**#How to compile**

mkdir build
cd build
cmake ..
make -j$(nproc)

**#How to run the project**

./HoleDetection input.ply
