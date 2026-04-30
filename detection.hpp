#pragma once
#include "util.hpp"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <queue>
#include <cmath>
#include <pcl/PolygonMesh.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using Edge = std::pair<int, int>;
using EdgeKey = uint64_t;

inline EdgeKey makeEdgeKey(int v1, int v2) {
    if (v1 > v2) std::swap(v1, v2);
    return (EdgeKey(v1) << 32) | EdgeKey(v2);
}

//Extract Open edges
inline std::vector<Edge>
extractMeshBoundaries(const pcl::PolygonMesh& mesh) {

    const size_t F = mesh.polygons.size();

#ifdef _OPENMP
    int T = omp_get_max_threads();
    std::vector<std::unordered_map<EdgeKey, uint32_t>> localMaps(T);
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        localMaps[tid].reserve(F * 2 / T);

#pragma omp for schedule(static)
        for (size_t f = 0; f < F; ++f) {
            const auto& poly = mesh.polygons[f];
            const size_t n = poly.vertices.size();
            for (size_t i = 0; i < n; ++i)
                ++localMaps[tid][makeEdgeKey(poly.vertices[i],
                                             poly.vertices[(i + 1) % n])];
        }
    }

    std::unordered_map<EdgeKey, uint32_t> edgeCount;
    edgeCount.reserve(F * 3);
    for (auto& m : localMaps)
        for (auto& [k, v] : m)
            edgeCount[k] += v;
#else
    std::unordered_map<EdgeKey, uint32_t> edgeCount;
    edgeCount.reserve(F * 3);
    for (const auto& poly : mesh.polygons)
        for (size_t i = 0; i < poly.vertices.size(); ++i)
            ++edgeCount[makeEdgeKey(poly.vertices[i],
                                    poly.vertices[(i + 1) % poly.vertices.size()])];
#endif

    std::vector<Edge> boundaryEdges;
    boundaryEdges.reserve(edgeCount.size());
    for (const auto& [k, c] : edgeCount)
        if (c == 1)
            boundaryEdges.emplace_back(int(k >> 32), int(k & 0xffffffff));

    return boundaryEdges;
}

inline std::vector<Edge>
extractOuterBoundary(const std::vector<Edge>& boundaries) {
    if (boundaries.empty()) return {};

    std::unordered_map<int, std::vector<int>> adj;
    for (auto& [a, b] : boundaries) {
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    std::unordered_set<int> visited;
    std::vector<Edge> largest;
    size_t maxSize = 0;

    for (auto& [start, _] : adj) {
        if (visited.count(start)) continue;

        std::vector<Edge> comp;
        std::queue<int> q;
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int nb : adj[v]) {
                comp.emplace_back(v, nb);
                if (!visited.count(nb)) {
                    visited.insert(nb);
                    q.push(nb);
                }
            }
        }

        if (comp.size() > maxSize) {
            maxSize = comp.size();
            largest.swap(comp);
        }
    }
    return largest;
}

inline std::unordered_set<int>
extractBoundaryVerticesOneRing(const pcl::PolygonMesh& mesh) {

    const size_t V = mesh.cloud.width;
    std::vector<int> faceCnt(V, 0), edgeCnt(V, 0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < mesh.polygons.size(); ++i) {
        const auto& poly = mesh.polygons[i];
        const size_t n = poly.vertices.size();
        for (int v : poly.vertices)
#pragma omp atomic
            ++faceCnt[v];
        for (size_t j = 0; j < n; ++j) {
#pragma omp atomic
            ++edgeCnt[poly.vertices[j]];
#pragma omp atomic
            ++edgeCnt[poly.vertices[(j + 1) % n]];
        }
    }

    std::unordered_set<int> boundaryVertices;
    for (size_t v = 0; v < V; ++v)
        if (edgeCnt[v] != faceCnt[v])
            boundaryVertices.insert(int(v));

    return boundaryVertices;
}

inline std::vector<Edge>
extractBoundaryLoops(const std::vector<Edge>& edges,
                     const std::unordered_set<int>& bVerts) {

    std::unordered_map<int, std::vector<int>> graph;
    for (auto& [a, b] : edges)
        if (bVerts.count(a) && bVerts.count(b)) {
            graph[a].push_back(b);
            graph[b].push_back(a);
        }

    std::unordered_set<EdgeKey> uniqueEdges;
    std::unordered_set<int> visited;

    for (auto& [start, _] : graph) {
        if (visited.count(start)) continue;

        std::queue<int> q;
        q.push(start);
        visited.insert(start);

        while (!q.empty()) {
            int v = q.front(); q.pop();
            for (int nb : graph[v]) {
                uniqueEdges.insert(makeEdgeKey(v, nb));
                if (!visited.count(nb)) {
                    visited.insert(nb);
                    q.push(nb);
                }
            }
        }
    }

    std::vector<Edge> loopEdges;
    for (auto k : uniqueEdges)
        loopEdges.emplace_back(int(k >> 32), int(k & 0xffffffff));
    return loopEdges;
}

//Computation of Gaussian curvature
inline std::unordered_map<int, double>
computeMeshCurvature(const pcl::PolygonMesh& mesh, double& adaptiveThreshold) {

    pcl::PointCloud<PointT> cloud;
    pcl::fromPCLPointCloud2(mesh.cloud, cloud);

    const size_t V = cloud.size();
    std::vector<double> angleSum(V, 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < mesh.polygons.size(); ++i) {
        const auto& f = mesh.polygons[i];
        if (f.vertices.size() != 3) continue;

        auto ang = [&](int a, int b, int c) {
            Eigen::Vector3d u(cloud[a].x - cloud[b].x,
                              cloud[a].y - cloud[b].y,
                              cloud[a].z - cloud[b].z);
            Eigen::Vector3d v(cloud[c].x - cloud[b].x,
                              cloud[c].y - cloud[b].y,
                              cloud[c].z - cloud[b].z);
            double cs = u.dot(v) / (u.norm() * v.norm());
            return std::acos(std::clamp(cs, -1.0, 1.0));
        };

#pragma omp atomic
        angleSum[f.vertices[0]] += ang(f.vertices[1], f.vertices[0], f.vertices[2]);
#pragma omp atomic
        angleSum[f.vertices[1]] += ang(f.vertices[0], f.vertices[1], f.vertices[2]);
#pragma omp atomic
        angleSum[f.vertices[2]] += ang(f.vertices[0], f.vertices[2], f.vertices[1]);
    }

    std::unordered_map<int, double> K;
    std::vector<double> vals;
    for (size_t v = 0; v < V; ++v) {
        double k = std::abs(2.0 * M_PI - angleSum[v]);
        if (k > 0.0) {
            K[int(v)] = k;
            vals.push_back(k);
        }
    }

    double mean = std::accumulate(vals.begin(), vals.end(), 0.0) / vals.size();
    double var = 0.0;
    for (double v : vals) var += (v - mean) * (v - mean);

    adaptiveThreshold = mean + 1.2 * std::sqrt(var / vals.size());
    return K;
}

//Computation of high curvature edges
inline std::vector<Edge>
extractHighCurvatureEdges(const pcl::PolygonMesh& mesh) {

    double T = 0.0;
    auto K = computeMeshCurvature(mesh, T);

#ifdef _OPENMP
    int Tn = omp_get_max_threads();
    std::vector<std::vector<Edge>> local(Tn);

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
#pragma omp for schedule(static)
        for (size_t i = 0; i < mesh.polygons.size(); ++i) {
            const auto& f = mesh.polygons[i];
            for (size_t j = 0; j < f.vertices.size(); ++j) {
                int a = f.vertices[j];
                int b = f.vertices[(j + 1) % f.vertices.size()];
                if ((K[a] + K[b]) * 0.5 > T)
                    local[tid].emplace_back(a, b);
            }
        }
    }

    std::vector<Edge> edges;
    for (auto& v : local)
        edges.insert(edges.end(), v.begin(), v.end());
#else
    std::vector<Edge> edges;
    for (const auto& f : mesh.polygons)
        for (size_t j = 0; j < f.vertices.size(); ++j) {
            int a = f.vertices[j];
            int b = f.vertices[(j + 1) % f.vertices.size()];
            if ((K[a] + K[b]) * 0.5 > T)
                edges.emplace_back(a, b);
        }
#endif

    return edges;
}

inline std::vector<Edge>
findHolesFromMesh(const pcl::PolygonMesh& mesh, int minClusterSize) {

    const auto boundaries = extractMeshBoundaries(mesh);
    //cout << "Extracted " << boundaries.size()
        // << " boundary edges." << endl;
    if (boundaries.empty()) return {};

    const auto bVerts = extractBoundaryVerticesOneRing(mesh);
    const auto loops  = extractBoundaryLoops(boundaries, bVerts);
    if (loops.empty()) return {};

    const auto highCurv = extractHighCurvatureEdges(mesh);
    cout << "Extracted " << highCurv.size()
         << " high-curvature edges." << endl;

    if (highCurv.empty()) {
        if (loops.size() < static_cast<size_t>(minClusterSize))
            return {};
        return loops;
    }

    std::unordered_set<EdgeKey> loopSet;
    loopSet.reserve(loops.size());
    for (const auto& e : loops)
        loopSet.insert(makeEdgeKey(e.first, e.second));

    std::vector<Edge> holeEdges;

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<Edge> local;
#pragma omp for nowait schedule(static)
        for (size_t i = 0; i < highCurv.size(); ++i) {
            const auto& e = highCurv[i];
            if (loopSet.count(makeEdgeKey(e.first, e.second)))
                local.push_back(e);
        }
#pragma omp critical
        holeEdges.insert(holeEdges.end(), local.begin(), local.end());
    }
#else
    for (const auto& e : highCurv)
        if (loopSet.count(makeEdgeKey(e.first, e.second)))
            holeEdges.push_back(e);
#endif

    if (holeEdges.size() < static_cast<size_t>(minClusterSize))
        return {};

    return holeEdges;
}

