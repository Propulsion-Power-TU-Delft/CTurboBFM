#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstddef>
#include <cstdint>
#include "types.hpp"

/**
 * @brief High-performance, self-contained K-D Tree for nearest-neighbor spatial searches in 2D and 3D.
 * 
 * Used for coarse-to-fine simulation restarts in CTurboBFM.
 * 
 * @tparam Dim Dimension of the search space (typically 2 or 3).
 */
template <size_t Dim>
class KDTree {
public:
    using Point = std::array<FloatType, Dim>;

    KDTree() = default;

    /**
     * @brief Construct KDTree from a collection of points.
     * @param points Input points. Indices in the tree correspond to indices in this vector.
     */
    explicit KDTree(const std::vector<Point>& points) {
        build(points);
    }

    /**
     * @brief Build the tree from a collection of points.
     * @param points Input points.
     */
    void build(const std::vector<Point>& points) {
        _nodes.clear();
        _root = -1;
        if (points.empty()) return;

        struct Item {
            Point pt;
            size_t origIdx;
        };

        std::vector<Item> items(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            items[i] = {points[i], i};
        }

        _nodes.reserve(points.size());
        _root = buildRecursive(items, 0, items.size(), 0);
    }

    /**
     * @brief Find the index of the nearest point to query.
     * @param query Query coordinates.
     * @return Original index of the nearest point in the points vector.
     */
    size_t findNearest(const Point& query) const {
        if (_nodes.empty() || _root < 0) {
            throw std::runtime_error("KDTree::findNearest called on empty tree");
        }
        size_t bestIdx = _nodes[_root].origIdx;
        FloatType bestDistSq = std::numeric_limits<FloatType>::max();
        searchRecursive(_root, query, bestIdx, bestDistSq);
        return bestIdx;
    }

    /**
     * @brief Batch query for multiple points, parallelized with OpenMP.
     * @param queries Vector of query coordinates.
     * @return Vector of indices of nearest points.
     */
    std::vector<size_t> findNearestBatch(const std::vector<Point>& queries) const {
        std::vector<size_t> results(queries.size());
        #pragma omp parallel for
        for (int64_t i = 0; i < static_cast<int64_t>(queries.size()); ++i) {
            results[i] = findNearest(queries[i]);
        }
        return results;
    }

    /**
     * @brief Check if tree is empty.
     */
    bool empty() const {
        return _nodes.empty();
    }

    /**
     * @brief Number of points stored in tree.
     */
    size_t size() const {
        return _nodes.size();
    }

private:
    struct Node {
        Point pt;
        size_t origIdx;
        size_t axis;
        int left;
        int right;
    };

    std::vector<Node> _nodes;
    int _root {-1};

    template <typename Item>
    int buildRecursive(std::vector<Item>& items, size_t start, size_t end, size_t depth) {
        if (start >= end) return -1;

        size_t axis = depth % Dim;
        size_t mid = start + (end - start) / 2;

        std::nth_element(
            items.begin() + start,
            items.begin() + mid,
            items.begin() + end,
            [axis](const Item& a, const Item& b) {
                return a.pt[axis] < b.pt[axis];
            }
        );

        int nodeIdx = static_cast<int>(_nodes.size());
        _nodes.push_back({items[mid].pt, items[mid].origIdx, axis, -1, -1});

        int leftChild = buildRecursive(items, start, mid, depth + 1);
        int rightChild = buildRecursive(items, mid + 1, end, depth + 1);

        _nodes[nodeIdx].left = leftChild;
        _nodes[nodeIdx].right = rightChild;

        return nodeIdx;
    }

    static FloatType distanceSq(const Point& a, const Point& b) {
        FloatType d2 = 0.0;
        for (size_t d = 0; d < Dim; ++d) {
            FloatType diff = a[d] - b[d];
            d2 += diff * diff;
        }
        return d2;
    }

    void searchRecursive(int nodeIdx, const Point& query, size_t& bestIdx, FloatType& bestDistSq) const {
        if (nodeIdx < 0) return;

        const Node& node = _nodes[nodeIdx];
        FloatType d2 = distanceSq(query, node.pt);
        if (d2 < bestDistSq) {
            bestDistSq = d2;
            bestIdx = node.origIdx;
        }

        size_t axis = node.axis;
        FloatType diff = query[axis] - node.pt[axis];
        int nearChild = (diff <= 0.0) ? node.left : node.right;
        int farChild  = (diff <= 0.0) ? node.right : node.left;

        searchRecursive(nearChild, query, bestIdx, bestDistSq);

        if (diff * diff < bestDistSq) {
            searchRecursive(farChild, query, bestIdx, bestDistSq);
        }
    }
};

/** @brief Helper to convert Vector3D to 3D KDTree Point */
inline std::array<FloatType, 3> toPoint3D(const Vector3D& v) {
    return {v.x(), v.y(), v.z()};
}

/** @brief Helper to convert Vector3D to 2D (x, y) KDTree Point */
inline std::array<FloatType, 2> toPointXY(const Vector3D& v) {
    return {v.x(), v.y()};
}

/** @brief Helper to convert Vector3D to 2D cylindrical (x, r) KDTree Point */
inline std::array<FloatType, 2> toPointXR(const Vector3D& v) {
    FloatType r = std::sqrt(v.y() * v.y() + v.z() * v.z());
    return {v.x(), r};
}
