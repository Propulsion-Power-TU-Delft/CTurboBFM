#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <cmath>
#include "kdtree.hpp"
#include "gtest/gtest.h"

TEST(KDTreeTest, TestBasic3D) {
    std::vector<std::array<FloatType, 3>> points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 1.0}
    };

    KDTree<3> tree(points);
    EXPECT_EQ(tree.size(), 5);

    // Exact matches
    EXPECT_EQ(tree.findNearest({0.0, 0.0, 0.0}), 0);
    EXPECT_EQ(tree.findNearest({1.0, 0.0, 0.0}), 1);
    EXPECT_EQ(tree.findNearest({0.0, 1.0, 0.0}), 2);
    EXPECT_EQ(tree.findNearest({0.0, 0.0, 1.0}), 3);
    EXPECT_EQ(tree.findNearest({1.0, 1.0, 1.0}), 4);

    // Nearby points
    EXPECT_EQ(tree.findNearest({0.1, 0.1, 0.1}), 0);
    EXPECT_EQ(tree.findNearest({0.9, 0.9, 0.8}), 4);
    EXPECT_EQ(tree.findNearest({0.8, 0.1, 0.0}), 1);
}

TEST(KDTreeTest, TestBasic2D) {
    std::vector<std::array<FloatType, 2>> points = {
        {0.0, 0.0},
        {2.0, 0.0},
        {0.0, 2.0},
        {2.0, 2.0}
    };

    KDTree<2> tree(points);
    EXPECT_EQ(tree.size(), 4);

    EXPECT_EQ(tree.findNearest({0.1, 0.2}), 0);
    EXPECT_EQ(tree.findNearest({1.9, 0.1}), 1);
    EXPECT_EQ(tree.findNearest({0.1, 1.9}), 2);
    EXPECT_EQ(tree.findNearest({1.8, 1.8}), 3);
}

TEST(KDTreeTest, TestAgainstBruteForce3D) {
    // Generate 500 pseudo-random points
    std::mt19937 rng(42);
    std::uniform_real_distribution<FloatType> dist(-10.0, 10.0);

    const size_t N = 500;
    std::vector<std::array<FloatType, 3>> points(N);
    for (size_t i = 0; i < N; ++i) {
        points[i] = {dist(rng), dist(rng), dist(rng)};
    }

    KDTree<3> tree(points);

    // Query 100 random points and verify KDTree matches brute force
    const size_t numQueries = 100;
    for (size_t q = 0; q < numQueries; ++q) {
        std::array<FloatType, 3> query = {dist(rng), dist(rng), dist(rng)};

        // Brute force search
        FloatType minD2 = std::numeric_limits<FloatType>::max();
        size_t minIdx = 0;
        for (size_t i = 0; i < N; ++i) {
            FloatType dx = query[0] - points[i][0];
            FloatType dy = query[1] - points[i][1];
            FloatType dz = query[2] - points[i][2];
            FloatType d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < minD2) {
                minD2 = d2;
                minIdx = i;
            }
        }

        size_t kdIdx = tree.findNearest(query);

        FloatType dx_kd = query[0] - points[kdIdx][0];
        FloatType dy_kd = query[1] - points[kdIdx][1];
        FloatType dz_kd = query[2] - points[kdIdx][2];
        FloatType kdD2 = dx_kd * dx_kd + dy_kd * dy_kd + dz_kd * dz_kd;

        // Distances must match exactly
        EXPECT_NEAR(kdD2, minD2, 1e-12);
        (void)minIdx;
    }
}

TEST(KDTreeTest, TestBatchQuery) {
    std::vector<std::array<FloatType, 3>> points = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    };

    KDTree<3> tree(points);

    std::vector<std::array<FloatType, 3>> queries = {
        {0.1, 0.0, 0.0},
        {0.9, 0.0, 0.0},
        {0.0, 0.8, 0.0}
    };

    auto batchResults = tree.findNearestBatch(queries);
    ASSERT_EQ(batchResults.size(), 3);
    EXPECT_EQ(batchResults[0], 0);
    EXPECT_EQ(batchResults[1], 1);
    EXPECT_EQ(batchResults[2], 2);
}

TEST(KDTreeTest, TestEdgeCases) {
    // Single point
    std::vector<std::array<FloatType, 3>> singlePoint = {{5.0, -3.0, 2.0}};
    KDTree<3> treeSingle(singlePoint);
    EXPECT_EQ(treeSingle.findNearest({100.0, 0.0, -50.0}), 0);

    // Empty tree
    KDTree<3> treeEmpty;
    EXPECT_THROW(treeEmpty.findNearest({0.0, 0.0, 0.0}), std::runtime_error);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
