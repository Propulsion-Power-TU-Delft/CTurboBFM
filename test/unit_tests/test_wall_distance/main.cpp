#include <iostream>
#include "mesh.hpp"
#include "config.hpp"
#include "gtest/gtest.h"
#include "solver_base.hpp"
#include "solver_euler.hpp"


TEST(WallDistanceTest, Along_X) {
    Config config("input_walls_along_x.ini");
    Mesh mesh(config);

    std::unique_ptr<SolverBase> solver;
    solver = std::make_unique<SolverEuler>(config, mesh);
    
    Matrix3D<Vector3D> vertices = solver->getVertices();
    Matrix3D<FloatType> wallDistance = solver->getWallDistance();

    size_t ni = vertices.sizeI();
    size_t nj = vertices.sizeJ();
    size_t nk = vertices.sizeK();

    std::cout << "Testing wall distance along x-direction..." << std::endl; 
    std::cout << "Mesh dimensions: " << ni << " x " << nj << " x " << nk << std::endl;

    FloatType LX = 9.0; // from the mesh file
    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                Vector3D vertex = vertices(i,j,k);
                FloatType expectedDistance = std::min(std::abs(vertex.x()), LX - std::abs(vertex.x()));
                ASSERT_DOUBLE_EQ(wallDistance(i,j,k), expectedDistance);
                std::cout << "Vertex (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ") - "
                          << "Computed distance: " << wallDistance(i,j,k) << ", "
                          << "Expected distance: " << expectedDistance << std::endl;
            }
        }
    }
}



TEST(WallDistanceTest, Along_Y) {
    Config config("input_walls_along_y.ini");
    Mesh mesh(config);

    std::unique_ptr<SolverBase> solver;
    solver = std::make_unique<SolverEuler>(config, mesh);
    
    Matrix3D<Vector3D> vertices = solver->getVertices();
    Matrix3D<FloatType> wallDistance = solver->getWallDistance();

    size_t ni = vertices.sizeI();
    size_t nj = vertices.sizeJ();
    size_t nk = vertices.sizeK();

    std::cout << "Testing wall distance along y-direction..." << std::endl; 
    std::cout << "Mesh dimensions: " << ni << " x " << nj << " x " << nk << std::endl;

    FloatType LY = 7.0; // from the mesh file
    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                Vector3D vertex = vertices(i,j,k);
                FloatType expectedDistance = std::min(std::abs(vertex.y()), LY - std::abs(vertex.y()));
                ASSERT_DOUBLE_EQ(wallDistance(i,j,k), expectedDistance);
                std::cout << "Vertex (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ") - "
                          << "Computed distance: " << wallDistance(i,j,k) << ", "
                          << "Expected distance: " << expectedDistance << std::endl;
            }
        }
    }
}



TEST(WallDistanceTest, Along_Z) {
    Config config("input_walls_along_z.ini");
    Mesh mesh(config);

    std::unique_ptr<SolverBase> solver;
    solver = std::make_unique<SolverEuler>(config, mesh);
    
    Matrix3D<Vector3D> vertices = solver->getVertices();
    Matrix3D<FloatType> wallDistance = solver->getWallDistance();

    size_t ni = vertices.sizeI();
    size_t nj = vertices.sizeJ();
    size_t nk = vertices.sizeK();

    std::cout << "Testing wall distance along z-direction..." << std::endl; 
    std::cout << "Mesh dimensions: " << ni << " x " << nj << " x " << nk << std::endl;

    FloatType LZ = 5.0; // from the mesh file
    for (size_t i = 0; i < ni; ++i) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t k = 0; k < nk; ++k) {
                Vector3D vertex = vertices(i,j,k);
                FloatType expectedDistance = std::min(std::abs(vertex.z()), LZ - std::abs(vertex.z()));
                ASSERT_DOUBLE_EQ(wallDistance(i,j,k), expectedDistance);
                std::cout << "Vertex (" << vertex.x() << ", " << vertex.y() << ", " << vertex.z() << ") - "
                          << "Computed distance: " << wallDistance(i,j,k) << ", "
                          << "Expected distance: " << expectedDistance << std::endl;
            }
        }
    }
}
