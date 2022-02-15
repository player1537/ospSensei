// Copyright 2014 Mark Harris (https://github.com/harrism)
// SPDX-License-Identifier: Apache-2.0
//
// Modifications (2022) by Tanner Hobson (https://github.com/player1537)

// stdlib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// SENSEI
#include <ConfigurableAnalysis.h>
#include <VTKDataAdaptor.h>

// VTK
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#define SOFTENING 1e-9f


typedef struct {
  float x, y, z;
} Point;


typedef struct {
  Point *pos;
  Point *vel;
} Bodies;

void randomize(Point *p, size_t n) {
  for (size_t i = 0; i < n; i++) {
    p[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    p[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    p[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void bodyForce(const Bodies &bodies, float dt, size_t n) {
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < n; i++) { 
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;

    for (size_t j = 0; j < n; j++) {
      float dx = bodies.pos[j].x - bodies.pos[i].x;
      float dy = bodies.pos[j].y - bodies.pos[i].y;
      float dz = bodies.pos[j].z - bodies.pos[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3;
      Fy += dy * invDist3;
      Fz += dz * invDist3;
    }

    bodies.vel[i].x += dt * Fx;
    bodies.vel[i].y += dt * Fy;
    bodies.vel[i].z += dt * Fz;
  }
}

// class NBodyDataAdaptor : public sensei::DataAdaptor {
// public:  // VTK API
//   static NBodyDataAdaptor *New();
//   senseiTypeMacro(NBodyDataAdaptor, sensei::DataAdaptor);

// protected:  // C++ API
//   NBodyDataAdaptor();
//   ~NBodyDataAdaptor() = default;
//   NBodyDataAdaptor(const NBodyDataAdaptor &) = delete;
//   NBodyDataAdaptor &operator=(const NBodyDataAdaptor &) = delete;

// public:  // NBody API
//   int Initialize(size_t nBodies, size_t nIters, float dt);

//   int SetBodies(Body *bodies);

// private:  // NBody member variables
//   size_t m_nBodies;
//   size_t m_nIters;
//   float m_dt;
//   Body *m_bodies{nullptr};

// public:  // SENSEI API

//   int GetNumberOfMeshes(unsigned int &numMeshes) override;

//   int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;

//   int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override;

//   int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;

//   int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;

//   int ReleaseData() override;
// };


// senseiNewMacro(NBodyDataAdaptor);

// NBodyDataAdaptor::NBodyDataAdaptor()
// {
// }


// int NBodyDataAdaptor::Initialize(int nBodies, int nIters, float dt) {
//   m_nBodies = nBodies;
//   m_nIters = nIters;
//   m_dt = dt;

//   return 0;
// }

// int NBodyDataAdaptor::SetBodies(Body *bodies) {
//   m_bodies = bodies;

//   return 0;
// }


// int NBodyDataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes) {
//   numMeshes = 1;

//   return 0;
// }

// int NBodyDataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) {
//   if (id != 0) {
//     SENSEI_ERROR("invalid mesh id " << id) // no semicolon
//     return -1;
//   }

//   metadata->MeshName = "mesh";
//   metadata->MeshType = VTK_POLY_DATA;
//   metadata->BlockType = VTK_POINT_SET;

//   metadata->NumGhostCells = 0;
//   metadata->NumGhostNodes = 0;

//   metadata->NumArrays = 1;
//   metadata->ArrayName.emplace_back("nbody");
//   metadata->ArrayCentering.emplace_back(vtkDataObject::POINT);
//   metadata->ArrayType.emplace_back(VTK_FLOAT);
//   metadata->ArrayComponents.emplace_back(6);

//   metadata->NumBlocks = 1;
//   metadata->NumBlocksLocal.emplace_back(1);

//   if (metadata->Flags.BlockDecompSet()) {
//     SENSEI_ERROR("BlockDecompSet support not implemented")// no semicolon
//     return -1;
//   }

//   if (metadata->Flags.BlockSizeSet()) {
//     SENSEI_ERROR("BlockSizeSet support not implemented")// no semicolon
//     return -1;
//   }

//   if (metadata->Flags.BlockExtentsSet()) {
//     SENSEI_ERROR("BlockExtentsSet support not implemented")// no semicolon
//     return -1;
//   }

//   if (metadata->Flags.BlockBoundsSet()) {
//     SENSEI_ERROR("BlockBoundsSet support not implemented")// no semicolon
//     return -1;
//   }

//   return 0;
// }

// int NBodyDataAdaptor::GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) {
//   if (meshName != "mesh") {
//     SENSEI_ERROR("only the mesh 'mesh' is supported")// no semicolon
//     return -1;
//   }

//   sensei::MeshMetadataFlags flags{0};
//   sensei::MeshMetadataPtr mmd = sensei::MeshMetadata::New(flags);
//   GetMeshMetadata(0, mmd);

//   vtkSmartPointer<vtkPolyData> polyData;
//   polyData = vtkSmartPointer<vtkPolyData>::New();
//   polyData->Initialize();

//   if (!structureOnly) {
//     vtkSmartPointer<vtkFloatArray> floatArrayPos{};
//     floatArrayPos = vtkSmartPointer<vtkFloatArray>::New();
//     floatArrayPos->Initialize();
//     floatArrayPos->SetName("position");
//     floatArrayPos->SetNumberOfComponents(3);
//     floatArrayPos->SetArray((float *)m_bodies->pos, 3 * m_nBodies, /*save=*/1);

//     vtkSmartPoint<vtkPoints> points{};
//     points = vtkSmartPointer<vtkPoints>::New();
//     points->Initialize();
//     points->SetData(floatArrayPos);

//     polyData->SetPoints(points);
//   }

//   mesh = polyData;

//   return 0;
// }

// int NBodyDataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) {
//   if (meshName != "mesh") {
//     SENSEI_ERROR("only the mesh 'mesh' is supported")// no semicolon
//     return -1;
//   }

//   // association == vtkDataObject::{Point,Cell,Field}
//   if (association != vtkDataObject::POINT) {
//     SENSEI_ERROR("only the point association is supported")// no semicolon
//     return -1;
//   }

//   if (arrayName != "velocity") {
//     SENSEI_ERROR("only the array 'velocity' is supported")// no semicolon
//     return -1;
//   }

//   vtkSmartPointer<vtkFloatArray> floatArrayVel{};
//   floatArrayVel = vtkSmartPointer<vtkFloatArray>::New();
//   floatArrayVel->Initialize();
//   floatArrayVel->SetName("velocity");
//   floatArrayVel->SetNumberOfComponents(3);
//   floatArrayVel->SetArray((float *)m_bodies->vel, 3 * m_nBodies, /*save=*/1);

//   vtkPolyData *polyData = dynamic_cast<vtkPolyData *>(mesh);
//   polyData->GetPointData()->AddArray(floatArrayVel);

//   return 0;
// }

// int NBodyDataAdaptor::AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) {
//   return 0;
// }

// int NBodyDataAdaptor::ReleaseData() {
//   if (m_bodies) {
//     m_bodies = nullptr;
//   }

//   return 0;
// }



int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  size_t nBodies = 30000;
  size_t nIters = 10;
  float dt = 0.01f; // time step
  std::string configFilename{"nbody.xml"};
  if (argc > 1 && argv[1][0] != '\0') nBodies = strtoul(argv[1], NULL, 10);
  if (argc > 2 && argv[2][0] != '\0') nIters = strtoul(argv[2], NULL, 10);
  if (argc > 3 && argv[3][0] != '\0') dt = atof(argv[3]);
  if (argc > 4 && argv[4][0] != '\0') configFilename = argv[4];

  Bodies bodies;
  bodies.pos = (Point *)malloc(nBodies * sizeof(*bodies.pos));
  bodies.vel = (Point *)malloc(nBodies * sizeof(*bodies.vel));

  randomize(bodies.pos, nBodies);
  randomize(bodies.vel, nBodies);

  vtkNew<sensei::ConfigurableAnalysis> analysisAdaptor;
  analysisAdaptor->SetCommunicator(MPI_COMM_WORLD);
  analysisAdaptor->Initialize(configFilename);

  for (int iter = 0; iter <= nIters; ++iter) {
    if (iter > 0) {
      bodyForce(bodies, dt, nBodies); // compute interbody forces

      for (int i = 0 ; i < nBodies; i++) { // integrate position
        bodies.pos[i].x += bodies.vel[i].x * dt;
        bodies.pos[i].y += bodies.vel[i].y * dt;
        bodies.pos[i].z += bodies.vel[i].z * dt;
      }
    }

    vtkNew<vtkFloatArray> floatArrayPos;
    floatArrayPos->Initialize();
    floatArrayPos->SetName("position");
    floatArrayPos->SetNumberOfComponents(3);
    floatArrayPos->SetArray((float *)bodies.pos, 3 * nBodies, /*save=*/0);

    vtkNew<vtkPoints> points;
    points->Initialize();
    points->SetData(floatArrayPos);

    vtkNew<vtkFloatArray> floatArrayVel;
    floatArrayVel->Initialize();
    floatArrayVel->SetName("velocity");
    floatArrayVel->SetNumberOfComponents(3);
    floatArrayVel->SetArray((float *)bodies.vel, 3 * nBodies, /*save=*/0);

    vtkNew<vtkPolyData> polyData;
    polyData->Initialize();
    polyData->SetPoints(points);
    polyData->GetPointData()->AddArray(floatArrayVel);

    vtkNew<sensei::VTKDataAdaptor> vtkDataAdaptor;
    vtkDataAdaptor->SetDataTime(iter * dt);
    vtkDataAdaptor->SetDataTimeStep(iter);
    vtkDataAdaptor->SetDataObject("bodies", polyData);

    analysisAdaptor->Execute(vtkDataAdaptor);

    vtkDataAdaptor->ReleaseData();
  }

  analysisAdaptor->Finalize();

  free(bodies.vel);
  free(bodies.pos);

  MPI_Finalize();

  return 0;
}
