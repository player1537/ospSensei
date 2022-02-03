// Copyright 2014 Mark Harris (https://github.com/harrism)
// SPDX-License-Identifier: Apache-2.0

// Modifications (2022) by Tanner Hobson (https://github.com/player1537)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SOFTENING 1e-9f

typedef struct {
  float x, y, z;
  float vx, vy, vz;
} Body;

void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void bodyForce(Body *p, float dt, int n) {
  #pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n; i++) { 
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }

    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

class NBodyDataAdapter : public sensei::DataAdapter {
public:  // VTK API
  static NBodyDataAdapter *New();
  senseiTypeMacro(NBodyDataAdapter, sensei::DataAdapter);

protected:  // C++ API
  NBodyDataAdapter();
  ~NBodyDataAdapter() = default;
  NBodyDataAdapter(const NBodyDataAdapter &) = delete;
  NBodyDataAdapter &operator=(const NBodyDataAdapter &) = delete;

public:  // NBody API
  int Initialize(int nBodies, int nIters, float dt);

  int SetBodies(Body *bodies);

private:  // NBody member variables
  int m_nBodies;
  int m_nIters;
  float m_dt;
  Body *m_bodies{nullptr};

public:  // SENSEI API

  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) override;

  int GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) override;

  int AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) override;

  int AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) override;

  int ReleaseData() override;
};


senseiNewMacro(NBodyDataAdapter);

NBodyDataAdapter::NBodyDataAdapter()
{
}


int NBodyDataAdapter::Initialize(int nBodies, int nIters, float dt) {
  m_nBodies = nBodies
  m_nIters = nIters;
  m_dt = dt;

  return 0;
}

int NBodyDataAdapter::SetBodies(Body *bodies) {
  m_bodies = bodies;

  return 0;
}


int NBodyDataAdapter::GetNumberOfMeshes(unsigned int &numMeshes) {
  numMeshes = 1;

  return 0;
}

int NBodyDataAdapter::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) {
  if (id != 0) {
    SENSEI_ERROR("invalid mesh id " << id) // no semicolon
    return -1;
  }

  metadata->MeshName = "mesh";
  metadata->MeshType = VTK_POLY_DATA;
  metadata->BlockType = VTK_POINT_SET;

  metadata->NumGhostCells = 0;
  metadata->NumGhostBlocks = 0;

  metadata->NumArrays = 1;
  metadata->ArrayName.emplace_back("nbody");
  metadata->ArrayCentering.emplace_back(vtkDataObject::POINT);
  metadata->ArrayType.emplace_back(VTK_FLOAT);
  metadata->ArrayComponents.emplace_back(6);

  metadata->NumBlocks = 1;
  metadata->NumBlocksLocal.emplace_back(1);

  if (metadata->Flags.BlockDecompSet()) {
    SENSEI_ERROR("BlockDecompSet support not implemented")// no semicolon
    return -1;
  }

  if (metadata->Flags.BlockSizeSet()) {
    SENSEI_ERROR("BlockSizeSet support not implemented")// no semicolon
    return -1;
  }

  if (metadata->Flags.BlockExtentsSet()) {
    SENSEI_ERROR("BlockExtentsSet support not implemented")// no semicolon
    return -1;
  }

  if (metadata->Flags.BlockBoundsSet()) {
    SENSEI_ERROR("BlockBoundsSet support not implemented")// no semicolon
    return -1;
  }

  return 0;
}

int NBodyDataAdapter::GetMesh(const std::string &meshName, bool structureOnly, vtkDataObject *&mesh) {
  if (meshName != "mesh") {
    SENSEI_ERROR("only the mesh 'mesh' is supported")// no semicolon
    return -1;
  }

  sensei::MeshMetadataFlags flags{0};
  sensei::MeshMetadataPtr mmd = sensei::MeshMetadataPtr::New(flags);
  GetMeshMetadata(0, mmd);

  vtkSmartPointer<vtkPolyData> polyData;
  polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->Initialize();

  return 0;
}

int NBodyDataAdapter::AddArray(vtkDataObject* mesh, const std::string &meshName, int association, const std::string &arrayName) {
  return 0;
}

int NBodyDataAdapter::AddGhostCellsArray(vtkDataObject* mesh, const std::string &meshName) {
  return 0;
}

int NBodyDataAdapter::ReleaseData() {
  if (m_bodies) {
    m_bodies = nullptr;
  }

  return 0;
}


int main(const int argc, const char** argv) {
  MPI_Init(&argc, &argv);

  int nBodies = 30000;
  int nIters = 10;
  float dt = 0.01f; // time step
  std::string configFilename{"nbody.xml"};
  if (argc > 1) nBodies = atoi(argv[1]);
  if (argc > 2) nIters = atoi(argv[2]);
  if (argc > 3) dt = atof(argv[3]);
  if (argc > 4) configFilename = argv[4];

  int bytes = nBodies*sizeof(Body);
  float *buf = (float*)malloc(bytes);
  Body *p = (Body*)buf;

  vtkSmartPointer<NBodyDataAdapter> dataAdapter;
  dataAdapter = vtkSmartPointer<NBodyDataAdapter>::New();
  dataAdapter->Initialize(nBodies, nIters, dt);
  dataAdapter->SetDataTimestep(-1);

  vtkSmartPointer<sensei::ConfigurableAnalysis> analysisAdapter;
  analysisAdapter = vtkSmartPointer<sensei::ConfigurableAnalysis>::New();
  analysisAdapter->SetCommunicator(MPI_COMM_WORLD);
  analysisAdapter->Initialize(configFilename)

  randomizeBodies(buf, 6*nBodies); // Init pos / vel data

  dataAdapter->SetBodies(p);
  dataAdapter->SetDataTime(0.0);
  dataAdapter->SetDataTimeStep(0);
  analysisAdapter->Execute(dataAdapter);

  dataAdapter->ReleaseData();

  for (int iter = 1; iter <= nIters; iter++) {
    bodyForce(p, dt, nBodies); // compute interbody forces

    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }

    dataAdapter->SetBodies(p);
    dataAdapter->SetDataTime(iter * dt);
    dataAdapter->SetDataTimeStep(iter);
    analysisAdapter->Execute(dataAdapter);

    dataAdapter->ReleaseData();
  }

  analysisAdapter->Finalize();

  analysisAdapter = nullptr;
  dataAdapter = nullptr;

  free(buf);

  MPI_Finalize();

  return 0;
}
