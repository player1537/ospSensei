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
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

// ospSensei
#include <ospSensei/OSPRayParticleVisualization.h>

#define SOFTENING 1e-9f


typedef struct {
  float x, y, z;
} Point;

struct Bodies {
  Bodies(size_t size_)
    : count(0)
    , size(size_)
    , pos(new Point[size_])
    , vel(new Point[size_])
  {};

  ~Bodies() {
    delete[] vel;
    delete[] pos;
  }

  size_t count;
  size_t size;
  Point *pos;
  Point *vel;

  void extents(float ret[6]);
  void reduce(float bounds[6], Bodies &ret);
  void randomize(size_t count_);
  void broadcast(const MPI_Comm &comm, int root);
  void step(float dt);
  vtkDataObject *vtk();
};

void Bodies::extents(float ret[6]) {
  if (count <= 0) {
    fprintf(stderr, "Error: cannot compute extents of empty object\n");
    std::exit(1);
  }

  ret[3*0+0] = ret[3*1+0] = pos[0].x;
  ret[3*0+1] = ret[3*1+1] = pos[0].y;
  ret[3*0+2] = ret[3*1+2] = pos[0].z;
  
  for (size_t i=1; i<count; ++i) {
    if (pos[i].x < ret[3*0+0]) ret[3*0+0] = pos[i].x;
    if (pos[i].y < ret[3*0+1]) ret[3*0+1] = pos[i].y;
    if (pos[i].z < ret[3*0+2]) ret[3*0+2] = pos[i].z;
    if (pos[i].x > ret[3*1+0]) ret[3*1+0] = pos[i].x;
    if (pos[i].y > ret[3*1+1]) ret[3*1+1] = pos[i].y;
    if (pos[i].z > ret[3*1+2]) ret[3*1+2] = pos[i].z;
  }
}

void Bodies::reduce(float bounds[6], Bodies &ret) {
  ret.count = 0;

  if (count == 0) {
    return;
  }

  float dataextents[6];
  extents(dataextents);

  for (size_t i=0; i<count; ++i) {
    float x = (pos[i].x - dataextents[3*0+0]) / (dataextents[3*1+0] - dataextents[3*0+0] + 1.0f);
    float y = (pos[i].y - dataextents[3*0+1]) / (dataextents[3*1+1] - dataextents[3*0+1] + 1.0f);
    float z = (pos[i].z - dataextents[3*0+2]) / (dataextents[3*1+2] - dataextents[3*0+2] + 1.0f);
    bool xwithin = bounds[3*0+0] <= x && x < bounds[3*1+0];
    bool ywithin = bounds[3*0+1] <= y && y < bounds[3*1+1];
    bool zwithin = bounds[3*0+2] <= z && z < bounds[3*1+2];

    if (xwithin && ywithin && zwithin) {
      ret.pos[ret.count] = pos[i];
      ret.vel[ret.count] = vel[i];
      ++ret.count;
    }
  }
}

void Bodies::randomize(size_t count_) {
  if (size < count_) {
    std::fprintf(stderr, "Error: Bodies size too small (%zu < requested %zu)\n", size, count_);
    std::exit(1);
  }

  count = count_;

  for (size_t i = 0; i < count; i++) {
    pos[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    pos[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    pos[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
  for (size_t i = 0; i < count; i++) {
    vel[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    vel[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    vel[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void Bodies::broadcast(const MPI_Comm &comm, int root) {
  MPI_Bcast(&count, 1, MPI_COUNT, root, comm);
  MPI_Bcast(pos, 3 * count, MPI_FLOAT, root, comm);
  MPI_Bcast(vel, 3 * count, MPI_FLOAT, root, comm);
}

void Bodies::step(float dt) {
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < count; i++) { 
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;

    for (size_t j = 0; j < count; j++) {
      float dx = pos[j].x - pos[i].x;
      float dy = pos[j].y - pos[i].y;
      float dz = pos[j].z - pos[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3;
      Fy += dy * invDist3;
      Fz += dz * invDist3;
    }

    vel[i].x += dt * Fx;
    vel[i].y += dt * Fy;
    vel[i].z += dt * Fz;
  }

  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0 ; i < count; i++) { // integrate position
    pos[i].x += vel[i].x * dt;
    pos[i].y += vel[i].y * dt;
    pos[i].z += vel[i].z * dt;
  }
}

vtkDataObject *Bodies::vtk() {
  vtkNew<vtkCellArray::ArrayType32> offsets;
  offsets->Initialize();
  offsets->SetNumberOfComponents(1);
  offsets->SetNumberOfTuples(count + 1);
  for (size_t i=0; i<count; ++i) {
    offsets->SetTypedComponent(i, 0, i);
  }
  offsets->SetTypedComponent(count, 0, count);

  vtkNew<vtkCellArray::ArrayType32> connectivity;
  connectivity->Initialize();
  connectivity->SetNumberOfComponents(1);
  connectivity->SetNumberOfTuples(count);
  for (size_t i=0; i<count; ++i) {
    connectivity->SetTypedComponent(i, 0, i);
  }

  vtkNew<vtkFloatArray> floatArrayPos;
  floatArrayPos->Initialize();
  floatArrayPos->SetName("position");
  floatArrayPos->SetNumberOfComponents(3);
  floatArrayPos->SetArray((float *)pos, 3 * count, /*save=*/1);

  vtkNew<vtkPoints> points;
  points->Initialize();
  points->SetData(floatArrayPos);

  vtkNew<vtkFloatArray> floatArrayVel;
  floatArrayVel->Initialize();
  floatArrayVel->SetName("velocity");
  floatArrayVel->SetNumberOfComponents(3);
  floatArrayVel->SetArray((float *)vel, 3 * count, /*save=*/1);

  vtkNew<vtkCellArray> cellArrayVerts;
  cellArrayVerts->Initialize();
  cellArrayVerts->SetData(offsets, connectivity);

  vtkPolyData *polyData = vtkPolyData::New();
  polyData->Initialize();
  polyData->SetPoints(points);
  polyData->GetPointData()->AddArray(floatArrayPos);
  polyData->GetPointData()->AddArray(floatArrayVel);
  polyData->SetVerts(cellArrayVerts);

  return static_cast<vtkDataObject *>(polyData);
}


int main(int argc, char** argv) {
  int provided;
  int success = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (success != MPI_SUCCESS) {
    SENSEI_ERROR("Error while initializing MPI")//no semicolon
    return 1;
  }

  if (provided != MPI_THREAD_MULTIPLE) {
    SENSEI_ERROR("MPI provided the wrong level of thread support")//no semicolon
    return 1;
  }

  size_t nBodies = 30000;
  size_t nIters = 10;
  float dt = 0.01f; // time step
  std::string configFilename{"nbody.xml"};
  if (argc > 1 && argv[1][0] != '\0') nBodies = strtoul(argv[1], NULL, 10);
  if (argc > 2 && argv[2][0] != '\0') nIters = strtoul(argv[2], NULL, 10);
  if (argc > 3 && argv[3][0] != '\0') dt = atof(argv[3]);
  if (argc > 4 && argv[4][0] != '\0') configFilename = argv[4];

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  size_t length = (nBodies + size - 1) / size;
  size_t start = rank * length;
  if (start + length > nBodies) {
    length = nBodies - start;
  }

  Bodies bodies(nBodies);
  Bodies filtered(nBodies);

  if (rank == 0) {
    bodies.randomize(nBodies);
  }

  bodies.broadcast(MPI_COMM_WORLD, 0);

  float bounds[6] = {
    0.0f, 0.0f, 0.0f,
    1.0f, 1.0f, 1.0f,
  };
  bounds[3*0+1] = (float)(rank + 0) / (float)size;
  bounds[3*1+1] = (float)(rank + 1) / (float)size;

  vtkNew<ospSensei::OSPRayParticleVisualization> analysisAdaptor;
  analysisAdaptor->SetCommunicator(MPI_COMM_WORLD);
  analysisAdaptor->SetWidth(512);
  analysisAdaptor->SetHeight(512);
  analysisAdaptor->Initialize();
  analysisAdaptor->SetMeshName("bodies");
  analysisAdaptor->SetArrayName("position");

  for (int iter = 0; iter <= nIters; ++iter) {
    SENSEI_STATUS("On iteration " << iter)//no semicolon
    if (iter > 0) {
      bodies.step(dt);
    }

    bodies.reduce(bounds, filtered);

    vtkSmartPointer<vtkDataObject> dataObject = filtered.vtk();

    vtkNew<vtkMultiBlockDataSet> multiBlockDataSet;
    multiBlockDataSet->Initialize();
    multiBlockDataSet->SetNumberOfBlocks(size);
    for (size_t i=0; i<size; ++i) {
      multiBlockDataSet->SetBlock(i, nullptr);
    }
    multiBlockDataSet->SetBlock(rank, dataObject);

    vtkNew<sensei::VTKDataAdaptor> vtkDataAdaptor;
    vtkDataAdaptor->SetDataTime(iter * dt);
    vtkDataAdaptor->SetDataTimeStep(iter);
    vtkDataAdaptor->SetDataObject("bodies", multiBlockDataSet);

    analysisAdaptor->Execute(vtkDataAdaptor);

    vtkDataAdaptor->ReleaseData();
  }

  analysisAdaptor->Finalize();
  analysisAdaptor.Reset();

  MPI_Finalize();

  return 0;
}
// vim: ts=2:sts=2
