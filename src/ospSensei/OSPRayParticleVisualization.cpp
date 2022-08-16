/**
 *
 */

// self
#include <ospSensei/OSPRayParticleVisualization.h>

// stdlib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// OSPRay
#include <ospray/ospray.h>
#include <ospray/ospray_util.h>

// SENSEI
#include <DataAdaptor.h>
#include <Error.h>

// VTK
#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUniformGrid.h>
#include <vtkUniformGridAMR.h>
#include <vtkUnsignedShortArray.h>


//---

namespace ospSensei {


//---

// helper function to write the rendered image as PPM file
static void writePPM(const char *fileName, int size_x, int size_y, const uint32_t *pixel) {
  using namespace std;

  FILE *file = fopen(fileName, "wb");
  if (!file) {
    fprintf(stderr, "fopen('%s', 'wb') failed: %d", fileName, errno);
    return;
  }
  fprintf(file, "P6\n%i %i\n255\n", size_x, size_y);
  unsigned char *out = (unsigned char *)alloca(3 * size_x);
  for (int y = 0; y < size_y; y++) {
    const unsigned char *in =
        (const unsigned char *)&pixel[(size_y - 1 - y) * size_x];
    for (int x = 0; x < size_x; x++) {
      out[3 * x + 0] = in[4 * x + 0];
      out[3 * x + 1] = in[4 * x + 1];
      out[3 * x + 2] = in[4 * x + 2];
    }
    fwrite(out, 3 * size_x, sizeof(char), file);
  }
  fprintf(file, "\n");
  fclose(file);
}


//---

struct OSPRayParticleVisualization::InternalsType : public vtkObject {
  virtual InternalsType *New();
  vtkSetMacro(Communicator, MPI_Comm);
  vtkSetMacro(Width, int);
  vtkSetMacro(Height, int);
  char const *Initialize();
  vtkSetMacro(NumberOfPoints, size_t);
  vtkSetMacro(PointPositions, float *);
  char const *Execute();
  char const *Finalize();

private:
  MPI_Comm Communicator;
  int Width;
  int Height;
  size_t NumberOfPoints;
  float *PointPositions;

  int FrameNumber{0};

  OSPDevice Device{nullptr};
  float GeometryRadius{0.1f};
  OSPData GeometrySpherePositionData{nullptr};
  OSPGeometry Geometry{nullptr};
  OSPMaterial Material{nullptr};
  OSPGeometricModel GeometricModel{nullptr};
  OSPGroup Group{nullptr};
  OSPInstance Instance{nullptr};
  OSPLight Light{nullptr};
  std::vector<float> WorldRegion;
  OSPData WorldRegionData{nullptr};
  OSPWorld World{nullptr};
  OSPCamera Camera{nullptr};
  OSPRenderer Renderer{nullptr};
  OSPFrameBuffer FrameBuffer{nullptr};
  OSPFuture Future{nullptr};
};

vtkStandardNewMacro(OSPRayParticleVisualization::InternalsType);

char const *OSPRayParticleVisualization::InternalsType::Initialize() {
  int rank, size;
  MPI_Comm_rank(Communicator, &rank);
  MPI_Comm_size(Communicator, &size);

  ospLoadModule("mpi");

  Device = ospNewDevice("mpiDistributed");
  ospDeviceCommit(Device);
  ospSetCurrentDevice(Device);

  GeometrySpherePositionData = nullptr;

  Geometry = ospNewGeometry("sphere");
  ospSetObject(Geometry, "sphere.position", GeometrySpherePositionData);
  ospSetFloat(Geometry, "radius", GeometryRadius);
  ospCommit(Geometry);

  Material = ospNewMaterial(nullptr, "obj");
  ospSetVec3f(Material, "kd", (rank % 4 == 0 ? 1.0f : 0.0f), (rank % 4 == 1 ? 1.0f : 0.0f), (rank % 4 == 2 ? 1.0f : 0.0f));
  ospCommit(Material);

  GeometricModel = ospNewGeometricModel(nullptr);
  ospSetObject(GeometricModel, "geometry", Geometry);
  ospSetObject(GeometricModel, "material", Material);
  ospCommit(GeometricModel);

  Group = ospNewGroup();
  ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, GeometricModel);
  ospCommit(Group);

  Instance = ospNewInstance(nullptr);
  ospSetObject(Instance, "group", Group);
  ospCommit(Instance);

  Light = ospNewLight("ambient");
  ospCommit(Light);

  WorldRegionData = nullptr;

  World = ospNewWorld();
  ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
  ospSetObjectAsData(World, "light", OSP_LIGHT, Light);
  ospSetObject(World, "region", WorldRegionData);
  ospCommit(World);

  Camera = ospNewCamera("perspective");
  ospSetFloat(Camera, "aspect", (float)Width / (float)Height);
  ospSetVec3f(Camera, "position", 0.0f, 0.0f, 0.75f);
  ospSetVec3f(Camera, "direction", 0.0f, 0.0f, -1.0f);
  ospSetVec3f(Camera, "up", 0.0f, 1.0f, 0.0f);
  ospCommit(Camera);

  Renderer = ospNewRenderer("mpiRaycast");
  ospSetInt(Renderer, "pixelSamples", 16);
  ospSetVec3f(Renderer, "backgroundColor", 0.0f, 0.0f, 0.0f);
  ospCommit(Renderer);

  FrameBuffer = ospNewFrameBuffer(Width, Height, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);

  Future = nullptr;

  return /*error=*/nullptr;
}

char const *OSPRayParticleVisualization::InternalsType::Execute() {
  int rank, size;
  MPI_Comm_rank(Communicator, &rank);
  MPI_Comm_size(Communicator, &size);

  float lo[3], hi[3];

  std::fprintf(stderr, "NumberOfPoints = %lu\n", NumberOfPoints);
  if (NumberOfPoints > 0) {
    std::fprintf(stderr, "PointPositions[0] = %+0.2f\n", PointPositions[0]);
    std::fprintf(stderr, "PointPositions[1] = %+0.2f\n", PointPositions[1]);
    std::fprintf(stderr, "PointPositions[2] = %+0.2f\n", PointPositions[2]);

    lo[0] = hi[0] = PointPositions[0*3+0];
    lo[1] = hi[1] = PointPositions[0*3+1];
    lo[2] = hi[2] = PointPositions[0*3+2];
    for (size_t i=1; i<NumberOfPoints; ++i) {
      if (PointPositions[i*3+0] < lo[0]) lo[0] = PointPositions[i*3+0];
      if (PointPositions[i*3+1] < lo[1]) lo[1] = PointPositions[i*3+1];
      if (PointPositions[i*3+2] < lo[2]) lo[2] = PointPositions[i*3+2];
      if (PointPositions[i*3+0] > hi[0]) hi[0] = PointPositions[i*3+0];
      if (PointPositions[i*3+1] > hi[1]) hi[1] = PointPositions[i*3+1];
      if (PointPositions[i*3+2] > hi[2]) hi[2] = PointPositions[i*3+2];
    }

    std::fprintf(stderr, "min = %+0.2f, %+0.2f, %+0.2f\n", lo[0], lo[1], lo[2]);
    std::fprintf(stderr, "max = %+0.2f, %+0.2f, %+0.2f\n", hi[0], hi[1], hi[2]);
  }

  if (GeometrySpherePositionData) {
    ospRelease(GeometrySpherePositionData);
    GeometrySpherePositionData = nullptr;
  }

  GeometrySpherePositionData = ospNewSharedData(PointPositions, OSP_VEC3F, NumberOfPoints);
  ospCommit(GeometrySpherePositionData);
  
  ospSetObject(Geometry, "sphere.position", GeometrySpherePositionData);
  ospCommit(Geometry);

  ospCommit(GeometricModel);

  ospCommit(Group);

  ospCommit(Instance);

  WorldRegion.clear();
  WorldRegionData = nullptr;

  if (NumberOfPoints > 0) {
    WorldRegion.insert(WorldRegion.end(), {
      lo[0]-GeometryRadius, lo[1]-GeometryRadius, lo[2]-GeometryRadius,
      hi[0]+GeometryRadius, hi[1]+GeometryRadius, hi[2]+GeometryRadius,
    });

    WorldRegionData = ospNewSharedData(WorldRegion.data(), OSP_BOX3F, WorldRegion.size() / 6);
    ospCommit(WorldRegionData);
  }

  ospSetObject(World, "region", WorldRegionData);
  ospCommit(World);

  ospResetAccumulation(FrameBuffer);
  Future = ospRenderFrame(FrameBuffer, Renderer, Camera, World);
  ospWait(Future, OSP_TASK_FINISHED);
  ospRelease(Future);
  Future = nullptr;

  if (rank == 0) {
    std::string filename = std::string("ospSensei.") + std::to_string(FrameNumber) + std::string(".ppm");
    const void *fb = ospMapFrameBuffer(FrameBuffer, OSP_FB_COLOR);
    writePPM(filename.c_str(), Width, Height, static_cast<const uint32_t *>(fb));
    ospUnmapFrameBuffer(fb, FrameBuffer);
  }

  ++FrameNumber;

  return /*error=*/nullptr;
}

char const *OSPRayParticleVisualization::InternalsType::Finalize() {
  ospRelease(FrameBuffer);
  FrameBuffer = nullptr;

  ospRelease(Renderer);
  Renderer = nullptr;

  ospRelease(Camera);
  Camera = nullptr;

  ospRelease(World);
  World = nullptr;

  ospRelease(Light);
  Light = nullptr;

  ospRelease(Instance);
  Instance = nullptr;

  ospRelease(Group);
  Group = nullptr;

  ospRelease(GeometricModel);
  GeometricModel = nullptr;

  ospRelease(Material);
  Material = nullptr;

  ospRelease(Geometry);
  Geometry = nullptr;

  if (GeometrySpherePositionData) {
    ospRelease(GeometrySpherePositionData);
    GeometrySpherePositionData = nullptr;
  }

  return /*error=*/nullptr;
}


//---

OSPRayParticleVisualization *OSPRayParticleVisualization::New() {
  auto result = new OSPRayParticleVisualization;
  result->InitializeObjectBase();
  return result;
}

int OSPRayParticleVisualization::Initialize() {
  this->Internals = new InternalsType;
  this->Internals->SetCommunicator(this->GetCommunicator());
  this->Internals->SetWidth(this->Width);
  this->Internals->SetHeight(this->Height);

  char const *error = this->Internals->Initialize();
  if (error != nullptr) {
    std::fprintf(stderr, "Error: %s\n", error);
    return /*failure=*/1;
  }

  return /*failure=*/0;
}

bool OSPRayParticleVisualization::Execute(sensei::DataAdaptor *dataAdaptor) {
  int rank, size;
  MPI_Comm_rank(GetCommunicator(), &rank);
  MPI_Comm_size(GetCommunicator(), &size);

  vtkDataObject *mesh;
  if (dataAdaptor->GetMesh(MeshName.c_str(), /*structureOnly=*/false, mesh)) {
    SENSEI_ERROR("Failed to get mesh '" << MeshName << "'")//no semicolon
    return /*success=*/false;
  }

  if (dataAdaptor->AddArray(mesh, MeshName.c_str(), vtkDataObject::POINT, ArrayName.c_str())) {
    SENSEI_ERROR("Failed to get array '" << ArrayName << "' from mesh '" << MeshName << "' block " << rank)//no semicolon
    return /*success=*/false;
  }

  vtkMultiBlockDataSet *multiBlockDataSet;
  multiBlockDataSet = dynamic_cast<vtkMultiBlockDataSet *>(mesh);
  if (!multiBlockDataSet) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' to be a vtkMultiBlockDataSet")//no semicolon
    return /*success=*/false;
  }

  if (multiBlockDataSet->GetNumberOfBlocks() != size) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' to have only " << size << " blocks")//no semicolon
    return /*success=*/false;
  }

  vtkDataObject *dataObject;
  dataObject = multiBlockDataSet->GetBlock(rank);
  if (dataObject == nullptr) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' block " << rank << " to be non-null")//no semicolon
    return /*success=*/false;
  }

  vtkPolyData *polyData;
  polyData = dynamic_cast<vtkPolyData *>(dataObject);
  if (!polyData) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' block " << rank << " to be a vtkPolyData")//no semicolon
    return /*success=*/false;
  }

  vtkPoints *points;
  points = polyData->GetPoints();
  if (!points) {
    SENSEI_ERROR("Expected points from mesh '" << MeshName << "' block " << rank << " to exist")//no semicolon
    return /*success=*/false;
  }

  vtkDataArray *dataArray;
  dataArray = points->GetData();
  if (!dataArray) {
    SENSEI_ERROR("Expected points from mesh '" << MeshName << "' block " << rank << " to have dataAdaptor")//no semicolon
    return /*success=*/false;
  }

  vtkFloatArray *floatArray;
  if (!(floatArray = vtkFloatArray::SafeDownCast(dataArray))) {
    SENSEI_ERROR("Expected array '" << ArrayName << "' from mesh '" << MeshName << "' block " << rank << " to be a vtkFloatArray")//no semicolon
    return /*success=*/false;
  }

  size_t nComponents = floatArray->GetNumberOfComponents();
  if (nComponents != 3) {
    SENSEI_ERROR("Expected array '" << ArrayName << "' from mesh '" << MeshName << "' block " << rank << " to be a vtkFloatArray with 3 components")//no semicolon
    return /*success=*/false;
  }

  if (!floatArray->HasStandardMemoryLayout()) {
    SENSEI_ERROR("Expected array '" << ArrayName << "' from mesh '" << MeshName << "' block " << rank << " to be a vtkFloatArray with 3 components and standard memory layout")//no semicolon
    return /*success=*/false;
  }

  size_t nPoints = floatArray->GetNumberOfTuples();
  float *positions = floatArray->GetPointer(0);

  this->Internals->SetCommunicator(GetCommunicator());
  this->Internals->SetNumberOfPoints(nPoints);
  this->Internals->SetPointPositions(positions);

  char const *error = this->Internals->Execute();
  if (error != nullptr) {
    std::fprintf(stderr, "Error: %s\n", error);
    return /*success=*/false;
  }

  return /*success=*/true;
}

int OSPRayParticleVisualization::Finalize() {
  char const *error = this->Internals->Finalize();
  if (error != nullptr) {
    std::fprintf(stderr, "Error: %s\n", error);
    return /*failure=*/1;
  }

  delete this->Internals;

  return /*failure=*/0;
}


//---

} /* namespace ospSensei */
