/**
 *
 */

// self
#include <ospSensei/OSPRayVisualization.h>

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


namespace ospSensei {

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


struct OSPRayVisualization::InternalsType {
  void Initialize(MPI_Comm comm, OSPRayVisualization::Type);
  bool Execute(MPI_Comm comm, OSPRayVisualization::Type, size_t nPoints, float *positions);
  bool Execute(MPI_Comm comm, OSPRayVisualization::Type, size_t nx, size_t ny, size_t nz, uint16_t *data);
  void Finalize(OSPRayVisualization::Type);

  // config
  OSPRayVisualization::Type Mode;
  int Width{256};
  int Height{256};

  int FrameNumber{0};
  int CommRank{-1};
  int CommSize{-1};

  struct {
    // particle config
    float GeometryRadius{0.1f};

    // particle member
    OSPData GeometrySpherePositionData{nullptr};
    OSPGeometry Geometry{nullptr};
    OSPMaterial Material{nullptr};
    OSPGeometricModel GeometricModel{nullptr};
  } Particle;

  struct {
    // volume member
    OSPData VolumeDataData{nullptr};
    OSPVolume Volume{nullptr};
    std::vector<float> TransferFunctionColor;
    OSPData TransferFunctionColorData{nullptr};
    std::vector<float> TransferFunctionOpacity;
    OSPData TransferFunctionOpacityData{nullptr};
    OSPTransferFunction TransferFunction{nullptr};
    OSPVolumetricModel VolumetricModel{nullptr};
  } Volume;

  struct {
    // isosurface config
    float GeometryIsovalue{64};

    // isosurface member
    OSPData VolumeDataData{nullptr};
    OSPVolume Volume{nullptr};
    OSPGeometry Geometry{nullptr};
    OSPMaterial Material{nullptr};
    OSPGeometricModel GeometricModel{nullptr};
  } Isosurface;

  OSPDevice Device{nullptr};
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


senseiNewMacro(OSPRayVisualization);


int OSPRayVisualization::Initialize() {
  this->Internals = new InternalsType;
  this->Internals->Initialize(GetCommunicator(), Mode);
  return /*failure=*/0;
}

void OSPRayVisualization::InternalsType::Initialize(MPI_Comm comm, OSPRayVisualization::Type mode) {
  MPI_Comm_rank(comm, &CommRank);
  MPI_Comm_size(comm, &CommSize);

  ospLoadModule("mpi");

  Device = ospNewDevice("mpiDistributed");
  ospDeviceCommit(Device);
  ospSetCurrentDevice(Device);

  if (mode == OSPRayVisualization::PARTICLE) {
    Particle.GeometrySpherePositionData = nullptr;

    Particle.Geometry = ospNewGeometry("sphere");
    ospSetObject(Particle.Geometry, "sphere.position", Particle.GeometrySpherePositionData);
    ospSetFloat(Particle.Geometry, "radius", Particle.GeometryRadius);
    ospCommit(Particle.Geometry);

    Particle.Material = ospNewMaterial(nullptr, "obj");
    ospSetVec3f(Particle.Material, "kd", (CommRank % 4 == 0 ? 1.0f : 0.0f), (CommRank % 4 == 1 ? 1.0f : 0.0f), (CommRank % 4 == 2 ? 1.0f : 0.0f));
    ospCommit(Particle.Material);

    Particle.GeometricModel = ospNewGeometricModel(nullptr);
    ospSetObject(Particle.GeometricModel, "geometry", Particle.Geometry);
    ospSetObject(Particle.GeometricModel, "material", Particle.Material);
    ospCommit(Particle.GeometricModel);
  
  } else if (mode == OSPRayVisualization::VOLUME) {
    Volume.VolumetricModel = nullptr;

  } else if (mode == OSPRayVisualization::ISOSURFACE) {
    Isosurface.GeometricModel = nullptr;
  }

  Group = ospNewGroup();
  if (mode == OSPRayVisualization::PARTICLE) {
    ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, Particle.GeometricModel);

  } else if (mode == OSPRayVisualization::VOLUME) {
    //ospSetObjectAsData(Group, "volume", OSP_VOLUMETRIC_MODEL, Volume.VolumetricModel);

  } else if (mode == OSPRayVisualization::ISOSURFACE) {
    // ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, Isosurface.GeometricModel);

  }
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
}

bool OSPRayVisualization::Execute(sensei::DataAdaptor *data) {
  MPI_Comm comm;
  comm = GetCommunicator();

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  vtkDataObject *mesh;
  if (data->GetMesh(MeshName.c_str(), /*structureOnly=*/false, mesh)) {
    SENSEI_ERROR("Failed to get mesh '" << MeshName << "'")//no semicolon
    return /*success=*/false;
  }

  int association;
  if (Mode == PARTICLE) {
    association = vtkDataObject::POINT;
  } else if (Mode == VOLUME) {
    association = vtkDataObject::CELL;
  } else if (Mode == ISOSURFACE) {
    association = vtkDataObject::CELL;
  } else {
    SENSEI_ERROR("Unexpected mode: " << Mode)//no semicolon
    return /*success=*/false;
  }

  if (data->AddArray(mesh, MeshName.c_str(), association, ArrayName.c_str())) {
    SENSEI_ERROR("Failed to get array '" << ArrayName << "' with association '" << association << "' from mesh '" << MeshName << "' block " << rank)//no semicolon
    return /*success=*/false;
  }

  if (Mode == PARTICLE) {
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
      SENSEI_ERROR("Expected points from mesh '" << MeshName << "' block " << rank << " to have data")//no semicolon
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

    return /*success=*/this->Internals->Execute(comm, Mode, nPoints, positions);
  
  } else if (Mode == VOLUME || Mode == ISOSURFACE) {
    vtkUniformGridAMR *uniformGridAMR;
    uniformGridAMR = dynamic_cast<vtkUniformGridAMR *>(mesh);
    if (uniformGridAMR == nullptr) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' to be a vtkUniformGridAMR")//no semicolon
      return /*success=*/false;
    }

    if (uniformGridAMR->GetNumberOfLevels() != 1) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' to have 1 level")//no semicolon
      return /*success=*/false;
    }

    if (uniformGridAMR->GetNumberOfDataSets(0) != size) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 to have '" << size << "' datasets")//no semicolon
      return /*success=*/false;
    }

    vtkUniformGrid *uniformGrid;
    uniformGrid = uniformGridAMR->GetDataSet(0, rank);
    if (uniformGrid == nullptr) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 index '" << rank << "' to be non-null")//no semicolon
      return /*success=*/false;
    }

    int dims[3];
    uniformGrid->GetDimensions(dims);

    vtkDataArray *dataArray;
    dataArray = uniformGrid->GetCellData()->GetArray(ArrayName.c_str());
    if (dataArray == nullptr) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 index '" << rank << "' cell data array '" << ArrayName << "' to be non-null")//no semicolon
      return /*success=*/false;
    }
    
    vtkUnsignedShortArray *unsignedShortArray;
    unsignedShortArray = vtkUnsignedShortArray::SafeDownCast(dataArray);
    if (unsignedShortArray == nullptr) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 index '" << rank << "' cell data array '" << ArrayName << "' to be a vtkUnsignedShortArray")//no semicolon
      return /*success=*/false;
    }
    
    if (unsignedShortArray->GetNumberOfComponents() != 1) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 index '" << rank << "' cell data array '" << ArrayName << "' to have 1 component")//no semicolon
      return /*success=*/false;
    }
    
    if (!unsignedShortArray->HasStandardMemoryLayout()) {
      SENSEI_ERROR("Expected mesh '" << MeshName << "' level 0 index '" << rank << "' cell data array '" << ArrayName << "' to have standard memory layout")//no semicolon
      return /*success=*/false;
    }

    return this->Internals->Execute(
      comm,
      Mode,
      dims[0], dims[1], dims[2],
      unsignedShortArray->GetPointer(0));
  }

  return /*success=*/false;
}

bool OSPRayVisualization::InternalsType::Execute(MPI_Comm comm, OSPRayVisualization::Type mode, size_t nx, size_t ny, size_t nz, uint16_t *volume) {
  uint16_t lo, hi;
  lo = hi = volume[0];
  for (size_t i=1; i<nx*ny*nz; ++i) {
    if (volume[i] < lo) lo = volume[i];
    if (volume[i] > hi) hi = volume[i];
  }
  std::fprintf(stderr, "lo = %u, hi = %u\n", lo, hi);

  if (mode == OSPRayVisualization::VOLUME) {
    if (Volume.VolumeDataData) {
      ospRelease(Volume.VolumeDataData);
      Volume.VolumeDataData = nullptr;
    }

    Volume.VolumeDataData = ospNewSharedData(volume, OSP_USHORT, nx, 0, ny, 0, nz, 0);
    ospCommit(Volume.VolumeDataData);

    Volume.Volume = ospNewVolume("structuredRegular");
    ospSetVec3f(Volume.Volume, "gridOrigin", -0.5f, -0.5f, -0.5f);
    ospSetVec3f(Volume.Volume, "gridSpacing", 1.0f/(float)nx, 1.0f/(float)ny, 1.0f/(float)nz);
    ospSetObject(Volume.Volume, "data", Volume.VolumeDataData);
    ospSetBool(Volume.Volume, "cellCentered", 1);
    ospCommit(Volume.Volume);

    if (Volume.TransferFunctionColorData) {
      ospRelease(Volume.TransferFunctionColorData);
      Volume.TransferFunctionColorData = nullptr;
    }

    Volume.TransferFunctionColor.clear();
    Volume.TransferFunctionColor.insert(Volume.TransferFunctionColor.end(), {
      (CommRank % 3 == 0 ? 1.0f : 0.0f),
      (CommRank % 3 == 1 ? 1.0f : 0.0f),
      (CommRank % 3 == 2 ? 1.0f : 0.0f),
      (CommRank % 3 == 0 ? 1.0f : 0.0f),
      (CommRank % 3 == 1 ? 1.0f : 0.0f),
      (CommRank % 3 == 2 ? 1.0f : 0.0f),
    });

    Volume.TransferFunctionColorData = ospNewSharedData(Volume.TransferFunctionColor.data(), OSP_VEC3F, Volume.TransferFunctionColor.size() / 3);
    ospCommit(Volume.TransferFunctionColorData);

    if (Volume.TransferFunctionOpacityData) {
      ospRelease(Volume.TransferFunctionOpacityData);
      Volume.TransferFunctionOpacityData = nullptr;
    }

    Volume.TransferFunctionOpacity.clear();
    Volume.TransferFunctionOpacity.insert(Volume.TransferFunctionOpacity.end(), {
      0.0f,
      1.0f,
    });

    Volume.TransferFunctionOpacityData = ospNewSharedData(Volume.TransferFunctionOpacity.data(), OSP_FLOAT, Volume.TransferFunctionOpacity.size() / 1);
    ospCommit(Volume.TransferFunctionOpacityData);

    Volume.TransferFunction = ospNewTransferFunction("piecewiseLinear");
    ospSetObject(Volume.TransferFunction, "color", Volume.TransferFunctionColorData);
    ospSetObject(Volume.TransferFunction, "opacity", Volume.TransferFunctionOpacityData);
    ospSetVec2f(Volume.TransferFunction, "valueRange", (float)0.0f, (float)hi);
    ospCommit(Volume.TransferFunction);

    Volume.VolumetricModel = ospNewVolumetricModel(nullptr);
    ospSetObject(Volume.VolumetricModel, "volume", Volume.Volume);
    ospSetObject(Volume.VolumetricModel, "transferFunction", Volume.TransferFunction);
    ospCommit(Volume.VolumetricModel);

    ospSetObjectAsData(Group, "volume", OSP_VOLUMETRIC_MODEL, Volume.VolumetricModel);
    ospCommit(Group);
  
  } else if (mode == OSPRayVisualization::ISOSURFACE) {
    if (Isosurface.VolumeDataData) {
      ospRelease(Isosurface.VolumeDataData);
      Isosurface.VolumeDataData = nullptr;
    }

    Isosurface.VolumeDataData = ospNewSharedData(volume, OSP_USHORT, nx, 0, ny, 0, nz, 0);
    ospCommit(Isosurface.VolumeDataData);

    Isosurface.Volume = ospNewVolume("structuredRegular");
    ospSetVec3f(Isosurface.Volume, "gridOrigin", -0.5f, -0.5f, -0.5f);
    ospSetVec3f(Isosurface.Volume, "gridSpacing", 1.0f/(float)nx, 1.0f/(float)ny, 1.0f/(float)nz);
    ospSetObject(Isosurface.Volume, "data", Isosurface.VolumeDataData);
    ospCommit(Isosurface.Volume);

    Isosurface.Geometry = ospNewGeometry("isosurface");
    ospSetFloat(Isosurface.Geometry, "isovalue", Isosurface.GeometryIsovalue);
    ospSetObject(Isosurface.Geometry, "volume", Isosurface.Volume);
    ospCommit(Isosurface.Geometry);

    Isosurface.Material = ospNewMaterial(nullptr, "obj");
    ospSetVec3f(Isosurface.Material, "kd", (CommRank % 4 == 0 ? 1.0f : 0.0f), (CommRank % 4 == 1 ? 1.0f : 0.0f), (CommRank % 4 == 2 ? 1.0f : 0.0f));
    ospCommit(Isosurface.Material);

    Isosurface.GeometricModel = ospNewGeometricModel(nullptr);
    ospSetObject(Isosurface.GeometricModel, "geometry", Isosurface.Geometry);
    ospSetObject(Isosurface.GeometricModel, "material", Isosurface.Material);
    ospCommit(Isosurface.GeometricModel);

    ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, Isosurface.GeometricModel);
    ospCommit(Group);

  }

  ospSetObject(Instance, "group", Group);
  ospCommit(Instance);

  ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
  ospSetObject(World, "region", WorldRegionData);
  ospCommit(World);

  ospResetAccumulation(FrameBuffer);
  Future = ospRenderFrame(FrameBuffer, Renderer, Camera, World);
  ospWait(Future, OSP_TASK_FINISHED);
  ospRelease(Future);
  Future = nullptr;

  if (CommRank == 0) {
    std::string filename = std::string("ospSensei.") + std::to_string(FrameNumber) + std::string(".ppm");
    const void *fb = ospMapFrameBuffer(FrameBuffer, OSP_FB_COLOR);
    writePPM(filename.c_str(), Width, Height, static_cast<const uint32_t *>(fb));
    ospUnmapFrameBuffer(fb, FrameBuffer);
  }

  ++FrameNumber;

  return /*success=*/true;
}

bool OSPRayVisualization::InternalsType::Execute(MPI_Comm comm, OSPRayVisualization::Type mode, size_t nPoints, float *positions) {
  float lo[3], hi[3];

  std::fprintf(stderr, "nPoints = %lu\n", nPoints);
  if (nPoints > 0) {
    std::fprintf(stderr, "positions[0] = %+0.2f\n", positions[0]);
    std::fprintf(stderr, "positions[1] = %+0.2f\n", positions[1]);
    std::fprintf(stderr, "positions[2] = %+0.2f\n", positions[2]);

    lo[0] = hi[0] = positions[0*3+0];
    lo[1] = hi[1] = positions[0*3+1];
    lo[2] = hi[2] = positions[0*3+2];
    for (size_t i=1; i<nPoints; ++i) {
      if (positions[i*3+0] < lo[0]) lo[0] = positions[i*3+0];
      if (positions[i*3+1] < lo[1]) lo[1] = positions[i*3+1];
      if (positions[i*3+2] < lo[2]) lo[2] = positions[i*3+2];
      if (positions[i*3+0] > hi[0]) hi[0] = positions[i*3+0];
      if (positions[i*3+1] > hi[1]) hi[1] = positions[i*3+1];
      if (positions[i*3+2] > hi[2]) hi[2] = positions[i*3+2];
    }

    std::fprintf(stderr, "min = %+0.2f, %+0.2f, %+0.2f\n", lo[0], lo[1], lo[2]);
    std::fprintf(stderr, "max = %+0.2f, %+0.2f, %+0.2f\n", hi[0], hi[1], hi[2]);
  }

  if (Particle.GeometrySpherePositionData) {
    ospRelease(Particle.GeometrySpherePositionData);
    Particle.GeometrySpherePositionData = nullptr;
  }

  Particle.GeometrySpherePositionData = ospNewSharedData(positions, OSP_VEC3F, nPoints);
  ospCommit(Particle.GeometrySpherePositionData);
  
  ospSetObject(Particle.Geometry, "sphere.position", Particle.GeometrySpherePositionData);
  ospCommit(Particle.Geometry);

  ospSetObject(Particle.GeometricModel, "geometry", Particle.Geometry);
  ospCommit(Particle.GeometricModel);

  ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, Particle.GeometricModel);
  ospCommit(Group);

  ospSetObject(Instance, "group", Group);
  ospCommit(Instance);

  WorldRegion.clear();
  WorldRegionData = nullptr;

  if (nPoints > 0) {
    WorldRegion.insert(WorldRegion.end(), {
      lo[0]-Particle.GeometryRadius, lo[1]-Particle.GeometryRadius, lo[2]-Particle.GeometryRadius,
      hi[0]+Particle.GeometryRadius, hi[1]+Particle.GeometryRadius, hi[2]+Particle.GeometryRadius,
    });

    WorldRegionData = ospNewSharedData(WorldRegion.data(), OSP_BOX3F, WorldRegion.size() / 6);
    ospCommit(WorldRegionData);
  }

  ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
  ospSetObject(World, "region", WorldRegionData);
  ospCommit(World);

  ospResetAccumulation(FrameBuffer);
  Future = ospRenderFrame(FrameBuffer, Renderer, Camera, World);
  ospWait(Future, OSP_TASK_FINISHED);
  ospRelease(Future);
  Future = nullptr;

  if (CommRank == 0) {
    std::string filename = std::string("ospSensei.") + std::to_string(FrameNumber) + std::string(".ppm");
    const void *fb = ospMapFrameBuffer(FrameBuffer, OSP_FB_COLOR);
    writePPM(filename.c_str(), Width, Height, static_cast<const uint32_t *>(fb));
    ospUnmapFrameBuffer(fb, FrameBuffer);
  }

  ++FrameNumber;

  return /*success=*/true;
}

void OSPRayVisualization::InternalsType::Finalize(OSPRayVisualization::Type mode) {
  ospRelease(FrameBuffer);
  ospRelease(Renderer);
  ospRelease(Camera);
  ospRelease(World);
  ospRelease(Light);
  ospRelease(Instance);
  ospRelease(Group);
  if (mode == OSPRayVisualization::PARTICLE) {
    ospRelease(Particle.GeometricModel);
    ospRelease(Particle.Material);
    ospRelease(Particle.Geometry);
    if (Particle.GeometrySpherePositionData) ospRelease(Particle.GeometrySpherePositionData);

  } else if (mode == OSPRayVisualization::VOLUME) {
    ospRelease(Volume.VolumetricModel);
    ospRelease(Volume.TransferFunction);
    ospRelease(Volume.TransferFunctionOpacityData);
    ospRelease(Volume.TransferFunctionColorData);
    ospRelease(Volume.Volume);
    ospRelease(Volume.VolumeDataData);
  }
  ospDeviceRelease(Device);

  ospShutdown();
}

int OSPRayVisualization::Finalize() {
  Internals->Finalize(Mode);
  delete this->Internals;

  return 0; // no error
}


} /* namespace ospSensei */
// vim: ts=2:sts=2
