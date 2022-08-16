/**
 *
 */

// self
#include <ospSensei/OSPRayUnstructuredVolumeVisualization.h>

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
#include <vtkCompositeDataIterator.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkMPI.h>
#include <vtkMPICommunicator.h>
#include <vtkMPIController.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPDistributedDataFilter.h>
#include <vtkPKdTree.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkUnstructuredGrid.h>


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

struct OSPRayUnstructuredVolumeVisualization::InternalsType : public vtkObject {
  static InternalsType *New();


  //--- Before Initialize()...

private:
  MPI_Comm Communicator;
  int Width;
  int Height;

public:
  vtkSetMacro(Communicator, MPI_Comm);
  vtkSetMacro(Width, int);
  vtkSetMacro(Height, int);


  //--- Initialize()...

  char const *Initialize();


  //--- Before Execute()...

private:
  float Bounds[6];
  size_t CellTypeSize;
  uint8_t *CellTypeData;
  size_t CellIndexSize;
  uint32_t *CellIndexData;
  size_t VertexPositionSize;
  float *VertexPositionData;
  size_t CellDataSize;
  float *CellDataData;
  size_t IndexSize;
  uint32_t *IndexData;

public:
  vtkSetVector6Macro(Bounds, float);

  void SetCellType(size_t _arg1, uint8_t *_arg2) {
    this->CellTypeSize = _arg1;
    this->CellTypeData = _arg2;
    this->Modified();
  }

  void SetCellIndex(size_t _arg1, uint32_t *_arg2) {
    this->CellIndexSize = _arg1;
    this->CellIndexData = _arg2;
    this->Modified();
  }

  void SetVertexPosition(size_t _arg1, float *_arg2) {
    this->VertexPositionSize = _arg1;
    this->VertexPositionData = _arg2;
    this->Modified();
  }

  void SetCellData(size_t _arg1, float *_arg2) {
    this->CellDataSize = _arg1;
    this->CellDataData = _arg2;
    this->Modified();
  }

  void SetIndex(size_t _arg1, uint32_t *_arg2) {
    this->IndexSize = _arg1;
    this->IndexData = _arg2;
    this->Modified();
  }


  //--- Execute()...
private:
  int FrameNumber{0};
  OSPDevice Device{nullptr};
  OSPData VolumeVertexPositionData{nullptr};
  OSPData VolumeIndexData{nullptr};
  OSPData VolumeCellIndexData{nullptr};
  OSPData VolumeCellDataData{nullptr};
  OSPData VolumeCellTypeData{nullptr};
  OSPVolume Volume{nullptr};
  std::vector<float> TransferFunctionColor;
  OSPData TransferFunctionColorData{nullptr};
  std::vector<float> TransferFunctionOpacity;
  OSPData TransferFunctionOpacityData{nullptr};
  OSPTransferFunction TransferFunction{nullptr};
  OSPVolumetricModel VolumetricModel{nullptr};
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

public:
  char const *Execute();


  //--- For Finalize()...

  char const *Finalize();
};

vtkStandardNewMacro(OSPRayUnstructuredVolumeVisualization::InternalsType);

char const *OSPRayUnstructuredVolumeVisualization::InternalsType::Initialize() {
  int rank, size;
  MPI_Comm_rank(Communicator, &rank);
  MPI_Comm_size(Communicator, &size);

  ospLoadModule("mpi");

  Device = ospNewDevice("mpiDistributed");
  ospDeviceCommit(Device);
  ospSetCurrentDevice(Device);

  Group = ospNewGroup();
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
  ospSetVec3f(Camera, "position", 0.0f, 0.0f, 10.0f);
  ospSetVec3f(Camera, "direction", 0.0f, 0.0f, -1.0f);
  ospSetVec3f(Camera, "up", 0.0f, 1.0f, 0.0f);
  ospCommit(Camera);

  Renderer = ospNewRenderer("mpiRaycast");
  ospSetInt(Renderer, "pixelSamples", 16);
  ospSetVec3f(Renderer, "backgroundColor", 1.0f, 0.0f, 1.0f);
  ospCommit(Renderer);

  FrameBuffer = ospNewFrameBuffer(Width, Height, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);

  Future = nullptr;

  return /*error=*/nullptr;
}

char const *OSPRayUnstructuredVolumeVisualization::InternalsType::Execute() {
  int rank, size;
  MPI_Comm_rank(Communicator, &rank);
  MPI_Comm_size(Communicator, &size);

  float lo, hi;
  lo = hi = CellDataData[0];
  for (size_t i=1; i<CellDataSize; ++i) {
    if (CellDataData[i] < lo) lo = CellDataData[i];
    if (CellDataData[i] > hi) hi = CellDataData[i];
  }
  std::fprintf(stderr, "lo = %f, hi = %f\n", lo, hi);

# define RESET(m_Data)                                                         \
  do {                                                                         \
    if (m_Data != nullptr) {                                                   \
      ospRelease(m_Data);                                                      \
      m_Data = nullptr;                                                        \
    }                                                                          \
  } while (0)

  RESET(VolumeCellTypeData);
  VolumeCellTypeData = ospNewSharedData(
    CellTypeData, OSP_UCHAR, CellTypeSize, 0, 1, 0, 1, 0);
  ospCommit(VolumeCellTypeData);

  RESET(VolumeCellIndexData);
  VolumeCellIndexData = ospNewSharedData(
    CellIndexData, OSP_UINT, CellIndexSize, 0, 1, 0, 1, 0);
  ospCommit(VolumeCellIndexData);

  RESET(VolumeVertexPositionData);
  VolumeVertexPositionData = ospNewSharedData(
    VertexPositionData, OSP_VEC3F, VertexPositionSize, 0, 1, 0, 1, 0);
  ospCommit(VolumeVertexPositionData);

  RESET(VolumeCellDataData);
  VolumeCellDataData = ospNewSharedData(
    CellDataData, OSP_FLOAT, CellDataSize, 0, 1, 0, 1, 0);
  ospCommit(VolumeCellDataData);

  RESET(VolumeIndexData);
  VolumeIndexData = ospNewSharedData(
    IndexData, OSP_UINT, IndexSize, 0, 1, 0, 1, 0);
  ospCommit(VolumeIndexData);

# undef RESET

  Volume = ospNewVolume("unstructured");
  ospSetObject(Volume, "vertex.position", VolumeVertexPositionData);
  ospSetObject(Volume, "index", VolumeIndexData);
  ospSetBool(Volume, "indexPrefixed", false);
  ospSetObject(Volume, "cell.index", VolumeCellIndexData);
  ospSetObject(Volume, "cell.data", VolumeCellDataData);
  ospSetObject(Volume, "cell.type", VolumeCellTypeData);
  ospSetFloat(Volume, "background", 0.0f);
  ospCommit(Volume);

  if (TransferFunctionColorData) {
    ospRelease(TransferFunctionColorData);
    TransferFunctionColorData = nullptr;
  }

  TransferFunctionColor.clear();
  TransferFunctionColor.insert(TransferFunctionColor.end(), {
    (rank % 3 == 0 ? 1.0f : 0.0f),
    (rank % 3 == 1 ? 1.0f : 0.0f),
    (rank % 3 == 2 ? 1.0f : 0.0f),
    (rank % 3 == 0 ? 1.0f : 0.0f),
    (rank % 3 == 1 ? 1.0f : 0.0f),
    (rank % 3 == 2 ? 1.0f : 0.0f),
  });

  TransferFunctionColorData = ospNewSharedData(TransferFunctionColor.data(), OSP_VEC3F, TransferFunctionColor.size() / 3);
  ospCommit(TransferFunctionColorData);

  if (TransferFunctionOpacityData) {
    ospRelease(TransferFunctionOpacityData);
    TransferFunctionOpacityData = nullptr;
  }

  TransferFunctionOpacity.clear();
  TransferFunctionOpacity.insert(TransferFunctionOpacity.end(), {
    0.0f,
    1.0f,
  });

  TransferFunctionOpacityData = ospNewSharedData(TransferFunctionOpacity.data(), OSP_FLOAT, TransferFunctionOpacity.size() / 1);
  ospCommit(TransferFunctionOpacityData);

  TransferFunction = ospNewTransferFunction("piecewiseLinear");
  ospSetObject(TransferFunction, "color", TransferFunctionColorData);
  ospSetObject(TransferFunction, "opacity", TransferFunctionOpacityData);
  ospSetVec2f(TransferFunction, "valueRange", (float)lo, (float)hi);
  ospCommit(TransferFunction);

  VolumetricModel = ospNewVolumetricModel(nullptr);
  ospSetObject(VolumetricModel, "volume", Volume);
  ospSetObject(VolumetricModel, "transferFunction", TransferFunction);
  ospCommit(VolumetricModel);

  ospSetObjectAsData(Group, "volume", OSP_VOLUMETRIC_MODEL, VolumetricModel);
  ospCommit(Group);

  ospCommit(Instance);

  WorldRegion.clear();
  WorldRegionData = nullptr;

  WorldRegion.insert(WorldRegion.end(), {
    Bounds[0], Bounds[1], Bounds[2], // xmin, ymin, zmin
    Bounds[3], Bounds[4], Bounds[5], // xmax, ymax, zmax
  });

  WorldRegionData = ospNewSharedData(WorldRegion.data(), OSP_BOX3F, WorldRegion.size() / 6);
  ospCommit(WorldRegionData);

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

# define RESET(m_Size, m_Data)                                                 \
  do {                                                                         \
    m_Size = 0;                                                                \
    m_Data = nullptr;                                                          \
  } while (0)

  RESET(CellTypeSize, CellTypeData);
  RESET(CellIndexSize, CellIndexData);
  RESET(VertexPositionSize, VertexPositionData);
  RESET(CellDataSize, CellDataData);
  RESET(IndexSize, IndexData);

# undef RESET

  return /*error=*/nullptr;
}

char const *OSPRayUnstructuredVolumeVisualization::InternalsType::Finalize() {
# define RESET(m_Name)                                                         \
  do {                                                                         \
    ospRelease(m_Name);                                                        \
    m_Name = nullptr;                                                          \
  } while (0)

  RESET(FrameBuffer);
  RESET(Renderer);
  RESET(Camera);
  RESET(World);
  RESET(Light);
  RESET(Instance);
  RESET(Group);
  RESET(VolumetricModel);
  RESET(TransferFunction);
  RESET(TransferFunctionOpacityData);
  RESET(TransferFunctionColorData);
  RESET(Volume);
  RESET(VolumeVertexPositionData);
  RESET(VolumeIndexData);
  RESET(VolumeCellIndexData);
  RESET(VolumeCellDataData);
  RESET(VolumeCellTypeData);

# undef RESET

  return /*error=*/nullptr;
}


//---

OSPRayUnstructuredVolumeVisualization *OSPRayUnstructuredVolumeVisualization::New() {
  auto result = new OSPRayUnstructuredVolumeVisualization;
  result->InitializeObjectBase();
  return result;
}

int OSPRayUnstructuredVolumeVisualization::Initialize() {
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

bool OSPRayUnstructuredVolumeVisualization::Execute(sensei::DataAdaptor *dataAdaptor) {
  int rank, size;
  MPI_Comm_rank(GetCommunicator(), &rank);
  MPI_Comm_size(GetCommunicator(), &size);

  vtkDataObject *mesh;
  if (dataAdaptor->GetMesh(MeshName.c_str(), /*structureOnly=*/false, mesh)) {
    SENSEI_ERROR("Failed to get mesh '" << MeshName << "'")//no semicolon
    return /*success=*/false;
  }

  if (dataAdaptor->AddArray(mesh, MeshName.c_str(), vtkDataObject::CELL, ArrayName.c_str())) {
    SENSEI_ERROR("Failed to get array '" << ArrayName << "' from mesh '" << MeshName << "' block " << rank)//no semicolon
    return /*success=*/false;
  }

  using DataSet = vtkMultiBlockDataSet;
  DataSet *dataSet;
  if ((dataSet = DataSet::SafeDownCast(mesh)) == nullptr) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' to be a vtkMultiBlockDataSet")//no semicolon
    return /*success=*/false;
  }

  using DataObject = vtkDataObject;
  DataObject *dataObject;

  using Iterator = vtkCompositeDataIterator;
  Iterator *iterator = dataSet->NewIterator();
  iterator->SkipEmptyNodesOn();
  size_t numberOfItemsInIterator = 0;
  for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal(); iterator->GoToNextItem()) {
    dataObject = iterator->GetCurrentDataObject();
    ++numberOfItemsInIterator;
  }

  if (numberOfItemsInIterator != 1) {
    SENSEI_ERROR("Expected mesh '" << MeshName << "' to only have 1 item (actual = " << numberOfItemsInIterator << ")")//no semicolon
    return /*success=*/false;
  }

  using UnstructuredGrid = vtkUnstructuredGrid;
  UnstructuredGrid *unstructuredGrid;
  if ((unstructuredGrid = UnstructuredGrid::SafeDownCast(dataObject)) == nullptr) {
    SENSEI_ERROR("Expected data object from mesh '" << MeshName << "' to be a vtkUnstructuredGrid")//no semicolon
    return /*success=*/false;
  }

  unstructuredGrid->Register(this);

  if (UseD3) {
    MPI_Comm comm = GetCommunicator();

    using OpaqueComm = vtkMPICommunicatorOpaqueComm;
    OpaqueComm *opaque = new OpaqueComm(&comm);

    using Communicator = vtkMPICommunicator;
    vtkNew<Communicator> communicator;
    communicator->InitializeExternal(opaque);

    delete opaque;
    opaque = nullptr;

    using Controller = vtkMPIController;
    vtkNew<Controller> controller;
    controller->Initialize();
    controller->SetCommunicator(communicator);

    using Filter = vtkPDistributedDataFilter;
    vtkNew<Filter> filter;
    filter->GetKdtree()->AssignRegionsRoundRobin();
    filter->SetInputData(unstructuredGrid);
    filter->SetBoundaryMode(0);
    filter->SetUseMinimalMemory(1);
    filter->SetMinimumGhostLevel(0);
    filter->RetainKdtreeOn();

    filter->Update();

    unstructuredGrid->UnRegister(this);
    unstructuredGrid = UnstructuredGrid::SafeDownCast(filter->GetOutput());
    unstructuredGrid->Register(this);
  }

  this->Internals->SetCommunicator(GetCommunicator());

  {
    unstructuredGrid->GetPoints()->Modified();
    double boundsDouble[6] = { 0.373737f, 0.373737f, 0.373737f, 0.373737f, 0.373737f, 0.373737f }; // xmin, xmax, ymin, ymax, zmin, zmax
    unstructuredGrid->GetBounds(boundsDouble);
    float boundsFloat[6] = { 0.37373f, 0.373737f, 0.373737f, 0.373737f, 0.373737f, 0.373737f };
    boundsFloat[0] = static_cast<float>(boundsDouble[0]); // xmin
    boundsFloat[1] = static_cast<float>(boundsDouble[2]); // ymin
    boundsFloat[2] = static_cast<float>(boundsDouble[4]); // zmin
    boundsFloat[3] = static_cast<float>(boundsDouble[1]); // xmax
    boundsFloat[4] = static_cast<float>(boundsDouble[3]); // ymax
    boundsFloat[5] = static_cast<float>(boundsDouble[5]); // zmax
    this->Internals->SetBounds(boundsFloat);
  }

  std::vector<uint8_t> volumeCellType{};
  {
    using Array = vtkUnsignedCharArray;
    Array *array = unstructuredGrid->GetCellTypesArray();
    if (array == nullptr) {
      SENSEI_ERROR("Expected CellTypesArray from ugrid from mesh '" << MeshName << "' to be a vtkUnsignedCharArray")//no semicolon
      return /*success=*/false;
    }
    volumeCellType.resize(array->GetNumberOfValues(), 0);
    for (size_t i=0; i<array->GetNumberOfValues(); ++i) {
      volumeCellType[i] = array->GetValue(i);
    }
    this->Internals->SetCellType(
      array->GetNumberOfTuples(), volumeCellType.data());
  }

  std::vector<uint32_t> volumeCellIndex{};
  {
    using Array = vtkIdTypeArray;
    Array *array = unstructuredGrid->GetCellLocationsArray();
    if (array == nullptr) {
      SENSEI_ERROR("Expected CellLocationsArray from ugrid from mesh '" << MeshName << "' to be a vtkIdTypeArray")//no semicolon
      return /*success=*/false;
    }
    volumeCellIndex.resize(array->GetNumberOfValues(), 0);
    for (size_t i=0; i<array->GetNumberOfValues(); ++i) {
      volumeCellIndex[i] = array->GetValue(i);
    }
    this->Internals->SetCellIndex(
      array->GetNumberOfTuples(), volumeCellIndex.data());
  }

  std::vector<float> volumeVertexPosition{};
  {
    using Array = vtkDoubleArray;
    Array *array = Array::SafeDownCast(unstructuredGrid->GetPoints()->GetData());
    if (array == nullptr) {
      SENSEI_ERROR("Expected GetPoints->GetData() from ugrid from mesh '" << MeshName << "' to be a vtkDoubleArray")//no semicolon
      return /*success=*/false;
    }
    volumeVertexPosition.resize(array->GetNumberOfValues(), 0.0f);
    for (size_t i=0; i<array->GetNumberOfValues(); ++i) {
      volumeVertexPosition[i] = array->GetValue(i);
    }
    this->Internals->SetVertexPosition(
      array->GetNumberOfTuples(), volumeVertexPosition.data());
  }

  std::vector<float> volumeCellData{};
  {
    auto array1 = unstructuredGrid->GetCellData()->GetAbstractArray(ArrayName.c_str());
    if (array1 == nullptr) {
      SENSEI_ERROR("Expected CellData->GetScalars() from ugrid from mesh '" << MeshName << "' to exist")//no semicolon
      return /*success=*/false;
    }
    using Array = vtkUnsignedShortArray;
    Array *array = Array::SafeDownCast(array1);
    if (array == nullptr) {
      SENSEI_ERROR("Expected CellData->GetScalars() from ugrid from mesh '" << MeshName << "' to be a vtkUnsignedShortArray")//no semicolon
      return /*success=*/false;
    }
    volumeCellData.resize(array->GetNumberOfValues(), 0.0f);
    for (size_t i=0; i<array->GetNumberOfValues(); ++i) {
      volumeCellData[i] = array->GetValue(i);
    }
    this->Internals->SetCellData(
      array->GetNumberOfTuples(), volumeCellData.data());
  }

  std::vector<uint32_t> volumeIndex{};
  {
    using Array = vtkTypeInt64Array;
    Array *array = Array::SafeDownCast(unstructuredGrid->GetCells()->GetConnectivityArray());
    if (array == nullptr) {
      SENSEI_ERROR("Expected Cells->GetConnectivityArray() from ugrid from mesh '" << MeshName << "' to be a vtkTypeInt64Array")//no semicolon
      return /*success=*/false;
    }
    volumeIndex.resize(array->GetNumberOfValues(), 0.0f);
    for (size_t i=0; i<array->GetNumberOfValues(); ++i) {
      volumeIndex[i] = array->GetValue(i);
    }
    this->Internals->SetIndex(
      array->GetNumberOfTuples(), volumeIndex.data());
  }

  char const *error = this->Internals->Execute();
  if (error != nullptr) {
    std::fprintf(stderr, "Error: %s\n", error);
    return /*success=*/false;
  }

  unstructuredGrid->UnRegister(this);

  return /*success=*/true;
}

int OSPRayUnstructuredVolumeVisualization::Finalize() {
  char const *error = this->Internals->Finalize();
  if (error != nullptr) {
    std::fprintf(stderr, "Error: %s\n", error);
    return /*failure=*/1;
  }

  this->Internals->Delete();
  this->Internals = nullptr;

  return /*failure=*/0;
}


//---

} /* namespace ospSensei */
