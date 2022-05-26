/**
 *
 */

#ifndef __OSPSENSEI_OSPRAYVISUALIZATION_H__
#define __OSPSENSEI_OSPRAYVISUALIZATION_H__


// sensei
#include <AnalysisAdaptor.h>


namespace ospSensei {

class OSPRayVisualization : public sensei::AnalysisAdaptor {
public:
  enum Type {
    PARTICLE,
    VOLUME,
  };

  static OSPRayVisualization *New();

  senseiTypeMacro(OSPRayVisualization, AnalysisAdaptor);

  void SetMode(Type _arg) {
    Mode = _arg;
  }

  int Initialize();

  void SetMeshName(const char *_arg) {
    MeshName = _arg;
  }

  void SetArrayName(const char *_arg) {
    ArrayName = _arg;
  }

  bool Execute(sensei::DataAdaptor *) override;
  int Finalize() override;

private:
  Type Mode;
  std::string MeshName;
  std::string ArrayName;

  struct InternalsType;
  InternalsType *Internals;
};


} /* namespace ospSensei */


#endif /* __OSPSENSEI_OSPRAYVISUALIZATION_H__ */
// vim: ts=2:sts=2
