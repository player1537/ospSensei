/**
 *
 */

#ifndef __PLUGIN_SENSEI_OSPRAYSTUDIOVISUALIZATION_H__
#define __PLUGIN_SENSEI_OSPRAYSTUDIOVISUALIZATION_H__


// SENSEI
#include <AnalysisAdaptor.h>

// OSPRay Studio
#include "sg/Node.h"
#include "sg/Scheduler.h"


namespace ospray {
namespace plugin_sensei {


class OSPRayStudioVisualization : public sensei::AnalysisAdaptor {
public:
  static OSPRayStudioVisualization *New();

  senseiTypeMacro(OSPRayStudioVisualization, AnalysisAdaptor);

  int Initialize();

  void SetRootNode(const ospray::sg::NodePtr &root);
  void SetScheduler(const ospray::sg::SchedulerPtr &scheduler);

  bool Execute(sensei::DataAdaptor *) override;
  int Finalize() override;

private:
  bool Execute();

  struct InternalsType;
  InternalsType *Internals;
};


} /* namespace plugin_sensei */
} /* namespace ospray */


#endif /* __PLUGIN_SENSEI_OSPRAYSTUDIOVISUALIZATION_H__ */
// vim: ts=2:sts=2
