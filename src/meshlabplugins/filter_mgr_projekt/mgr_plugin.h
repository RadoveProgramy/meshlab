#ifndef FILTER_PLUGIN_H
#define FILTER_PLUGIN_H

#include <common/plugins/interfaces/filter_plugin.h>
#include <vcg/complex/algorithms/stat.h>
#include <QObject>

using namespace vcg;
using namespace std;

class Mgr_plugin : public QObject, public FilterPlugin
{
	Q_OBJECT
	MESHLAB_PLUGIN_IID_EXPORTER(FILTER_PLUGIN_IID)
	Q_INTERFACES(FilterPlugin)

private:

	static float globalMinZ;
	static float  globalMaxZ;

	static float CURVED_REGION_START;
	static float  CURVED_REGION_END;
	static float  CURVED_MAX_Z;  

public:
	enum { FP_FIRST, FP_SECOND, FP_THIRD };

	Mgr_plugin();
	QString pluginName() const;
	virtual QString filterName(ActionIDType filter) const;
	QString pythonFilterName(ActionIDType f) const;
	virtual QString filterInfo(ActionIDType filter) const;
	virtual RichParameterList initParameterList(const QAction* a, const MeshDocument& md);
	std::map<std::string, QVariant> applyFilter(
		const QAction* action,
		const RichParameterList& parameters,
		MeshDocument& md,
		unsigned int& postConditionMask,
		vcg::CallBackPos* cb);

	FilterArity         filterArity(const QAction* act) const { return FilterPlugin::NONE; }
	virtual FilterClass getClass(const QAction* a) const;
	QString             filterScriptFunctionName(ActionIDType filterID);
	virtual int         getRequirements(const QAction* filter);

	virtual int getPreConditions(const QAction* filter) const;
	virtual int postCondition(const QAction* filter) const;
};

#endif // FILTER_PLUGIN_H
