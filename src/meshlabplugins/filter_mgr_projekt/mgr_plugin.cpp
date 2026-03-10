#include "mgr_plugin.h"

#define VCG_USE_VERTEX_VFADJ
#define VCG_USE_FACE_VFADJ
#define VCG_USE_VERTEX_CURVATURE
#define VCG_USE_VERTEX_CURVATURE_DIR

#include <QFileDialog>
#include <QImage>
#include <QPainter>
#include <cmath>
#include <common/plugins/interfaces/filter_plugin.h>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <math.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/complex.h>
#include <vector>
#include <random>

using namespace std;
using namespace tri;
using namespace vcg;

// Staticke premenne na ulozenie minimalnej a maximalnej vysky modelu
float Mgr_plugin::globalMinZ = 0.f;
float Mgr_plugin::globalMaxZ = 1.f;

// Zakladne hodnoty pre zakriveny filter
float Mgr_plugin::CURVED_REGION_START = -0.5f;
float Mgr_plugin::CURVED_REGION_END   = 0.5f;
float Mgr_plugin::CURVED_MAX_Z        = 1.f;


// Struktura pre ulozenie shape spectrum deskriptora
struct ShapeSpectrumDescriptor
{
	std::vector<float> spectrumBins;
	float              planarAreaRatio   = 0.0f;
	float              singularAreaRatio = 0.0f;
};


// Pomocna funkcia na obmedzenie hodnoty v rozsahu [minVal, maxVal]
template<typename T>
T myClamp(T val, T minVal, T maxVal)
{
	return std::max(minVal, std::min(val, maxVal));
}

// Vypocet plochy trojuholnika
float computeFaceArea(const CMeshO::FaceType& f)
{
	return DoubleArea(f) / 2.0f;
}

// Vypocet priemernej zakrivenosti tvare na zaklade uhlov medzi normalami
float computeFaceMeanCurvature(
	const CMeshO&           mesh,
	const CMeshO::FaceType& f,
	FILE* debugFile = nullptr)
{
	Point3f n0       = f.N();
	float   accAngle = 0.0f;
	int     count    = 0;
	for (int i = 0; i < 3; ++i) {
		auto* adj = f.FFp(i);
		if (adj && adj != &f && !adj->IsD()) {
			accAngle += Angle(n0, adj->N());
			++count;
		}
	}
	float result = (count > 0) ? accAngle / count : 0.f;
	if (debugFile)
		fprintf(debugFile, "Mean curvature: %.6f\n", result);
	return result;
}



// Výpočet 3D Shape Spectrum deskriptora (v1) podľa MPEG-7.
ShapeSpectrumDescriptor computeShapeSpectrumDescriptor(const CMeshO& mesh, int numBins = 90)
{
	ShapeSpectrumDescriptor descriptor;
	descriptor.spectrumBins.resize(numBins, 0.0f);  // Inicializácia histogramu

	float       totalArea = 0.f, planarArea = 0.f, singularArea = 0.f;
	const float curvatureThreshold = 0.1f;

	// Aktualizácia topológie a výpočtu normál
	UpdateTopology<CMeshO>::FaceFace((CMeshO&) mesh);
	UpdateNormal<CMeshO>::PerFaceNormalized((CMeshO&) mesh);

	// Nepovinný debug výpis
	//FILE* debugFile = fopen("spectrum_debug.txt", "w");
	FILE* debugFile = NULL;
	if (debugFile)
		fprintf(debugFile, "Shape Spectrum Descriptor Debug Output\n");

	
	// Iterácia cez všetky tváre v modeli	
	for (const auto& f : mesh.face) {
		if (f.IsD())
			continue; // Preskočenie odstránených tvárí

		float area = computeFaceArea(f);
		totalArea += area;

		//zistenie zakrivenosti
		float curvature = computeFaceMeanCurvature(mesh, f, debugFile);
		if (debugFile)
			fprintf(debugFile, "Area: %.6f\n", area);

		if (curvature < curvatureThreshold) {
			planarArea += area;
			continue;
		}

		int borderCount = 0;
		for (int i = 0; i < 3; ++i) {
			auto* adj = f.FFp(i);
			if (!adj || adj == &f || adj->IsD())
				++borderCount;
		}
		if (borderCount > 0) {
			singularArea += area;
			continue;
		}

		// Výpočet shape indexu z curvature:
		float shapeIdx = 0.5f - (1.0f / float(M_PI)) * atan(curvature);
		shapeIdx       = myClamp(shapeIdx, 0.0f, 1.0f - 1e-6f);
		int binIdx     = std::min(int(shapeIdx * numBins), numBins - 1);
		descriptor.spectrumBins[binIdx] += area;
	}

	// Normalizácia histogramu
	if (totalArea > 0.f) {
		for (auto& val : descriptor.spectrumBins)
			val /= totalArea;
		descriptor.planarAreaRatio   = planarArea / totalArea;
		descriptor.singularAreaRatio = singularArea / totalArea;
	}

	// Výpis pre ladenie
	if (debugFile) {
		fprintf(debugFile, "\nNormalized Spectrum:\n");
		for (size_t i = 0; i < descriptor.spectrumBins.size(); ++i)
			fprintf(debugFile, "Bin %2zu: %.6f\n", i, descriptor.spectrumBins[i]);
		fprintf(
			debugFile,
			"Planar Ratio: %.4f\nSingular Ratio: %.4f\n",
			descriptor.planarAreaRatio,
			descriptor.singularAreaRatio);
		fclose(debugFile);
	}

	return descriptor;
}


// Výpočet 3D Shape Spectrum deskriptora (v2) podľa MPEG-7.
ShapeSpectrumDescriptor computeShapeSpectrumFromCurvature(CMeshO& mesh, int numBins = 90)
{
	ShapeSpectrumDescriptor descriptor;
	descriptor.spectrumBins.resize(numBins, 0.0f); // Inicializácia histogramu

	float totalArea    = 0.f;
	float planarArea   = 0.f;
	float singularArea = 0.f;

	const float curvatureThreshold = 0.05f;
	const int   minBorderEdges     = 2;

	// Nepovinný debug výstup do súboru
	//FILE* debug = fopen("spectrum_debug_v2.txt", "w");
	FILE* debug = NULL;
	if (debug)
		fprintf(debug, "== MPEG-7 Shape Spectrum (v2) Debug ==\n");


	// Iterácia cez všetky tváre
	for (const auto& f : mesh.face) {
		if (f.IsD())
			continue;

		float area = computeFaceArea(f);
		totalArea += area;

		// Spočítanie počtu hraničných/degenerovaných susedov
		int borderCount = 0;
		for (int i = 0; i < 3; ++i) {
			auto* adj = f.FFp(i);
			if (!adj || adj == &f || adj->IsD())
				++borderCount;
		}

		if (borderCount >= minBorderEdges) {
			singularArea += area;
			if (debug)
				fprintf(debug, "Skipped face (border), area=%.6f\n", area);
			continue;
		}

		// Vypočítaj  hodnoty K1 a K2 z vrcholov tváre
		float k1 = 0.f, k2 = 0.f;
		int   valid = 0;
		for (int i = 0; i < 3; ++i) {
			const auto* v = f.V(i);
			if (!v->IsD() && std::isfinite(v->K1()) && std::isfinite(v->K2())) {
				k1 += v->K1();
				k2 += v->K2();
				++valid;
			}
		}

		if (valid == 0) {
			if (debug)
				fprintf(debug, "Skipped face (no valid curvatures)\n");
			continue;
		}

		k1 /= valid;
		k2 /= valid;

		// Kombinovaná zakrivenosť (pre určenie "plochosti")
		float ka = std::sqrt(k1 * k1 + k2 * k2);
		if (ka < curvatureThreshold) {
			planarArea += area;
			if (debug)
				fprintf(debug, "Skipped face (planar), ka=%.6f, area=%.6f\n", ka, area);
			continue;
		}

		// Výpočet shape indexu podľa MPEG-7
		float shapeIdx = 0.5f - (1.0f / float(M_PI)) * std::atan((k1 + k2) / (k1 - k2 + 1e-6f));
		shapeIdx       = myClamp(shapeIdx, 0.0f, 1.0f - 1e-6f);

		// Pridanie plochy do príslušného histogramového koša
		int binIdx     = std::min(int(shapeIdx * numBins), numBins - 1);
		descriptor.spectrumBins[binIdx] += area;

		if (debug)
			fprintf(
				debug,
				"face: k1=%.6f, k2=%.6f, idx=%.3f, bin=%d, area=%.6f\n",
				k1,
				k2,
				shapeIdx,
				binIdx,
				area);
	}

	// Normalizácia histogramu
	if (totalArea > 0.f) {
		for (auto& val : descriptor.spectrumBins)
			val /= totalArea;
		descriptor.planarAreaRatio   = planarArea / totalArea;
		descriptor.singularAreaRatio = singularArea / totalArea;
	}

	// Debug výpis normalizovaného histogramu
	if (debug) {
		fprintf(debug, "\n== Normalized Histogram ==\n");
		for (size_t i = 0; i < descriptor.spectrumBins.size(); ++i)
			fprintf(debug, "bin %02zu: %.6f\n", i, descriptor.spectrumBins[i]);
		fprintf(
			debug,
			"Planar: %.4f, Singular: %.4f\n",
			descriptor.planarAreaRatio,
			descriptor.singularAreaRatio);
		fclose(debug);
	}

	return descriptor;
}



float roundUpToNice(float val)
{
	float scale = pow(10.0f, floor(log10(val)));
	float n     = val / scale;
	if (n <= 1.0f)
		return 1.0f * scale;
	else if (n <= 2.0f)
		return 2.0f * scale;
	else if (n <= 5.0f)
		return 5.0f * scale;
	else
		return 10.0f * scale;
}



// Funkcia na vykreslenie a uloženie histogramu (napr. shape spectrum) ako obrázka.
// Výstupom je .bmp alebo .png obrázok s osami, popismi a krivkou histogramu.
void saveSpectrumHistogramImage(const std::vector<float>& data, const std::string& filename)
{
	const int width  = 800;  // šírka obrázka v pixeloch
	const int height = 400;  // výška obrázka v pixeloch
	const int margin = 50;   // okraj pre osi a popisy

	QImage img(width, height, QImage::Format_ARGB32);
	img.fill(Qt::white);  // výplň pozadia na bielo

	QPainter painter(&img);
	painter.setRenderHint(QPainter::Antialiasing);
	painter.setPen(QPen(Qt::black, 2));
	QFont font = painter.font();
	font.setPointSize(8);
	painter.setFont(font);

	// Určenie maximálnej hodnoty v histograme (pre škálovanie Y osi)
	float peakVal = *std::max_element(data.begin(), data.end());
	float maxVal  = roundUpToNice(peakVal * 1.5f);

	// počet binov (stĺpcov histogramu)
	const int   bins   = static_cast<int>(data.size());
	const float xScale = float(width - 2 * margin) / bins;
	const float yScale = float(height - 2 * margin) / maxVal;

	// Kreslenie osí 
	painter.drawLine(margin, height - margin, margin, margin);                  // Y
	painter.drawLine(margin, height - margin, width - margin, height - margin); // X

	// Kreslenie pomocných čiar a popisov Y osi
	const int yTicks = 5;
	for (int i = 0; i <= yTicks; ++i) {
		float yValue = i * maxVal / yTicks;
		int   y      = height - margin - int(yValue * yScale);
		painter.drawLine(margin - 5, y, margin + 5, y);
		QString label = QString::number(yValue, 'f', 2);
		painter.drawText(margin - 40, y + 4, label);
	}

	// Vykreslenie histogramovej krivky
	painter.setPen(QPen(Qt::blue, 2));
	for (int i = 1; i < bins; ++i) {
		int x1 = margin + int((i - 1) * xScale);
		int y1 = height - margin - int(data[i - 1] * yScale);
		int x2 = margin + int(i * xScale);
		int y2 = height - margin - int(data[i] * yScale);
		painter.drawLine(x1, y1, x2, y2);
	}

	// Popisy X osi (indexy binov) 
	const int xLabelStep = std::max(1, bins / 10);
	painter.setPen(Qt::black);
	for (int i = 0; i < bins; i += xLabelStep) {
		int x = margin + int(i * xScale);
		painter.drawLine(x, height - margin - 5, x, height - margin + 5);
		QString label = QString::number(i);
		painter.drawText(x - 10, height - margin + 20, label);
	}

	// Uloženie obrázka 
	bool ok = img.save(QString::fromStdString(filename));


	//FILE* debugFile = fopen("save_histogram_debug.txt", "w");
	FILE* debugFile = NULL;
	if (debugFile) {
		fprintf(debugFile, "Saved to: %s\nSuccess: %s\n", filename.c_str(), ok ? "true" : "false");
		fclose(debugFile);
	}

}


// Vypočítanie D2 histogram tvaru podľa Osada et al. (2002).
std::vector<float> computeD2Histogram(const CMeshO& mesh, int numSamples = 100000, int numBins = 90)
{
	std::vector<float> hist(numBins, 0.0f);
	if (mesh.vert.empty())
		return hist; // ak mesh nemá žiadne vertexy, vrátime prázdny histogram


	// Získanie platných (neodstránených) vertexov
	std::vector<const CVertexO*> validVerts;
	for (const auto& v : mesh.vert)
		if (!v.IsD())
			validVerts.push_back(&v);


	if (validVerts.size() < 2)
		return hist; // málo bodov na porovnanie


	// Príprava na sampling
	std::random_device rd;
	std::mt19937 rng(rd());
	std::uniform_int_distribution<size_t> dist(0, validVerts.size() - 1);


	float              maxDist = 0.f;
	std::vector<float> distances;
	distances.reserve(numSamples);

	// Vzorkovanie náhodných dvojíc bodov a meranie vzdialeností 
	for (int i = 0; i < numSamples; ++i) {
		const auto* v1 = validVerts[dist(rng)];
		const auto* v2 = validVerts[dist(rng)];

		if (v1 == v2) 
			continue; // preskočnie identického bodu

		float d = Distance(v1->cP(), v2->cP());
		distances.push_back(d);
		maxDist = std::max(maxDist, d);
	}

	// naplnenie histogramu
	for (float d : distances) {
		int bin = int((d / maxDist) * numBins);
		if (bin >= numBins)
			bin = numBins - 1;
		hist[bin] += 1.f;
	}

	// Normalizácia
	for (float& val : hist)
		val /= distances.size();

	return hist;
}



Mgr_plugin::Mgr_plugin()
{
	typeList = {FP_FIRST, FP_SECOND, FP_THIRD};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}



QString Mgr_plugin::pluginName() const
{
	return "Mgr project plugin";
}



QString Mgr_plugin::filterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_FIRST: return QString("Mgr plugin - create mesh");
	case FP_SECOND: return QString("Mgr plugin - filter");
	case FP_THIRD: return QString("Mgr plugin - descriptor");
	default: assert(0); return QString();
	}
}



QString Mgr_plugin::pythonFilterName(ActionIDType f) const
{
	switch (f) {
	case FP_FIRST: return QString("");
	case FP_SECOND: return QString("");
	case FP_THIRD: return QString("");
	default: assert(0); return QString();
	}
}



QString Mgr_plugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_FIRST:
		return tr(
			"This function allows users to create 3d model from depthmap."
			"<br>");
	case FP_SECOND:
		return tr("Filter / push vertices based on height from flat or curved reference.");
	case FP_THIRD: return tr("Export descriptor to file.");
	default: assert(0);
	}
	return "";
}



RichParameterList Mgr_plugin::initParameterList(const QAction* a, const MeshDocument& md)
{
	RichParameterList parlst;
	switch (ID(a)) {
	case FP_FIRST:
		parlst.addParam(RichFileOpen(
			"depthMap",
			"",
			QStringList() << "*.png " << "*.jpg",
			"Load depth map",
			"Choose a depth map image"));
		parlst.addParam(RichDynamicFloat(
			"scaleXY", 0.01f, 0.000001f, 0.3f, "scaleXY", "Scaling for X/Y coordinates"));
		parlst.addParam(
			RichDynamicFloat("scaleZ", 1.0f, 0.01f, 5.0f, "scaleZ", "Scaling for Z (depth)"));
		break;
	case FP_SECOND: {
		float avgZ = (Mgr_plugin::globalMinZ + Mgr_plugin::globalMaxZ) / 2.0f;
		parlst.addParam(RichDynamicFloat(
			"thresholdZ",
			avgZ,
			globalMinZ,
			globalMaxZ,
			"Min height",
			"Minimum relative Z to keep"));

		QStringList modes;
		modes << "flat" << "curved" << "flatten edges";
		parlst.addParam(RichEnum("referencePlane", 0, modes, "Reference", "Base reference shape"));
		break;
	}
	case FP_THIRD: {
		QStringList options;
		options << "Shape Spectrum - MPEG 7 v1 " << "D2 Histogram";
		//		<< "Shape Spectrum - MPEG 7 v2";
		parlst.addParam(
			RichEnum("descriptorType", 0, options, "Descriptor", "Select descriptor to export"));
		parlst.addParam(
			RichFileSave("outputFile", "", "*.bmp", "Save descriptor as"));
		break;
	}

	default: break;
	}
	return parlst;
}



FilterPlugin::FilterClass Mgr_plugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_FIRST: return FilterPlugin::Other;
	case FP_SECOND: return FilterPlugin::Other;
	case FP_THIRD: return FilterPlugin::Other;
	default: assert(0); return FilterPlugin::Generic;
	}
}



QString Mgr_plugin::filterScriptFunctionName(ActionIDType filterID)
{
	switch (filterID) {
	case FP_FIRST: return QString("First option");
	case FP_SECOND: return QString("Second option");
	case FP_THIRD: return QString("Third option");
	default: assert(0);
	}
	return NULL;
}



std::map<std::string, QVariant> Mgr_plugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	CallBackPos*)
{
	if (ID(filter) == FP_FIRST) {

		// Načítanie cesty k obrázku z parametrov
		QString fileName = par.getString("depthMap");
		if (fileName.isEmpty())
			return {};
		
		// Získanie mierky pre XY a Z os
		float scaleXY = par.getDynamicFloat("scaleXY");
		float scaleZ  = par.getDynamicFloat("scaleZ");

		// Načítanie obrázku / hĺbkovej mapy
		QImage img(fileName);
		if (img.isNull())
			return {};

		int width  = img.width();
		int height = img.height();

		// Vytvorenie nového mesh modelu
		MeshModel* mm   = md.addNewMesh("MeshFromMap", "Generated from depth");
		CMeshO&    mesh = mm->cm;

		// Maticová štruktúra na ukladanie indexov vertexov
		std::vector<std::vector<int>> grid(height, std::vector<int>(width));
		int index = 0;

		// Generovanie vertexov z obrázka (jas → výška Z)
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				QColor color = img.pixelColor(x, y);
				double gray  = 0.299 * color.red() + 0.587 * color.green() + 0.114 * color.blue();
				float  z     = float(gray / 255.0f) * scaleZ;

				CVertexO v;
				v.P().X() = x * scaleXY;
				v.P().Y() = y * scaleXY;
				v.P().Z() = z;

				mesh.vert.push_back(v);
				grid[y][x] = index++;
			}
		}
		mesh.vn = mesh.vert.size(); // aktualizácia počtu vertexov

		// Výpočet výšok
		float minZ_ = std::numeric_limits<float>::max();
		float maxZ_ = std::numeric_limits<float>::lowest();

		// Výpočet min a max výšky (Z)
		for (const auto& v : mesh.vert) {
			minZ_ = std::min(minZ_, v.P().Z());
			maxZ_ = std::max(maxZ_, v.P().Z());
		}

		// Uloženie do statických premenných
		Mgr_plugin::globalMinZ = minZ_;
		Mgr_plugin::globalMaxZ = maxZ_;

		// Centrovanie na stred (X, Y)
		float centerX = (width - 1) * 0.5f * scaleXY;
		float centerY = (height - 1) * 0.5f * scaleXY;
		for (auto& v : mesh.vert) {
			v.P().X() -= centerX;
			v.P().Y() -= centerY;
		}

		// Vytvorenie trojuholníkov z gridu – výber diagonály podľa vzdialenosti
		for (int y = 0; y < height - 1; ++y) {
			for (int x = 0; x < width - 1; ++x) {
				int i0 = grid[y][x];
				int i1 = grid[y][x + 1];
				int i2 = grid[y + 1][x];
				int i3 = grid[y + 1][x + 1];

				float d1 = Distance(mesh.vert[i0].P(), mesh.vert[i3].P()); // v0–v3
				float d2 = Distance(mesh.vert[i1].P(), mesh.vert[i2].P()); // v1–v2


				if (d1 < d2) {
					// diagonála v0–v3
					CFaceO f1, f2;
					f1.V(0) = &mesh.vert[i0];
					f1.V(1) = &mesh.vert[i1];
					f1.V(2) = &mesh.vert[i3];
					mesh.face.push_back(f1);
					mesh.fn++;

					f2.V(0) = &mesh.vert[i0];
					f2.V(1) = &mesh.vert[i3];
					f2.V(2) = &mesh.vert[i2];
					mesh.face.push_back(f2);
					mesh.fn++;
				}
				else {
					// diagonála v1–v2
					CFaceO f1, f2;
					f1.V(0) = &mesh.vert[i0];
					f1.V(1) = &mesh.vert[i1];
					f1.V(2) = &mesh.vert[i2];
					mesh.face.push_back(f1);
					mesh.fn++;

					f2.V(0) = &mesh.vert[i2];
					f2.V(1) = &mesh.vert[i1];
					f2.V(2) = &mesh.vert[i3];
					mesh.face.push_back(f2);
					mesh.fn++;
				}
			}
		}

		// Vytvorenie spodku modelu
		int topLeft     = grid[0][0];
		int topRight    = grid[0][width - 1];
		int bottomLeft  = grid[height - 1][0];
		int bottomRight = grid[height - 1][width - 1];

		float minZ = std::numeric_limits<float>::max();
		for (const auto& v : mesh.vert)
			minZ = std::min(minZ, v.P().Z());

		// Vytvorenie 4 nových vertexov so Z = minZ
		int i0 = mesh.vert.size();
		mesh.vert.emplace_back(mesh.vert[topLeft]);
		mesh.vert.back().P().Z() = minZ;
		int i1                   = mesh.vert.size();
		mesh.vert.emplace_back(mesh.vert[topRight]);
		mesh.vert.back().P().Z() = minZ;
		int i2                   = mesh.vert.size();
		mesh.vert.emplace_back(mesh.vert[bottomLeft]);
		mesh.vert.back().P().Z() = minZ;
		int i3                   = mesh.vert.size();
		mesh.vert.emplace_back(mesh.vert[bottomRight]);
		mesh.vert.back().P().Z() = minZ;
		mesh.vn                  = mesh.vert.size();

		// Dva trojuholníky na základni
		CFaceO fb1, fb2;
		fb1.V(0) = &mesh.vert[i0];
		fb1.V(1) = &mesh.vert[i1];
		fb1.V(2) = &mesh.vert[i2];
		mesh.face.push_back(fb1);

		fb2.V(0) = &mesh.vert[i2];
		fb2.V(1) = &mesh.vert[i1];
		fb2.V(2) = &mesh.vert[i3];
		mesh.face.push_back(fb2);
		mesh.fn += 2;

		// Nastavenie curved pásma a maxZ podľa centrovaného bboxu
		UpdateBounding<CMeshO>::Box(mesh);
		float minY__        = mesh.bbox.min.Y();
		float maxY__        = mesh.bbox.max.Y();
		float minZ__        = mesh.bbox.min.Z();
		float maxZ__        = mesh.bbox.max.Z();
		CURVED_REGION_START = minY__ + 0.105f * (maxY__ - minY__); 
		CURVED_REGION_END   = minY__ + 0.775f * (maxY__ - minY__); 
		CURVED_MAX_Z        = maxZ__;

		// Normály a box
		UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
		UpdateNormal<CMeshO>::PerVertexNormalized(mesh);
		// UpdateNormal<CMeshO>::PerFace(mesh);
		// UpdateNormal<CMeshO>::PerVertex(mesh);
		// UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh);
	}

	else if (ID(filter) == FP_SECOND) {
		// Získanie aktuálneho modelu a meshu
		MeshModel* mm   = md.mm();
		CMeshO&    mesh = mm->cm;

		// Načítanie výškového prahu (Z)
		float thresholdZ = par.getDynamicFloat("thresholdZ");

		// Zistenie režimu filtrácie
		QStringList modes;
		modes << "flat" << "curved" << "flatten edges";
		QString mode = modes[par.getEnum("referencePlane")];

		// Aktualizácia bounding boxu
		UpdateBounding<CMeshO>::Box(mesh);

		int   movedCount = 0;  // Počet posunutých vertexov
		float maxShift   = 0.f; // Maximálna zmena výšky

		//FILE* logFile = fopen("debug_output.txt", "w");
		FILE* logFile = NULL;
		if (!logFile) {
			qWarning("Unable to open debug_output.txt");
		}
		else {
			fprintf(logFile, "=== FP_SECOND DEBUG LOG ===\n");
			fprintf(logFile, "Mode: %s\n", mode.toUtf8().constData());
			fprintf(logFile, "ThresholdZ: %.4f\n", thresholdZ);
			fprintf(
				logFile,
				"Bounding box X range: %.2f .. %.2f (DimX=%.2f)\n",
				mesh.bbox.min.X(),
				mesh.bbox.max.X(),
				mesh.bbox.DimX());
		}

		// Režim 1: flat
		if (mode == "flat") {
			// Všetky vertexy pod prahom sa posunú na Z = 0
			for (auto& v : mesh.vert) {
				if (v.P().Z() < thresholdZ) {
					if (logFile)
						fprintf(logFile, "Vertex moved (flat): Z=%.4f -> Z=0.0\n", v.P().Z());
					maxShift  = std::max(maxShift, thresholdZ - v.P().Z());
					v.P().Z() = 0.f;
					movedCount++;
				}
			}
		}

		// Režim 2: curved
		else if (mode == "curved") {
			float flatStartY  = CURVED_REGION_START;
			float flatEndY    = CURVED_REGION_END;
			float minY        = mesh.bbox.min.Y();
			float maxY        = mesh.bbox.max.Y();
			float flatReturnY = minY + 0.835f * (maxY - minY);

			int movedCount = 0;

			for (size_t i = 0; i < mesh.vert.size(); ++i) {
				auto& v = mesh.vert[i];
				if (v.IsD())
					continue;

				// Preskočenie spodných 4 vertexov (súčasť základne)
				if (i >= mesh.vert.size() - 4)
					continue;

				float y     = v.P().Y();
				float baseZ = 0.f;

				// Oblasti mimo zakrivenia sa posunú na Z=0
				if (y > flatReturnY) {
					if (v.P().Z() != 0.f) {
						v.P().Z() = 0.f;
						movedCount++;
					}
					continue;
				}

				// Výpočet zakrivenej výšky základne pre vertex podľa jeho pozície v osi Y
				if (y < flatStartY) {
					float dy = (y - minY) / (flatStartY - minY + 1e-6f);
					baseZ    = CURVED_MAX_Z * (1.f - dy) * (1.f - dy);
				}
				else if (y > flatEndY) {
					float dy = (y - flatEndY) / (flatReturnY - flatEndY + 1e-6f);
					baseZ    = CURVED_MAX_Z * dy * dy;
				}
				else {
					baseZ = 0.f;
				}

				// Posunutie vertexu, ak je pod základňou + threshold
				if ((v.P().Z() - baseZ) < thresholdZ) {
					v.P().Z() = baseZ;
					movedCount++;
				}
			}

			// Prepočet normál a bounding boxu po úpravách
			UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
			UpdateNormal<CMeshO>::PerVertexNormalized(mesh);
			UpdateBounding<CMeshO>::Box(mesh);

			qDebug("FP_SECOND (curved): moved %d vertices", movedCount);
		}

		// Režim 3: flatten edges
		else if (mode == "flatten edges") {
			// Okraje mimo zakriveného pásma sa stiahnu na Z = 0
			for (auto& v : mesh.vert) {
				float y = v.P().Y();
				if (y < CURVED_REGION_START || y > CURVED_REGION_END) {
					if (logFile) {
						fprintf(
							logFile,
							"Flattening edge vertex: X=%.2f Z=%.4f -> Z=0.0\n",
							y,
							v.P().Z());
					}
					maxShift  = std::max(maxShift, std::abs(v.P().Z()));
					v.P().Z() = 0.f;
					movedCount++;
				}
			}
		}

		// Uzavretie log súboru
		if (logFile) {
			fprintf(logFile, "Total moved vertices: %d\n", movedCount);
			fprintf(logFile, "Max shift: %.4f\n", maxShift);
			fclose(logFile);
		}

		// Aktualizácia normál a bounding boxu po všetkých zmenách
		UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
		UpdateNormal<CMeshO>::PerVertexNormalized(mesh);
		UpdateBounding<CMeshO>::Box(mesh);


		qDebug("FP_SECOND: moved %d vertices (max shift %.4f). Log saved.", movedCount, maxShift);

	}

	else if (ID(filter) == FP_THIRD) {
		// Získaj aktuálny mesh
		CMeshO& mesh = md.mm()->cm;

		// Otvorenie log súboru (momentálne deaktivované)
		//FILE* logFile1 = fopen("debug_output1.txt", "w");
		FILE* logFile1 = NULL;
		if (!logFile1) {
			qWarning("Unable to open debug_output1.txt");
		}

		// Získanie cieľového súboru pre uloženie histogramu
		QString outputFile = par.getString("outputFile");

		// Výber typu deskriptora podľa vstupného parametra
		QStringList modes;
		modes << "Shape Spectrum - MPEG 7 v1" << "D2 Histogram" << "Shape Spectrum - MPEG 7 v2";
		QString descType = modes[par.getEnum("descriptorType")];

		// Logovanie výberu 
		if (logFile1) {
			fprintf(logFile1, "%s\n", outputFile.toUtf8().constData());
			fprintf(logFile1, "%s\n", descType.toUtf8().constData());
		}

		// Vypočítanie a uloženie D2 histogramu
		if (descType == "D2 Histogram") {
			auto hist = computeD2Histogram(mesh);
			saveSpectrumHistogramImage(hist, outputFile.toStdString());
			qDebug("D2 Histogram saved to %s", qUtf8Printable(outputFile));
		}

		// Vypočítanie a uloženie Shape Spectrum podľa MPEG-7 v1
		else if (descType == "Shape Spectrum - MPEG 7 v1") {
			auto spectrum = computeShapeSpectrumDescriptor(mesh); 
			saveSpectrumHistogramImage(
				spectrum.spectrumBins, outputFile.toStdString());
			qDebug("Shape Spectrum saved to %s", qUtf8Printable(outputFile));
		}

		// Vypočítanie a uloženie Shape Spectrum podľa MPEG-7 v2
		else if (descType == "Shape Spectrum - MPEG 7 v2") {

			UpdateTopology<CMeshO>::FaceFace(mesh);
			UpdateTopology<CMeshO>::VertexFace(mesh);
			UpdateNormal<CMeshO>::PerFaceNormalized(mesh);
			UpdateNormal<CMeshO>::PerVertexNormalized(mesh);
			UpdateCurvature<CMeshO>::PrincipalDirections(mesh);

			// Výpočet deskriptora
			auto spectrum = computeShapeSpectrumFromCurvature(mesh);
			

			saveSpectrumHistogramImage(spectrum.spectrumBins, outputFile.toStdString());
			qDebug("Shape Spectrum (v2) saved to %s", qUtf8Printable(outputFile));
		}

		if (logFile1) {
			fclose(logFile1);
		}

	}

	else {
		wrongActionCalled(filter);
	}
	return {};
}



int Mgr_plugin::postCondition(const QAction* filter) const
{
	switch (ID(filter)) {
	case FP_FIRST: return MeshModel::MM_NONE;
	case FP_SECOND: return MeshModel::MM_VERTCOORD;
	case FP_THIRD: return MeshModel::MM_NONE;
	default: assert(0);
	}
	return MeshModel::MM_NONE;
}



int Mgr_plugin::getRequirements(const QAction* filter){
	switch (ID(filter)) {
	case FP_FIRST: return MeshModel::MM_VERTCOORD;
	case FP_SECOND:
		return MeshModel::MM_VERTCOORD | MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACENORMAL |
			   MeshModel::MM_FACEFLAG;
	case FP_THIRD:
		return MeshModel::MM_VERTCOORD | MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACENORMAL |
			   MeshModel::MM_FACEFLAG | MeshModel::MM_VERTNORMAL | MeshModel::MM_WEDGNORMAL |
			   MeshModel::MM_VERTFLAG | MeshModel::MM_VERTQUALITY | MeshModel::MM_VERTCURV |
			   MeshModel::MM_VERTCURVDIR;
	default: assert(0); return MeshModel::MM_NONE;
	}
}



int Mgr_plugin::getPreConditions(const QAction* filter) const{
	switch (ID(filter)) {
    case FP_FIRST:
	case FP_SECOND:
	case FP_THIRD: return MeshModel::MM_FACEFACETOPO | MeshModel::MM_FACENORMAL;
	default: return MeshModel::MM_NONE;
	}
}



MESHLAB_PLUGIN_NAME_EXPORTER(Mgr_plugin)
