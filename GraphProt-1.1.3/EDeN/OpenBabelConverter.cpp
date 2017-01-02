#include "Utility.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"
#include "OpenBabelConverter.h"

using namespace std;

#ifdef USEOBABEL
void OpenBabelConverter::Convert(OpenBabel::OBMol& aObMol, GraphClass& oG) {
	//status
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	vector<bool> edge_status(3, false);
	edge_status[0] = false; //edge dead
	edge_status[1] = false; //edge abstraction_of
	edge_status[2] = false; //edge part_of

	map<string, int> index_map_nominal_to_real;
	for (OpenBabel::OBMolAtomIter a(aObMol); a; ++a) {
		string label = stream_cast<string>(a->GetAtomicNum());
		string nominal_vertex_index = stream_cast<string>(a->GetIdx());
		unsigned real_vertex_index = oG.InsertVertex();
		index_map_nominal_to_real[nominal_vertex_index] = real_vertex_index;
		vector<string> vertex_symbolic_attribute_list;
		vertex_symbolic_attribute_list.push_back(label);
		oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
		oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
	}
	for (OpenBabel::OBMolBondIter b(aObMol); b; ++b) {
		string nominal_src_index = stream_cast<string>(b->GetBeginAtomIdx());
		string nominal_dest_index = stream_cast<string>(b->GetEndAtomIdx());
		string label;
		if (b->IsAromatic()) {
			label = "a";
		} else if (b->IsSingle()) {
			label = "s";
		} else if (b->IsDouble()) {
			label = "d";
		} else if (b->IsTriple()) {
			label = "t";
		}
		assert(index_map_nominal_to_real.count(nominal_src_index)>0);
		assert(index_map_nominal_to_real.count(nominal_dest_index)>0);
		vector<string> edge_symbolic_attribute_list;
		edge_symbolic_attribute_list.push_back(label);
		unsigned real_src_index = index_map_nominal_to_real[nominal_src_index];
		unsigned real_dest_index = index_map_nominal_to_real[nominal_dest_index];
		unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
		oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
		oG.SetEdgeStatusAttributeList(edge_index, edge_status);
		unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
		oG.SetEdgeSymbolicAttributeList(reverse_edge_index, oG.GetEdgeSymbolicAttributeList(edge_index));
		oG.SetEdgeStatusAttributeList(reverse_edge_index, oG.GetEdgeStatusAttributeList(edge_index));
	}
}

void OpenBabelConverter::ConvertOpenBabelFormatToGraph(istream* pIn, GraphClass& oG, string aFormat) {
	OpenBabel::OBMol obMol;
	obconversion.SetInFormat(aFormat.c_str());
	bool success = obconversion.Read(&obMol, pIn);
	if (!success) return;
	Convert(obMol, oG);
}

void OpenBabelConverter::ConvertOpenBabelFormatToGraph(string aFileName, GraphClass& oG, string aFormat) {
	OpenBabel::OBMol obMol;
	if (mInitStatus==false){
		mInitStatus=true;
		obconversion.SetInFormat(aFormat.c_str());
		bool success = obconversion.ReadFile(&obMol, aFileName);
		if (!success) return;
		Convert(obMol, oG);
	} else {
		bool success = obconversion.Read(&obMol);
		if (!success) return;
		Convert(obMol, oG);
	}
}
#endif
