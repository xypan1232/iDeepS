/* -*- mode:c++ -*- */
#ifndef GRAPHCLASS_H
#define GRAPHCLASS_H

#include "Utility.h"
#include "BaseGraphClass.h"

//CONSTANTS
/**
 Definition of offsets in the STATUS vector for vertices and their
 semantic meaning.
 */
const unsigned LABEL_VERTEX_ATTRIBUTE_ID = 0;
const unsigned KERNEL_POINT_ID = 0;
const unsigned KIND_ID = 1;
const unsigned VIEWPOINT_ID = 2;
const unsigned DEAD_ID = 3;
const unsigned ABSTRACTION_ID = 4;

/**
 Definition of offsets in the STATUS vector for edges and their
 semantic meaning.
 */
const unsigned EDGE_ATTRIBUTE_ID = 0;
const unsigned EDGE_DEAD_ID = 0;
const unsigned EDGE_ABSTRACTIONOF_ID = 1;
const unsigned EDGE_PARTOF_ID = 2;

///Adds semantics to the raw interface of BaseGraphClass.
///
///These concepts are introduced: ViewPoint, KernelPoint, Alive, Abstraction.
///
///Offers functions to compute pairwise distance information and to extract subgraphs.
class GraphClass: public BaseGraphClass {
	friend ostream& operator<<(ostream& out, const GraphClass& aG);
private:
public:
	GraphClass(const string& aGid = "temp") :
			mGraphID(aGid), mSliceIdList(1) {
	}
	string GetGraphID() const;
	void SetGraphID(string aGraphID);
	vector<unsigned> GetVertexAdjacentList(unsigned aID) const;
	bool IsSliced(void) const {
		return mSliceIdList.size() > 1;
	}
	ostream& Output(ostream& out) const;
	void ExportGraph(const string& aFilename, const string& aFormat) const;

	void SetSliceID(unsigned aID, const std::string& aSliceID);
	string GetSliceID(unsigned aID) const;
	void SetVertexKernelPoint(unsigned aID, bool state);
	bool GetVertexKernelPoint(unsigned aID) const;
	void SetVertexKind(unsigned aID, bool kind);
	bool GetVertexKind(unsigned aID) const;
	void SetVertexViewPoint(unsigned aID, bool state);
	bool GetVertexViewPoint(unsigned aID) const;
	void SetVertexAlive(unsigned aID, bool state);
	bool GetVertexAlive(unsigned aID) const;
	void SetVertexAbstraction(unsigned aID, bool state);
	bool GetVertexAbstraction(unsigned aID) const;
	void KillVertices(std::string s);
	void SetVertexDead(unsigned aID, bool state);
	bool GetVertexDead(unsigned aID) const;
	void SetVertexLabel(unsigned aID, const string& aLabel);
	string GetVertexLabel(unsigned aID) const;
	string GetVertexLabelConcatenated(unsigned aID) const;
	void SetEdgeLabel(unsigned aSrcID, unsigned aDestID, const string& aLabel);
	string GetEdgeLabel(unsigned aSrcID, unsigned aDestID) const;
	string GetEdgeLabel(unsigned aEdgeID) const;
	string GetEdgeLabelConcatenated(unsigned aSrcID, unsigned aDestID) const;
	string GetEdgeLabelConcatenated(unsigned aEdgeID) const;
	bool GetEdgeAbstractionOf(unsigned aSrcID, unsigned aDestID) const;
	bool GetEdgeAbstractionOf(unsigned aEdgeID) const;
	void SetEdgeAbstractionOf(unsigned aSrcID, unsigned aDestID, bool aStatus);
	void SetEdgeAbstractionOf(unsigned aEdgeID, bool aStatus);
	bool GetEdgePartOf(unsigned aSrcID, unsigned aDestID) const;
	bool GetEdgePartOf(unsigned aEdgeId) const;
	void SetEdgePartOf(unsigned aSrcID, unsigned aDestID, bool aStatus);
	void SetEdgePartOf(unsigned aEdgeID, bool aStatus);
	bool Check() const;
	void ComputePairwiseDistanceInformation(int aHorizon=-1, vector<unsigned> aViewPointList = vector<unsigned>()) const;
	vector<unsigned> GetFixedDistanceVertexIDList(unsigned aSrcID, int aDistance) const;
	int PairwiseDistance(unsigned aSrcID, unsigned aDestID) const;
	set<unsigned> GetUnionShortestPathsVertexIDList(unsigned aSrcID, unsigned aDestID, unsigned aDistance) const;
	set<unsigned> GetUnionShortestPathsVertexIDList(unsigned aSrcID, unsigned aDestID, unsigned aDistance, unsigned aRadius)const;
	set<unsigned> GetUnionThickShortestPathsVertexIDSet(unsigned aSrcID, unsigned aDestID, unsigned aDistance, unsigned aThickness) const;
	vector<unsigned> GetNeighborhoodVertexIDList(unsigned aSrcID, unsigned aDistance) const;

	void SaveAsDotFile(const string& aFilename) const;
	void SaveAsGMLFile(const string& aFilename) const;
	void SaveAsGDLFile(const string& aFilename) const;
	void SaveAsCSVFile(const string& aFilename) const;
	void SaveAsGspanFile(const string& aFilename) const;
	void SaveAsMatlabFile(const string& aFilename, vector<unsigned>& aViewPointList)const;
	//needed for  DDK kernel
    void SingleVertexBoundedBreadthFirstVisitTree(unsigned aRootVertexIndex, int aRadius, map<pair<unsigned, unsigned>, int>& oSrcDestMaptoDistance, map<unsigned, vector<unsigned> >& MapNodetoParents, map<unsigned, vector<unsigned> >& MapNodetoChildren, vector<unsigned>& TopologicalSort, int & maxLevel) const;

protected:
	void SingleVertexBoundedBreadthFirstVisit(unsigned aRootVertexIndex, int aRadius, map<pair<unsigned, unsigned>, int>& oSrcDestMaptoDistance) const;
	string vertex_label_serialize(unsigned v, string separator = "\n") const;
	string edge_label_serialize(unsigned edge_id) const;

protected:
	string mGraphID;
	vector<string> mSliceIdList;

	mutable unsigned mHorizon;
	mutable vector< vector< vector<unsigned> > > mSrcDistanceMaptoDestList; //semantic: <src_id,distance> \mapsto list of dest_id
	mutable map<pair<unsigned, unsigned>, int> mSrcDestMaptoDistance; //semantic: <src_id,dest_id> \mapsto distance
//	mutable map<pair<unsigned, int>, vector<unsigned> > mSrcDistanceMaptoDestListDebug; //semantic: <src_id,distance> \mapsto list of dest_id
};
#endif

