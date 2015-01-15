
/* *************************************************************************************************

http://www.voronoi.com/wiki/index.php?title=Spatial_Data_Structures

The quad-edge graph is a structure for representing the topology of a well-formed manifold mesh. No
actual geometry is represented, only the connectivity of faces and vertices. Faces with any number
of sides are supported naturally.

Every edge of a well-formed, closed mesh is adjacent to two faces and two vertices. The quad-edge
structure expresses this connectivity as a group of four directed edges, two opposed directed edges
between vertices and two between faces (an edge of the dual mesh). Additionally, each directed edge
has a pointer to the next quad-edge around its origin, whether that is a vertex or a face.

Stated another way, the quad-edge graph represents the set of directed loops of faces around each
vertex, and the directed loops of vertices around each face. For each edge of the mesh, two vertex
loops and two face loops intersect in a quad-edge that expresses the connectivity of the vertices
and faces and allows for organized, efficient traversal across the mesh.

There is no distinction in how faces and vertices are treated, the two are interchangeable as far as
the quadedge structure is concerned (it encodes the primal and dual graph of the mesh connectivity
equally). The assignment of one set of indices to faces and another to vertices is arbitrary. The
original paper describing quad-edges defined an algebra on them. This implementation takes a
different approach, with the concept of quad-edges used to organize navigation of a graph, rather
than being mathematical entities that operations are performed on.

Mesh requirements: (for the data structure, specific algorithms may have stricter requirements)

* Face and vertex indices must range from 0 to N-1, where N is the number of faces or vertices.
* Faces may have any number of vertices.
* Vertices may belong to any number of faces.
* Mesh must be manifold and closed: (The data structure can technically represent meshes with holes,
but closed meshes are assumed for efficiency and simplicity. Dummy faces can be used to represent
holes.)
* Each edge is adjacent to exactly 2 faces.
* Faces sharing a vertex form a closed fan.
* Face vertices must be in a consistent order. This is assumed to be CCW order, though this is only
a convention.


Implementation details:
Conceptually, a quad edge is composed of four directed edges, each directed edge having an origin,
destination, and "next" pointer. A naive implementation would entail a great deal of redundancy.
In this implementation, quad-edges don't exist as a distinct type, a quad-edge is instead a group of
4 directed edges packed in the order: Origin, Right, Destination, Left, with the first dir-edge
aligned to a multiple of the size of a full quad-edge.
The address of any directed edge is a valid means for accessing the quad-edge, and any other edge
of the quad-edge can be accessed with some simple arithmetic.

In the following descriptions, the "forward" direction is considered to be counter-clockwise with
respect to the left face or around the origin vertex.
 
  lnext<-    <-dnext       dprev->    ->rprev   
         dest                     dest          
      left  right              left  right      
         orig                     orig          
  onext<-    <-rnext       lprev->    ->oprev   

lnext: next around left in CCW order
rnext: next around right in CCW order
dnext: next around destination in CCW order
onext: next around origin in CCW order

lprev: previous around left in CCW order
rprev: previous around right in CCW order
dprev: previous around destination in CCW order
oprev: previous around origin in CCW order


null_idx is provided as a null value allowing use of 0 as an index within the structure. The value
of null_idx is all-bits-set.

Total memory usage of the QuadEdgeGraph structure with 32-bit indices is:
8 indices (32 bytes) per quadedge (8 bytes per directed edge)
1 index (4 bytes) per vertex
1 index (4 bytes) per face
Plus an unknown amount of standard library overhead.


Construction:
Use AddEdge() calls to generate quad-edges with connectivity information for vertices and faces.
Call SortAllEdgeLists() to correct edge lists.

Use NewEdge() calls to generate quad-edges.
Assign diredge origins and next pointers manually.

Alternatively, call GenerateFromFaces().


************************************************************************************************* */

#ifndef QUADEDGE_H
#define QUADEDGE_H

#include <vector>
#include <map>
#include <cstdint>
#include <memory>

namespace qeg {

// A simple vector type that makes some guarantees and adds some functionality not in the standard
// library vectors. In particular, it has more predictable memory consumption and has some alignment
// features.
// Stored types are assumed to be plain-old-data types that can be handled as raw data.
// Elements are aligned such that begin()%ALIGN == 0. ALIGN-1 additional bytes are allocated to
// ensure the allocation contains a contiguous block with the desired alignment. No guarantees are
// made as to the validity of memory accesses with the resulting alignment.
template<typename T, size_t ALIGN = 1>
struct aligned_vector {
  private:
    uint8_t * allocation;
    T * _contents;
    size_t _size;
    size_t _capacity;
    
  public:
    aligned_vector(size_t s = 0):
        allocation{nullptr},
        _contents{nullptr},
        _size{0},
        _capacity{0}
    {
        reserve(s);
        _size = s;
    }
    ~aligned_vector() {delete[] allocation;}
    
    auto reserve(size_t newCapacity) -> void {
        auto oldAlloc = allocation;
        auto oldCont = _contents;
        
        auto contentsBytes = sizeof(T)*newCapacity;
        auto allocBytes = contentsBytes + ALIGN - 1;
        auto tmpalloc = new uint8_t[allocBytes];
        auto tmpaligned = reinterpret_cast<void*>(tmpalloc);
        if(!std::align(ALIGN, contentsBytes, tmpaligned, allocBytes))
        {
            delete[] tmpalloc;
            throw std::runtime_error("Could not allocate aligned_vector storage");
        }
        allocation = tmpalloc;
        _contents = reinterpret_cast<T*>(tmpaligned);
        _capacity = newCapacity;
        
        if(_size != 0)
        {
            _size = std::min(_size, _capacity);
            memcpy(_contents, oldCont, sizeof(T)*_size);
        }
        
        delete[] oldAlloc;
    }
    
    // Resize the vector, reallocating and copying data if necessary for tightest fit.
    auto resize(size_t newSize) -> void {
        reserve(newSize);
        _size = newSize;
    }
    
    // Grow the vector, reallocating if necessary, and return a pointer to the new entries.
    auto grow(size_t num) -> T * {
        if(_size + num > _capacity)
        {
            reserve(std::max(_size + num, _capacity*2));
        }
        
        auto items = _contents + _size;
        _size += num;
        return items;
    }
    
    // Truncate the vector, dropping entries from the end.
    auto shrink(size_t num) -> void {_size -= num;}
    
    // Reallocate the vector to minimize unused space.
    auto compact() -> void {reserve(_size);}
    
    // Append a single entry to the vector.
    auto push_back(T val) -> void {
        auto item = grow(1);
        *item = val;
    }
    
    // Non-order-preserving delete operation.
    // Overwrite values in the middle of the vector with values from the end, and then truncate the
    // vector.
    auto quick_delete(T * range_begin, T * range_end) -> void {
        auto numToDelete = range_end - range_begin;
        auto numToMove = std::min(end() - range_end, numToDelete);
        memcpy(range_begin, end() - numToMove, sizeof(T)*numToMove);
        shrink(numToDelete);
    }
    
    auto begin() -> T * {return _contents;}
    auto end() -> T * {return _contents + _size;}
    
    auto alignment() const -> size_t {return ALIGN;}
    auto size() const -> size_t {return _size;}
    auto capacity() const -> size_t {return _capacity;}
    
    // TODO: remove bounds checking or make optional for testing.
    auto operator[](size_t idx) -> T & {
        if(idx >= _size) {throw std::runtime_error("Access out of bounds");}
        return _contents[idx];}
    auto operator[](size_t idx) const -> const T & {
        if(idx >= _size) {throw std::runtime_error("Access out of bounds");}
        return _contents[idx];}
};


// *************************************************************************************************
class QuadEdgeGraph {
  public:
    using idx_t = uint32_t;
    using eidx_t = uint32_t;
    using eoff_t = int32_t;
    enum null_idx_t: uint32_t {null_idx = ~idx_t{0}};// ODR relaxation isn't working!
    
    enum dir_t {
        kOrig = 0,
        kRight = 1,
        kDest = 2,
        kLeft = 3
    };
    
    // Note: DirEdge instances should not be created and manipulated outside of a QuadEdgeGraph, as
    // they rely on the relative locations in memory of other DirEdge instances.
    struct DirEdge {
        idx_t origin;// face/vertex index for origin of edge.
        eoff_t next;// Offset of next edge in CCW order.
        
        // Return a dir-edge of the same quad-edge, relative to this one.
        // Always yields a valid iterator.
        auto FromRight() -> DirEdge * {return modadd_ptr(this, 1);}
        auto FromRight() const -> const DirEdge * {return modadd_ptr(this, 1);}
        auto FromDest() -> DirEdge * {return modadd_ptr(this, 2);}
        auto FromDest() const -> const DirEdge * {return modadd_ptr(this, 2);}
        auto FromLeft() -> DirEdge * {return modadd_ptr(this, 3);}
        auto FromLeft() const -> const DirEdge * {return modadd_ptr(this, 3);}
        
        // Return true if next pointer is valid. Should always be false for completed data
        // structures.
        auto HasNext() const -> bool {return next != 0;}
        auto SetNext(DirEdge * de) -> void {next = de - this;}
        
        // Return the next/previous dir-edge for this diredge.
        // Note: it is safe to use consecutive next calls without checks for validity, the result
        // will simply be this. This can simplify some null pointer checks. A valid data structure
        // will not contain any null pointers, however.
        auto Next() -> DirEdge * {return this + next;}
        auto Next() const -> const DirEdge * {return this + next;}
        auto Prev() -> DirEdge * {return FromRight()->Next()->FromRight();}
        auto Prev() const -> const DirEdge * {return FromRight()->Next()->FromRight();}
        
        // Return a face/vector index of the current quad-edge, relative to this dir-edge.
        auto Origin() const -> const idx_t & {return origin;}
        auto Origin() -> idx_t & {return origin;}
        auto Right() const -> const idx_t & {return FromRight()->origin;}
        auto Right() -> idx_t & {return FromRight()->origin;}
        auto Dest() const -> const idx_t & {return FromDest()->origin;}
        auto Dest() -> idx_t & {return FromDest()->origin;}
        auto Left() const -> const idx_t & {return FromLeft()->origin;}
        auto Left() -> const idx_t & {return FromLeft()->origin;}
    };
    
  protected:
    // The *Edges tables map face/vertex indices to indices of directed edges that are members of
    // the edge loops around those faces/vertices.
    std::vector<eidx_t> faceEdges;// face idx to edge ptr
    std::vector<eidx_t> vertexEdges;// vertex idx to edge ptr
    
    // quadedges array must be aligned to multiple of one full quad-edge for address calculations.
    // This is 32 bytes, a stricter requirement than any likely machine alignment requirements
    aligned_vector<DirEdge, 4*sizeof(DirEdge)> quadedges;
    
    // Compute base + i%4 for base aligned to multiple of 4.
    // static auto modadd(eidx_t idx, unsigned int n) -> eidx_t {
    //     return (idx & (~eidx_t{3})) | ((idx + n) & eidx_t{3});
    // }
    static auto modadd_ptr(DirEdge * ptr, unsigned int n) -> DirEdge * {
        auto addr = reinterpret_cast<uintptr_t>(ptr)/sizeof(DirEdge);
        addr = (addr & (~uintptr_t{3})) | ((addr + n) & uintptr_t{3});
        return reinterpret_cast<DirEdge *>(addr*sizeof(DirEdge));
    }
    static auto modadd_ptr(const DirEdge * ptr, unsigned int n) -> const DirEdge * {
        auto addr = reinterpret_cast<uintptr_t>(ptr)/sizeof(DirEdge);
        addr = (addr & (~uintptr_t{3})) | ((addr + n) & uintptr_t{3});
        return reinterpret_cast<DirEdge *>(addr*sizeof(DirEdge));
    }
    
    auto InitEdge(DirEdge * de, idx_t origin, eidx_t & listHead) -> void;
    
  public:
    QuadEdgeGraph() {}
    
    // Construct the quad-edge graph from a number of faces, number of vertices, and a function or
    // lambda of the form:
    // getFace(idx_t, ind & a, ind & b, ind & c)
    template<typename GetNVertsFn, typename GetVertFn>
    auto GenerateFromFaces(size_t numFaces, size_t numVertices,
                           GetNVertsFn getNumVertices,
                           GetVertFn getVertex) -> void;
    
    // Allocate a new quad-edge and return a pointer to the first dir-edge.
    auto NewEdge() -> DirEdge * {return quadedges.grow(4);};
    
    // Connect two faces (left and right) and two vertices (orig and dest) with a quad-edge, and
    // compute next pointers.
    // A pointer to the new quad-edge is returned.
    auto AddEdge(idx_t orig, idx_t right, idx_t dest, idx_t left) -> DirEdge *;
    
    // Reorder the entries of edge list so they are in order of traversal around the face/vertex.
    // If edge list is incomplete, throw an error.
    auto SortEdgeList(DirEdge * de) -> void;
    
    // Sort edge list for every vertex and face.
    auto SortAllEdgeLists() -> void;
    
    // Return a unique ID for the quad-edge. The ID is independent of the perspective the
    // quad-edge is referenced from (it is the same for all 4 directed edges).
    auto DirEdgeIdx(const DirEdge * de) const -> eidx_t {return de - &quadedges[0];}
    auto QuadEdgeID(const DirEdge * de) const -> uint32_t {return DirEdgeIdx(de)/4;}
    
    auto NumQuadEdges() const -> size_t {return quadedges.size()/4;}
    
    // Get begin/end iterators for iterating entire set of quad-edges in arbitrary order.
    auto begin() -> DirEdge * {return EdgeForIndex(0);}
    static auto begin(QuadEdgeGraph & qeg) -> DirEdge * {return qeg.begin();}
    auto end() -> DirEdge * {return EdgeForIndex(quadedges.size());}
    static auto end(QuadEdgeGraph & qeg) -> DirEdge * {return qeg.end();}
    
    auto begin() const -> const DirEdge * {return EdgeForIndex(0);}
    static auto begin(const QuadEdgeGraph & qeg) -> const DirEdge * {return qeg.begin();}
    auto end() const -> const DirEdge * {return EdgeForIndex(quadedges.size());}
    static auto end(const QuadEdgeGraph & qeg) -> const DirEdge * {return qeg.end();}
    
    // Get quad-edge pointer by face index.
    auto EdgeForFace(idx_t idx) -> DirEdge * {return EdgeForIndex(faceEdges[idx]);}
    auto EdgeForFace(idx_t idx) const -> const DirEdge * {return EdgeForIndex(faceEdges[idx]);}
    
    // Get quad-edge pointer by vertex index.
    auto EdgeForVertex(idx_t idx) -> DirEdge * {return EdgeForIndex(vertexEdges[idx]);}
    auto EdgeForVertex(idx_t idx) const -> const DirEdge * {return EdgeForIndex(vertexEdges[idx]);}
    
    // Get dir-edge pointer by quad-edge index.
    auto EdgeForIndex(eidx_t idx) -> DirEdge * {return &quadedges[idx];}
    auto EdgeForIndex(eidx_t idx) const -> const DirEdge * {return &quadedges[idx];}
    
    auto GetFaceEdges() const -> const std::vector<eidx_t> & {return faceEdges;}
    auto GetVertexEdges() const -> const std::vector<eidx_t> & {return vertexEdges;}
    
    // Get face/vertex index for direction.
    auto operator[](eidx_t ptr) const -> const DirEdge & {return quadedges[ptr];}
    auto operator[](eidx_t ptr) -> DirEdge & {return quadedges[ptr];}
    
    // Iterate over each quad-edge of the graph.
    template<typename Fn>
    auto EachQuadEdge(Fn fn) -> void;
    
    // Count number of edges for given face/vertex.
    auto CountEdgesInLoop(DirEdge * de) -> size_t {
        auto count = size_t{0};
        EachEdgeInLoop(de, [&](DirEdge *){++count;});
        return count;
    }
    
    // Iterate over each edge of given face.
    template<typename Fn>
    auto EachEdgeInLoop(DirEdge * de, Fn fn) -> void;
    
    // Iterate over each strip of edges matching a predicate, starting with given edge.
    // Assumes triangular faces and a predicate that is always true for either 0 or 2 edges of a
    // face.
    template<typename Pred, typename Fn>
    auto EachTriEdgeWithPred(DirEdge * startEdge, Pred pred, Fn fn) -> void;
};


// *************************************************************************************************
inline auto operator<<(std::ostream & ostrm, const QuadEdgeGraph::DirEdge & val) -> std::ostream &
{
    ostrm << "(" << val.origin << " " << val.next << ") ";
    return ostrm;
}


// *************************************************************************************************
inline auto QuadEdgeGraph::
    SortEdgeList(DirEdge * startDE) -> void
{
    // Take start as beginning of sorted list. Compute destination index.
    // There is no general comparison operator, only a test for whether a node should follow
    // directly from the previous one. In addition, the number of edges is typically going to be
    // small, only 3 in the case of triangle meshes.
    // Therefore, we simply perform a selection sort by repeatedly doing linear searches of the
    // unsorted list and appending elements to the sorted one.
    
    // sortedTail is the last entry of the sorted portion of the list.
    // sortedTail->next is the first entry of the unsorted portion of the list.
    // We're looking for node X where quadedges[X->next].origin == Dir(sortedTail->origin, kDest).
    auto sortedTail = startDE;
    auto desiredOrigin = sortedTail->Dest();
    
    while(sortedTail->HasNext() && sortedTail->Next() != startDE)
    {
        // Skip over any portion that is already sorted.
        auto thisNode = sortedTail->Next();
        if(thisNode->Origin() == desiredOrigin)
        {
            sortedTail = thisNode;
            desiredOrigin = sortedTail->Dest();
            continue;
        }
        
        // We have stepped over any sorted portion and need to search for an entry with origin ==
        // desiredOrigin.
        // sortedTail->next points to the head of the unsorted list, which is guaranteed to not be
        // the desired entry, the list head, or null.
        auto prevNode = sortedTail->Next();
        thisNode = prevNode->Next();
        while(thisNode != startDE && thisNode->Origin() != desiredOrigin)
        {
            if(!thisNode->HasNext())
            {
                // Did not find the desired entry.
                throw std::runtime_error("QuadEdgeGraph::SortEdgeList(): edge list is incomplete.");
            }
            prevNode = thisNode;
            thisNode = thisNode->Next();
        }
        
        // Found the right entry. Extract it from the unsorted list...
        if(thisNode->HasNext())
        {
            prevNode->SetNext(thisNode->Next());
        }
        else
        {
            prevNode->next = 0;
        }
        
        // ...and insert at the tail of the sorted list.
        thisNode->SetNext(sortedTail->Next());// sortedTail always has a next at this point
        sortedTail->SetNext(thisNode);
        
        // Update new sorted list tail and desired origin.
        sortedTail = thisNode;
        desiredOrigin = sortedTail->Dest();
    }
    
    // Close the loop if open.
    if(!sortedTail->HasNext())
    {
        sortedTail->SetNext(startDE);
    }
}


// *************************************************************************************************
inline auto QuadEdgeGraph::
    SortAllEdgeLists() -> void
{
    for(auto edgeIdx: vertexEdges)
    {
        if(edgeIdx != null_idx)
        {
            SortEdgeList(EdgeForIndex(edgeIdx));
        }
    }
    
    for(auto edgeIdx: faceEdges)
    {
        if(edgeIdx != null_idx)
        {
            SortEdgeList(EdgeForIndex(edgeIdx));
        }
    }
}


// *************************************************************************************************
inline auto QuadEdgeGraph::
    InitEdge(DirEdge * de, idx_t origin, eidx_t & listHead) -> void
{
    auto deIdx = DirEdgeIdx(de);
    de->origin = origin;
    de->next = (listHead == null_idx)? 0 : eoff_t(listHead - deIdx);
    listHead = deIdx;
}


// *************************************************************************************************
// Add an edge to the quad-edge graph. Adjacency information for the faces and vertices is stored,
// and the directed edges are linked into a preliminary unordered list. Once all quad-edges are
// created, call LinkEdges() to compute the correct next pointers.
// diredges around a face/vertex may be allocated in any order, and the diredge structure does not
// contain enough information to recover the proper order. However, it can still function to list
// the diredges for each face, the order being corrected in a second pass, so just insert each
// diredge as the head of the edge list for each face/vertex.
inline auto QuadEdgeGraph::
    AddEdge(idx_t orig, idx_t right, idx_t dest, idx_t left) -> DirEdge *
{
    auto maxIdx = std::max(orig, dest);
    if(maxIdx >= vertexEdges.size())
        vertexEdges.resize(maxIdx + 1, null_idx);
    
    maxIdx = std::max(left, right);
    if(maxIdx >= faceEdges.size())
        faceEdges.resize(maxIdx + 1, null_idx);
    
    auto qeIdx = quadedges.size();
    auto qe = quadedges.grow(4);
    
    // Add the diredge from the right to the list around the origin, the diredge from the dest to
    // the list around the right, etc.
    InitEdge(&qe[0], orig, faceEdges[left]);
    InitEdge(&qe[1], right, vertexEdges[orig]);
    InitEdge(&qe[2], dest, faceEdges[right]);
    InitEdge(&qe[3], left, vertexEdges[dest]);
    
    return qe;
}


// *************************************************************************************************
// Compute an integer ID code for an directed edge defined by two indices.
// This ID may then be used to determine adjacent faces, using the fact that the corresponding side
// of a neighboring face has the ID with the opposite order.
// 
//      C/E
// A    /    D
//    B/F
// 
inline auto EncodeEdgeID(uint64_t a, uint64_t b) -> uint64_t {return (a << 32) | b;}


// *************************************************************************************************
template<typename GetNVertsFn, typename GetVertFn>
auto QuadEdgeGraph::
    GenerateFromFaces(size_t numFaces, size_t numVertices,
                      GetNVertsFn getNumVertices,
                      GetVertFn getVertex) -> void
{
    // Build tables mapping all edges to faces, for adjacency tests.
    auto faceMap = std::map<uint64_t, idx_t>();
    for(auto faceIdx = idx_t{0}; faceIdx < numFaces; ++faceIdx)
    {
        auto numVertices = getNumVertices(faceIdx);
        auto vertA_Idx = idx_t{0};
        auto vertB_Idx = getVertex(faceIdx, numVertices - 1);
        for(auto vertID = 0; vertID < numVertices; ++vertID)
        {
            vertA_Idx = vertB_Idx;
            vertB_Idx = getVertex(faceIdx, vertID);
            faceMap[EncodeEdgeID(vertA_Idx, vertB_Idx)] = faceIdx;
        }
    }
    
    auto visitedFaces = std::vector<bool>(numFaces, false);
    
    // Adjust size of vertex and face edge tables.
    // This is merely an optimization to avoid reallocating them as edges are added.
    vertexEdges.resize(numVertices, null_idx);
    faceEdges.resize(numFaces, null_idx);
    
    // For each edge of each face, if the neighboring face hasn't been visited, add a quad-edge for
    // that face edge. Mark each face as visited when done.
    for(auto faceIdx = idx_t{0}; faceIdx < numFaces; ++faceIdx)
    {
        auto numVertices = getNumVertices(faceIdx);
        auto vertA_Idx = idx_t{0};
        auto vertB_Idx = getVertex(faceIdx, numVertices - 1);
        
        for(auto endingVert = 0; endingVert < numVertices; ++endingVert)
        {
            vertA_Idx = vertB_Idx;
            vertB_Idx = getVertex(faceIdx, endingVert);
            
            auto neighbor = faceMap.find(EncodeEdgeID(vertB_Idx, vertA_Idx));
            if(neighbor == faceMap.end())
            {
                throw std::runtime_error("No neighboring face exists, mesh contains a hole.");
            }
                
            if(!visitedFaces[neighbor->second])
            {
                // Add this edge.
                // We are circling the center face in CCW order, so the center face is left.
                AddEdge(vertA_Idx, neighbor->second, vertB_Idx, faceIdx);
            }
            // Else, we already did all edges of the neighboring face, so this edge has already been
            // added.
        }
        
        visitedFaces[faceIdx] = true;
    }
    
    SortAllEdgeLists();
}


// *************************************************************************************************
template<typename Fn>
auto QuadEdgeGraph::
    EachQuadEdge(Fn fn) -> void
{
    for(auto qedge = begin(); qedge != end(); qedge += 4)
    {
        fn(qedge);
    }
}


// *************************************************************************************************
template<typename Fn>
auto QuadEdgeGraph::
    EachEdgeInLoop(DirEdge * startEdge, Fn fn) -> void
{
    auto currentEdge = startEdge;
    do {
        fn(currentEdge);
        currentEdge = currentEdge->Next();
    } while(currentEdge != startEdge);
} 


// *************************************************************************************************
template<typename Pred, typename Fn>
auto QuadEdgeGraph::
    EachTriEdgeWithPred(DirEdge * startEdge, Pred predicate, Fn fn) -> void
{
    // Traverse edges for which predicate returns true, starting with startEdge.
    // Assume triangles, with predicate true for 0 or 2 edges.
    // If pred returns false for one edge, assume it is true for the remaining edge.
    // Halt immediately if initial edge fails the predicate.
    if(!predicate(startEdge))
        return;
    
    auto currentEdge = startEdge;
    do {
        fn(currentEdge);
        
        auto nbrEdge = currentEdge->Next();
        if(!predicate(nbrEdge))
        {
            // Edge failed predicate, choose remaining edge.
            nbrEdge = nbrEdge->Next();
        }
        // Continue from neighboring face by flipping to opposing dir-edge.
        currentEdge = nbrEdge->FromDest();
    }
    while(currentEdge != startEdge);
}


// *************************************************************************************************
// QEG_Builder
// *************************************************************************************************
// This class manages some additional data of use while constructing or performing more complex
// modifications of quad-edge graphs, but is unnecessary for general use.
// template<typename qeg_t>
// class QEG_Builder {
//     using idx_t = qeg_t::idx_t;
//     using eidx_t = qeg_t::eidx_t;
    
//     qeg_t & qeg;
//     std::vector<eidx_t> edgeFreelist;
//     std::vector<idx_t> faceIdxFreelist;
//     std::vector<idx_t> vertexIdxFreelist;
    
//   public:
//     QEG_Builder(qeg_t & g): qeg{g} {}
    
//     // SplitEdge()
//     // SplitTriangles()
//     // DeleteEdge()
    
} // namespace qeg
#endif // QUADEDGE_H
