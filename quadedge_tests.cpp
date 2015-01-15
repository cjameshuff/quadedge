
#include <iostream>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <vector>

#include "lambdatest.h"

#include "quadedge.h"

using namespace std;
using namespace lambdatest;

LAMBDATEST_DECL

using QEG = qeg::QuadEdgeGraph;



auto PrintQEG(QEG & graph) -> void
{
    std::clog << "Graph has " << graph.NumQuadEdges() << " quad-edges." << std::endl;
    
    std::clog << "Face edges: (" << graph.GetFaceEdges().size() << ")" << std::endl;
    auto ctr = 0;
    for(auto fedge: graph.GetFaceEdges())
    {
        std::clog << fedge << " ";
        if(++ctr % 8 == 0)
            std::clog << std::endl;
    }
    std::clog << std::endl;
    
    std::clog << "Vertex edges: (" << graph.GetVertexEdges().size() << ")" << std::endl;
    ctr = 0;
    for(auto vedge: graph.GetVertexEdges())
    {
        std::clog << vedge << " ";
        if(++ctr % 8 == 0)
            std::clog << std::endl;
    }
    std::clog << std::endl;
    
    std::clog << "Quad edges: (" << graph.NumQuadEdges() << ")" << std::endl;
    for(auto qedgeID = 0; qedgeID < graph.NumQuadEdges(); ++qedgeID)
    {
        for(auto dedgeID = 0; dedgeID < 4; ++dedgeID)
        {
            auto dedgeIdx = qedgeID*4 + dedgeID;
            auto dedge = graph.EdgeForIndex(dedgeIdx);
            std::clog << dedgeIdx << *dedge;
        }
        std::clog << std::endl;
    }
    std::clog << std::endl;
}


// Generate a test graph for a tetrahedron.
// To avoid confusion, vertex indices are offset by 100, face indices by 200.
// There are 6 quad-edges, 24 directed edges.
// Vertices are A (0), B (1), C (2), D (3).
// Edges are AB, BC, CA, AD, BD, CD.
// Faces are ABC (0), ADB (1), BDC (2), CDA (3)
//
//         D
//        / \
//     3 / 2 \ 1
//      / BDC \
//     C_______B
//    / \ ABC / \
// 2 / 3 \ 0 / 1 \ 2
//  / CDA \ / ADB \
// D-------A-------D
//     1       3
auto QEG_Tetrahedron(QEG & graph) -> void
{
    graph.AddEdge(100, 201, 101, 200);// AB
    graph.AddEdge(101, 202, 102, 200);// BC
    graph.AddEdge(102, 203, 100, 200);// CA
    graph.AddEdge(100, 203, 103, 201);// AD
    graph.AddEdge(101, 201, 103, 202);// BD
    graph.AddEdge(102, 202, 103, 203);// CD
    graph.SortAllEdgeLists();
}


// Generate a test graph for a cube composed of 4 quads.
// To avoid confusion, vertex indices are offset by 100, face indices by 200.
// There are 12 quad-edges, 48 directed edges.
// Faces are L (0), R (1), T (2), B (3), N (4), F(5) (Left, Right, Top, Bottom, Near, Far)
// Vertices are LBN (0), LBF (1), LTN (2), LTF (3), RBN (4), RBF (5), RTN (6), RTF (7).
// Edges are .
//
//     F     F     F
//  1-----3-----7-----5
//  |     |     |     |
// B|  L  |  T  |  R  |B
//  |     |     |     |
//  0-----2-----6-----4
//     N  |     |  N
//       L|  N  |R
//        |     |
//        0-----4
//        |     |
//       L|  B  |R
//        |     |
//        1-----5
//        |     |
//       L|  F  |R
//        |     |
//        3-----7
//           T
//
auto QEG_Cube(QEG & graph) -> void
{
    enum {LBN, LBF, LTN, LTF, RBN, RBF, RTN, RTF};
    int faces[6][4] = {
        {LBF, LBN, LTN, LTF},// Left
        {RTF, RTN, RBN, RBF},// Right
        {LTN, RTN, RTF, LTF},// Top
        {LBF, RBF, RBN, LBN},// Bottom
        {LBN, RBN, RTN, LTN},// Near
        {LTF, RTF, RBF, LBF} // Far
    };
    
    graph.GenerateFromFaces(6, 8,
        [](QEG::idx_t faceIdx){return 4;},
        [&](QEG::idx_t faceIdx, QEG::idx_t vertIdx){return faces[faceIdx][vertIdx];});
}


auto main(int argc, const char * argv[]) -> int
{
    auto vec = (qeg::aligned_vector<int, 1>*){nullptr};
    
    TestGroup("aligned_vector<>", [&]{
    
    // Setup for default constructor
    Setup([&]{
        vec = new qeg::aligned_vector<int, 1>();
    });
    
    // Teardown for all tests
    Teardown([&]{delete vec;});
    
    Test("default constructor", [&]{
        ShouldEq("vec->size() == 0", vec->size(), 0);
        ShouldEq("vec->capacity() == 0", vec->capacity(), 0);
        ShouldEq("vec->begin() == vec->end()", vec->begin(), vec->end());
    });
    
    // Setup with 16 values
    Setup([&]{
        vec = new qeg::aligned_vector<int, 1>(16);
        for(int j = 0; j < 16; ++j)
        {
            (*vec)[j] = 100 + j;
        }
    });
    
    Test("allocation", [&]{
        ShouldEq("vec->size() == 16", vec->size(), 16);
        ShouldEq("vec->capacity() == 16", vec->capacity(), 16);
        ShouldEq("vec->begin() + 16 == vec->end()", vec->begin() + 16, vec->end());
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109, 110, 111,
            112, 113, 114, 115
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->reserve()", [&]{
        auto oldBegin = vec->begin();
        
        vec->reserve(20);
        
        ShouldNEq("vec->begin() != oldBegin", vec->begin(), oldBegin);
        ShouldEq("vec->size() == 16", vec->size(), 16);
        ShouldEq("vec->capacity() == 20", vec->capacity(), 20);
        ShouldEq("vec->begin() + 16 == vec->end()", vec->begin() + 16, vec->end());
        
        oldBegin = vec->begin();
        vec->reserve(10);
        
        ShouldNEq("vec->begin() != oldBegin", vec->begin(), oldBegin);
        ShouldEq("vec->size() == 10", vec->size(), 10);
        ShouldEq("vec->capacity() == 10", vec->capacity(), 10);
        ShouldEq("vec->begin() + 10 == vec->end()", vec->begin() + 10, vec->end());
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->grow()", [&]{
        auto newEntries = vec->grow(4);
        newEntries[0] = 1001;
        newEntries[1] = 1002;
        newEntries[2] = 1003;
        newEntries[3] = 1004;
        
        ShouldEq("vec->size() == 20", vec->size(), 20);
        ShouldEq("vec->capacity() == 32", vec->capacity(), 32);
        ShouldEq("vec->begin() + 20 == vec->end()", vec->begin() + 20, vec->end());
        ShouldEq("newEntries == vec->begin() + 16", newEntries, vec->begin() + 16);
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109, 110, 111,
            112, 113, 114, 115,
            1001, 1002, 1003, 1004
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->shrink()", [&]{
        vec->shrink(4);
        
        ShouldEq("vec->size() == 12", vec->size(), 12);
        ShouldEq("vec->capacity() == 16", vec->capacity(), 16);
        ShouldEq("vec->begin() + 12 == vec->end()", vec->begin() + 12, vec->end());
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109, 110, 111
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->compact()", [&]{
        vec->reserve(20);
        
        ShouldEq("vec->size() == 16", vec->size(), 16);
        ShouldEq("vec->capacity() == 20", vec->capacity(), 20);
        ShouldEq("vec->begin() + 16 == vec->end()", vec->begin() + 16, vec->end());
        
        vec->compact();
        
        ShouldEq("vec->size() == 16", vec->size(), 16);
        ShouldEq("vec->capacity() == 16", vec->capacity(), 16);
        ShouldEq("vec->begin() + 16 == vec->end()", vec->begin() + 16, vec->end());
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109, 110, 111,
            112, 113, 114, 115
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->resize()", [&]{
        // Resize adjusts the capacity to be an exact fit.
        vec->resize(20);
        
        ShouldEq("vec->size() == 20", vec->size(), 20);
        ShouldEq("vec->capacity() == 20", vec->capacity(), 20);
        ShouldEq("vec->begin() + 20 == vec->end()", vec->begin() + 20, vec->end());
        
        vec->resize(10);
        
        ShouldEq("vec->size() == 10", vec->size(), 10);
        ShouldEq("vec->capacity() == 10", vec->capacity(), 10);
        ShouldEq("vec->begin() + 10 == vec->end()", vec->begin() + 10, vec->end());
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->push_back()", [&]{
        vec->push_back(777);
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            104, 105, 106, 107,
            108, 109, 110, 111,
            112, 113, 114, 115,
            777
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    Test("ec->quick_delete()", [&]{
        vec->quick_delete(vec->begin() + 4, vec->begin() + 8);
        
        auto expected = std::vector<int>{{
            100, 101, 102, 103,
            112, 113, 114, 115,
            108, 109, 110, 111
        }};
        MShouldEq("contents", vec->begin(), &expected[0], expected.size());
    });
    
    });// TestGroup("aligned_vector")
    
    
    auto qegraph = (QEG*){nullptr};
    
    TestGroup("QuadEdge", [&]{
    
    // Setup for default constructor
    Setup([&]{
        qegraph = new QEG();
    });
    Teardown([&]{delete qegraph;});
    
    Test("default constructor", [&]{
        ShouldEq("NumQuadEdges()", qegraph->NumQuadEdges(), 0);
    });
    
    Test("AddEdge, DirEdgeIdx, FromLeft/Right/Dest", [&]{
        auto qe = qegraph->AddEdge(10, 11, 12, 13);
        ShouldEq("qe->origin", qe->origin, 10);
        ShouldEq("qe->next", qe->next, 0);
        
        ShouldEq("qe->DirEdgeIdx()", qegraph->DirEdgeIdx(qe), 0);
        ShouldEq("qe->FromRight().DirEdgeIdx()", qegraph->DirEdgeIdx(qe->FromRight()), 1);
        ShouldEq("qe->FromDest().DirEdgeIdx()", qegraph->DirEdgeIdx(qe->FromDest()), 2);
        ShouldEq("qe->FromLeft().DirEdgeIdx()", qegraph->DirEdgeIdx(qe->FromLeft()), 3);
        ShouldEq("qe->FromLeft().FromRight().DirEdgeIdx()", qegraph->DirEdgeIdx(qe->FromLeft()->FromRight()), 0);
        
        ShouldEq("qe->origin", qe->origin, 10);
        ShouldEq("qe->FromRight()->origin", qe->FromRight()->origin, 11);
        ShouldEq("qe->FromDest()->origin", qe->FromDest()->origin, 12);
        ShouldEq("qe->FromLeft()->origin", qe->FromLeft()->origin, 13);
        ShouldEq("qe->next", qe->next, 0);
        ShouldEq("qe->FromRight()->next", qe->FromRight()->next, 0);
        ShouldEq("qe->FromDest()->next", qe->FromDest()->next, 0);
        ShouldEq("qe->FromLeft()->next", qe->FromLeft()->next, 0);
        
        // Face and vertex lists should point to appropriate dir-edges
        auto & qeg_fedges = qegraph->GetFaceEdges();
        auto & qeg_vedges = qegraph->GetVertexEdges();
        // Edge around the origin comes from the right.
        ShouldEq("qeg_fedges[10] == qe->FromRight()", qeg_vedges[10], qegraph->DirEdgeIdx(qe->FromRight()));
        // Edge around the right comes from the destination.
        ShouldEq("qeg_vedges[11] == qe->FromDest()", qeg_fedges[11], qegraph->DirEdgeIdx(qe->FromDest()));
        // Edge around the destination comes from the left.
        ShouldEq("qeg_fedges[12] == qe->FromLeft()", qeg_vedges[12], qegraph->DirEdgeIdx(qe->FromLeft()));
        // Edge around the left comes from the origin.
        ShouldEq("qeg_vedges[13] == qe", qeg_fedges[13], qegraph->DirEdgeIdx(qe));
    });
    
    Test("AddEdge, multiple edges around triangle", [&]{
        // Add edges for triangle 10, with vertices 0, 1, and 2, and neighbored by faces 11, 12, and
        // 13.
        auto qeAB = qegraph->AddEdge(0, 11, 1, 10);
        auto qeBC = qegraph->AddEdge(1, 12, 2, 10);
        auto qeCA = qegraph->AddEdge(2, 13, 0, 10);
        
        // PrintQEG(*qegraph);
        
        // Now, there should be dir-edges forming a complete loop around face 10, in a
        // null-terminated linked list in arbitrary order.
        auto edge0 = qegraph->EdgeForFace(10);
        auto edge1 = edge0->Next();
        auto edge2 = edge1->Next();
        
        Should("edge0 has a next", edge0->HasNext());
        Should("edge1 has a next", edge1->HasNext());
        Should("edge2 does not have a next", !edge2->HasNext());
        
        // TODO: test in more detail.
    });
    
    
    Test("SortEdgeList", [&]{
        // Add edges for triangle 10, with vertices 0, 1, and 2, and neighbored by faces 11, 12, and
        // 13.
        auto qeAB = qegraph->AddEdge(0, 11, 1, 10);
        auto qeBC = qegraph->AddEdge(1, 12, 2, 10);
        auto qeCA = qegraph->AddEdge(2, 13, 0, 10);
        
        PrintQEG(*qegraph);
        
        std::clog << "Edge list for face 10:" << std::endl;
        auto edge = qegraph->EdgeForFace(10);
        while(1)
        {
            std::clog << qegraph->DirEdgeIdx(edge) << "(" << edge->origin << " " << edge->next << ") ";
            if(edge->HasNext())
                edge = edge->Next();
            else break;
        }
        std::clog << std::endl;
        
        std::clog << "Sorting edge list..." << std::endl;
        qegraph->SortEdgeList(qegraph->EdgeForFace(10));
        
        PrintQEG(*qegraph);
        
        // Now, there should be dir-edges forming a complete loop around face 10, in a
        // closed linked list in sequential order (dest of one == orig of next).
        auto edge0 = qegraph->EdgeForFace(10);
        auto edge1 = edge0->Next();
        auto edge2 = edge1->Next();
        
        std::clog << qegraph->DirEdgeIdx(edge0) << *edge0 << ") ";
        std::clog << qegraph->DirEdgeIdx(edge1) << *edge1 << ") ";
        std::clog << qegraph->DirEdgeIdx(edge2) << *edge2 << ") ";
        std::clog << std::endl;
        
        ShouldEq("edge2->Next() -> edge0", edge2->Next(), edge0);
        
        ShouldEq("edge0 -> edge1", edge0->Dest(), edge1->Origin());
        ShouldEq("edge1 -> edge2", edge1->Dest(), edge2->Origin());
        ShouldEq("edge2 -> edge0", edge2->Dest(), edge0->Origin());
        
        // TODO: add some more assertions to automate checks that the correct edges are in the list
    });
    
    
    Test("Traversal test", [&]{
        // Test navigation of a full mesh structure.
        QEG_Tetrahedron(*qegraph);
        PrintQEG(*qegraph);
        
        auto edge0 = qegraph->EdgeForFace(200);
        auto edge1 = edge0->Next();
        auto edge2 = edge1->Next();
        
        // We should have 3 edges of triangle 200, in proper order.
        ShouldEq("edge0 -> edge1", edge0->Dest(), edge1->Origin());
        ShouldEq("edge1 -> edge2", edge1->Dest(), edge2->Origin());
        ShouldEq("edge2 -> edge0", edge2->Dest(), edge0->Origin());
        
        edge0 = qegraph->EdgeForVertex(100);
        edge1 = edge0->Next();
        edge2 = edge1->Next();
        
        // We should have 3 edges of triangle 200, in proper order.
        ShouldEq("edge0 -> edge1", edge0->Dest(), edge1->Origin());
        ShouldEq("edge1 -> edge2", edge1->Dest(), edge2->Origin());
        ShouldEq("edge2 -> edge0", edge2->Dest(), edge0->Origin());
        
        ShouldEq("edge0 == edge1->Prev()", edge0, edge1->Prev());
        ShouldEq("edge1 == edge2->Prev()", edge1, edge2->Prev());
        ShouldEq("edge2 == edge0->Prev()", edge2, edge0->Prev());
    });
    
    
    Test("GenerateFromFaces", [&]{
        // Test generation from faces.
        QEG_Cube(*qegraph);
        PrintQEG(*qegraph);
        
        for(auto faceIdx = QEG::idx_t{0}; faceIdx < 6; ++faceIdx)
        {
            auto count = qegraph->CountEdgesInLoop(qegraph->EdgeForFace(faceIdx));
            ShouldEq("face edges == 4", count, 4);
        }
        
        for(auto vertIdx = QEG::idx_t{0}; vertIdx < 8; ++vertIdx)
        {
            auto count = qegraph->CountEdgesInLoop(qegraph->EdgeForVertex(vertIdx));
            ShouldEq("vertex edges == 3", count, 3);
        }
    });
    
    
    });// TestGroup("QuadEdge")
    
    return EXIT_SUCCESS;
}
