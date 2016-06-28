#include <algorithm>
#include <iostream>
#include <fstream>
#include <climits>
#include <limits>
#include <vector>
#include <cmath>
#include <queue>
#include <list>
#include <set>
#include <map>
#include "lib_SimpleObject.h"

using namespace std;

class PairInfo {
public:
  PairInfo() : cost(0), vertexPair(pair<int, int>(0, 0)), joinedVertexPosition(vector<double>(3, 0)) {}
  PairInfo(double newCost, pair<int, int> newVertexPair, const vector<double>& newJoinedVertexPosition) : cost(newCost), vertexPair(newVertexPair), joinedVertexPosition(newJoinedVertexPosition) {}
  PairInfo(const PairInfo& newPairInfo) : cost(newPairInfo.cost), vertexPair(newPairInfo.vertexPair), joinedVertexPosition(newPairInfo.joinedVertexPosition) {}
  
  double cost;
  pair<int, int> vertexPair;
  vector<double> joinedVertexPosition;
};

bool operator<(const PairInfo& left, const PairInfo& right) {
  return left.cost < right.cost;
}

template <typename T>
vector<T> sort3(const T& v1, const T& v2, const T& v3) {
  vector<T> result;
  if (v1 < v2) {
    if (v2 < v3) {
      result.push_back(v1);
      result.push_back(v2);
      result.push_back(v3);
    }
    else if (v1 < v3) {
      result.push_back(v1);
      result.push_back(v3);
      result.push_back(v2);
    }
    else {
      result.push_back(v3);
      result.push_back(v1);
      result.push_back(v2);
    }
  }
  else {
    if (v2 > v3) {
      result.push_back(v3);
      result.push_back(v2);
      result.push_back(v1);
    }
    else if (v1 > v3) {
      result.push_back(v2);
      result.push_back(v3);
      result.push_back(v1);
    }
    else {
      result.push_back(v2);
      result.push_back(v1);
      result.push_back(v3);
    }
  }
  return result;
}

vector<double> getEquationOfTriangle(const vector<vector<double> >& vertices, const vector<vector<int> >& triangles, int i) {
  vector<double> result(4, 0);
  double a = (vertices[triangles[i][1]][1] - vertices[triangles[i][0]][1]) * 
    (vertices[triangles[i][2]][2] - vertices[triangles[i][0]][2]) -
    (vertices[triangles[i][2]][1] - vertices[triangles[i][0]][1]) * 
    (vertices[triangles[i][1]][2] - vertices[triangles[i][0]][2]);
  double b = (vertices[triangles[i][1]][2] - vertices[triangles[i][0]][2]) * 
    (vertices[triangles[i][2]][0] - vertices[triangles[i][0]][0]) -
    (vertices[triangles[i][2]][2] - vertices[triangles[i][0]][2]) * 
    (vertices[triangles[i][1]][0] - vertices[triangles[i][0]][0]);
  double c = (vertices[triangles[i][1]][0] - vertices[triangles[i][0]][0]) * 
    (vertices[triangles[i][2]][1] - vertices[triangles[i][0]][1]) -
    (vertices[triangles[i][2]][0] - vertices[triangles[i][0]][0]) * 
    (vertices[triangles[i][1]][1] - vertices[triangles[i][0]][1]);
  double length = sqrt(a * a + b * b + c * c);
  a /= length;
  b /= length;
  c /= length;
  double d = - a * vertices[triangles[i][0]][0] 
    - b * vertices[triangles[i][0]][1]
    - c * vertices[triangles[i][0]][2];
  
  result[0] = a;
  result[1] = b;
  result[2] = c;
  result[3] = d;
  return result;
}

vector<double> getEquationOfTriangle(const SimpleOBJ::Vec3f* m_pVertexList, const SimpleOBJ::Array<int,3>* m_pTriangleList, int i) {
  /*
  cout << "Point 1: " << m_pVertexList[m_pTriangleList[i][0] - 1].x << " " << m_pVertexList[m_pTriangleList[i][0] - 1].y << " " << m_pVertexList[m_pTriangleList[i][0] - 1].z << endl;
  cout << "Point 2: " << m_pVertexList[m_pTriangleList[i][1] - 1].x << " " << m_pVertexList[m_pTriangleList[i][1] - 1].y << " " << m_pVertexList[m_pTriangleList[i][1] - 1].z << endl;
  cout << "Point 3: " << m_pVertexList[m_pTriangleList[i][2] - 1].x << " " << m_pVertexList[m_pTriangleList[i][2] - 1].y << " " << m_pVertexList[m_pTriangleList[i][2] - 1].z << endl;
  */
  vector<double> result(4, 0);
  double a = (m_pVertexList[m_pTriangleList[i][1] - 1].y - m_pVertexList[m_pTriangleList[i][0] - 1].y) * 
    (m_pVertexList[m_pTriangleList[i][2] - 1].z - m_pVertexList[m_pTriangleList[i][0] - 1].z) -
    (m_pVertexList[m_pTriangleList[i][2] - 1].y - m_pVertexList[m_pTriangleList[i][0] - 1].y) * 
    (m_pVertexList[m_pTriangleList[i][1] - 1].z - m_pVertexList[m_pTriangleList[i][0] - 1].z);
  double b = (m_pVertexList[m_pTriangleList[i][1] - 1].z - m_pVertexList[m_pTriangleList[i][0] - 1].z) * 
    (m_pVertexList[m_pTriangleList[i][2] - 1].x - m_pVertexList[m_pTriangleList[i][0] - 1].x) -
    (m_pVertexList[m_pTriangleList[i][2] - 1].z - m_pVertexList[m_pTriangleList[i][0] - 1].z) * 
    (m_pVertexList[m_pTriangleList[i][1] - 1].x - m_pVertexList[m_pTriangleList[i][0] - 1].x);
  double c = (m_pVertexList[m_pTriangleList[i][1] - 1].x - m_pVertexList[m_pTriangleList[i][0] - 1].x) * 
    (m_pVertexList[m_pTriangleList[i][2] - 1].y - m_pVertexList[m_pTriangleList[i][0] - 1].y) -
    (m_pVertexList[m_pTriangleList[i][2] - 1].x - m_pVertexList[m_pTriangleList[i][0] - 1].x) * 
    (m_pVertexList[m_pTriangleList[i][1] - 1].y - m_pVertexList[m_pTriangleList[i][0] - 1].y);
  double length = sqrt(a * a + b * b + c * c);
  a /= length;
  b /= length;
  c /= length;
  double d = - a * m_pVertexList[m_pTriangleList[i][0] - 1].x 
    - b * m_pVertexList[m_pTriangleList[i][0] - 1].y
    - c * m_pVertexList[m_pTriangleList[i][0] - 1].z;

  result[0] = a;
  result[1] = b;
  result[2] = c;
  result[3] = d;
  //cout << a << " * x + " << b << " * y + " << c << " * z + " << d << " = 0" << endl;
  //cout << a << "*(" << m_pVertexList[m_pTriangleList[i][0] - 1].x << ")+(" << b << ")*(" << m_pVertexList[m_pTriangleList[i][0] - 1].y << ")+(" << c << ")*(" << m_pVertexList[m_pTriangleList[i][0] - 1].z
  //  << ")+(" << d << ")" << endl;
  /*
  cout << "Point 1: " << m_pVertexList[m_pTriangleList[i][0] - 1].x << " " << m_pVertexList[m_pTriangleList[i][0] - 1].y << " " << m_pVertexList[m_pTriangleList[i][0] - 1].z << endl;
  cout << "Point 2: " << m_pVertexList[m_pTriangleList[i][1] - 1].x << " " << m_pVertexList[m_pTriangleList[i][1] - 1].y << " " << m_pVertexList[m_pTriangleList[i][1] - 1].z << endl;
  cout << "Point 3: " << m_pVertexList[m_pTriangleList[i][2] - 1].x << " " << m_pVertexList[m_pTriangleList[i][2] - 1].y << " " << m_pVertexList[m_pTriangleList[i][2] - 1].z << endl;
  */
  //cout << endl;
  return result;
}

vector<double> getMidPoint(const pair<int, int>& vertexPair, const vector<vector<double> >& Q, const vector<vector<double> >& vertices) {
  vector<double> midPoint(4, 0);
  double delta = Q[0][0] * Q[1][1] * Q[2][2] - Q[0][0] * Q[1][2] * Q[1][2] - Q[0][1] * Q[0][1] * Q[2][2] 
    + Q[0][1] * Q[1][2] * Q[0][2] + Q[0][2] * Q[0][1] * Q[1][2] - Q[0][2] * Q[1][1] * Q[0][2];
  if (delta > 1e-7) {
    double deltaX = Q[0][3] * Q[1][1] * Q[2][2] - Q[0][3] * Q[1][2] * Q[1][2] - Q[1][3] * Q[0][1] * Q[2][2] 
      + Q[1][3] * Q[1][2] * Q[0][2] + Q[2][3] * Q[0][1] * Q[1][2] - Q[2][3] * Q[1][1] * Q[0][2];
    double deltaY = Q[0][0] * Q[1][3] * Q[2][2] - Q[0][0] * Q[2][3] * Q[1][2] - Q[0][1] * Q[0][3] * Q[2][2] 
      + Q[0][1] * Q[2][3] * Q[0][2] + Q[0][2] * Q[0][3] * Q[1][2] - Q[0][2] * Q[1][3] * Q[0][2];
    double deltaZ =  Q[0][0] * Q[1][1] * Q[2][3] - Q[0][0] * Q[1][2] * Q[1][3] - Q[0][1] * Q[0][1] * Q[2][3] 
      + Q[0][1] * Q[1][2] * Q[0][3] + Q[0][2] * Q[0][1] * Q[1][3] - Q[0][2] * Q[1][1] * Q[0][3];
    midPoint[0] = -deltaX / delta;
    midPoint[1] = -deltaY / delta;
    midPoint[2] = -deltaZ / delta;
  }
  else {
    midPoint[0] = (vertices[vertexPair.first][0] + vertices[vertexPair.second][0]) / 2.0;
    midPoint[1] = (vertices[vertexPair.first][1] + vertices[vertexPair.second][1]) / 2.0;
    midPoint[2] = (vertices[vertexPair.first][2] + vertices[vertexPair.second][2]) / 2.0;
  }
  midPoint[3] = 1;
  /* 
  cout << "Point 1:" << vertices[vertexPair.first][0] << ", " << vertices[vertexPair.first][1] << ", " << vertices[vertexPair.first][2] << endl;
  cout << "Point 2:" << vertices[vertexPair.second][0] << ", " << vertices[vertexPair.second][1] << ", " << vertices[vertexPair.second][2] << endl;
  cout << "midPoint:" << midPoint[0] << ", " << midPoint[1] << ", " << midPoint[2] << endl;
  */ 
  return midPoint;
}

PairInfo getPairInfoFromVertexPair(const pair<int, int>& vertexPair, const vector<vector<double> >& vertices, const vector<vector<vector<double> > >& errorMatrixOfVertex) {
  vector<vector<double> > Q(4, vector<double>(4, 0));
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      Q[i][j] = errorMatrixOfVertex[vertexPair.first][i][j] + errorMatrixOfVertex[vertexPair.second][i][j];
    }
  }
  vector<double> midPoint = getMidPoint(vertexPair, Q, vertices);
  vector<double> tempVec(4, 0); 
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      tempVec[i] += midPoint[j] * Q[j][i];
    }
  }
  double deltaV = 0;
  for (int i = 0; i < 4; ++i) {
    deltaV += tempVec[i] * midPoint[i];
  }
  return PairInfo(deltaV, vertexPair, midPoint);
}

void addEdgesOfTriangleToSet(const vector<int>& triangle, set<PairInfo>& s, const vector<vector<double> >& vertices, const vector<vector<vector<double> > >& errorMatrixOfVertex) {
  s.insert(triangle[0] < triangle[1] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[1]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[0]), vertices, errorMatrixOfVertex));
  s.insert(triangle[0] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[2]), vertices, errorMatrixOfVertex)
  : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[0]), vertices, errorMatrixOfVertex));
  s.insert(triangle[1] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[2]), vertices, errorMatrixOfVertex)
  : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[1]), vertices, errorMatrixOfVertex));
}

void addAdjacentEdgesOfTriangleToSet(int adjacentVertex, const vector<int>& triangle, set<PairInfo>& s, const vector<vector<double> >& vertices, const vector<vector<vector<double> > >& errorMatrixOfVertex) {
  if (triangle[0] == adjacentVertex) {
    s.insert(triangle[0] < triangle[1] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[1]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[0]), vertices, errorMatrixOfVertex));
    s.insert(triangle[0] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[2]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[0]), vertices, errorMatrixOfVertex));
  }
  else if (triangle[1] == adjacentVertex) {
    s.insert(triangle[0] < triangle[1] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[1]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[0]), vertices, errorMatrixOfVertex));
    s.insert(triangle[1] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[2]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[1]), vertices, errorMatrixOfVertex));
  }
  else if (triangle[2] == adjacentVertex) {
    s.insert(triangle[0] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[0], triangle[2]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[0]), vertices, errorMatrixOfVertex));
    s.insert(triangle[1] < triangle[2] ? getPairInfoFromVertexPair(pair<int, int>(triangle[1], triangle[2]), vertices, errorMatrixOfVertex)
    : getPairInfoFromVertexPair(pair<int, int>(triangle[2], triangle[1]), vertices, errorMatrixOfVertex));
  }
}

int main(int argc, char **argv) {
  // Process command line arguments
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <input obj file> <output obj file> <vertex number>" << endl;
    return 0;
  }
  
  unsigned int targetVertexNum;
  sscanf(argv[3], "%u", &targetVertexNum);

  // Define variables
  SimpleOBJ::CSimpleObject objFile;
  objFile.LoadFromObj(argv[1]);
  vector<vector<double> > vertices(objFile.m_nVertices + 1, vector<double>(3, 0));
  vector<vector<int> > triangles(objFile.m_nTriangles + 1, vector<int>(3, 0)); 
  vector<set<int> > trianglesOfVertex(objFile.m_nVertices + 1, set<int>());
  vector<vector<double> > equationOfTriangle(objFile.m_nTriangles + 1, vector<double>(4, 0));
  vector<vector<vector<double> > > errorMatrixOfVertex(objFile.m_nVertices + 1, vector<vector<double> >(4, vector<double>(4, 0)));
  vector<pair<int, int> > vertexPairs;
  set<PairInfo> costOfVertexPairs;
  vector<bool> triangleDisabled(objFile.m_nTriangles + 1, false);
  vector<bool> vertexDisabled(objFile.m_nVertices + 1, false);
  
  // Get trianglesOfVertex, equationOfTriangle, vertexPairs from m_pTriangleList, m_pVertexList
  cout << "Calculating trianglesOfVertex, equationOfTriangle, vertexPairs from m_pTriangleList, m_pVertexList..." << endl;
  for (int i = 0; i < objFile.m_nTriangles; ++i) {
    trianglesOfVertex[objFile.m_pTriangleList[i][0]].insert(i + 1);
    trianglesOfVertex[objFile.m_pTriangleList[i][1]].insert(i + 1);
    trianglesOfVertex[objFile.m_pTriangleList[i][2]].insert(i + 1);
    triangles[i + 1][0] = objFile.m_pTriangleList[i][0];
    triangles[i + 1][1] = objFile.m_pTriangleList[i][1];
    triangles[i + 1][2] = objFile.m_pTriangleList[i][2];
    
    equationOfTriangle[i + 1] = getEquationOfTriangle(objFile.m_pVertexList, objFile.m_pTriangleList, i);
    vertexPairs.push_back(pair<int, int>(objFile.m_pTriangleList[i][0], objFile.m_pTriangleList[i][1]));
    vertexPairs.push_back(pair<int, int>(objFile.m_pTriangleList[i][1], objFile.m_pTriangleList[i][2]));
    vertexPairs.push_back(pair<int, int>(objFile.m_pTriangleList[i][2], objFile.m_pTriangleList[i][0]));
  }
  
  // Get vertices, errorMatrixOfVertex from m_pVertexList, trianglesOfVertex, equationOfTriangle
  cout << "Calculating vertices, errorMatrixOfVertex from m_pVertexList, trianglesOfVertex, equationOfTriangle..." << endl;
  for (int i = 0; i < objFile.m_nVertices; ++i) {
    vertices[i + 1][0] = objFile.m_pVertexList[i].x;
    vertices[i + 1][1] = objFile.m_pVertexList[i].y;
    vertices[i + 1][2] = objFile.m_pVertexList[i].z;
    for (set<int>::iterator it = trianglesOfVertex[i + 1].begin(); it != trianglesOfVertex[i + 1].end(); ++it) {
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
          errorMatrixOfVertex[i + 1][j][k] += equationOfTriangle[*it][j] * equationOfTriangle[*it][k];
        }
      }
    }
  }
  
  // Get costOfVertexPairs from vertexPairs, errorMatrixOfVertex
  cout << "Calculating costOfVertexPairs from vertexPairs, errorMatrixOfVertex..." << endl;
  for (vector<pair<int, int> >::iterator it = vertexPairs.begin(); it != vertexPairs.end(); ++it) {
    costOfVertexPairs.insert(getPairInfoFromVertexPair(*it, vertices, errorMatrixOfVertex));
  }

  int removedVertices = 0;
  while (vertices.size() - 1 - removedVertices > targetVertexNum) {
    set<pair<int, int> > vertexPairsConnectedWithFirstVertex;
    set<PairInfo> addedEdges;
    set<PairInfo> removedEdges;
    set<int> outerVertices;
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].end(); ++it) {
      addEdgesOfTriangleToSet(triangles[*it], removedEdges, vertices, errorMatrixOfVertex);
      for (int i = 0; i < 3; ++i) {
        if (triangles[*it][i] != (*costOfVertexPairs.begin()).vertexPair.first && triangles[*it][i] != (*costOfVertexPairs.begin()).vertexPair.second) {
          for (set<int>::iterator it2 = trianglesOfVertex[triangles[*it][i]].begin(); it2 != trianglesOfVertex[triangles[*it][i]].end(); ++it2) {
            addAdjacentEdgesOfTriangleToSet(triangles[*it][i], triangles[*it2], removedEdges, vertices, errorMatrixOfVertex);
          }
        }
      }
    }
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].end(); ++it) {
      addEdgesOfTriangleToSet(triangles[*it], removedEdges, vertices, errorMatrixOfVertex);
      for (int i = 0; i < 3; ++i) {
        if (triangles[*it][i] != (*costOfVertexPairs.begin()).vertexPair.first && triangles[*it][i] != (*costOfVertexPairs.begin()).vertexPair.second) {
          for (set<int>::iterator it2 = trianglesOfVertex[triangles[*it][i]].begin(); it2 != trianglesOfVertex[triangles[*it][i]].end(); ++it2) {
            addAdjacentEdgesOfTriangleToSet(triangles[*it][i], triangles[*it2], removedEdges, vertices, errorMatrixOfVertex);
          }
        }
      }
    }
    // Update vertices, triangles, trianglesOfVertex, vertexDisabled, triangleDisabled
    vector<set<int>::iterator> elementsToRemove;
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].end(); ++it) {
      int thirdVertex = -1;
      if (triangles[*it][0] == (*costOfVertexPairs.begin()).vertexPair.first) {
        if (triangles[*it][1] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 2;
        }
        else if (triangles[*it][2] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 1;
        }
        else {
          vertexPairsConnectedWithFirstVertex.insert(triangles[*it][1] < triangles[*it][2] ? pair<int, int>(triangles[*it][1], triangles[*it][2]) : pair<int, int>(triangles[*it][2], triangles[*it][1]));
        }
      }
      else if (triangles[*it][1] == (*costOfVertexPairs.begin()).vertexPair.first) {
        if (triangles[*it][0] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 2;
        }
        else if (triangles[*it][2] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 0;
        }
        else {
          vertexPairsConnectedWithFirstVertex.insert(triangles[*it][0] < triangles[*it][2] ? pair<int, int>(triangles[*it][0], triangles[*it][2]) : pair<int, int>(triangles[*it][2], triangles[*it][0]));
        }
      }
      else {
        if (triangles[*it][0] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 1;
        }
        else if (triangles[*it][1] == (*costOfVertexPairs.begin()).vertexPair.second) {
          thirdVertex = 0;
        }
        else {
          vertexPairsConnectedWithFirstVertex.insert(triangles[*it][0] < triangles[*it][1] ? pair<int, int>(triangles[*it][0], triangles[*it][1]) : pair<int, int>(triangles[*it][1], triangles[*it][0]));
        }
      }
      if (thirdVertex != -1) {
        trianglesOfVertex[triangles[*it][thirdVertex]].erase(trianglesOfVertex[triangles[*it][thirdVertex]].find(*it));
        trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].erase(trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].find(*it));
        triangleDisabled[*it] = true;
        elementsToRemove.push_back(it);
      }
    }
    for (vector<set<int>::iterator>::iterator it = elementsToRemove.begin(); it != elementsToRemove.end(); ++it) {
      trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].erase(*it);
    }
    elementsToRemove.clear();
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].end(); ++it) {
      if (triangles[*it][0] != (*costOfVertexPairs.begin()).vertexPair.second && triangles[*it][1] != (*costOfVertexPairs.begin()).vertexPair.second && triangles[*it][2] != (*costOfVertexPairs.begin()).vertexPair.second) {
        cout << "Exception: invalid triangle." << endl;
      }
      pair<int, int> vertexPair;
      bool found = false;
      if (triangles[*it][0] == (*costOfVertexPairs.begin()).vertexPair.second) {
        vertexPair = triangles[*it][1] < triangles[*it][2] ? pair<int, int>(triangles[*it][1], triangles[*it][2]) : pair<int, int>(triangles[*it][2], triangles[*it][1]);
        if (vertexPairsConnectedWithFirstVertex.count(vertexPair) == 0) {
          triangles[*it][0] = (*costOfVertexPairs.begin()).vertexPair.first;
        }
        else {
          found = true;
        }
      }
      else if (triangles[*it][1] == (*costOfVertexPairs.begin()).vertexPair.second) {
        vertexPair = triangles[*it][0] < triangles[*it][2] ? pair<int, int>(triangles[*it][0], triangles[*it][2]) : pair<int, int>(triangles[*it][2], triangles[*it][0]);
        if (vertexPairsConnectedWithFirstVertex.count(vertexPair) == 0) {
          triangles[*it][1] = (*costOfVertexPairs.begin()).vertexPair.first;
        }
        else {
          found = true;
        }
      }
      else {
        vertexPair = triangles[*it][0] < triangles[*it][1] ? pair<int, int>(triangles[*it][0], triangles[*it][1]) : pair<int, int>(triangles[*it][1], triangles[*it][0]);
        if (vertexPairsConnectedWithFirstVertex.count(vertexPair) == 0) {
          triangles[*it][2] = (*costOfVertexPairs.begin()).vertexPair.first;
        }
        else {
          found = true;
        }
      }
      if (found) {
        trianglesOfVertex[vertexPair.first].erase(trianglesOfVertex[vertexPair.first].find(*it));
        trianglesOfVertex[vertexPair.second].erase(trianglesOfVertex[vertexPair.second].find(*it));
        triangleDisabled[*it] = true;
        elementsToRemove.push_back(it);
      }
      else {
        trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].insert(*it);
      }
    }
    for (vector<set<int>::iterator>::iterator it = elementsToRemove.begin(); it != elementsToRemove.end(); ++it) {
      trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.second].erase(*it);
    }
    vertices[(*costOfVertexPairs.begin()).vertexPair.first] = (*costOfVertexPairs.begin()).joinedVertexPosition;
    vertexDisabled[(*costOfVertexPairs.begin()).vertexPair.second] = true;
    
    set<int> adjacentVertices;
    adjacentVertices.insert((*costOfVertexPairs.begin()).vertexPair.first);
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].end(); ++it) {
      equationOfTriangle[*it] = getEquationOfTriangle(vertices, triangles, *it);
      for (int i = 0; i < 3; ++i) {
        adjacentVertices.insert(triangles[*it][i]);
      }
    }

    for (set<int>::iterator it = adjacentVertices.begin(); it != adjacentVertices.end(); ++it) {
      for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
          errorMatrixOfVertex[*it][j][k] = 0;
        }
      } 
      for (set<int>::iterator it2 = trianglesOfVertex[*it].begin(); it2 != trianglesOfVertex[*it].end(); ++it2) {
        for (int j = 0; j < 4; ++j) {
          for (int k = 0; k < 4; ++k) {
            errorMatrixOfVertex[*it][j][k] += equationOfTriangle[*it2][j] * equationOfTriangle[*it2][k];
          }
        }
      }
    }
    
    for (set<int>::iterator it = trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].begin(); it != trianglesOfVertex[(*costOfVertexPairs.begin()).vertexPair.first].end(); ++it) {
      addEdgesOfTriangleToSet(triangles[*it], addedEdges, vertices, errorMatrixOfVertex);
      for (int i = 0; i < 3; ++i) {
        if (triangles[*it][i] != (*costOfVertexPairs.begin()).vertexPair.first) {
          for (set<int>::iterator it2 = trianglesOfVertex[triangles[*it][i]].begin(); it2 != trianglesOfVertex[triangles[*it][i]].end(); ++it2) {
            addAdjacentEdgesOfTriangleToSet(triangles[*it][i], triangles[*it2], addedEdges, vertices, errorMatrixOfVertex);
          }
        }
      }
    }
        
    // TODO:
    // 有些顶点被删除了，有些顶点被修改了，有些顶点是新增加的，但removedEdges和updatedEdges总数为O(1)的
    // 把costOfVertexPairs由vector改为set 
    // 之前已经把需要删除和更新的边记录为pair<int, int>的形式（顶点编号对）
    // 由此直接得到PairInfo的vertexPair变量
    // PairInfo的cost变量可以由errorMatrixOfVertex重构出来 
    // PairInfo的joinedVertexPosition也可以由errorMatrixOfVertex重构出来 
    // 删除时可以用O(log(n))的时间，先find再erase即可 
    // 修改时也可以用O(log(n))的时间，先find再erase再insert即可 
    
    for (set<PairInfo>::iterator it = removedEdges.begin(); it != removedEdges.end(); ++it) {
      set<PairInfo>::iterator position = costOfVertexPairs.find(*it);
      if (position != costOfVertexPairs.end()) {
        costOfVertexPairs.erase(position);
      }
      /*
      else {
        cout << "Exception: pair to remove not found." << endl;
        cout << "pair to remove: (" << it->vertexPair.first << ", " << it->vertexPair.second << ")" << endl;
        cout << "removedVertices = " << removedVertices << endl;
      }
      */
    }
    
    for (set<PairInfo>::iterator it = addedEdges.begin(); it != addedEdges.end(); ++it) {
      costOfVertexPairs.insert(*it);
    }

    if (++removedVertices % 1000 == 0) {
      cout << removedVertices << " vertices removed" << endl;
    }
  }
  
  ofstream fout(argv[2]);
  int vertexId = 0;
  vector<int> decrementOfId;
  int currentDecrement = 0;
  decrementOfId.push_back(currentDecrement);
  for (vector<vector<double> >::iterator it = vertices.begin() + 1; it != vertices.end(); ++it) {
    if (!vertexDisabled[++vertexId]) {
      fout << "v " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
    }
    else {
      ++currentDecrement; 
    }
    decrementOfId.push_back(currentDecrement);
  }
  
  int triangleId = 0;
  for (vector<vector<int> >::iterator it = triangles.begin() + 1; it != triangles.end(); ++it) {
    if (!triangleDisabled[++triangleId]) {
      fout << "f " << (*it)[0] - decrementOfId[(*it)[0]] << " " << (*it)[1] - decrementOfId[(*it)[1]]  << " " << (*it)[2] - decrementOfId[(*it)[2]] << endl;
    }
  }
  
  fout.close();
  return 0;
}
