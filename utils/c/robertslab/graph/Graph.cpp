/*
 * Copyright 2017 Johns Hopkins University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
 *
 * Author(s): Elijah Roberts
 */

/*
 * BioNetGen resources:
 *
 * https://www.csb.pitt.edu/Faculty/Faeder/Publications/Reprints/Faeder_2009.pdf
 */

#include <list>
#include <set>
#include <sstream>

#include "robertslab/graph/Graph.h"

using std::list;
using std::set;
using std::string;

namespace robertslab {
namespace graph {

GraphMapping::GraphMapping()
{
}

GraphMapping::GraphMapping(Graph* sourceGraph, Graph* targetGraph)
:sourceGraph(sourceGraph),targetGraph(targetGraph)
{
}

void GraphMapping::appendMappings(const GraphMapping& additionalMapping)
{
    for (int i=0; i<additionalMapping.sourceVertices.size(); i++)
    {
        sourceVertices.push_back(additionalMapping.sourceVertices[i]);
        targetVertices.push_back(additionalMapping.targetVertices[i]);
        forwardMapping[additionalMapping.sourceVertices[i]] = additionalMapping.targetVertices[i];
        reverseMapping[additionalMapping.targetVertices[i]] = additionalMapping.sourceVertices[i];
    }
    for (int i=0; i<additionalMapping.unusedSourceVertices.size(); i++)
    {
        unusedSourceVertices.push_back(additionalMapping.unusedSourceVertices[i]);
    }
    for (int i=0; i<additionalMapping.unusedTargetVertices.size(); i++)
    {
        unusedTargetVertices.push_back(additionalMapping.unusedTargetVertices[i]);
    }
}


void GraphMapping::setGraphs(Graph* sourceGraph, Graph* targetGraph)
{
    this->sourceGraph = sourceGraph;
    this->targetGraph = targetGraph;
}

void GraphMapping::addMapping(Vertex* sourceVertex, Vertex* targetVertex)
{
    sourceVertices.push_back(sourceVertex);
    targetVertices.push_back(targetVertex);
    forwardMapping[sourceVertex] = targetVertex;
    reverseMapping[targetVertex] = sourceVertex;
}

bool GraphMapping::containsSourceVertex(Vertex* sourceVertex)
{
    return (forwardMapping.count(sourceVertex) > 0);
}

int GraphMapping::getNumberMatches()
{
    return sourceVertices.size();
}

Vertex* GraphMapping::getSourceVertex(int index)
{
    return sourceVertices[index];
}

Vertex* GraphMapping::getTargetVertex(int index)
{
    return targetVertices[index];
}

Vertex* GraphMapping::getSourceVertex(Vertex* targetVertex)
{
    return reverseMapping[targetVertex];
}

Vertex* GraphMapping::getTargetVertex(Vertex* sourceVertex)
{
    return forwardMapping[sourceVertex];
}

void GraphMapping::addUnusedSourceVertex(Vertex* vertex)
{
    unusedSourceVertices.push_back(vertex);
}

void GraphMapping::addUnusedTargetVertex(Vertex* vertex)
{
    unusedTargetVertices.push_back(vertex);
}

int GraphMapping::getNumberUnusedSourceVertices()
{
    return unusedSourceVertices.size();
}

int GraphMapping::getNumberUnusedTargetVertices()
{
    return unusedTargetVertices.size();
}

Vertex* GraphMapping::getUnusedSourceVertex(int index)
{
    return unusedSourceVertices[index];
}

Vertex* GraphMapping::getUnusedTargetVertex(int index)
{
    return unusedTargetVertices[index];
}

bool GraphMapping::isSourceVertexUsed(Vertex* vertex)
{
    for (int i=0; i<unusedSourceVertices.size(); i++)
    {
        if (unusedSourceVertices[i] == vertex) return false;
    }
    return true;
}

bool GraphMapping::isTargetVertexUsed(Vertex* vertex)
{
    for (int i=0; i<unusedTargetVertices.size(); i++)
    {
        if (unusedTargetVertices[i] == vertex) return false;
    }
    return true;
}

string GraphMapping::getString(bool includeGraphs)
{
    if (sourceGraph == NULL || targetGraph == NULL) return "NULL";

    std::stringstream ss;
    bool firstLine = true;
    for (int i=0; i<sourceVertices.size(); i++)
    {
        if (firstLine) firstLine = false;
        else ss << "\n";
        if (includeGraphs) ss << "<" << sourceGraph->getString() << "> ";
        ss << sourceVertices[i]->getString().c_str();
        ss << " -> ";
        if (includeGraphs) ss << "<" << targetGraph->getString() << "> ";
        ss << targetVertices[i]->getString().c_str();
    }
    for (int i=0; i<unusedSourceVertices.size(); i++)
    {
        if (firstLine) firstLine = false;
        else ss << "\n";
        if (includeGraphs) ss << "<" << sourceGraph->getString() << "> ";
        ss << unusedSourceVertices[i]->getString().c_str();
        ss << " -> ";
        if (includeGraphs) ss << "<" << targetGraph->getString() << "> ";
        ss << "0";
    }
    for (int i=0; i<unusedTargetVertices.size(); i++)
    {
        if (firstLine) firstLine = false;
        else ss << "\n";
        if (includeGraphs) ss << "<" << sourceGraph->getString() << "> ";
        ss << "0";
        ss << " -> ";
        if (includeGraphs) ss << "<" << targetGraph->getString() << "> ";
        ss << unusedTargetVertices[i]->getString().c_str();
    }
    return ss.str();
}

void GraphMapping::reverse()
{
    Graph* tmp1 = sourceGraph;
    sourceGraph = targetGraph;
    targetGraph = tmp1;

    vector<Vertex*> tmp2 = sourceVertices;
    sourceVertices = targetVertices;
    targetVertices = tmp2;

    forwardMapping.clear();
    reverseMapping.clear();
    for (int i=0; i<sourceVertices.size(); i++)
    {
        forwardMapping[sourceVertices[i]] = targetVertices[i];
        reverseMapping[targetVertices[i]] = sourceVertices[i];
    }

    vector<Vertex*> tmp3 = unusedSourceVertices;
    unusedSourceVertices = unusedTargetVertices;
    unusedTargetVertices = tmp3;
}

Vertex::Vertex()
:nextClearedMark(0)
{

}

void Vertex::setMark(string mark)
{
    marks.insert(mark);
}

void Vertex::setMarkOnConnected(string mark)
{
    if (marks.count(mark) == 0)
    {
        marks.insert(mark);
        for (int i=0; i<getMaxNumberEdges(); i++)
        {
            Vertex* v = getEdge(i);
            if (v != NULL) v->setMarkOnConnected(mark);
        }
    }
}

bool Vertex::hasMark(string mark)
{
    return (marks.count(mark) == 1);
}

void Vertex::clearMark(string mark)
{
    marks.erase(mark);
}

void Vertex::clearMarkOnConnected(string mark, string clearedMark)
{
    // If this is the root, generate a unqie mark to use for marking vertices that have been cleared.
    if (clearedMark == "") clearedMark = "clearMarkOnConnected"+std::to_string((unsigned long long)this)+std::to_string(nextClearedMark++);

    // If this vertex hasn't been cleared yet.
    if (marks.count(clearedMark) == 0)
    {
        // Clear the mark and mark that this vertex has been cleared.
        marks.erase(mark);
        marks.insert(clearedMark);

        // Clear any connected vertices.
        for (int i=0; i<getMaxNumberEdges(); i++)
        {
            Vertex* v = getEdge(i);
            if (v != NULL) v->clearMarkOnConnected(mark, clearedMark);
        }
    }
}

void Vertex::clearAllMarks()
{
    marks.clear();
}

set<Vertex*> Vertex::getConnectedVertices(bool root)
{
    if (root) clearMarkOnConnected("getConnectedVertices");

    // Add this vertex to the set.
    set<Vertex*> ret;
    ret.insert(this);
    setMark("getConnectedVertices");

    // Add any vertices from the edges.
    for (int i=0; i<getMaxNumberEdges(); i++)
    {
        Vertex* v2 = getEdge(i);
        if (v2 != NULL && !v2->hasMark("getConnectedVertices"))
        {
            set<Vertex*> ret2=getEdge(i)->getConnectedVertices(false);
            ret.insert(ret2.begin(), ret2.end());
        }
    }
    return ret;
}

int Vertex::findEdgeLeadingTo(Vertex* destination)
{
    for (int i=0; i<getMaxNumberEdges(); i++)
        if (getEdge(i) == destination)
            return i;
    return -1;
}

bool Graph::addEdge(Vertex* v1, int index1, Vertex* v2, int index2)
{
    // Make sure we can find both vertices.
    int vertexCount=0;
    for (int i=0; i<getNumberVertices(); i++)
    {
        if (getVertex(i) == v1) vertexCount++;
        if (getVertex(i) == v2) vertexCount++;
    }
    if (vertexCount != 2) return false;

    // Make sure that both verticies can be connected at the specified indicies.
    if (index1 >= v1->getMaxNumberEdges() || v1->getEdge(index1) != NULL) return false;
    if (index2 >= v2->getMaxNumberEdges() || v2->getEdge(index2) != NULL) return false;

    // Add the edges.
    v1->addEdge(index1, v2, index2);
    v2->addEdge(index2, v1, index1);

    return true;
}

bool Graph::removeEdge(Vertex* v1, Vertex* v2)
{
    // Make sure we can find both vertices.
    int vertexCount=0;
    for (int i=0; i<getNumberVertices(); i++)
    {
        if (getVertex(i) == v1) vertexCount++;
        if (getVertex(i) == v2) vertexCount++;
    }
    if (vertexCount != 2) return false;

    // Make sure the vertices are linked by an edge.
    int l1=-1;
    for (int i=0; i<v1->getMaxNumberEdges(); i++)
    {
        if (v1->getEdge(i) == v2)
        {
            l1=i;
            break;
        }
    }
    int l2=-1;
    for (int i=0; i<v2->getMaxNumberEdges(); i++)
    {
        if (v2->getEdge(i) == v1)
        {
            l2=i;
            break;
        }
    }
    if (l1 == -1 || l2 == -1) return false;

    // Remove the edges.
    v1->removeEdge(l1);
    v2->removeEdge(l2);
    return true;
}

bool Graph::hasMarkOnAny(string mark)
{
    for (int i=0; i<getNumberVertices(); i++)
        if (getVertex(i)->hasMark(mark))
            return true;
    return false;
}

void Graph::clearMark(string mark)
{
    for (int i=0; i<getNumberVertices(); i++)
        getVertex(i)->clearMark(mark);
}

void Graph::clearAllMarks()
{
    for (int i=0; i<getNumberVertices(); i++)
        getVertex(i)->clearAllMarks();
}

bool Graph::isIsomorphic(Graph* target)
{
    // If both graphs have no vertices, they are isomorphic.
    if (getNumberVertices() == 0 && target->getNumberVertices() == 0) return true;

    // Otherwise, if the target has no vertices, they cannot be isomorphic.
    if (target->getNumberVertices() == 0) return false;

    // Go through all of the vertices in the graph and search for a match to the first vertex in the target.
    for (int i=0; i<getNumberVertices(); i++)
    {
        // See if the vertices match.
        if (isIsomorphicFromVertices(target, getVertex(i), target->getVertex(0)))
        {
            return true;
        }
    }
    return false;
}

bool Graph::isIsomorphicFromVertices(Graph* target, Vertex* sourceVertex, Vertex* targetVertex, bool permitSubgraph, GraphMapping* mapping, Vertex* sourceVertexOrigin, Vertex* targetVertexOrigin)
{
    //printf("Checking for isomorphic subgraph %s and %s, vertices %s and %s\n",getString().c_str(), target->getString().c_str(), sourceVertex->getString().c_str(), targetVertex->getString().c_str()); fflush(stdout);

    // If this is the first call, clear the marks for graphs.
    bool isFirstVertex = (sourceVertexOrigin == NULL && targetVertexOrigin == NULL);
    if (isFirstVertex) target->clearMark("isIsomorphicFromVertices");

    // Mark that we have visited this vertex.
    targetVertex->setMark("isIsomorphicFromVertices");

    // If the vertices don't match, return false.
    if (!sourceVertex->matches(targetVertex)) return false;

    // Make sure that the same edges take us back to the vertices we came from.
    if (!isFirstVertex && sourceVertex->findEdgeLeadingTo(sourceVertexOrigin) != targetVertex->findEdgeLeadingTo(targetVertexOrigin)) return false;

    // Check each child edge.
    for (int i=0; i<targetVertex->getMaxNumberEdges(); i++)
    {
        // If the target has an edge.
        if (targetVertex->getEdge(i) != NULL)
        {
            // If the source is missing an edge, return false.
            if (sourceVertex->getEdge(i) == NULL) return false;

            // Make sure that we have not already processed the other end of this edge.
             if (!targetVertex->getEdge(i)->hasMark("isIsomorphicFromVertices"))
             {
                // Recursively follow any edges.
                if (!isIsomorphicFromVertices(target, sourceVertex->getEdge(i), targetVertex->getEdge(i), permitSubgraph, mapping, sourceVertex, targetVertex))
                    return false;
             }
        }

        // If we are not permitting subgraphs, make sure if the target doesn't have an edge neither does the source.
        else if (!permitSubgraph && targetVertex->getEdge(i) == NULL)
        {
            // If the source has an edge, return false.
            if (sourceVertex->getEdge(i) != NULL) return false;
        }
    }

    //printf("Yes, is an isomorphic subgraph %X and %X\n",sourceVertex,targetVertex); fflush(stdout);
    if (mapping != NULL) mapping->addMapping(sourceVertex, targetVertex);
    return true;
}

list<Vertex*> Graph::findConnectedSubgraphs()
{
    // Clear the mark from the whole graph.
    clearMark("findConnectedSubgraphs");

    // Go through each vertex.
    list<Vertex*> anchors;
    for (int i=0; i<getNumberVertices(); i++)
    {
        // If the vertex is not marked.
        Vertex* v = getVertex(i);
        if (!v->hasMark("findConnectedSubgraphs"))
        {
            // Add it to the list.
            anchors.push_back(v);

            // Mark it and anything connected to it.
            v->setMarkOnConnected("findConnectedSubgraphs");
        }
    }

    return anchors;
}

GraphMapping Graph::findGraphMapping(Graph* target)
{
    // Find all of the subgraphs in common.
    list<GraphMapping> componentSubgraphs = findCommonSubgraphs(target);

    // Create a new mapping that is the union of all the common subgraphs.
    GraphMapping mapping(this, target);
    for (auto it=componentSubgraphs.begin(); it != componentSubgraphs.end(); it++)
    {
        GraphMapping subgraphMapping = *it;
        for (int i=0; i<subgraphMapping.getNumberMatches(); i++)
            mapping.addMapping(subgraphMapping.getSourceVertex(i), subgraphMapping.getTargetVertex(i));
    }

    // Find any source vertices that were unused.
    for (int i=0; i<getNumberVertices(); i++)
    {
        // See if we should ignore this vertex.
        if (!getVertex(i)->hasMark("findCommonSubgraphs"))
        {
            mapping.addUnusedSourceVertex(getVertex(i));
            getVertex(i)->setMark("findCommonSubgraphs");
        }
    }

    // Find any target vertices that were unused.
    for (int i=0; i<target->getNumberVertices(); i++)
    {
        // See if we should ignore this vertex.
        if (!target->getVertex(i)->hasMark("findCommonSubgraphs"))
        {
            mapping.addUnusedTargetVertex(target->getVertex(i));
            target->getVertex(i)->setMark("findCommonSubgraphs");
        }
    }

    return mapping;
}

list<GraphMapping> Graph::findCommonSubgraphs(Graph* target)
{
    list<GraphMapping> ret;

    // Clear any marks on the graphs.
    clearMark("findCommonSubgraphs");
    target->clearMark("findCommonSubgraphs");

    // Loop until we can't find any more subgraphs in common.
    while (true)
    {
        GraphMapping nextCommonSubgraph = findLargestCommonSubgraph(target, "findCommonSubgraphs");

        // If the subgraph was empty, we are done.
        if (nextCommonSubgraph.getNumberMatches() == 0) break;

        // Mark the vertices as in use.
        for (int i=0; i<nextCommonSubgraph.getNumberMatches(); i++)
        {
            nextCommonSubgraph.getSourceVertex(i)->setMark("findCommonSubgraphs");
            nextCommonSubgraph.getTargetVertex(i)->setMark("findCommonSubgraphs");
        }

        // Add the subgraph to the list.
        ret.push_back(nextCommonSubgraph);
    }

    return ret;
}

GraphMapping Graph::findLargestCommonSubgraph(Graph* target, string ignoreMark)
{
    GraphMapping maxSubgraph;

    //printf("Finding largest common subgraph between [%s] and [%s]\n",getString().c_str(),target->getString().c_str()); fflush(stdout);

    // Go through all of the vertices in the source.
    for (int i=0; i<getNumberVertices(); i++)
    {
        // See if we should ignore this vertex.
        if (ignoreMark != "" && getVertex(i)->hasMark(ignoreMark)) continue;

        // Go through all of the vertices in the target.
        for (int j=0; j<target->getNumberVertices(); j++)
        {
            // See if we should ignore this vertex.
            if (ignoreMark != "" && target->getVertex(j)->hasMark(ignoreMark)) continue;

            //printf("Finding largest common subgraph between [%s] anchor=[%s] and [%s] anchor=[%s]\n",getString().c_str(),getVertex(i)->getString().c_str(),target->getString().c_str(),target->getVertex(j)->getString().c_str()); fflush(stdout);

            // Check for a mapping between these two vertices.
            GraphMapping subgraph(this, target);
            isIsomorphicFromVertices(target, getVertex(i), target->getVertex(j), true, &subgraph);
            if (subgraph.getNumberMatches() > maxSubgraph.getNumberMatches())
                maxSubgraph = subgraph;

            // Check for a mapping between these two vertices.
            GraphMapping reverseSubgraph(this, target);
            target->isIsomorphicFromVertices(this, target->getVertex(j), getVertex(i), true, &reverseSubgraph);
            if (reverseSubgraph.getNumberMatches() > maxSubgraph.getNumberMatches())
            {
                reverseSubgraph.reverse();
                maxSubgraph = reverseSubgraph;
            }
        }
    }

    //printf("Largest common subgraph was:\n%s\n", maxSubgraph.getString().c_str());

    return maxSubgraph;
}

list<GraphMapping> Graph::findAllIsomorphicSubgraphs(Graph* subgraph)
{
    //printf("Finding all isomorphic subgraphs between %s and %s\n",getString().c_str(),subgraph->getString().c_str()); fflush(stdout);

    list<GraphMapping> ret;

    // Go through all of the vertices in the graph and search for a match to the first vertex in the target.
    for (int i=0; i<getNumberVertices(); i++)
    {
        GraphMapping mapping(this, subgraph);

        // See if the vertices match.
        if (isIsomorphicFromVertices(subgraph, getVertex(i), subgraph->getVertex(0), true, &mapping))
        {
            // TODO: check to ensure that this mapping is not a duplicate of a previous mapping.

            ret.push_back(mapping);
        }
    }

    //printf("Found all isomorphic subgraphs between %s and %s: %ld\n",getString().c_str(),subgraph->getString().c_str(), ret.size()); fflush(stdout);

    return ret;
}


}
}

