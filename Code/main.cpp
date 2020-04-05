/*
 * File:   main.cpp
 * Author: angel
 *
 * Created on April 5, 2014, 11:55 PM
 */

// C / C++ program for Dijkstra's shortest path algorithm for adjacency
// list representation of graph

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

// A structure to represent a node in adjacency list

struct AdjListNode {
    int dest;
    double length;
    int capacity;
    bool srcORsink;
    struct AdjListNode* next;
    //    struct AdjListNode* prev;
};

// A structure to represent an adjacency liat

struct AdjList {
    int prev;
    struct AdjListNode *head; // pointer to head node of list
};

// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)

struct Graph {
    int V;
    struct AdjList* array;
};

// A utility function to create a new adjacency list node

 struct AdjListNode* newAdjListNode(int dest, int weight) {
    struct AdjListNode* newNode =
            (struct AdjListNode*) malloc(sizeof (struct AdjListNode));
    double epsilon;
    newNode->dest = dest;
    newNode->length = 1.0;
    newNode->capacity = weight;
    newNode->next = NULL;
    //    newNode->prev = NULL;
    return newNode;
}

// A utility function that creates a graph of V vertices

struct Graph* createGraph(int V) {
    struct Graph* graph = (struct Graph*) malloc(sizeof (struct Graph));
    graph->V = V;

    // Create an array of adjacency lists.  Size of array will be V
    graph->array = (struct AdjList*) malloc(V * sizeof (struct AdjList));

    // Initialize each adjacency list as empty by making head as NULL
    for (int i = 0; i < V; ++i)
        graph->array[i].head = NULL;

    return graph;
}

double getEdgeLength(struct Graph* graph, int src, int dest) {
    AdjListNode* iteration = graph->array[src].head;
    while (iteration != NULL) {
        if (iteration->dest == dest) {

            return iteration->length;
        }
        iteration = iteration->next;
    }
    return -1;

}

int getEdgeCapacity(struct Graph* graph, int src, int dest) {
    AdjListNode* iteration = graph->array[src].head;

    while (iteration != NULL) {
        if (iteration->dest == dest) {

            return iteration->capacity;
        }
        iteration = iteration->next;
    }
    return -1;

}

// Adds an edge to an undirected graph

void addEdge(struct Graph* graph, int src, int dest, int weight) {
    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    if (-1 == getEdgeLength(graph, src, dest)) {
        struct AdjListNode* newNode = newAdjListNode(dest, weight);
        newNode->next = graph->array[src].head;
        graph->array[src].head = newNode;


    }
    else {
        AdjListNode* iteration = graph->array[src].head;
        while (iteration != NULL) {
            if (iteration->dest == dest) {

                iteration->capacity = weight;
                break;
            }
            iteration = iteration->next;
        }

    }



}

void updateEdgeLength(struct Graph* graph, int src, int dest, double newLength) {
    AdjListNode* iteration = graph->array[src].head;

    while (iteration != NULL) {
        if (iteration->dest == dest) {

            iteration->length = newLength;
            break;
        }
        iteration = iteration->next;
    }
}



// Structure to represent a min heap node

struct MinHeapNode {
    int v;
    double dist;
};

// Structure to represent a min heap

struct MinHeap {
    int size; // Number of heap nodes present currently
    int capacity; // Capacity of min heap
    int *pos; // This is needed for decreaseKey()
    struct MinHeapNode **array;
};

// A utility function to create a new Min Heap Node

struct MinHeapNode* newMinHeapNode(int v, double dist) {
    struct MinHeapNode* minHeapNode =
            (struct MinHeapNode*) malloc(sizeof (struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->dist = dist;
    return minHeapNode;
}

// A utility function to create a Min Heap

struct MinHeap* createMinHeap(int capacity) {
    struct MinHeap* minHeap =
            (struct MinHeap*) malloc(sizeof (struct MinHeap));
    minHeap->pos = (int *) malloc(capacity * sizeof (int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array =
            (struct MinHeapNode**) malloc(capacity * sizeof (struct MinHeapNode*));
    return minHeap;
}

// A utility function to swap two nodes of min heap. Needed for min heapify

void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b) {
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()

void minHeapify(struct MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size &&
            minHeap->array[left]->dist < minHeap->array[smallest]->dist)
        smallest = left;

    if (right < minHeap->size &&
            minHeap->array[right]->dist < minHeap->array[smallest]->dist)
        smallest = right;

    if (smallest != idx) {
        // The nodes to be swapped in min heap
        MinHeapNode *smallestNode = minHeap->array[smallest];
        MinHeapNode *idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is ampty or not

int isEmpty(struct MinHeap* minHeap) {
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap

struct MinHeapNode* extractMin(struct MinHeap* minHeap) {
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decreasy dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap

void decreaseKey(struct MinHeap* minHeap, int v, double dist) {
    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its dist value
    minHeap->array[i]->dist = dist;

    // Travel up while the complete tree is not heapified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not

bool isInMinHeap(struct MinHeap *minHeap, int v) {
    if (minHeap->pos[v] < minHeap->size)
        return true;
    return false;
}

// A utility function used to print the solution

void printArr(double dist[], int n) {
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; ++i)
        printf("%d \t\t %f\n", i, dist[i]);
}
//  int prev[100];
// The main function that calulates distances of shortest paths from src to all
// vertices. It is a O(ELogV) function

double dijkstra(struct Graph* graph, int src, int dest) {
    int V = graph->V; // Get the number of vertices in graph
    double dist[V]; // dist values used to pick minimum weight edge in cut
    //    int prev[V];
    // minHeap represents set E

    struct MinHeap* minHeap = createMinHeap(V);

    // Initialize min heap with all vertices. dist value of all vertices
    for (int v = 0; v < V; ++v) {
        dist[v] = INT_MAX;
        graph->array[v].prev = -1;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }

    // Make dist value of src vertex as 0 so that it is extracted first
    minHeap->array[src] = newMinHeapNode(src, dist[src]);
    minHeap->pos[src] = src;
    dist[src] = 0;
    decreaseKey(minHeap, src, dist[src]);

    // Initially size of min heap is equal to V
    minHeap->size = V;
    //alind
    //    graph->array[src].head->prev= NULL;
    // In the followin loop, min heap contains all nodes
    // whose shortest distance is not yet finalized.
    while (!isEmpty(minHeap)) {
        // Extract the vertex with minimum distance value
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v; // Store the extracted vertex number

        // Traverse through all adjacent vertices of u (the extracted
        // vertex) and update their distance values
        struct AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl != NULL) {
            int v = pCrawl->dest;

            // If shortest distance to v is not finalized yet, and distance to v
            // through u is less than its previously calculated distance
            if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX &&
                    pCrawl->length + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->length;
                graph->array[v].prev = u;
                //alind update the previous node
                //                 graph->array[v].head->prev= graph->array[u].head;

                // update distance value in min heap also
                decreaseKey(minHeap, v, dist[v]);
            }

            pCrawl = pCrawl->next;
        }
    }

    // print the calculated shortest distances
//    printArr(dist, V);

    return dist[dest];
}

void updatePath(struct Graph* graph, int sink, int source,double epsilon) {
    int p = sink;
    int minc = INT_MAX;

    do {

        if (getEdgeCapacity(graph, graph->array[p].prev, p) < minc) {

            minc = getEdgeCapacity(graph, graph->array[p].prev, p);
        }

        p = graph->array[p].prev;
    } while (p != source);
    p = sink;

    do {
        updateEdgeLength(graph, graph->array[p].prev, p, getEdgeLength(graph, graph->array[p].prev, p)*(1+epsilon*(minc/getEdgeCapacity(graph, graph->array[p].prev, p))));

        p = graph->array[p].prev;
    } while (p != source);
}

void MMF(struct Graph* graph, int* source, int* sink, int k) {
    double distr[k];
    double alpha = .99;
    int min = -1;
    double delta;
    double L=0;
    double epsilon=.09;
    double D=0;
    for (int i = 0; i < k; i++) {
            distr[i] = dijkstra(graph, source[i], sink[i]);
            if (distr[i] > L) {
                L = distr[i];
//                min = i;
            }
        }
    L = 3;

    for (int i = 0; i < graph->V; i++)
        for (int j = 0; j < graph->V; j++)
            if(getEdgeCapacity(graph,i,j)!= -1)
            updateEdgeLength(graph,i,j,(1+epsilon)*pow((1+epsilon)*L,-1/epsilon));




    while (alpha < 1) {
        for (int i = 0; i < k; i++) {
            distr[i] = dijkstra(graph, source[i], sink[i]);

            if (distr[i] < alpha || i == 0) {

                alpha = distr[i];
                min = i;
            }
        }
        dijkstra(graph, source[min], sink[min]);

        updatePath(graph, sink[min], source[min], epsilon);
    }
    fflush(stdout);

fflush(stdout);
  for (int i = 0; i < graph->V; i++)
       for (int j = 0; j < graph->V; j++)
           if(getEdgeCapacity(graph,i,j)!= -1)
               D +=  getEdgeCapacity(graph,i,j)*getEdgeLength(graph,i,j);
       printf("Alpha = %f",D);
}
// Driver program to test above functions

int main() {
    // create the graph given in above fugure
    int V = 14;
    int a[] = {0,1,2};
    int b[] = {11,12,13};
    struct Graph* graph = createGraph(V);
//    addEdge(graph, 0, 1, 4);
//    addEdge(graph, 0, 7, 8);
//    addEdge(graph, 1, 2, 8);
//    addEdge(graph, 1, 7, 11);
//    addEdge(graph, 2, 3, 7);
//    addEdge(graph, 2, 8, 2);
//    addEdge(graph, 2, 5, 4);
//    addEdge(graph, 3, 4, 9);
//    addEdge(graph, 3, 5, 14);
//    addEdge(graph, 4, 5, 10);
//    addEdge(graph, 5, 6, 2);
//    addEdge(graph, 6, 7, 1);
//    addEdge(graph, 6, 8, 6);
//    addEdge(graph, 7, 8, 7);
//
//        addEdge(graph, 0, 1, 5);
//    addEdge(graph, 0, 2, 30);
//    addEdge(graph, 4, 2, 30);
//    addEdge(graph, 4, 5, 30);
//    addEdge(graph, 2, 3, 10);
//    addEdge(graph, 3,5,30);
//    addEdge(graph, 3,1,30);

    addEdge(graph,0,3,5);
    addEdge(graph,0,5,5);
    addEdge(graph,1,4,5);
    addEdge(graph,1,6,5);
    addEdge(graph,2,5,5);
    addEdge(graph,2,6,15);
    addEdge(graph,3,7,5);
    addEdge(graph,4,7,5);
    addEdge(graph,4,8,5);
    addEdge(graph,4,9,5);
    addEdge(graph,5,8,5);
    addEdge(graph,5,9,5);
    addEdge(graph,6,9,10);
    addEdge(graph,6,10,10);
    addEdge(graph,7,11,10);
    addEdge(graph,8,11,10);
    addEdge(graph,8,12,10);
    addEdge(graph,9,12,10);
    addEdge(graph,9,13,10);
    addEdge(graph,10,12,10);
    addEdge(graph,10,13,10);

    MMF(graph, a, b, 3);
    // addEdge(graph, 7, 8, 5);
    //    addEdge(graph, 2, 8, 24);

    //    dijkstra(graph, 0);

    //    printf("\n\n%d",graph->array[4].head->prev->dest);
    //    printf("\n\n%d",graph->array[5].head->dest);
    return 0;
}
