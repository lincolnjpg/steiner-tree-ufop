#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstdio>
#include <climits>
#include <chrono>

/*Mensages*/
const char MSG_PARAMS_ERROR[] =
  "Numero de parametros esta incorreto (informe apenas 1 arquivo de entrada)";
const char MSG_READFILE_ERROR[] =
  "Operacao de leitura do arquivo falhou. Codigo do erro:";
const char MSG_OPEN_FILE_ERROR[] =
  "Nao foi possivel abrir o arquivo especificado";

using namespace std;
using namespace std::chrono;

/*Structs*/
typedef struct tShortestPathInfo
{
  unsigned int weight;
  vector<int> path;
} ShortestPathInfo;

typedef struct tVertex
{
  int id = -1;
  vector<ShortestPathInfo> shortestPathInfo;
} Vertex;

typedef struct tEdge
{
  int vertexId;
  int vertex2Id;
  unsigned int weight;
  vector<int> realPath;
} Edge;

typedef struct tDijkstraVertex : Vertex
{
  unsigned int distanceFromSource;
  int previousVertexId;
} DijkstraVertex;

typedef struct tPointerVertex
{
  int id;
  DijkstraVertex* realVertex;
} PointerVertex;

typedef struct tAdjacencyInfo
{
  int id;
  unsigned int weight;
} AdjacencyInfo;

typedef struct tAdjacencies
{
  Vertex vertex;
  vector<AdjacencyInfo> adjacencies;
} Adjacencies;

typedef struct tGraph
{
  vector<Adjacencies>* adjacencyList;
  vector<Adjacencies*>* terminalList;
} Graph;

typedef struct tSteinerTree : Graph
{
  unsigned int totalWeight;
} SteinerTree;

/*Typedefs*/
typedef vector<vector<ShortestPathInfo>> ShortestPaths;

/*Declaracao de funcoes*/
unsigned int readFile(FILE*, Graph**);
unsigned int readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);
bool compare(DijkstraVertex*&, DijkstraVertex*&);
ShortestPaths* dijkstra(Graph*);
Graph* createCompleteGraph(ShortestPaths*);
void insertMST(Graph*, ShortestPaths*, const int&, const unsigned int&,
               const int&, const int&, vector<int>);
void removeMST(Graph*, const unsigned int&, const unsigned int&);
bool compareKruskal(Edge&, Edge&);
Graph* generateMST(Graph*, map<int, int>*);
bool hasCycle(Graph*, set<int>*, map<int, int>*, int, int);
void mstEdgeAdjustment(Graph*, ShortestPaths*, map<int, int>*);
void insertSteiner(SteinerTree*, const int&, const unsigned int&, const int&,
                   const unsigned int&);
SteinerTree* generateSteinerTree(Graph*, Graph*);
void printResults(SteinerTree*);
//test
void output_graphviz(SteinerTree*, Graph*);

/*Funcao principal*/
int main(int argc, char* argv[])
{
  unsigned int returnedValue;
  Graph* inputGraph = nullptr;
  Graph* terminalsCompleteGraph = nullptr;
  Graph* mst = nullptr;
  map<int, int> distinctVertices;
  SteinerTree* steinerTree = nullptr;

  ShortestPaths* dijkstraResult;
  FILE* file = nullptr;

  if (argc == 2)
  {
    file = fopen(argv[1], "r");
  }
  else
  {
    printf("%s\n", MSG_PARAMS_ERROR);

    return -1;
  }

  if (!ferror(file))
  {
    returnedValue = readFile(file, &inputGraph);

    if (returnedValue)
      printf("%s %u\n", MSG_READFILE_ERROR, returnedValue);

    //trecho medicao tempo (1) - inicio
    duration<double> time_span2;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    //trecho medicao tempo (1) - fim

    dijkstraResult = dijkstra(inputGraph);

    terminalsCompleteGraph = createCompleteGraph(dijkstraResult);
    //Algoritmo de Kruskal
    mst = generateMST(terminalsCompleteGraph, &distinctVertices);
    //Ajusta arestas da MST, de acordo com o caminho mais curto
    //mstEdgeAdjustment(mst, dijkstraResult, &distinctVertices);
    //Gera arvore de Steiner
    steinerTree = generateSteinerTree(mst, inputGraph);

    //trecho medicao tempo (2) - inicio
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    //trecho medicao tempo (2) - fim

    cout << time_span.count() << endl;

    printResults(steinerTree);

    //test
    output_graphviz(steinerTree, inputGraph);

    fclose(file);

    delete inputGraph;
    delete terminalsCompleteGraph;

    return returnedValue;
  }
  else
  {
    printf("%s\n", MSG_OPEN_FILE_ERROR);

    return -2;
  }
}

/*Funcao principal de leitura do arquivo*/
unsigned int readFile(FILE* file, Graph** inputGraph)
{
  char aux[100];
  unsigned int numVertices;
  unsigned int numEdges;
  unsigned int numTerminals;
  unsigned int returnedValue;
  bool error;

  //Le linha de inicio da secao "Grafo"
  fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);

  if (ferror(file))
    return 1;

  //Le numero de vertices
  fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numVertices);

  if (ferror(file))
    return 2;

  //Le numero de arestas
  fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numEdges);

  if (ferror(file))
    return 3;

  //Cria o grafo, de acordo com numero de vertices lido
  *inputGraph = new Graph;
  (*inputGraph)->adjacencyList = new vector<Adjacencies>(numVertices);
  (*inputGraph)->terminalList = nullptr;

  //Chama funcao que realiza a leitura das arestas
  returnedValue = readEdges(file, *inputGraph, numEdges);

  if (returnedValue)
    return returnedValue;

  //Descarta "End"
  fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);

  if (ferror(file))
    return 4;

  //Le linha de inicio da secao "Terminais"
  fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);

  if (ferror(file))
    return 5;

  //Le numero de terminais
  fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numTerminals);

  if (ferror(file))
    return 6;

  //Chama funcao que realiza a leitura dos terminais
  error = readTerminals(file, *inputGraph, numTerminals);

  if (error)
    return 7;

  //Sem erro na manipulacao do arquivo
  return 0;
}

/*Funcao responsavel pela leitura das arestas*/
unsigned int readEdges(FILE* file, Graph* inputGraph,
                       const unsigned int& numEdges)
{
  char aux[100];
  int vertex, vertex2;
  unsigned int weight;

  AdjacencyInfo auxAdjacencyInfo;

  for (unsigned int i = 0; i < numEdges; i++)
  {
    fscanf(file, "%c %d %d %u %*[\r] %*[\n]", aux, &vertex, &vertex2, &weight);

    if (ferror(file))
      return 1;

    //vertice u -> vertice v
    auxAdjacencyInfo.id = vertex2;
    auxAdjacencyInfo.weight = weight;

    if (inputGraph->adjacencyList->at(vertex - 1).adjacencies.size() == 0)
      inputGraph->adjacencyList->at(vertex - 1).vertex.id = vertex;

    inputGraph->adjacencyList->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);

    //vertice v -> vertice u
    auxAdjacencyInfo.id = vertex;
    auxAdjacencyInfo.weight = weight;

    if (inputGraph->adjacencyList->at(vertex2 - 1).adjacencies.size() == 0)
      inputGraph->adjacencyList->at(vertex2 - 1).vertex.id = vertex2;

    inputGraph->adjacencyList->at(vertex2 - 1).adjacencies.push_back(auxAdjacencyInfo);
  }

  /*Verifica se grafo eh conexo*/
  for (unsigned int i = 0; i < inputGraph->adjacencyList->size(); i++)
    if (inputGraph->adjacencyList->at(i).vertex.id == -1)
      return 2;

  return 0;
}

/*Funcao responsavel pela leitura dos terminais*/
bool readTerminals(FILE* file, Graph* inputGraph, const unsigned int& numTerminals)
{
  char aux[100];
  int vertex;

  inputGraph->terminalList = new vector<Adjacencies*>(numTerminals);

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    fscanf(file, "%c %d %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      return true;

    inputGraph->terminalList->at(i) = &(inputGraph->adjacencyList->at(vertex - 1));
  }

  return false;
}

/*Funcao auxiliar que usada na criacao/ajuste do heap binario
para o algoritmo de Dijkstra*/
bool compare(DijkstraVertex*& ptr, DijkstraVertex*& ptr2)
{
  return ptr->distanceFromSource > ptr2->distanceFromSource;
}

/*Funcao que executa o algoritmo de Dijkstra*/
ShortestPaths* dijkstra(Graph* inputGraph)
{
  unsigned int newDistance;
  DijkstraVertex auxVertex;
  AdjacencyInfo neighbourInfo;
  vector<DijkstraVertex*> fakeVertices(inputGraph->adjacencyList->size());
  vector<DijkstraVertex> openVertices(inputGraph->adjacencyList->size());
  ShortestPaths* result = new ShortestPaths(inputGraph->terminalList->size());

  /*Inicializacao do conjunto que representa os vertices abertos*/
  for (unsigned int i = 0; i < inputGraph->adjacencyList->size(); i++)
  {
    /*Inicializa ID do vertice*/
    auxVertex.id = i + 1;    
    /*Insere vertice corrente na lista de vertices abertos*/
    //openVertices.push_back(auxVertex); //test
    openVertices.at(i) = auxVertex;
    /*Vertice aberto i passa a apontar para o vertice Dijkstra i*/
    fakeVertices.at(i) = &openVertices.at(i);
  } 

  /*Executa algoritmo de Dijkstra para cada um dos vertices terminais*/
  for (unsigned int i = 0; i < inputGraph->terminalList->size(); i++)
  {
    /*Armazena ID do vertice (terminal) atual*/
    int currentTerminalId = inputGraph->terminalList->at(i)->vertex.id;

    for (unsigned int i = 0; i < inputGraph->adjacencyList->size(); i++)
    {
      /*Indica que a distancia entre o vertice corrente e o vertice
      de origem ainda nao foi definida*/
      openVertices.at(i).distanceFromSource = UINT_MAX;
      /*Indica que o ID do vertice anterior ao vertice corrente
      ainda nao foi definido (incializado com o valor do proprio ID)*/
      openVertices.at(i).previousVertexId = i + 1;
    }

    /*Fara com que o vertice (terminal) atual seja usado como vertice de origem*/
    //fakeVertices.at(currentTerminalId - 1)->distanceFromSource = 0;
    openVertices.at(currentTerminalId - 1).distanceFromSource = 0;
    /*Variavel que ajuda no controle do tamanho do conjunto de vertices abertos*/
    unsigned int j = 0;

    while (fakeVertices.size() - j > 0)
    {
      /*Organiza lista de forma que ela comportesse como um heap binario minimo*/
      make_heap(fakeVertices.begin(), fakeVertices.end() - j, compare);
      /*Sempre considera o primeiro vertice do heap*/
      auxVertex = *(fakeVertices.at(0)); //teste - pegar o primeiro terminal, ao inves desse

      for (unsigned int k = 0;
           k < inputGraph->adjacencyList->at(auxVertex.id - 1).adjacencies.size();
           k++)
      {
        neighbourInfo = inputGraph->adjacencyList->at(auxVertex.id - 1).
          adjacencies.at(k);

        /*Evita calcular distancia duas vezes*/
        if (neighbourInfo.id != auxVertex.previousVertexId)
          newDistance = auxVertex.distanceFromSource + neighbourInfo.weight;

        if (newDistance < openVertices.at(neighbourInfo.id - 1).distanceFromSource)
        {          
          openVertices.at(neighbourInfo.id - 1).distanceFromSource = newDistance;
          openVertices.at(neighbourInfo.id - 1).previousVertexId = auxVertex.id;
        }
      }

      /*"Remove" primeiro objeto do heap*/
      pop_heap(fakeVertices.begin(), fakeVertices.end() - j);
      /*Essa variavel limita o numero de vertices abertos (pop_heap nao remove - apenas
      joga para o final da lista*/
      j++;
    }

    ShortestPathInfo auxShortestPathInfo;

    /*Laco que identifica rota e seu custo minimo*/
    for (unsigned int j = 0; j < inputGraph->terminalList->size(); j++)
    {
      /*Variavel que indica o ID do terminal para o qual deseja-se
      calcular a distancia menos custosa*/
      int nextId = inputGraph->terminalList->at(j)->vertex.id;
      /*Nao ha necessidade de calcular rota para si proprio*/
      if (currentTerminalId != nextId)
      {
        /*Obtem custo da distancia entre o terminal corrente e seu adjacente (outro terminal)*/
        auxShortestPathInfo.weight = openVertices.at(nextId - 1).distanceFromSource;
        /*Certifica-se de que a rota 'i' nao se confunda com a rota 'j'*/
        auxShortestPathInfo.path.clear();

        //while (nextId != graph->terminalList->at(i + 1)->vertex.id)
        while (nextId != openVertices.at(nextId - 1).previousVertexId)
        {   
          /*Insere ID no vertice na rota*/
          auxShortestPathInfo.path.push_back(nextId);
          /*Atualiza*/
          nextId = openVertices.at(nextId - 1).previousVertexId;
        }

        /*Insere ID do vertice de origem*/
        auxShortestPathInfo.path.push_back(nextId);
        /*Insere informacao da rota na lista de terminais*/
        result->at(i).push_back(auxShortestPathInfo);
      }
    }
  }

  return result;
}

/*Funcao que cria um grafo completo, em que os vertices sao os terminais e o
peso das arestas que os une eh a distancia de menor custo entre cada um dos
terminais (valor calculado via Dijkstra)*/
Graph* createCompleteGraph(ShortestPaths* dijkstraResult)
{
  Graph* completeGraph = new Graph;
  completeGraph->adjacencyList = new vector<Adjacencies>(dijkstraResult->size());
  AdjacencyInfo adjacencyInfo;
  ShortestPathInfo pathInfo;
  unsigned int numTerminals = dijkstraResult->size();
  unsigned int numAdjacencies;
  unsigned int firstVertex;
  unsigned int lastVertex;

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    numAdjacencies = dijkstraResult->at(i).size();

    for (unsigned int j = 0; j < numAdjacencies; j++)
    {
      firstVertex = dijkstraResult->at(i).at(j).path.at(dijkstraResult->at(i).
        at(j).path.size() - 1);
      lastVertex = dijkstraResult->at(i).at(j).path.at(0);
      completeGraph->adjacencyList->at(i).vertex.id = firstVertex;
      adjacencyInfo.id = lastVertex;
      adjacencyInfo.weight = dijkstraResult->at(i).at(j).weight;
      completeGraph->adjacencyList->at(i).adjacencies.push_back(adjacencyInfo);
      pathInfo.path = dijkstraResult->at(i).at(j).path;
      pathInfo.weight = dijkstraResult->at(i).at(j).weight;
      completeGraph->adjacencyList->at(i).vertex.shortestPathInfo.push_back(
        pathInfo);
    }
  }

  return completeGraph;
}

/*Funcao que insere um vertice na MST*/
/*void insertMST(Graph* mst, const int& vertex, const unsigned int& index,
               const int& vertex2, const int& weight, vector<int> realPath)
*/
void insertMST(Graph* mst, Edge& edge, const unsigned int& index,
               const unsigned int& index2)
{
  AdjacencyInfo adjacencyInfo;
  ShortestPathInfo pathInfo;

  //u - v
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = edge.vertex2Id;
  adjacencyInfo.weight = edge.weight;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  mst->adjacencyList->at(index).vertex.id = edge.vertexId;
  pathInfo.path = edge.realPath;
  pathInfo.weight = edge.weight;
  mst->adjacencyList->at(index).vertex.shortestPathInfo.push_back(pathInfo);
  mst->adjacencyList->at(index).adjacencies.push_back(adjacencyInfo);

  //v - u
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = edge.vertexId;
  adjacencyInfo.weight = edge.weight;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  mst->adjacencyList->at(index2).vertex.id = edge.vertex2Id;
  reverse(edge.realPath.begin(), edge.realPath.end());
  pathInfo.path = edge.realPath;
  pathInfo.weight = edge.weight;
  mst->adjacencyList->at(index2).vertex.shortestPathInfo.push_back(pathInfo);
  mst->adjacencyList->at(index2).adjacencies.push_back(adjacencyInfo);
}

/*Funcao que remove um vertice da MST*/
void removeMST(Graph* mst, const unsigned int& index,
               const unsigned int& index2)
{
  if (mst->adjacencyList->at(index).adjacencies.size() > 0)
  {
    mst->adjacencyList->at(index).adjacencies.pop_back();

    if (mst->adjacencyList->at(index).adjacencies.size() == 0)
      mst->adjacencyList->at(index).vertex.id = -1;
  }

  if (mst->adjacencyList->at(index).vertex.shortestPathInfo.size() > 0)
    mst->adjacencyList->at(index).vertex.shortestPathInfo.pop_back();

  if (mst->adjacencyList->at(index2).adjacencies.size() > 0)
  {
    mst->adjacencyList->at(index2).adjacencies.pop_back();

    if (mst->adjacencyList->at(index2).adjacencies.size() == 0)
      mst->adjacencyList->at(index2).vertex.id = -1;
  }

  if (mst->adjacencyList->at(index2).vertex.shortestPathInfo.size() > 0)
    mst->adjacencyList->at(index2).vertex.shortestPathInfo.pop_back();
}

/*Funcao auxiliar que eh usada na criacao/ajuste do heap binario
para o algoritmo de Kruskal*/
bool compareKruskal(Edge& edge, Edge& edge2)
{
  return edge.weight < edge2.weight;
}

/*Função que verifica se um grafo possui algum ciclo, atraves
do uso do algoritmo DFS*/
bool hasCycle(Graph* mst, set<int>* visitedVertices,
              map<int, int>* distinctVertices, int sourceIndex,
              int parentId)
{
  unsigned int i = 0;
  int vertexId = mst->adjacencyList->at(sourceIndex).vertex.id;
  int neighbourId;
  int index;

  visitedVertices->insert(vertexId);

  while (i < mst->adjacencyList->at(sourceIndex).adjacencies.size())
  {
    neighbourId = mst->adjacencyList->at(sourceIndex).adjacencies.at(i).id;
    index = distinctVertices->find(neighbourId)->second;
    i++;

    if (visitedVertices->find(neighbourId) == visitedVertices->end())
    {
      /*Vizinho corrente ainda nao foi visitado*/

      /*Vizinho corrente eh inserido no conjunto de vertices visitados*/
      visitedVertices->insert(neighbourId);
      /*Chamada recursiva*/
      if (hasCycle(mst, visitedVertices, distinctVertices, index, vertexId))
        /*Interrompe busca caso seja detectado*/
        return true;
    }
    else if (neighbourId != parentId)
    {
      /*Vizinho corrente ja foi visitado e foi alcancado
      a partir de outro vertice: ciclo detectado*/
      return true;
    }
  }

  return false;
}

Graph* generateMST(Graph* completeGraph, map<int, int>* distinctVertices)
{
  //TODO: distinctVertices pode ser local... era pointeiro pq usava em outra funcao... usava isso pra preencher mst com shortest path (agora vem de completeGraph)
  /*Usar essa ideia de 'size' na funcao acima...*/
  unsigned int size = completeGraph->adjacencyList->size();
  vector<Edge>* edgeList = new vector<Edge>;
  Graph* mst = new Graph;
  mst->adjacencyList = new vector<Adjacencies>(size);
  Edge edge;

  set<int> visitedVertices;
  pair<map<int, int>::iterator, bool> resultPair;
  pair<map<int, int>::iterator, bool> resultPair2;
  unsigned int vertexIndex;
  unsigned int vertex2Index;

  /*Gera lista de arestas existentes no grafo completo, recebido como parametro*/
  for (unsigned int i = 0; i < size; i++)
  {
    for (unsigned int j = i + 1; j < size; j++)
    {
      edge.vertexId = completeGraph->adjacencyList->at(i).vertex.id;
      edge.vertex2Id = completeGraph->adjacencyList->at(i).adjacencies.at(j - 1).id;
      edge.weight = completeGraph->adjacencyList->at(i).adjacencies.at(j - 1).weight;
      edge.realPath = completeGraph->adjacencyList->at(i).vertex.
          shortestPathInfo.at(j - 1).path;
      edgeList->push_back(edge);
    }
  }

  /*Ordena, de forma crescente, a lista de arestas*/
  sort(edgeList->begin(), edgeList->end(), compareKruskal);

  unsigned int i = 0;
  unsigned int insertedEdges = 0;
  unsigned int insertionIndex = 0;

  /*Numa arvore, n = m - 1*/
  while (insertedEdges < size - 1)
  {
    resultPair = distinctVertices->insert(
      pair<int, int>(edgeList->at(i).vertexId, insertionIndex));
    vertexIndex = resultPair.first->second;

    if (resultPair.second)
      insertionIndex++;

    resultPair2 = distinctVertices->insert(
      pair<int, int>(edgeList->at(i).vertex2Id, insertionIndex));
    vertex2Index = resultPair2.first->second;

    if (resultPair2.second)
      insertionIndex++;

    /*Atualiza arvore geradora minima com nova aresta e vertice(s)*/
    insertMST(mst, edgeList->at(i), vertexIndex, vertex2Index);
    /*
    insertMST(mst, edgeList->at(i).vertexId, vertexIndex,
              edgeList->at(i).vertex2Id, edgeList->at(i).weight);
    insertMST(mst, edgeList->at(i).vertex2Id, vertex2Index,
              edgeList->at(i).vertexId, edgeList->at(i).weight);
    */

    visitedVertices.clear();

    if (hasCycle(mst, &visitedVertices, distinctVertices, vertexIndex, -1))
    {
      /*Ciclo detectado*/
      removeMST(mst, vertexIndex, vertex2Index);

      if (resultPair.second)
        /*Se vertex foi inserido, sera removido*/
        distinctVertices->erase(resultPair.first);

      if (resultPair2.second)
        /*Se vertex2 foi inserido, sera removido*/
        distinctVertices->erase(resultPair2.first);

      //insertedEdges--;
    }
    else
    {
      insertedEdges++;
    }

    i++;
  }

  return mst;
}

void insertSteiner(SteinerTree* steinerTree, const int& vertexId,
                   const unsigned int& index, const int& vertex2Id,
                   const unsigned int& index2)
{
  AdjacencyInfo adjacencyInfo;

  //u - v
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = vertex2Id;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  steinerTree->adjacencyList->at(index).vertex.id = vertexId;
  steinerTree->adjacencyList->at(index).adjacencies.push_back(adjacencyInfo);

  //v - u
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = vertexId;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  steinerTree->adjacencyList->at(index2).vertex.id = vertex2Id;
  steinerTree->adjacencyList->at(index2).adjacencies.push_back(adjacencyInfo);
}

SteinerTree* generateSteinerTree(Graph* mst, Graph* inputGraph)
{
  SteinerTree* steinerTree;
  Edge edge;
  set<int> visitedVertices;
  map<int, int> distinctVertices;
  pair<map<int, int>::iterator, bool> resultPair;
  pair<map<int, int>::iterator, bool> resultPair2;
  unsigned int numPaths;
  unsigned int pathSize;
  unsigned int weight;
  unsigned int vertexIndex;
  unsigned int vertex2Index;
  unsigned int insertionIndex = 0;
  int vertexId;
  int vertex2Id;

  steinerTree = new SteinerTree;
  steinerTree->adjacencyList = new vector<Adjacencies>(inputGraph->adjacencyList->size());
  steinerTree->totalWeight = 0;

  for (unsigned int i = 0; i < mst->adjacencyList->size(); i++)
  {
    numPaths = mst->adjacencyList->at(i).vertex.shortestPathInfo.size();

    for (unsigned int j = 0; j < numPaths; j++)
    {
      pathSize = mst->adjacencyList->at(i).vertex.shortestPathInfo.at(j).path.size();
      vector<int>& path = mst->adjacencyList->at(i).vertex.shortestPathInfo.at(j).path;

      for (unsigned int k = 0; k < pathSize - 1; k++)
      {
        resultPair = distinctVertices.insert(
          pair<int, int>(path.at(k), insertionIndex));
        vertexIndex = resultPair.first->second;

        if (resultPair.second)
          insertionIndex++;

        resultPair2 = distinctVertices.insert(
          pair<int, int>(path.at(k + 1), insertionIndex));
        vertex2Index = resultPair2.first->second;

        if (resultPair2.second)
          insertionIndex++;

        insertSteiner(steinerTree, path.at(k), vertexIndex, path.at(k + 1), vertex2Index);

        visitedVertices.clear();

        if (hasCycle(steinerTree, &visitedVertices, &distinctVertices, vertexIndex, -1))
        {
          /*Ciclo detectado*/
          removeMST(steinerTree, vertexIndex, vertex2Index);

          if (resultPair.second)
            /*Se vertex foi inserido, sera removido*/
            distinctVertices.erase(resultPair.first);

          if (resultPair2.second)
            /*Se vertex2 foi inserido, sera removido*/
            distinctVertices.erase(resultPair2.first);
        }
      }
    }
  }

  for (unsigned int i = 0; i < steinerTree->adjacencyList->size(); i++)
  {
    if (steinerTree->adjacencyList->at(i).vertex.id != -1)
    {
      vertexId = steinerTree->adjacencyList->at(i).vertex.id;

      for (unsigned int j = 0; j < steinerTree->adjacencyList->at(i).adjacencies.size(); j++)
      {
        vertex2Id = steinerTree->adjacencyList->at(i).adjacencies.at(j).id;

        for (unsigned int k = 0; k < inputGraph->adjacencyList->at(vertexId - 1).adjacencies.size(); k++)
        {
          if (inputGraph->adjacencyList->at(vertexId - 1).adjacencies.at(k).id == vertex2Id)
          {
            weight = inputGraph->adjacencyList->at(vertexId - 1).adjacencies.at(k).weight;
            steinerTree->totalWeight += weight;

            break;
          }
        }
      }
    }
  }

  /*O ajuste abaixo eh necessario, pois o algoritmo percorre todas as adjacencias,
  adicionando o peso de cada aresta duas vezes*/
  steinerTree->totalWeight /= 2;

  return steinerTree;
}

void printResults(SteinerTree* steinerTree)
{
  unsigned int i = 0;
  int vertexId;
  set<int> distinctVertices;

  while (steinerTree->adjacencyList->at(i).vertex.id != -1)
  {
    vertexId = steinerTree->adjacencyList->at(i).vertex.id;
    distinctVertices.insert(vertexId);

    for (unsigned int j = 0; j < steinerTree->adjacencyList->at(i).adjacencies.size(); j++)
    {
      vertexId = steinerTree->adjacencyList->at(i).adjacencies.at(j).id;
      distinctVertices.insert(vertexId);
    }

    i++;
  }

  printf("Vertices (%u): ", distinctVertices.size());

  for (set<int>::iterator it = distinctVertices.begin(); it != distinctVertices.end(); it++)
    printf("%d, ", *it);

  printf("Custo: %d\n", steinerTree->totalWeight);
}

void output_graphviz(SteinerTree* steinerTree, Graph* inputGraph)
{
  unsigned int i = 0;
  int vertexId;
  int vertex2Id;
  set<string> distinctPairs;
  FILE* file;
  string str, str2;

  file = fopen("output.dot", "w");

  while (steinerTree->adjacencyList->at(i).vertex.id != -1)
  {
    vertexId = steinerTree->adjacencyList->at(i).vertex.id;

    for (unsigned int j = 0; j < steinerTree->adjacencyList->at(i).adjacencies.size(); j++)
    {
      vertex2Id = steinerTree->adjacencyList->at(i).adjacencies.at(j).id;
      str = to_string(vertexId) + " -- " + to_string(vertex2Id);
      str2 = to_string(vertex2Id) + " -- " + to_string(vertexId);

      if (distinctPairs.find(str2) == distinctPairs.end())
        distinctPairs.insert(str);
    }

    i++;
  }

  unsigned int index;
  int v1, v2, weight;

  fprintf(file, "graph G\n{\nnode [shape = circle width=0.3 fixedsize=true fontsize=10]\n");

  for (set<string>::iterator it = distinctPairs.begin(); it != distinctPairs.end(); it++)
  {
     index= it->find(" -- ");
     v1 = stoi(it->substr(0, it->size() - (it->size() - (index + 1)) - 1));
     v2 = stoi(it->substr(index + 4, it->size() - (index + 3)));

     for (unsigned int i = 0; i < inputGraph->adjacencyList->at(v1 -1).adjacencies.size(); i++)
     {
       if (inputGraph->adjacencyList->at(v1 - 1).adjacencies.at(i).id == v2)
       {
         weight = inputGraph->adjacencyList->at(v1 - 1).adjacencies.at(i).weight;

         fprintf(file, "%d -- %d [label=%d width=0.3 fontsize=10]\n", v1, v2, weight);

         break;
       }
     }
  }

  fprintf(file, "}\n");

  fclose(file);
}
