#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdio>
#include <climits>

/*Mensages*/
const char MSG_PARAMS_ERROR[] =
  "Numero de parametros esta incorreto (informe apenas 1 arquivo de entrada)";
const char MSG_READFILE_ERROR[] =
  "Operacao de leitura do arquivo falhou. Codigo do erro:";
const char MSG_OPEN_FILE_ERROR[] =
  "Nao foi possivel abrir o arquivo especificado";

using namespace std;

/*Structs*/
typedef struct tAdjacencyInfo
{
  int id;
  unsigned int weight;
} AdjacencyInfo;

typedef struct tVertex
{
  int id = -1;
  bool isTerminal;
} Vertex;

typedef struct tEdge
{
  int vertexId;
  int vertex2Id;
  unsigned int weight;
} Edge;

typedef struct tMinimalRouteInfo
{
  unsigned int weight;
  vector<int> route;
} MinimalRouteInfo;

typedef struct tDijkstraVertex : Vertex
{
  unsigned int distanceFromSource;
  int previousVertexId;
  vector<MinimalRouteInfo> minimalRouteInfo;
} DijkstraVertex;

typedef struct tPointerVertex
{
  int id;
  DijkstraVertex* realVertex;
} PointerVertex;

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

/*Declaracao de funcoes*/
unsigned int readFile(FILE*, Graph**);
unsigned int readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);
bool compare(DijkstraVertex*&, DijkstraVertex*&);
vector<vector<MinimalRouteInfo>>* shortestPath(Graph*);
Graph* createCompleteGraph(vector<vector<MinimalRouteInfo>>*);
void generateMST(Graph*);
bool hasCycle(Graph*);

/*Funcao principal*/
int main(int argc, char* argv[])
{
  unsigned int returnedValue;
  Graph* inputGraph = NULL;
  Graph* terminalsCompleteGraph = NULL;

  vector<vector<MinimalRouteInfo>>* dijkstraResult;
  FILE* file = NULL;

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

    dijkstraResult = shortestPath(inputGraph);

    terminalsCompleteGraph = createCompleteGraph(dijkstraResult);
    //Algoritmo de Kruskal
    generateMST(terminalsCompleteGraph);

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

  //Test - begin
  /*
  for (unsigned int i = 0; i < numVertices; i++)
  {
    printf("vetor na posicao %u, %u:\n", i + 1, (*graph)->adjacencyList->at(i).vertex.isTerminal);

    for (unsigned int j = 0; j < (*graph)->adjacencyList->at(i).adjacencies.size(); j++)
      printf("%u, %u\n", (*graph)->adjacencyList->at(i).adjacencies.at(j).id,
             (*graph)->adjacencyList->at(i).adjacencies.at(j).weight);
  }
  */
  //Test - end

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
    fscanf(file, "%c %u %u %u %*[\r] %*[\n]", aux, &vertex, &vertex2, &weight);

    if (ferror(file))
      return 1;

    //vertice u -> vertice v
    auxAdjacencyInfo.id = vertex2;
    auxAdjacencyInfo.weight = weight;

    if (inputGraph->adjacencyList->at(vertex - 1).adjacencies.size() == 0)
    {
      inputGraph->adjacencyList->at(vertex - 1).vertex.id = vertex;
      inputGraph->adjacencyList->at(vertex - 1).vertex.isTerminal = false;
    }

    inputGraph->adjacencyList->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);

    //vertice v -> vertice u
    auxAdjacencyInfo.id = vertex;
    auxAdjacencyInfo.weight = weight;

    if (inputGraph->adjacencyList->at(vertex2 - 1).adjacencies.size() == 0)
    {
      inputGraph->adjacencyList->at(vertex2 - 1).vertex.id = vertex2;
      inputGraph->adjacencyList->at(vertex2 - 1).vertex.isTerminal = false;
    }

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
    fscanf(file, "%c %u %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      return true;

    inputGraph->adjacencyList->at(vertex - 1).vertex.isTerminal = true;
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
vector<vector<MinimalRouteInfo>>* shortestPath(Graph* inputGraph)
{
  unsigned int newDistance;
  DijkstraVertex auxVertex;
  AdjacencyInfo neighbourInfo;
  vector<DijkstraVertex*> fakeVertices(inputGraph->adjacencyList->size());
  vector<DijkstraVertex> openVertices(inputGraph->adjacencyList->size());
  vector<vector<MinimalRouteInfo>>* result =
    new vector<vector<MinimalRouteInfo>>(inputGraph->terminalList->size());

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
          //TODO:armazenar isto fora do no... quem sabe
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

    MinimalRouteInfo auxMinimalRouteInfo;    

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
        auxMinimalRouteInfo.weight = openVertices.at(nextId - 1).distanceFromSource;
        /*Certifica-se de que a rota 'i' nao se confunda com a rota 'j'*/
        auxMinimalRouteInfo.route.clear();

        //while (nextId != graph->terminalList->at(i + 1)->vertex.id)
        while (nextId != openVertices.at(nextId - 1).previousVertexId)
        {   
          /*Insere ID no vertice na rota*/
          auxMinimalRouteInfo.route.push_back(nextId);
          /*Atualiza*/
          nextId = openVertices.at(nextId - 1).previousVertexId;
        }

        /*Insere ID do vertice de origem*/
        auxMinimalRouteInfo.route.push_back(nextId);
        /*Insere informacao da rota na lista de terminais*/
        result->at(i).push_back(auxMinimalRouteInfo);
      }
    }
  }

  return result;
}

/*Funcao que cria um grafo completo, em que os vertices sao os terminais e o
peso das arestas que os une eh a distancia de menor custo entre cada um dos
terminais (valor calculado via Dijkstra)*/
/*
vector<vector<Edge>>* createCompleteGraph(vector<vector<MinimalRouteInfo>>* dijkstraResult)
{
  vector<vector<Edge>>* completeGraph =
    new vector<vector<Edge>>(dijkstraResult->size());
  Edge edge;
  unsigned int numVertices = dijkstraResult->size();
  //unsigned int numAdjacencies;
  unsigned int firstVertex;
  unsigned int lastVertex;
  unsigned int weight;

  for (unsigned int i = 0; i < numVertices; i++)
  {
    for (unsigned int j = 0; j < numVertices; j++)
    {
      firstVertex = dijkstraResult->at(i).at(j).route.at(dijkstraResult->at(i).
        at(j).route.size() - 1);
      lastVertex = dijkstraResult->at(i).at(j).route.at(0);
      weight = dijkstraResult->at(i).at(j).weight;
      edge.vertexId = firstVertex;
      edge.vertex2Id = lastVertex;
      edge.weight = weight;
      completeGraph->at(i).push_back(edge);
    }
  }

  return completeGraph;
}
*/
Graph* createCompleteGraph(vector<vector<MinimalRouteInfo>>* dijkstraResult)
{
  Graph* completeGraph = new Graph;
  completeGraph->adjacencyList = new vector<Adjacencies>(dijkstraResult->size());
  AdjacencyInfo adjacencyInfo;
  unsigned int numVertices = dijkstraResult->size();
  unsigned int numAdjacencies;
  unsigned int firstVertex;
  unsigned int lastVertex;

  for (unsigned int i = 0; i < numVertices; i++)
  {
    numAdjacencies = dijkstraResult->at(i).size();

    for (unsigned int j = 0; j < numAdjacencies; j++)
    {
      firstVertex = dijkstraResult->at(i).at(j).route.at(dijkstraResult->at(i).
        at(j).route.size() - 1);
      lastVertex = dijkstraResult->at(i).at(j).route.at(0);
      completeGraph->adjacencyList->at(i).vertex.id = firstVertex;
      adjacencyInfo.id = lastVertex;
      adjacencyInfo.weight = dijkstraResult->at(i).at(j).weight;
      completeGraph->adjacencyList->at(i).adjacencies.push_back(adjacencyInfo);
    }
  }

  return completeGraph;
}

/*Funcao auxiliar que usada na criacao/ajuste do heap binario
para o algoritmo de Kruskal*/
bool compareKruskal(Edge& edge, Edge& edge2)
{
  return edge.weight < edge2.weight;
}

/*Função que verifica se um grafo possui algum ciclo, atraves
do uso do algoritmo DFS*/
bool hasCycle(Graph* mst)
{

}

void generateMST(Graph* completeGraph)
{
  /*Usar essa ideia de 'size' na funcao acima...*/
  unsigned int size = completeGraph->adjacencyList->size();
  vector<Edge>* edgeList = new vector<Edge>;
  Graph mst;
  mst.adjacencyList = new vector<Adjacencies>(size);
  AdjacencyInfo adjacenyInfo;
  Edge edge;
  set<int> treeVertexSet;

  for (unsigned int i = 0; i < size; i++)
  {
    for (unsigned int j = i + 1, k = 0; j < size; j++, k++)
    {
      edge.vertexId = completeGraph->adjacencyList->at(i).vertex.id;
      edge.vertex2Id = completeGraph->adjacencyList->at(i).adjacencies.at(j - 1).id;
      edge.weight = completeGraph->adjacencyList->at(i).adjacencies.at(j - 1).weight;
      edgeList->push_back(edge);
    }
  }

  sort(edgeList->begin(), edgeList->end(), compareKruskal);
  treeVertexSet.insert(edgeList->at(0).vertexId);
  treeVertexSet.insert(edgeList->at(0).vertex2Id);
  adjacenyInfo.id = edgeList->at(0).vertex2Id;
  adjacenyInfo.weight = edgeList->at(0).weight;

  /*mst.adjacencyList->at(edgeList->at(0).vertexId - 1).adjacencies.at(0).id =
    ;
  mst.adjacencyList->at(edgeList->at(0).vertexId - 1).adjacencies.at(0).weight =
    edgeList->at(0).weight;*/
  mst.adjacencyList->at(edgeList->at(0).vertexId - 1).adjacencies.
    push_back(adjacenyInfo);

  unsigned int i = 1;

  while (treeVertexSet.size() < size)
  {
    if (!hasCycle(&mst))
    {
      adjacenyInfo.id = edgeList->at(i).vertex2Id;
      adjacenyInfo.weight = edgeList->at(i).weight;
      mst.adjacencyList->at(edgeList->at(i).vertexId - 1).adjacencies.
        push_back(adjacenyInfo);
      treeVertexSet.insert(edgeList->at(i).vertexId);
      treeVertexSet.insert(edgeList->at(i).vertex2Id);
    }

    i++;
  }
}
