#include <iostream>
#include <vector>
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
bool readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);
bool compare(DijkstraVertex*&, DijkstraVertex*&);
vector<vector<MinimalRouteInfo>>* shortestPath(Graph*);
Graph* createCompleteGraph(vector<vector<MinimalRouteInfo>>*);

/*Funcao principal*/
int main(int argc, char* argv[])
{
  unsigned int returnedValue;
  Graph* graph = NULL;
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
    returnedValue = readFile(file, &graph);

    if (returnedValue)
      printf("%s %u\n", MSG_READFILE_ERROR, returnedValue);

    dijkstraResult = shortestPath(graph);

    //test
    createCompleteGraph(dijkstraResult);

    fclose(file);
    delete graph;

    return returnedValue;
  }
  else
  {
    printf("%s\n", MSG_OPEN_FILE_ERROR);

    return -2;
  }
}

/*Funcao principal de leitura do arquivo*/
unsigned int readFile(FILE* file, Graph** graph)
{
  char aux[100];
  unsigned int numVertices, numEdges, numTerminals;
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
  *graph = new Graph;
  (*graph)->adjacencyList = new vector<Adjacencies>(numVertices);
  (*graph)->terminalList = nullptr;

  //Chama funcao que realiza a leitura das arestas
  error = readEdges(file, *graph, numEdges);

  if (error)
    return 4;

  //Descarta "End"
  fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);

  if (ferror(file))
    return 5;

  //Le linha de inicio da secao "Terminais"
  fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);

  if (ferror(file))
    return 6;

  //Le numero de terminais
  fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numTerminals);

  if (ferror(file))
    return 7;

  //Chama funcao que realiza a leitura dos terminais
  error = readTerminals(file, *graph, numTerminals);

  if (error)
    return 8;

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
bool readEdges(FILE* file, Graph* graph, const unsigned int& numEdges)
{
  char aux[100];
  int vertex, vertex2;
  unsigned int weight;

  AdjacencyInfo auxAdjacencyInfo;

  for (unsigned int i = 0; i < numEdges; i++)
  {
    fscanf(file, "%c %u %u %u %*[\r] %*[\n]", aux, &vertex, &vertex2, &weight);

    if (ferror(file))
      return true;

    //vertice u -> vertice v
    auxAdjacencyInfo.id = vertex2;
    auxAdjacencyInfo.weight = weight;

    if (graph->adjacencyList->at(vertex - 1).adjacencies.size() == 0)
    {
      graph->adjacencyList->at(vertex - 1).vertex.id = vertex;
      graph->adjacencyList->at(vertex - 1).vertex.isTerminal = false;
      graph->adjacencyList->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
    else
    {
      graph->adjacencyList->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);
    }

    //vertice v -> vertice u
    auxAdjacencyInfo.id = vertex;
    auxAdjacencyInfo.weight = weight;

    if (graph->adjacencyList->at(vertex2 - 1).adjacencies.size() == 0)
    {
      graph->adjacencyList->at(vertex2 - 1).vertex.id = vertex2;
      graph->adjacencyList->at(vertex2 - 1).vertex.isTerminal = false;
      graph->adjacencyList->at(vertex2 - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
    else
    {
      graph->adjacencyList->at(vertex2 - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
  }

  return false;
}

/*Funcao responsavel pela leitura dos terminais*/
bool readTerminals(FILE* file, Graph* graph, const unsigned int& numTerminals)
{
  char aux[100];
  int vertex;

  graph->terminalList = new vector<Adjacencies*>(numTerminals);

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    fscanf(file, "%c %u %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      return true;

    graph->adjacencyList->at(vertex - 1).vertex.isTerminal = true;
    graph->terminalList->at(i) = &(graph->adjacencyList->at(vertex - 1));
  }

  return false;
}

/*Funcao auxiliar (usada na criacao/ajuste do heap binario)*/
bool compare(DijkstraVertex*& ptr, DijkstraVertex*& ptr2)
{
  return ptr->distanceFromSource > ptr2->distanceFromSource;
}

/*Funcao que executa o algoritmo de Dijkstra*/
vector<vector<MinimalRouteInfo>>* shortestPath(Graph* graph)
{
  unsigned int newDistance;
  DijkstraVertex auxVertex;
  AdjacencyInfo neighbourInfo;
  vector<DijkstraVertex*> fakeVertices(graph->adjacencyList->size());
  vector<DijkstraVertex> openVertices(graph->adjacencyList->size());
  vector<vector<MinimalRouteInfo>>* result =
    new vector<vector<MinimalRouteInfo>>(graph->terminalList->size());

  /*Inicializacao do conjunto que representa os vertices abertos*/
  for (unsigned int i = 0; i < graph->adjacencyList->size(); i++)
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
  for (unsigned int i = 0; i < graph->terminalList->size(); i++)
  {
    /*Armazena ID do vertice (terminal) atual*/
    int currentTerminalId = graph->terminalList->at(i)->vertex.id;

    for (unsigned int i = 0; i < graph->adjacencyList->size(); i++)
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

      for (unsigned int k = 0; k < graph->adjacencyList->at(auxVertex.id - 1).adjacencies.size(); k++)
      {
        neighbourInfo = graph->adjacencyList->at(auxVertex.id - 1).adjacencies.at(k);

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
    for (unsigned int j = 0; j < graph->terminalList->size(); j++)
    {
      /*Variavel que indica o ID do terminal para o qual deseja-se
      calcular a distancia menos custosa*/
      int nextId = graph->terminalList->at(j)->vertex.id;
      /*Nao ha necessidade de calcular rota para si proprio*/
      if (currentTerminalId != nextId)
      {
        /*Obtem custo da distancia entre o terminal corrente e sua adjacencia corrente*/
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
}
