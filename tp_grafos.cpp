#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstdio>
#include <climits>
#include <chrono>
#include <cstring>

/*Mensages*/
const char MSG_PARAMS_ERROR[] =
  "Numero de parametros esta incorreto (Arquivo de entrada deve ser informado. Exibicao do tempo eh opcional [-t])";
const char MSG_READFILE_ERROR[] =
  "Operacao de leitura do arquivo falhou. Codigo do erro:";
const char MSG_OPEN_FILE_ERROR[] =
  "Nao foi possivel abrir o arquivo especificado";
const char MSG_TIME_ARG_EXPECTED[] =
  "Caso exista, terceiro argumento deve ser -t";

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
/*GraphViz*/
void output_graphviz(SteinerTree* steinerTree, Graph* inputGraph);

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

  if (argc == 3)
  {
    if (strcmp(argv[2], "-t"))
    {
      printf("%s\n", MSG_TIME_ARG_EXPECTED);

      return -1;
    }
  }
  else if (argc < 2)
  {
    /*O arquivo de entrada deve ser enviado como argumento para
    o programa principal*/
    printf("%s\n", MSG_PARAMS_ERROR);

    return -2;
  }

  /*Tenta abrir arquivo de entrada*/
  file = fopen(argv[1], "r");

  if (!ferror(file))
  {
    /*Chama funcao que realiza leitura no arquivo de entrada*/
    returnedValue = readFile(file, &inputGraph);

    if (returnedValue)
      printf("%s %u\n", MSG_READFILE_ERROR, returnedValue);

    /*Trecho de medicao tempo (1): inicio*/
    high_resolution_clock::time_point t1 = high_resolution_clock::now();    
    /*Trecho de medicao tempo (1): fim*/

    /*Cacula o caminho de de menor custo para todos os vertices
    do tipo terminal*/
    dijkstraResult = dijkstra(inputGraph);
    /*Monta o fecho metrico (grafo completo, cujos vertices sao todos os
    vertices do tipo terminal*/
    terminalsCompleteGraph = createCompleteGraph(dijkstraResult);
    /*Chama funcao que obtem a arvore geradora minima, usando como grafo
    de entrada, o grafo completo, obtido acima*/
    mst = generateMST(terminalsCompleteGraph, &distinctVertices);
    /*Chama funcao que gera a arvore de Steiner*/
    steinerTree = generateSteinerTree(mst, inputGraph);

    /*Trecho medicao tempo (2): inicio*/
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);        
    /*Trecho medicao tempo (2): fim*/
    /*Imprime resultados*/
    printResults(steinerTree);
    /*Caso tenha recebido -t como argumento, imprime tempo de execucao*/
    if (!strcmp(argv[2], "-t"))
      printf("Tempo de execucao: %lf segundos\n", time_span.count());
    /*GraphViz*/
    /*output_graphviz(steinerTree, inputGraph);*/

    /*Fecha arquivo de entrada*/
    fclose(file);

    return returnedValue;
  }
  else
  {
    /*Nao foi possivel abrir o arquivo de entrada*/
    printf("%s\n", MSG_OPEN_FILE_ERROR);

    return -3;
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
    /*Le os dados e uma aresta*/
    fscanf(file, "%c %d %d %u %*[\r] %*[\n]", aux, &vertex, &vertex2, &weight);

    if (ferror(file))
      return 1;

    /*Como o grafo nao eh direcionado, arestas deverao ser incluidas nos
    dois sentidos*/

    /*vertice u -> vertice v*/
    auxAdjacencyInfo.id = vertex2;
    auxAdjacencyInfo.weight = weight;

    if (inputGraph->adjacencyList->at(vertex - 1).adjacencies.size() == 0)
      inputGraph->adjacencyList->at(vertex - 1).vertex.id = vertex;

    /*Adiciona adjacencia a variavel que representa o grafo de entrada*/
    inputGraph->adjacencyList->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);

    /*vertice v -> vertice u*/
    auxAdjacencyInfo.id = vertex;
    auxAdjacencyInfo.weight = weight;


    if (inputGraph->adjacencyList->at(vertex2 - 1).adjacencies.size() == 0)
      inputGraph->adjacencyList->at(vertex2 - 1).vertex.id = vertex2;

    /*Adiciona adjacencia a variavel que representa o grafo de entrada*/
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

  /*Aloca memoria suficiente para acomodar um ponteiro para cada terminal
  (Essa lista de ponteiros faz com que nao seja necessario percorrer toda
  a lista de vertices, testando qual deles eh terminal*/
  inputGraph->terminalList = new vector<Adjacencies*>(numTerminals);

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    /*Le linha que representa um terminal*/
    fscanf(file, "%c %d %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      /*Indica que houve erro de leitura*/
      return true;

    /*Aponta o ponteiro i para seu respectivo terminal*/
    inputGraph->terminalList->at(i) = &(inputGraph->adjacencyList->at(vertex - 1));
  }

  /*Indica que nao houve erro de leitura*/
  return false;
}

/*Funcao auxiliar que usada na criacao/ajuste do heap minimo binario
para o algoritmo de Dijkstra*/
bool compare(DijkstraVertex*& ptr, DijkstraVertex*& ptr2)
{
  /*Compara duas distancias ate o vertice de origem*/
  return ptr->distanceFromSource > ptr2->distanceFromSource;
}

/*Funcao que executa o algoritmo de Dijkstra
Obs.: A lista de ponteiros para vertices do tipo Dijkstra eh necessaria para
evitar que atualizacoes dos custos/distancias sejam feitas em vertices errados.
Como os vertices abertos sao acessados atraves de seu indice na lista de
vertices abertos, o processo de criacao/ajuste do heap faria com que essa
atualizacao ficasse impossivel.
*/
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
    openVertices.at(i) = auxVertex;
    /*Ponteiro para um vertice aberto eh inicializado*/
    fakeVertices.at(i) = &openVertices.at(i);
  } 

  /*Executa algoritmo de Dijkstra para cada um dos vertices terminais*/
  for (unsigned int i = 0; i < inputGraph->terminalList->size(); i++)
  {
    /*Armazena ID do terminal corrente*/
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

    /*Define que o terminal corrente seja usado como vertice de origem*/
    openVertices.at(currentTerminalId - 1).distanceFromSource = 0;
    /*Variavel que ajuda no controle do tamanho do conjunto de vertices abertos*/
    unsigned int j = 0;

    while (fakeVertices.size() - j > 0)
    {
      /*Organiza lista de forma que ela comportesse como um heap binario minimo*/
      make_heap(fakeVertices.begin(), fakeVertices.end() - j, compare);
      /*Sempre considera o primeiro vertice do heap*/
      auxVertex = *(fakeVertices.at(0));

      for (unsigned int k = 0;
           k < inputGraph->adjacencyList->at(auxVertex.id - 1).adjacencies.size();
           k++)
      {
        /*Obtem <vertice anterior> e <distancia ate o vertice fonte> do
        do vizinho corrente*/
        neighbourInfo = inputGraph->adjacencyList->at(auxVertex.id - 1).
          adjacencies.at(k);

        /*Evita calcular distancia duas vezes*/
        if (neighbourInfo.id != auxVertex.previousVertexId)
          newDistance = auxVertex.distanceFromSource + neighbourInfo.weight;

        /*Atualiza distancia do vertice vizinho, caso necessario*/
        if (newDistance < openVertices.at(neighbourInfo.id - 1).distanceFromSource)
        {          
          openVertices.at(neighbourInfo.id - 1).distanceFromSource = newDistance;
          openVertices.at(neighbourInfo.id - 1).previousVertexId = auxVertex.id;
        }
      }

      /*"Remove" primeiro objeto do heap
      Obs.: Remocao logica. Objeto na primeira posicao eh trocado com objeto
      da ultima posicao. A variavel j auxilia a identificar quantos objetos
      ainda restam no estrutura de dados*/
      pop_heap(fakeVertices.begin(), fakeVertices.end() - j);
      /*Essa variavel limita o numero de vertices abertos (pop_heap nao remove - apenas
      joga para o final da lista*/
      j++;
    }

    /*Variavel auxiliar que armazenara informacoes sobre o caminho de
    menor custo*/
    ShortestPathInfo auxShortestPathInfo;

    /*Laco principal de identificacao do caminho e seu custo minimo*/
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

        /*Montagem do caminho propriamente dito, que ocorre de tras pra
        frente*/
        while (nextId != openVertices.at(nextId - 1).previousVertexId)
        {   
          /*Insere ID no vertice na rota*/
          auxShortestPathInfo.path.push_back(nextId);
          /*Atualiza ID do proximo vertice a ser verificado*/
          nextId = openVertices.at(nextId - 1).previousVertexId;
        }

        /*Insere ID do vertice de origem*/
        auxShortestPathInfo.path.push_back(nextId);
        /*Insere informacao da rota na lista de terminais*/
        result->at(i).push_back(auxShortestPathInfo);
      }
    }
  }

  /*A variavel result contem as ditancias menos custosas entre
  todos os terminais do grafo de entrada*/
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

  /*Para cada terminal*/
  for (unsigned int i = 0; i < numTerminals; i++)
  {
    numAdjacencies = dijkstraResult->at(i).size();

    /*Para cada, eventualmente, falsa adjacencia do terminal corrente*/
    for (unsigned int j = 0; j < numAdjacencies; j++)
    {
      /*Insere falsa adjacencia na variavel que representa o grafo completo*/
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
void insertMST(Graph* mst, Edge& edge, const unsigned int& index,
               const unsigned int& index2)
{
  AdjacencyInfo adjacencyInfo;
  ShortestPathInfo pathInfo;

  /*Como o grafo nao eh direcionado, a insercao devera ocorrer nos dois sentidos*/

  /*u -> v*/
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

  /*v -> u*/
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
void removeMST(Graph* mst, const unsigned int& index, const unsigned int& index2)
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

  /*Insere vertice de origem no conjunto de vertices ja visitados*/
  visitedVertices->insert(vertexId);

  /*Para cada vizinho do vertice de origem*/
  while (i < mst->adjacencyList->at(sourceIndex).adjacencies.size())
  {
    /*Identifica vizinho*/
    neighbourId = mst->adjacencyList->at(sourceIndex).adjacencies.at(i).id;
    /*Identica o indice onde esse vizinho encontra-se na arvore geradora minima*/
    index = distinctVertices->find(neighbourId)->second;
    i++;

    /*Verifica se vizinho corrente ainda nao foi visitado*/
    if (visitedVertices->find(neighbourId) == visitedVertices->end())
    {
      /*Vizinho corrente ainda nao foi visitado*/

      /*Vizinho corrente eh inserido no conjunto de vertices visitados*/
      visitedVertices->insert(neighbourId);
      /*Chamada recursiva*/
      if (hasCycle(mst, visitedVertices, distinctVertices, index, vertexId))
        /*Ciclo detectado: interrompe busca*/
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

/*Funcao responsavel por construir uma arvore geradora minima, a partir
do grafo completo (fecho metrico), usando o algoritmo de Kruskal*/
Graph* generateMST(Graph* completeGraph, map<int, int>* distinctVertices)
{
  unsigned int size = completeGraph->adjacencyList->size();
  vector<Edge>* edgeList = new vector<Edge>;
  Graph* mst = new Graph;
  mst->adjacencyList = new vector<Adjacencies>(size);
  Edge edge;

  /*Variavel utilizada pela funcao hasCycle(). Eh declarada aqui pq hasCycle()
  eh recursiva. Caso contrario nao funcionaria... variavel seria redeclarada a
  cada chamada*/
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
  /*Como durante a construcao da arvore gerado minima nao eh possivel prever quais
  vertices serao inseridos primeiro, ja que isso depende das arestas de menor peso,
  e como as listas sao acessadas via indice, eh necessario armazenar o indica de
  insercao de cada vertice em algum lugar. A variavel insertionIndex, em conjunto
  com a estrutura de dados distinctVertices realizam esse controle*/
  unsigned int insertionIndex = 0;

  /*Laco eh realizado ate que o numero de arestas inseridas seja igual a
  n - 1*/
  while (insertedEdges < size - 1)
  {
    /*map.insert() retorna um iterator para um par, em que o primeiro
    elemento eh objeto inserido na estrutura de dados e o segundo eh um
    valor booleano, que indica se a insercao ocorreu com sucesso*/
    resultPair = distinctVertices->insert(
      pair<int, int>(edgeList->at(i).vertexId, insertionIndex));
    vertexIndex = resultPair.first->second;

    if (resultPair.second)
      /*Caso um novo vertice tenha sido inserido, atualiza variavel
      que indica o indice onde a proxima insercao ocorrera*/
      insertionIndex++;

    resultPair2 = distinctVertices->insert(
      pair<int, int>(edgeList->at(i).vertex2Id, insertionIndex));
    vertex2Index = resultPair2.first->second;

    if (resultPair2.second)
      insertionIndex++;

    /*Atualiza arvore geradora minima com nova aresta e vertice(s)*/
    insertMST(mst, edgeList->at(i), vertexIndex, vertex2Index);

    /*Conjunto de vertices visitados pela funcao que executa o DFS deve
    ser reiniciado a cada chamada.*/
    visitedVertices.clear();

    /*Verifica se a aresta recem inserida adicionou um ciclo ao grafo*/
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
    }
    else
    {
      /*Contador de arestas inseridas com sucesso*/
      insertedEdges++;
    }

    /*Contador de quantidade de arestas existentes no grafo completo*/
    i++;
  }

  return mst;
}

/*Funcao que insere os vertices de uma aresta na arvore de Steiner*/
void insertSteiner(SteinerTree* steinerTree, const int& vertexId,
                   const unsigned int& index, const int& vertex2Id,
                   const unsigned int& index2)
{
  AdjacencyInfo adjacencyInfo;

  /*Como o grafo nao eh direcionado, a insercao devera ocorrer nos dois sentidos*/

  /*u -> v*/
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = vertex2Id;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  steinerTree->adjacencyList->at(index).vertex.id = vertexId;
  steinerTree->adjacencyList->at(index).adjacencies.push_back(adjacencyInfo);

  /*v -> u*/
  //Atualiza dados da adjacencia a ser inserida na MST
  adjacencyInfo.id = vertexId;

  /*Vertice vertex  eh inserido, juntamente com sua adjacencia, na
  arvore geradora minima*/
  steinerTree->adjacencyList->at(index2).vertex.id = vertex2Id;
  steinerTree->adjacencyList->at(index2).adjacencies.push_back(adjacencyInfo);
}

/*Funcao que gera a arvore de Steiner, a partir da arvore geradora minima e do
grafo de entrada (obter pesos).
Obs.: A ideia dessa funcao eh substituir cada aresta da arvore geradora minima
pelo caminho mais curto, determinado via Dijkstra. Porem, como isso nao garante
que o grafo resultante nao contera ciclos, eh preciso fazer essa verificacao.*/
SteinerTree* generateSteinerTree(Graph* mst, Graph* inputGraph)
{
  SteinerTree* steinerTree;
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

  /*Como nao eh possivel determinar a quantidade de exata de vertices que irao
  compor a arvore de Steiner, aloca memoria para n vertices*/
  steinerTree = new SteinerTree;
  steinerTree->adjacencyList = new vector<Adjacencies>(inputGraph->adjacencyList->size());
  steinerTree->totalWeight = 0;

  /*Para cada vertice da arvore geradora minima*/
  for (unsigned int i = 0; i < mst->adjacencyList->size(); i++)
  {
    numPaths = mst->adjacencyList->at(i).vertex.shortestPathInfo.size();

    /*Para cada adjacencia*/
    for (unsigned int j = 0; j < numPaths; j++)
    {
      /*Identifica comprimento do caminho mais curto*/
      pathSize = mst->adjacencyList->at(i).vertex.shortestPathInfo.at(j).path.size();
      vector<int>& path = mst->adjacencyList->at(i).vertex.shortestPathInfo.at(j).path;

      /*Para cada par de vertices do caminho mais curto, tenta inserir tal
      aresta na arvore de Steiner*/
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

        /*Chama funcao que realiza a insercao dos vertices na arvore de Steiner*/
        insertSteiner(steinerTree, path.at(k), vertexIndex, path.at(k + 1), vertex2Index);

        /*Conjunto de vertices visitados pela funcao que executa o DFS deve
        ser reiniciado a cada chamada.*/
        visitedVertices.clear();

        /*Verifica se a aresta recem inserida adicionou um ciclo ao grafo*/
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

  /*Identifica o peso de cada aresta que compoe a arvore de Steiner*/
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

/*Funcao que imprime os resultados*/
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

  printf("Vertices (%lu): ", distinctVertices.size());

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
