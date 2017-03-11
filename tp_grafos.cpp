#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <climits>

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

typedef struct tPointerVertex2
{
  int id;
  char c;
} PointerVertex2;

typedef struct tAdjacencies
{
  Vertex vertex;
  vector<AdjacencyInfo> adjacencies;
} Adjacencies;

typedef struct tGraphX //usar isso ao inves do Graph anterior
{
  Adjacencies adjacencyList;
  vector<Vertex> terminalList;
} GraphX;

/*Typedefs*/
typedef vector<Adjacencies> Graph;

/*Declaracao de funcoes*/
unsigned int readFile(FILE*, Graph**);
bool readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);
bool compare(DijkstraVertex*&, DijkstraVertex*&);
void shortestPath(Graph*);

/*Funcao principal*/
int main(int argc, char* argv[])
{
  unsigned int returnedValue;
  Graph* graph = NULL;
  FILE* file = NULL;

  if (argc == 2)
  {
    file = fopen(argv[1], "r");
  }
  else
  {
    printf("Numero de parametros esta incorreto (informe apenas 1 arquivo de entrada)\n");

    return -1;
  }

  if (!ferror(file))
  {
    returnedValue = readFile(file, &graph);

    if (returnedValue)
      printf("Erro (readFile): %u\n", returnedValue);

    shortestPath(graph);

    fclose(file);
    delete graph;

    return returnedValue;
  }
  else
  {
    printf("Nao foi possivel ler o arquivo\n");

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
  *graph = new Graph(numVertices);

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
  for (unsigned int i = 0; i < numVertices; i++)
  {
    printf("vetor na posicao %u, %u:\n", i + 1, (*graph)->at(i).vertex.isTerminal);

    for (unsigned int j = 0; j < (*graph)->at(i).adjacencies.size(); j++)
      printf("%u, %u\n", (*graph)->at(i).adjacencies.at(j).id, (*graph)->at(i).adjacencies.at(j).weight);
  }
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

    if (graph->at(vertex - 1).adjacencies.size() == 0)
    {
      graph->at(vertex - 1).vertex.id = vertex;
      graph->at(vertex - 1).vertex.isTerminal = false;
      graph->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
    else
    {
      graph->at(vertex - 1).adjacencies.push_back(auxAdjacencyInfo);
    }

    //vertice v -> vertice u
    auxAdjacencyInfo.id = vertex;
    auxAdjacencyInfo.weight = weight;

    if (graph->at(vertex2 - 1).adjacencies.size() == 0)
    {
      graph->at(vertex2 - 1).vertex.id = vertex2;
      graph->at(vertex2 - 1).vertex.isTerminal = false;
      graph->at(vertex2 - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
    else
    {
      graph->at(vertex2 - 1).adjacencies.push_back(auxAdjacencyInfo);
    }
  }

  return false;
}

/*Funcao responsavel pela leitura dos terminais*/
bool readTerminals(FILE* file, Graph* graph, const unsigned int& numTerminals)
{
  char aux[100];
  int vertex;

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    fscanf(file, "%c %u %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      return true;

    graph->at(vertex - 1).vertex.isTerminal = true;
  }

  return false;
}

/*Funcao auxiliar (usada na criacao/ajuste do heap binario)*/
bool compare(DijkstraVertex*& ptr, DijkstraVertex*& ptr2)
{
  return ptr->distanceFromSource > ptr2->distanceFromSource;
}

/*Funcao que executa o algoritmo de Dijkstra*/
void shortestPath(Graph* graph)
{
  unsigned int newDistance;
  DijkstraVertex auxVertex;
  AdjacencyInfo neighbourInfo;
  vector<DijkstraVertex*> openVertices(graph->size());
  vector<PointerVertex2> test(100);
  //vector<DijkstraVertex> openVertices;
  vector<DijkstraVertex> resultVertices(graph->size()); //rever nome...
  vector<DijkstraVertex> terminals;

  /*Inicializacao do conjunto que representa os vertices abertos*/
  for (unsigned int i = 0; i < graph->size(); i++)
  {
    /*Inicializa ID do vertice*/
    auxVertex.id = i + 1;
    /*Indica que a distancia entre o vertice corrente e o vertice
    de origem ainda nao foi definida*/
    auxVertex.distanceFromSource = UINT_MAX;
    /*Indica que o rotulo do vertice anterior ao vertice corrente
    ainda nao foi definido (incializado com o valor do proprio rotulo)*/
    auxVertex.previousVertexId = i + 1;
    /*Insere vertice corrente na lista de vertices abertos*/
    //resultVertices.push_back(auxVertex); //test
    resultVertices.at(i) = auxVertex;
    /*Caso o vertice i seja um terminal, ele eh adicionado a lista
    de terminais tambem*/
    /*Vertice aberto i passa a apontar para o vertice Dijkstra i*/
    openVertices.at(i) = &resultVertices.at(i);

    if (graph->at(i).vertex.isTerminal)
    {
      auxVertex.isTerminal = true;

      terminals.push_back(auxVertex);
    }
  } 

  //test - inicio
  //for (unsigned int i = 0; i < resultVertices.size(); i++)
    //resultVertices.at(i).distanceFromSource = 0;
  //test - end

  /*Executa algoritmo de Dijkstra para cada um dos vertices terminais*/
  for (unsigned int i = 0; i < terminals.size(); i++)
  {    
    //TODO: fazer loop no terminalList.size() e usar apenas os vertices de terminalList
    openVertices.at(i)->distanceFromSource = 0; //teste - inicializar o primeiro terminal ao inves desse

    unsigned int j = 0;

    while (openVertices.size() - j > 0)
    {
      /*Organiza lista de forma que ela comportesse como um heap binario minimo*/
      //teste - inicio
      //acho que nao precisarei de make_heap aqui, pois, inicialmente, apenas a primeira posicao eh inicializada com 0
      make_heap(openVertices.begin(), openVertices.end() - j, compare);
      //teste - fim
      /*Sempre considera o primeiro vertice do heap*/
      auxVertex = *(openVertices.at(0)); //teste - pegar o primeiro terminal, ao inves desse

      for (unsigned int k = 0; k < graph->at(auxVertex.id - 1).adjacencies.size(); k++)
      {
        neighbourInfo = graph->at(auxVertex.id - 1).adjacencies.at(k);

        /*Evita calcular distancia duas vezes*/
        if (neighbourInfo.id != auxVertex.previousVertexId)
          newDistance = auxVertex.distanceFromSource + neighbourInfo.weight;

        if (newDistance < resultVertices.at(neighbourInfo.id - 1).distanceFromSource)
        {
          //TODO:armazenar isto fora do no... quem sabe
          resultVertices.at(neighbourInfo.id - 1).distanceFromSource = newDistance;
          resultVertices.at(neighbourInfo.id - 1).previousVertexId = auxVertex.id;
        }
      }
      //teste - inicio
      //if (j == 0)
        //make_heap(openVertices.begin(), openVertices.end() - j, compare);
      //teste - fim

      /*Primeiro elemento do heap troca de lugar com o ultimo e heap eh reorganizado
      (Obs.: pop_heap so deve ser executado apos a lista openVertices ter sido atualizada.
      Caso contrario, os valores seriam atribuidos a posicoes diferentes)*/
      pop_heap(openVertices.begin(), openVertices.end() - j);
      /*Essa variavel limita o numero de vertices abertos (pop_heap nao remove - apenas
      joga para o final da lista*/
      j++;
    }

    /*Identifica rota e seu custo minimo*/
    MinimalRouteInfo auxMinimalRouteInfo;    

    //for (unsigned int j = 0; j < graph->at(terminals.at(i).id - 1).adjacencies.size(); j++)
    for (unsigned int j = 0; j < terminals.size(); j++)
    {
      /*Nao ha necessidade de calcular rota para si proprio*/
      if (i != j)
      {
        /*Identifica rotulo do vertice que faz ajdacencia com o terminal corrente*/
        unsigned int nextId = j + 1;
        /*Obtem custo da distancia entre o terminal corrente e sua adjacencia corrente*/
        auxMinimalRouteInfo.weight = resultVertices.at(nextId - 1).distanceFromSource;
        /**/
        auxMinimalRouteInfo.route.clear();

        while (nextId != i + 1)
        {
          //TODO: arrumar comentarios          
          /*Obtem rotulo do vertice anterior no caminho terminal -> adjacencia correntes*/
          auxMinimalRouteInfo.route.push_back(nextId);
          /*Atualiza rotulo*/
          nextId = resultVertices.at(nextId - 1).previousVertexId;
        }

        /*Insere rotulo do vertice de origem*/
        auxMinimalRouteInfo.route.push_back(i + 1);
        /*Insere informacao da rota na lista de terminais*/
        terminals.at(i).minimalRouteInfo.push_back(auxMinimalRouteInfo);
      }
    }
  }
}

/*TODO: distancias e rotulos dos vertices anteriores devem ser gravados em outra estrutura
de dados - semelhante a openVertices.*/
