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
} MinimalRoute;

typedef struct tDijkstraVertex : Vertex
{
  unsigned int distanceFromSource;
  int previousVertexId;
  vector<MinimalRoute> minimalRouteInfo;
} DijkstraVertex;

typedef struct tAdjacencies
{
  Vertex vertex;
  vector<AdjacencyInfo> adjacencies;
} Adjacencies;

/*Typedefs*/
typedef vector<Adjacencies> Graph;

/*Declaracao de funcoes*/
unsigned int readFile(FILE*, Graph**);
bool readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);
bool compare(DijkstraVertex&, DijkstraVertex&);
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
bool compare(DijkstraVertex& vertex, DijkstraVertex& vertex2)
{
  return vertex.distanceFromSource > vertex2.distanceFromSource;
}

/*Funcao que executa o algoritmo de Dijkstra*/
void shortestPath(Graph* graph)
{
  unsigned int newDistance;
  DijkstraVertex auxVertex;
  AdjacencyInfo neighbourInfo;
  vector<DijkstraVertex> openVertices;
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
    openVertices.push_back(auxVertex);
    /*Caso o vertice i seja um terminal, ele eh adicionado a lista
    de terminais tambem*/
    if (graph->at(i).vertex.isTerminal)
    {
      auxVertex.isTerminal = true;

      terminals.push_back(auxVertex);
    }
  }

  /*Executa algoritmo de Dijkstra para cada um dos vertices terminais*/
  for (unsigned int i = 0; i < terminals.size(); i++)
  {    
    openVertices.at(i).distanceFromSource = 0;

    unsigned int j = 0;

    while (openVertices.size() - j > 0)
    {
      make_heap(openVertices.begin(), openVertices.end() - j, compare);
      auxVertex = openVertices.at(0);
      pop_heap(openVertices.begin(), openVertices.end() - j, compare);

      for (unsigned int k = 0; k < graph->at(auxVertex.id - 1).adjacencies.size(); k++)
      {
        neighbourInfo = graph->at(auxVertex.id - 1).adjacencies.at(k);
        newDistance = auxVertex.distanceFromSource + neighbourInfo.weight;

        if (newDistance < openVertices.at(neighbourInfo.id - 1).distanceFromSource)
        {
          terminals.at(neighbourInfo.id - 1).distanceFromSource = newDistance;
          terminals.at(neighbourInfo.id - 1).previousVertexId = auxVertex.id;
        }
      }

      /*Essa variavel limita o numero de verices abertos (pop_heap nao remove - so
      joga para o final da lista*/
      j++;
    }
  }
}

/*TODO: distancias e rotulos dos vertices anteriores devem ser gravados em outra estrutura
de dados - semelhante a openVertices.*/
