#include <iostream>
#include <vector>
#include <cstdio>

using namespace std;

/*Structs*/
typedef struct tAdjacencyInfo
{
  unsigned int id;
  unsigned int weight;
} AdjacencyInfo;

typedef struct tVertex
{
  unsigned int id = -1;
  bool isTerminal;
} Vertex;

typedef struct tAdjacencies
{
  Vertex vertex;
  vector<AdjacencyInfo> adjacencies;
} Adjacencies;


/*Typedefs*/
typedef vector<Adjacencies> Graph;

/*Declaracao de funcoes*/
unsigned int readFile(FILE*);
bool readEdges(FILE*, Graph*, const unsigned int&);
bool readTerminals(FILE*, Graph*, const unsigned int&);

/*Funcao principal*/
int main(int argc, char* argv[])
{
  unsigned int error;
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
    error = readFile(file);

    if (error)
      printf("Erro (readFile): %u\n", error);

    fclose(file);

    return error;
  }
  else
  {
    printf("Nao foi possivel ler o arquivo\n");

    return -2;
  }

  return 0;
}

/*Funcao principal de leitura do arquivo*/
unsigned int readFile(FILE* file)
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
  Graph graph(numVertices);

  //Chama funcao que realiza a leitura das arestas
  error = readEdges(file, &graph, numEdges);

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
  error = readTerminals(file, &graph, numTerminals);

  if (error)
    return 8;

  //Test - begin
  for (unsigned int i = 0; i < numVertices; i++)
  {
    printf("vetor na posicao %u, %u:\n", i + 1, graph.at(i).vertex.isTerminal);

    for (unsigned int j = 0; j < graph.at(i).adjacencies.size(); j++)
    {
      printf("%u, %u\n", graph.at(i).adjacencies.at(j).id, graph.at(i).adjacencies.at(j).weight);
    }
  }
  //Test - end

  //Sem erro na manipulacao do arquivo
  return 0;
}

/*Funcao responsavel pela leitura das arestas*/
bool readEdges(FILE* file, Graph* graph, const unsigned int& numEdges)
{
  char aux[100];
  unsigned int vertex, vertex2, weight;

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
  unsigned int vertex;

  for (unsigned int i = 0; i < numTerminals; i++)
  {
    fscanf(file, "%c %u %*[\r] %*[\n]", aux, &vertex);

    if (ferror(file))
      return true;

    graph->at(vertex - 1).vertex.isTerminal = true;
  }

  return false;
}
