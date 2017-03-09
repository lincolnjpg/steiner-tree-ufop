#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

typedef struct tVertex
{
  bool isTerminal;
  unsigned int id;
  unsigned int weight;
} Vertex;

typedef vector<vector<Vertex>> VertexList;

typedef struct tGraph
{
  VertexList list;
} Graph;

void readFile()
{
  char aux[100];
  unsigned int numVertex, numEdges, numTerminals, vertex, vertex2, weight;
  FILE* file = NULL;

  file = fopen("/home/junin/workspace/qtcreator/graph_theory_tp/1.stp", "r");

  if (file != NULL)
  {    
    //Descarta primeira linha
    fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);
    printf("%s\n", aux);
    //Le numero de vertices
    fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numVertex);
    printf("%s, %u\n", aux, numVertex);
    //Le numero de arestas
    fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numEdges);
    printf("%s, %u\n", aux, numEdges);

    VertexList graph(numVertex);

    for (unsigned int i = 0; i < numEdges; i++)
    {
      fscanf(file, "%c %u %u %u %*[\r] %*[\n]", aux, &vertex, &vertex2, &weight);

      Vertex v;
      v.id = vertex2 - 1;
      v.weight = weight;

      graph.at(vertex - 1).push_back(v);

      v.id = vertex - 1;

      graph.at(vertex2 - 1).push_back(v);
    }

    //Test - begin
    for (unsigned int i = 0; i < numVertex; i++)
    {
      printf("vetor na posicao %u:\n", i + 1);

      for (unsigned int j = 0; j < graph.at(i).size(); j++)
      {
        printf("%u, %u\n", graph.at(i).at(j).id + 1, graph.at(i).at(j).weight);
      }
    }
    //Test - end

    //Descarta "end"
    fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);
    //Descarta "Section Terminals"
    fscanf(file, "%[A-Z a-z] %*[\r] %*[\n]", aux);
    //Le numero de terminais
    fscanf(file, "%[A-Z a-z] %u %*[\r] %*[\n]", aux, &numTerminals);

    for (unsigned int i = 0; i < numTerminals; i++)
    {
      fscanf(file, "%c %u %*[\r] %*[\n]", aux, &vertex);

      /*Como atualizar os vertices com a informacao "terminal",
      que ha vertices duplicados na lista de adjacencia? De repente,
      pode ser mesmo melhor deixar para incluir os vertices na lista
      apos definir quais deles sao terminais. Alem disso, estava pensando
      numa outra forma de representar o grafo no nosso programa (que tb
      usa lista de adjacencia*/
    }

    fclose(file);
  }  
}

int main()
{
  readFile();

  return 0;
}
