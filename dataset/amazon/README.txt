This is s a dataset extracted from the Amazon.
The folder contains the following files:
1. vertex.txt
 Each line represents a vertex in Amazon. Each line starts with the vertex_id, following by vertex type.

2.edge.txt 
 Each line represents an edge in Amazon. Each line starts with the edge_id, following by edge type.

3. graph.txt 
 Each line represents an adjacent array. Each line starts with the vertex_id, following by a list of neighbor_vertex_id and edge_id.

The vertex types and edge types are numbered as follow:


<Item> : 0;
<User> : 1;
<View> : 2;
<Brand> : 3;

<Item->User> 0;
<Item->View> 1;
<Item->Brand> 2;
<User->Item> : 3;
<View->Item> : 4;
<Brand->Item> : 5;

Type	  number
Item	  2,753
User	  6,170
View	  3,857
Brand	  334

Edge    90340