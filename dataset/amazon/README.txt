This is s a dataset extracted from the Amazon.
The folder contains the following files:
1. vertex.txt
 Each line represents a vertex in Amazon. Each line starts with the vertex_id, following by vertex type.

2. edge.txt
 Each line represents an edge in Amazon. Each line starts with the start_id, following by target_id.

2.edge.txt ***
 Each line represents an edge in Amazon. Each line starts with the edge_id, following by edge type.


3. graph.txt
 Each line represents an adjacent array. Each line starts with the vertex_id, following by a list of neighbor_vertex_id and -1.

3. graph.txt ***
 Each line represents an adjacent array. Each line starts with the vertex_id, following by a list of neighbor_vertex_id and edge_id.


The vertex types and edge types are numbered as follow:


<Paper> : 0;
<Author> : 1;
<Venue> : 2;
<Topic> : 3;

<Paper->Author> 0;
<Paper->Venue> 1;
<Paper->Topic> 2;
<Author->Paper> : 3;
<Venue->Paper> : 4;
<Topic->Paper> : 5;

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