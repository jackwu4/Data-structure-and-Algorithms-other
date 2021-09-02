/*
* A main function for you to build and run your
* own tests with.
* This file is not part of the marking, so you
* can do anything you want here.
*/
#include <iostream>

#include "directed_graph_algorithms.cpp"

int main() {
	directed_graph<char> g;
	g.add_vertex('a');
	g.add_vertex('b');
	g.add_vertex('c');
	g.add_vertex('d');
	g.add_vertex('e');
	g.add_edge('a','b');
	g.add_edge('b','c');
	g.add_edge('c','d');
	g.add_edge('d','e');
	g.add_edge('e','b');
	std::cout << g.num_edges() << std::endl;
	std::cout << g.num_vertices() << std::endl;
	for (auto vertex_iterator = g.begin(); vertex_iterator != g.end(); vertex_iterator++) {
		//std::cout << *vertex_iterator << " ";
	}
	directed_graph<char> u = undirected_graph(g);
	std::cout << u.num_edges() << std::endl;
	std::cout << u.num_vertices() << std::endl;
	std::vector<char> visited;
	depth_first_componentsHelper(u,*g.begin(),visited);
	//std::cout << is_dag(g) << std::endl;
	//directed_graph<std::string> s;
	//std::cout << is_dag(s) << std::endl;
	//std::cout << "ham "<<is_hamiltonian_dag(g) << std::endl;
	//std::list<char> sort = topological_sort(g);
  	//for(auto const &v:sort){
	//	std::cout << v << std::endl;
	//}
}
