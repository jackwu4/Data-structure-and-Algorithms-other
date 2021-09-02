/*
 * Notice that the list of included headers has
 * expanded a little. As before, you are not allowed
 * to add to this.
 */
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>
#include <cstddef>
#include <string>
#include <utility>
#include <algorithm>
#include <limits>
#include <optional>
#include <exception>
#include <stdexcept>

#include "directed_graph.hpp"

/*
 * Computes whether the input is a Directed Acyclic Graph (DAG).
 * A digraph is a DAG if there is no vertex that has a cycle.
 * A cycle is a non-empty set of [out-]edges that starts at one 
 * vertex, and returns to it.
 */
template <typename vertex>
bool is_dag(const directed_graph<vertex> & d) {
  	
	return topological_sort(d).size() == d.num_vertices();  //if topologicalsort's size is same as amount of vertices then its a dag
}                                                         //any graph that has a loop, the loop will not be in the topological sort

/*
 * Computes a topological ordering of the vertices.
 * For every vertex u in the order, and any of its
 * neighbours v, v appears later in the order than u.
 */
template <typename vertex>
std::list<vertex> topological_sort(const directed_graph<vertex> & d) {
  std::list<vertex> sorted;
  directed_graph<vertex> g = d;
  
  for(;;){                                                                       //a infinite loop until reach break
  	auto vertex_iterator = g.begin();                                            //initialise iterator with the first vertex
  	for (vertex_iterator = g.begin(); vertex_iterator != g.end(); vertex_iterator++){
		  if(g.in_degree(*vertex_iterator) == 0){                                    //loop through vertices and when found a leaf
			  break;                                                                   //break the current vertices loop
			}
	  }
	  if(vertex_iterator != g.end()){                                              //check if iterator reached the end
		  sorted.push_back(*vertex_iterator);                                        //if not, add vertex to the list 
		  g.remove_vertex(*vertex_iterator);                                         //remove vertex from the graph (including it's edges)
	   }else{
		  break;                                                                     //if reached the end, means all non dag value are found and added to the list
	   }                                                                           //break the infinite loop
  }
  return sorted;
}

/*
 * Given a DAG, computes whether there is a Hamiltonian path.
 * a Hamiltonian path is a path that visits every vertex
 * exactly once.
 */
template <typename vertex>
bool is_hamiltonian_dag(const directed_graph<vertex> & d) {
	bool result = false;                                                          //final result
  if(d.num_vertices()==0 || d.num_vertices() == 1)                              //if number of vertices is 1 or 0
    return true;                                                                //its a hamiltonian dag
  std::unordered_map<vertex, bool> visited;                                     //create a map which being passed around to updates whether if a vertex has been visited or not
  for(auto iter = d.begin(); iter!=d.end(); iter++){                            
    visited[*iter] = false;                                                     //add all vertex in the graph to the map and set it to false
    result = result || ham_help(d, visited, 0, *iter);                          //start the path from every vertex by running hamhelp method
  }
    
  return result;
}

template <typename vertex>
bool ham_help(const directed_graph<vertex> & d, std::unordered_map<vertex,bool> &visited, int currLength, const vertex &u) {
  bool result = false;
  visited[u] = true;
  if((currLength + 1) == d.num_vertices())                                      //check if the current path length is the same as all vertices
    return true;                                                                //if true it means all vertices has been reached with this path
  
  for(auto u_neighbour = d.nbegin(u); u_neighbour!=d.nend(u); u_neighbour++){   //iteration through all neighbours of u (breath first travels)
    if(!visited[*u_neighbour]){                                                 //if u not visited 
    result = result || ham_help(d, visited, currLength+1, *u_neighbour);        //go to them and from them find all its neighbour and go to them
    }                                                                           //everytime when a neighbour is found the path length increase by 1
  }                                                                             //when a path is found result will be true
  visited[u] = false;                                                           
  
  return result;
}

/*
 * Computes the weakly connected components of the graph.
 * A [weak] component is the smallest subset of the vertices
 * such that the in and out neighbourhood of each vertex in
 * the set is also contained in the set.
 */
template <typename vertex>
std::vector<std::vector<vertex>> components(const directed_graph<vertex> & d) {
  std::vector<std::vector<vertex>> componentsResult;
  directed_graph<vertex> undirected = undirected_graph(d);                            //turn graph in to undirected
  std::vector<vertex> unvisited;                                                      //a vector of vertex storing all havent visited vertices
  
  for(auto iter = d.begin(); iter != d.end(); iter++){
    unvisited.push_back(*iter);                                                       //add all vertices to the vector
  }
  
  while(!unvisited.empty()){
    std::vector<vertex> visited;                                                      //temp location storing all weak connected components
    visited.push_back(*unvisited.begin());
    depth_first_componentsHelper(undirected,*unvisited.begin(),visited);
    componentsResult.push_back(visited);                                              //after founding all weak connected components add the vector to the Result
    for(auto v : visited){
      unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), v) , unvisited.end());  //remove all visited vertices from the unvisited vector
    }
  }
  
  return componentsResult;
}

template <typename vertex>
directed_graph<vertex> undirected_graph(const directed_graph<vertex> & d) {
  directed_graph<vertex> g(d);
  
  for(auto iter = d.begin(); iter!= d.end(); iter++)
    for(auto neighbour = d.nbegin(*iter); neighbour != d.nend(*iter); neighbour++)
      g.add_edge(*neighbour,*iter);                                                   //adding reverse edge to make the graph undirected
  

  return g;
}

template <typename vertex>
void depth_first_componentsHelper(const directed_graph<vertex> & d, vertex u, std::vector<vertex> & visited){
  std::stack<vertex> unprocessed;
  unprocessed.push(u);
  while(!unprocessed.empty()){
    vertex current = unprocessed.top();
    unprocessed.pop();
    //visited.push_back(current);
    //std::cout << current << std::endl;
    for(auto iter = d.nbegin(current); iter != d.nend(current); iter++){              //iteration through all neighbour of current vertices
      if(std::find(visited.begin(),visited.end(), *iter) == visited.end()){           //if the iterator's pointing vertex is not in visited vector
        visited.push_back(*iter);                                                     //add it in to visited
        unprocessed.push(*iter);                                                      //add it in to unprocessed
        //std::cout << *iter << std::endl;
        
      }
    }
  }
  //std::cout << visited.size() << std::endl;
}
/*
 * Computes the strongly connected components of the graph.
 * A strongly connected component is a subset of the vertices
 * such that for every pair u, v of vertices in the subset,
 * v is reachable from u and u is reachable from v.
 */
// template <typename vertex>
// std::vector<vertex> strong_component(const directed_graph<vertex> & d, const vertex & u){
//   std::vector<vertex> component;
//   std::unordered_map<vertex, bool> visited;
//   std::stack<vertex> stack;
//   stack.push(u);
//   while(!stack.empty()){
//     vertex curr = stack.top();
//     stack.pop();
//     visited[curr] = true;
//     component.push_back(curr);
//     for(auto iter = d.nbegin(curr); iter != d.nend(curr); iter++){
//        if(visited.find(*iter) == visited.end()){
//          visited[*iter] = true;
//          component.push_back(*iter);
//        }else{
//          return component;
//        }
//     }
//   }
//   return std::vector<vertex>();
// }

template <typename vertex>
std::vector<std::vector<vertex>> strongly_connected_components(const directed_graph<vertex> & d) {
  std::vector<std::vector<vertex>> result;
  
  if(d.num_vertices() == 0)                                            //Test Empty graph
    return std::vector<std::vector<vertex>>();
  else if(d.num_vertices() == 1){                                      //Test Single vertex
    std::vector<vertex> component;
    component.push_back(*d.begin());                                   //Save vertex in the temp component
    result.push_back(component);                                       //Add to result
    return result;
  }else if(d.num_edges() == 0){                                        //Test EdgeslessGraph
    for(auto iter = d.begin(); iter != d.end(); iter++){               //Iterate all vertex
      std::vector<vertex> temp;                                        //For each vertex it will be saved in its own vector                
      temp.push_back(*iter);                                           //Because each of them will be their own component
      result.push_back(temp);                                          //add them to the result
    }
    return result;
  }
  directed_graph<vertex> g(d);                                         //Test Single Dag or several Dag Strong components
  for(;;){                                                             //Infinte loop
  	auto iter = g.begin();                                             //a iterate that can be used outside the loop
  	for (iter = g.begin(); iter != g.end(); iter++){                   //iterating through all vertex
		  if(g.in_degree(*iter) == 0 || g.out_degree(*iter) == 0){         //gets the ones that is a leaf which is not possible for them to have a loop   
			  break;                                              
			  }
	  }
	if(iter != g.end()){                                                 //if a vertex meets the requirement above have been found
		std::vector<vertex> temp;                                          //save it in a temp vector and add it to the result as its own component
    temp.push_back(*iter);
    result.push_back(temp);                             
		g.remove_vertex(*iter);                                            //after it has been added, remove it then search again and see if any new leaf popped up due to the remove.
	}else{
		break;                                                             //after all a single component been added break the loop
	  }                                                       
  }
  
  
//   std::vector<vertex> strongTemp;
//   for(;;){                                                             //Infinte loop
//   	auto iter = g.begin();                                             //a iterate that can be used outside the loop
//   	for (iter = g.begin(); iter != g.end(); iter++){                   //iterating through all vertex
      
//       strongTemp = strong_component(g,*iter);
// 		  if(!strongTemp.empty()){                                             //gets the ones that is a leaf which is not possible for them to have a loop   
// 			  break;                                              
// 			  }
// 	  }
// 	if(iter != g.end()){                                                 //if a vertex meets the requirement above have been found
// 		                                                                   //save it in a temp vector and add it to the result as its own component
    
//     result.push_back(strongTemp);  
//     for(auto iter: strongTemp){
//       if(g.in_degree(iter) < 2 || g.out_degree(iter) < 2){
//         g.remove_vertex(iter);
//       }   
//     }
// 	}else{
// 		break;                                                             //after all a single component been added break the loop
// 	  }                                                       
//   }
  
  
  
  
  // std::unordered_map<vertex,bool> visited;
  // std::stack<vertex> stack;
  // //stack.push(*g.begin());
  // std::vector<vertex> strong_component;
  // while(!stack.empty()){
  //   vertex curr = stack.top();
  //   stack.pop();
  //   if(visited.find(curr) == visited.end()){  
  //     strong_component.push_back(curr);
  //     visited[curr] = true;
  //     for(auto iter = g.nbegin(curr); iter != g.nend(curr); iter++){
  //       stack.push(*iter);
  //       //strong_component.push_back(*iter);
  //     }
  //   }else{
  //     result.push_back(strong_component);
  //     for(vertex v: strong_component){
  //       if(g.in_degree(v) == 1 || g.out_degree(v) == 1){
  //         g.remove_vertex(v);
  //       }
  //     }
  //     //strong_component = std::vector<vertex>();
  //   }
  // }
  
  
  
  return result;
}





/*
 * Computes the shortest distance from u to every other vertex
 * in the graph d. The shortest distance is the smallest number
 * of edges in any path from u to the other vertex.
 * If there is no path from u to a vertex, set the distance to
 * be the number of vertices in d plus 1.
 */
template <typename vertex>
std::unordered_map<vertex, std::size_t> shortest_distances(const directed_graph<vertex> & d, const vertex & u) {
	std::unordered_map<vertex, std::size_t> distances;
  distances[u] = 0;
  
  std::queue<vertex> queue;                                        //bft
  queue.push(u);                                                   //add starting vertex to the queue
  while(!queue.empty()){                                           //run till queue empty
    vertex curr = queue.front();                                   //get first one in the queue
    queue.pop();                                                   //remove the first one in the queue
    for(auto iter = d.nbegin(curr); iter != d.nend(curr); iter++){ //find all curr's neighbour
      if(distances.find(*iter) == distances.end()){                //check if neighbour has already been added to the distances map
        queue.push(*iter);                                         //if not add it to the queue
        distances[*iter] = distances[curr] +1;                     //and add it to the map with its neighbour's distance +1
      }else{
        if(distances[*iter] > distances[curr] +1){                 //if neighbour is in the distance map
          distances[*iter] = distances[curr] +1;                   //check if current distance is short than the distance stored in the map
        }                                                          //if so update the shorter distance
      }
    }
  }                                                                //after checking all connected vertex from u
  if(distances.size() != d.num_vertices()){                        //check for if all vertices are added to the map
      for(auto iter = d.begin(); iter!= d.end(); iter++)           //not added ones are disconnected ones
        if(distances.find(*iter)== distances.end())                //loop through all vertex and check if its in the map
          distances[*iter] = d.num_vertices()+1;                   //if not in the map, add it to the map and set its distance to number of vertices +1
    }                                                              //just like the description says
  return distances;
}

















