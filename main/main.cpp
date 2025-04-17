#include "graph.h"
#include "generator.h"


int main() {
	//std::ifstream fo;
	//fo.open("graph_new300000.txt"); //114599 vertices 119666
	//std::ofstream fi;
	//fi.open("graph_withPparts300000.txt");
	//generator G;
	//G.createNewGraph(fo, fi,3);
	//fo.close();
	//fi.close();
	
	std::ifstream f;
	f.open("../../newTests/graph_new300000.txt");
	if (f.is_open()){ 
		graph G(f); 
		std::vector<double> result;
		auto start= std::chrono::steady_clock::now();
		result = G.findShortWays(7);
		auto end = std::chrono::steady_clock::now();
		std::cout << std::endl;
		std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end- start).count() << std::endl;
		
		/*std::cout << "\n";
		for (int i = 0; i<G.getVertices(); i++) {
			if (result[i] != -initializConstD) {
				std::cout << result[i] << " " << i << "\n";
			}
		}*/
	}
	f.close();
	return 0;
}