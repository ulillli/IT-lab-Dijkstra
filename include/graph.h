#pragma once
#include <fstream> 
#include <queue>
#include "priority_queue.h"
#include "priority_queue2.h"
#include <omp.h>
#include <chrono>

const int initializConstI = -1000000;
const double initializConstD = -10000000.;
class DijkstraAlgorithm {
public: 
	std::vector<double> Algorithm1(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) { //отправить вершину с индексом [0,adjList.size()-1]
		int countV(adjList.size());
		std::vector<double> path(countV, -initializConstD);

		std::vector<bool> visited(countV);
		path[start] = 0; //обработка стартовой вершины

		for (int i = 0; i < countV; i++) {
			int near = -1;
			for (int j = 0; j < countV; j++) { //поиск ближайшей вершины
				if (!visited[j] && ((near == -1) || (path[j] < path[near]))) {
					near = j;
				}
			}
			visited[near] = true;
			for (int v = 0; v < adjList[near].size(); v++) { // релаксация 
				int currV = adjList[near][v].first;
				double currW = adjList[near][v].second;
				if (path[near] + currW < path[currV]) {
					path[currV] = path[near] + currW;
				}
			}
			if (i % 10000 == 0) {
				std::cout << i << std::endl;
			}
		}
		return path;
	}
	std::vector<double> Algorithm2(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) { //отправить вершину с индексом [0,adjList.size()-1]
		int countV(adjList.size());
		std::vector<std::pair<double, double>> data = {{0,start}};

		std::vector<double> path(countV, -initializConstD);

		priority_queue<double, double> Q(data);
		path[start] = 0;

		for (int i = 0; i < countV; i++) {
			int near = Q.top().second;
			if (near != INFINIT) {
				Q.pop();
				for (int v = 0; v < adjList[near].size(); v++) { // релаксация 
					int currV = adjList[near][v].first;
					double currW = adjList[near][v].second;

					if (path[near] + currW < path[currV]) {
						path[currV] = path[near] + currW;
						Q.push(std::make_pair(path[currV], currV));
					}
				}
			}
			else break;
		}
		return path;
	}
	std::vector<double> Algorithm3(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) { //отправить вершину с индексом [0,adjList.size()-1]
		int countV(adjList.size());
		std::vector<std::pair<double, double>> data = { {0,start} };

		std::vector<double> path(countV, -initializConstD);
		std::priority_queue<std::pair<double, double>, std::vector<std::pair<double, double>>, CompareForStdQueue<std::pair<double, double>>> Q;
		Q.push({ 0,start });
		path[start] = 0;
		for (int i = 0; i < countV; i++) {
				int near = Q.top().second;
				if (near != INFINIT) {
					{
					Q.pop(); 
					}
					for (int v = 0; v < adjList[near].size(); v++) { // релаксация 
						int currV = adjList[near][v].first;
						double currW = adjList[near][v].second;
						{
							if (path[near] + currW < path[currV]) {
								path[currV] = path[near] + currW; 
								Q.push(std::make_pair(path[currV], currV)); 
							}
						}
				}
			}
		}
		return path;
	}
	// первым проходом bfs делим на области связности 
	// каждую области делим bfs

	std::vector<std::vector<int>> devision_with_bfs(const std::vector<std::vector<std::pair<int, double>>> adjList, int start, int P) {
		std::vector<std::vector<int>> V(P);
		std::queue<std::pair<int,double>> line;
		line.push({ start, 0});
		while (!line.empty()) {
			std::pair<int, double> tmp = line.front();
			line.pop();
			/*for (int i = 0; i < adjList[tmp.first].size(); i++) {
				if (distance[data[tmp][i]] == INF) {
					distance[data[tmp][i]] = distance[tmp] + 1;
					line.push(data[tmp][i]);
				}
			}*/
		}
		/*std::vector<int> distance(countV, INF);
		std::queue<int> line;

		line.push(start);
		distance[start] = 0;

		while (!line.empty()) {
			int tmp = line.front();
			line.pop();
			for (int i = 0; i < data[tmp].size(); i++) {
				if (distance[data[tmp][i]] == INF) {
					distance[data[tmp][i]] = distance[tmp] + 1;
					line.push(data[tmp][i]);
				}
			}
		}*/

		return V;
	}
	std::vector<double> Algorithm4(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) {
		//номер каждой вершины в графе, должен быть меньше countV, чтобы вообще построить граф (нумерация с 0 и до countV)
		int countV=adjList.size();

		int P=1; //count of streams

		std::vector<std::vector<int>> V(P); //разделение вершин на потоки
		int j = 0,i=0;
		while (j <countV) {
			if ((j % (countV/P+1)==0)&&(j!=0)) {
				i++;
			}
			V[i].push_back(j);
			j++;
		}
		std::vector<int> indexOfStream(countV);
		for (int i = 0; i < P; i++) {
			for (int j = 0; j < V[i].size(); j++) {
				indexOfStream[V[i][j]] = i;
			}
		}
		//creating queues 
		CompareForMinHeap<double> C;
		std::vector<priority_queue2<double, int>> Q_d(P);
		std::vector<priority_queue2<double, int>> Q_in(P);
		std::vector<priority_queue2<double, int>> Q_out(P);
		for (int i = 0; i < P; i++) {
			Q_d[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_in[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_out[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
		}
		std::vector<double> d(countV, -initializConstD);
		std::vector<double> in(countV, -initializConstD); //min входящее ребро в вершину
		std::vector<double> out(countV, -initializConstD); //min исходящее ребро из вершины 
		//инициализация массивов in, out
		for (int i = 0; i < countV; i++) {
			for (int j=0; j < adjList[i].size(); j++) {
				out[i] = min(out[i],adjList[i][j].second);
				in[adjList[i][j].first] = min(in[adjList[i][j].first], adjList[i][j].second);
			}
		}
		std::vector<std::pair<double, int>> S;

		d[start] = 0.;
		S.push_back({0.,start});
		Q_d[0].push({ 0.,start });
		Q_in[0].push({ 0.-in[start],start});
		Q_out[0].push({ 0.+out[start],start });
		omp_set_num_threads(2);
		std::vector<bool> D(countV,0); //показывает, вычислено ли расстояние до вершины или нет
		D[start] = 0;
		auto start_ = std::chrono::steady_clock::now();
		int iterations = 0;
		while (checkQueues(Q_d)) { 
			iterations++; // посчитать количество итераций 
			double L=-initializConstD, M=-initializConstD;
#pragma omp parallel for
			for (int i = 0; i < P; i++) {
				L = min(L, Q_out[i].top().first); 
				M = min(M, Q_d[i].top().first);
			}
			std::vector<std::vector<int>> R(P); //формируем множество R и удаляем вершины из очередей шаг 3,4
			auto start1 = std::chrono::steady_clock::now();
#pragma omp parallel for
			for (int i = 0; i < P; i++) { 
				while ((Q_d[i].top().first <= L)) { 
					R[i].push_back(Q_d[i].top().second);
					D[Q_d[i].top().second] = 1;
					Q_in[i].deleteKey(Q_d[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_d[i].top().second);
					Q_d[i].pop();
				}
				while ((Q_in[i].top().first <= M)) {
					R[i].push_back(Q_in[i].top().second);
					D[Q_in[i].top().second] = 1;
					Q_d[i].deleteKey(Q_in[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_in[i].top().second);
					Q_in[i].pop();
				}
			}
			auto end1 = std::chrono::steady_clock::now();
			std::vector<std::vector<std::pair<double, int>>> Z(P); // формируем множество пар {Xi,Zi} шаг 5
			auto start2 = std::chrono::steady_clock::now();
#pragma omp parallel for
			for (int i = 0; i < P; i++) { //выбираем номер потока 
#pragma omp critical
				for (int j = 0; j < R[i].size(); j++) { //выбрали вершину из i-го потока
					for (int k = 0; k < adjList[R[i][j]].size(); k++) { //выбираем смежную ей вершину из списка смежности 
						std::pair<int, double> vertex = adjList[R[i][j]][k];
						if (D[vertex.first]==0) { // проверяем смежную вершину на "необработанность"
							Z[i].push_back({ vertex.second,vertex.first}); // добавляем {ребро, смежная вершина}
							if (d[vertex.first] > d[R[i][j]] + vertex.second) {
								d[vertex.first] = d[R[i][j]] + vertex.second; 
								for (int p = 0; p < P; p++) { //обновление очередей 
									if (!(Q_d[p].decreaseKey(vertex.first, d[R[i][j]] + vertex.second)) && (indexOfStream[R[i][j]] == p))
									{
										Q_d[p].push({ d[R[i][j]] + vertex.second,vertex.first});
									} 
									if (!(Q_in[p].decreaseKey(vertex.first, d[R[i][j]] + vertex.second-in[vertex.first]))&&(indexOfStream[R[i][j]]==p))
									{Q_in[p].push({ d[R[i][j]] + vertex.second - in[vertex.first],vertex.first});
								}
									if (!(Q_out[p].decreaseKey(vertex.first, d[R[i][j]] + vertex.second+out[vertex.first]))&&(indexOfStream[R[i][j]]==p))
									{Q_out[p].push({d[R[i][j]] + vertex.second + out[vertex.first],vertex.first});
									}
								}
							}
						}
					}
				}
			}
			auto end2 = std::chrono::steady_clock::now();
			//std::cout << std::endl << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
			//std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << std::endl;
		}
		std::cout << "cout of iterations = " << iterations << std::endl;
		return d;
	}
	std::vector<double> Algorithm5(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) {
		//номер каждой вершины в графе, должен быть меньше countV, чтобы вообще построить граф (нумерация с 0 и до countV)
		int countV = adjList.size();

		int P = 1; //count of streams

		std::vector<std::vector<int>> V(P); //разделение вершин на потоки
		int j = 0, i = 0;
		while (j < countV) {
			if ((j % (countV / P + 1) == 0) && (j != 0)) {
				i++;
			}
			V[i].push_back(j);
			j++;
		}
		std::vector<int> indexOfStream(countV);
		for (int i = 0; i < P; i++) {
			for (int j = 0; j < V[i].size(); j++) {
				indexOfStream[V[i][j]] = i;
			}
		}
		//creating queues 
		CompareForMinHeap<double> C;
		std::vector<priority_queue2<double, int>> Q_d(P);
		std::vector<priority_queue2<double, int>> Q_in(P);
		std::vector<priority_queue2<double, int>> Q_out(P);
		for (int i = 0; i < P; i++) {
			Q_d[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_in[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_out[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
		}
		std::vector<double> d(countV, -initializConstD);
		std::vector<double> in(countV, -initializConstD); //min входящее ребро в вершину
		std::vector<double> out(countV, -initializConstD); //min исходящее ребро из вершины 
		//инициализация массивов in, out
		for (int i = 0; i < countV; i++) {
			for (int j = 0; j < adjList[i].size(); j++) {
				out[i] = min(out[i], adjList[i][j].second);
				in[adjList[i][j].first] = min(in[adjList[i][j].first], adjList[i][j].second);
			}
		}
		std::vector<std::pair<double, int>> S;

		d[start] = 0.;
		S.push_back({ 0.,start });
		Q_d[0].push({ 0.,start });
		Q_in[0].push({ 0. - in[start],start });
		Q_out[0].push({ 0. + out[start],start });
		omp_set_num_threads(1);
		std::vector<bool> D(countV, 0); //показывает, вычислено ли расстояние до вершины или нет
		D[start] = 0;
		int time1 = 0,time2=0, time3=0, time4=0;
		auto start_ = std::chrono::steady_clock::now();
		int iterations = 0;
		while (checkQueues(Q_d)) {
			iterations++; // посчитать количество итераций 
			double L = -initializConstD, M = -initializConstD;
//#pragma omp parallel for
			auto start1 = std::chrono::steady_clock::now();
			for (int i = 0; i < P; i++) {
				L = min(L, Q_out[i].top().first);
				M = min(M, Q_d[i].top().first);
			}
			auto end1 = std::chrono::steady_clock::now();
			std::vector<std::vector<int>> R(P); //формируем множество R и удаляем вершины из очередей шаг 3,4
			auto start2 = std::chrono::steady_clock::now();
//#pragma omp parallel for
			int count1 = 0, count2 = 0;
			for (int i = 0; i < P; i++) {
				while ((Q_d[i].top().first <= L)) {
					count1++;
					R[i].push_back(Q_d[i].top().second);
					D[Q_d[i].top().second] = 1;
					Q_in[i].deleteKey(Q_d[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_d[i].top().second);
					Q_d[i].pop();
				}
				while ((Q_in[i].top().first <= M)) {
					count2++;
					R[i].push_back(Q_in[i].top().second);
					D[Q_in[i].top().second] = 1;
					Q_d[i].deleteKey(Q_in[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_in[i].top().second);
					Q_in[i].pop();
				}
			}
			//std::cout << count1 << " " << count2 << std::endl;
			auto end2 = std::chrono::steady_clock::now();
			std::vector<std::vector<std::pair<double, int>>> Z(P); // формируем множество пар {Xi,Zi} шаг 5
			auto start3 = std::chrono::steady_clock::now();
//#pragma omp parallel for
			for (int i = 0; i < P; i++) { //выбираем номер потока 
				for (int j = 0; j < R[i].size(); j++) { //выбрали вершину из i-го потока
					for (int k = 0; k < adjList[R[i][j]].size(); k++) { //выбираем смежную ей вершину из списка смежности 
						std::pair<int, double> vertex = adjList[R[i][j]][k];
						if (D[vertex.first] == 0) { Z[i].push_back({ d[R[i][j]] + vertex.second,vertex.first }); } // добавляем {ребро+расстояние до R[i][j], смежная вершина}
					}
				}
			}
			auto end3 = std::chrono::steady_clock::now();
//#pragma omp parallel for
			auto start4 = std::chrono::steady_clock::now();
			for (int i = 0; i < P; i++) {
				for (int j = 0; j < Z[i].size(); j++) {
					std::pair<double, int> vertex = Z[i][j]; //{расстояние до вершины, сама вершина}
					if (d[vertex.second] > vertex.first) {
						d[vertex.second] = vertex.first;
						for (int p = 0; p < P; p++) {
//#pragma omp critical
							if (!Q_d[p].decreaseKey(vertex.second, vertex.first)) {
								Q_d[p].push({ vertex.first,vertex.second }); 
							}
//#pragma omp critical
							if (!Q_in[p].decreaseKey(vertex.second, vertex.first - in[vertex.second])) {
								Q_in[p].push({ vertex.first - in[vertex.second],vertex.second });
							}
//#pragma omp critical
							if (!Q_out[p].decreaseKey(vertex.second, vertex.first + out[vertex.second])) {
								Q_out[p].push({ vertex.first + out[vertex.second],vertex.second });
							}
						}
					}
				}
			}
			auto end4 = std::chrono::steady_clock::now();
			time1 += std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
			time2 += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
			time3 += std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();
			time4 += std::chrono::duration_cast<std::chrono::milliseconds>(end4 - start4).count();
			//std::cout << std::endl << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
			//std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << std::endl;
		}
		std::cout << "time of the first cycle = " << time1 << std::endl;
		std::cout << "time of the second cycle = " << time2 << std::endl;
		std::cout << "time of the third cycle = " << time3 << std::endl;
		std::cout << "time of the fourth cycle = " << time4 << std::endl;

		std::cout << "cout of iterations = " << iterations << std::endl;
		return d;
	}
	std::vector<double> Algorithm6(const std::vector<std::vector<std::pair<int, double>>> adjList, int start) {
		
		int countV = adjList.size();

		int P = 4; 

		std::vector<std::vector<int>> V(P); 
		int j = 0, i = 0;
		while (j < countV) {
			if ((j % (countV / P + 1) == 0) && (j != 0)) {
				i++;
			}
			V[i].push_back(j);
			j++;
		}
		std::vector<int> indexOfStream(countV);
		for (int i = 0; i < P; i++) {
			for (int j = 0; j < V[i].size(); j++) {
				indexOfStream[V[i][j]] = i;
			}
		}
		//creating queues 
		CompareForMinHeap<double> C;
		std::vector<priority_queue2<double, int>> Q_d(P);
		std::vector<priority_queue2<double, int>> Q_in(P);
		std::vector<priority_queue2<double, int>> Q_out(P);
		for (int i = 0; i < P; i++) {
			Q_d[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_in[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
			Q_out[i] = priority_queue2<double, int>(std::vector<std::pair<double, int>>(), &C, countV);
		}
		std::vector<double> d(countV, -initializConstD);
		std::vector<double> in(countV, -initializConstD); //min входящее ребро в вершину
		std::vector<double> out(countV, -initializConstD); //min исходящее ребро из вершины 
		
		//инициализация массивов in, out
		for (int i = 0; i < countV; i++) {
			for (int j = 0; j < adjList[i].size(); j++) {
				out[i] = min(out[i], adjList[i][j].second);
				in[adjList[i][j].first] = min(in[adjList[i][j].first], adjList[i][j].second);
			}
		}
		
		std::vector<std::pair<double, int>> S;

		d[start] = 0.;
		Q_d[0].push({ 0.,start });
		Q_in[0].push({ 0. - in[start],start });
		Q_out[0].push({ 0. + out[start],start });
		omp_set_num_threads(4);
		std::vector<bool> D(countV, 0); //показывает, вычислено ли расстояние до вершины или нет
		D[start] = 0;
		auto start_ = std::chrono::steady_clock::now();
		int iterations = 0;
		int time1 = 0, time2 = 0;
		int iterations_while1 = 0, iterations_while2 = 0,iterations_for=0;
		while (checkQueues(Q_d)) {
			iterations++; // посчитать количество итераций 
			double L = -initializConstD, M = -initializConstD;
#pragma omp parallel for
			for (int i = 0; i < P; i++) {
				L = min(L, Q_out[i].top().first);
				M = min(M, Q_d[i].top().first);
			}
			
			std::vector<std::vector<int>> R(P); //формируем множество R и удаляем вершины из очередей шаг 3,4
			auto start1 = std::chrono::steady_clock::now();
#pragma omp parallel for
			for (int i = 0; i < P; i++) {
				while ((Q_d[i].top().first <= L)) {
					iterations_while1++;
					R[i].push_back(Q_d[i].top().second);
					D[Q_d[i].top().second] = 1;
					Q_in[i].deleteKey(Q_d[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_d[i].top().second);
					Q_d[i].pop();
				}
				while ((Q_in[i].top().first <= M)) {
					iterations_while2++;
					R[i].push_back(Q_in[i].top().second);
					D[Q_in[i].top().second] = 1;
					Q_d[i].deleteKey(Q_in[i].top().second); //или тут можно просто делать pop
					Q_out[i].deleteKey(Q_in[i].top().second);
					Q_in[i].pop();
				}
			}
			auto end1 = std::chrono::steady_clock::now();
			auto start2 = std::chrono::steady_clock::now();
			for (int i = 0; i < P; i++) { //выбираем номер потока 
				for (int j = 0; j < R[i].size(); j++) { //выбрали вершину из i-го потока
					for (int k = 0; k < adjList[R[i][j]].size(); k++) { //выбираем смежную ей вершину из списка смежности 
						std::pair<int, double> vertex = adjList[R[i][j]][k];
						if (D[vertex.first] == 0) { 
							double x = d[R[i][j]] + vertex.second;
							int z = vertex.first;
							int p =indexOfStream[R[i][j]];
							if (d[vertex.first] > x) {
								d[vertex.first] = x;
								iterations_for++;
								if (!Q_d[p].decreaseKey(z, x)) {
									Q_d[p].push({x,z});
								}
								if (!Q_in[p].decreaseKey(z, x-in[vertex.first])) {
									Q_in[p].push({x-in[vertex.first],z});
								}
								if (!Q_out[p].decreaseKey(z,x+out[vertex.first])) {
									Q_out[p].push({x+out[vertex.first],z});
								}
							}
						}
					}
				}
			}
			auto end2 = std::chrono::steady_clock::now();
			time1 += std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
			time2 += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
			//std::cout << std::endl << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
			//std::cout << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count() << std::endl;
		}
		auto end_ = std::chrono::steady_clock::now();
		std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count() << std::endl;
		/*std::cout << std::endl<<"Generall count of iterations = " << iterations << std::endl;
		std::cout << "The count of iterations in the first while = " << iterations_while1 << std::endl;
		std::cout << "The count of iterations in the second while = " << iterations_while2 << std::endl;
		std::cout << "The count of iterations in for = " << iterations_for << std::endl;*/
		std::cout << "time1 = " << time1 << " time2 = " << time2 << std::endl;
		return d;
	}
	double min(double value1, double value2) {
		if (value1 <= value2) {
			return value1;
		}
		else return value2;
	}
	bool checkQueues(std::vector<priority_queue2<double, int>>& Q) {
		bool flag = 0;
#pragma omp parallel for
		for (int i = 0; i < Q.size(); i++) {
			if (Q[i].Size()>0) {
				flag = 1;
				break;
			}
		}
		return flag;
	}
}; 


class fileParsing {
public:
	void ParsingInMatrix(std::ifstream& file) {
		int countV = 0;
		file >> countV;
		int v1 = initializConstI, v2 = initializConstI;
		double weight = initializConstD;
		std::vector<std::vector<double>> matrix(countV);
		for (int i = 0; i < countV; i++) matrix[i] = std::vector<double>(countV);

		while (file >> v1 >> v2 >> weight) {
			matrix[v1][v2] = weight;
		}
	}
	std::vector<std::vector<std::pair<int, double>>> ParsingInAdjList(std::ifstream& file) {
		int countV = 0;
		file >> countV;
		int v1 = initializConstI, v2 = initializConstI;
		double weight = initializConstD;
		std::vector<std::vector<std::pair<int, double>>> data(countV);
		std::cout << countV << std::endl;
		while (file >> v1 >> v2 >> weight) {
			data[v1 - 1].push_back(std::make_pair(v2 - 1, weight));
			data[v2 - 1].push_back(std::make_pair(v1 - 1, weight));
		}
		return data;
	}
};
class graph {
	std::vector<std::vector<std::pair<int, double>>> adjList;
	int countV;
public:
	int getVertices() { return countV; }
	int max(int x, int y) {
		if (x > y) return x;
		return y;
	}
	graph(std::ifstream& file) {
		fileParsing Pars;
		if (file.is_open()) {
			adjList = Pars.ParsingInAdjList(file);
			countV = adjList.size();
		}
		int count_zero_vertex = 0;
		int maxDegreeOfVertex = -1000000000;
		int sum = 0;
		for (int i = 0; i < countV; i++) {
			maxDegreeOfVertex = max(adjList[i].size(), maxDegreeOfVertex); 
			sum += adjList[i].size();
			if (adjList[i].size() == 0) count_zero_vertex++;
		}
		/*std::cout << "Maximum degree = " << maxDegreeOfVertex << std::endl;
		std::cout << "Zero vertex = " << (double)sum / countV << std::endl;
		std::cout << "Number of vertex = ";
		for (int i = 0; i < countV; i++) {
			if (adjList[i].size() == maxDegreeOfVertex) {
				std::cout << i << std::endl;
			}
		}*/
	}
	void printAdjList() {
		for (int i = 0; i < countV; i++) {
			std::cout << i << " ";
			for (int j = 0; j < adjList[i].size(); j++) {
				std::cout << adjList[i][j].first << " ";
			}
			std::cout << std::endl;
		}
	}
	
	std::vector<double> findShortWays(int start) {
		DijkstraAlgorithm D;
		return D.Algorithm6(adjList, start);
	}
};
