#pragma once
#include <fstream> 
#include <string>
int START = 0;
int END = 100000;
int K1 = 1000000;
int K2 = K1%100000+1;
double K3 = K1 * 1.5;
class generator {
public:

	std::pair<int,int> parsString(const std::string s) {
		std::pair<int, int> result;
		int j = 0;
		std::string first, second;
		while (s[j] != 32 && !s.empty() && (j < s.size())) {
			first += s[j];
			j++;
			}
		result.first = stoi(first);
		if (j + 1 < s.size()) {
			j++;
			while (j < s.size()) {
				second += s[j];
				j++;
			}
			result.second = stoi(second);
		}
		else result.second = -1;
		return result;

	}
	void createNewGraph(std::ifstream& file_out, std::ofstream& file_in,int P) {
		std::string str;
		getline(file_out, str);
		int countV= (int)parsString(str).first;
		std::vector<std::vector<int>> V(P);
		int j = 0, i = 0;
		while (j < countV) {
			if ((j % (countV / P + 1) == 0) && (j != 0)) {
				i++;
			}
			V[i].push_back(j);
			j++;
		}
		file_in << countV << "\n";
		for (int p = 0; p < P; p++) {
			for (int k = 0; k < V[p].size()-1; k++) {
				int x = V[i][k];
				int y = V[i][k+1];
				double r= rand() % (END - START + 1) + x % END;
				file_in << x << " " << y << " " << r << "\n";
			}
		}
	}
	void addWeightForVertices(std::ifstream& file_out, std::ofstream& file_in) {
		if (file_out.is_open()&&file_in.is_open()) {
			int x=0, y=0;
			int count;
			std::string str;
			getline(file_out, str);
			count = (int)parsString(str).first; //считали число вершин
			file_in << K1 << "\n"; //записали число вершин
			int m = 0;
			srand(time(NULL));
			while (m<K3) {
				getline(file_out, str);
				m++;
				x = (int)parsString(str).first;
				y = (int)parsString(str).second;
				if (y >= K1) {
					y = (rand()+rand())%K1;
					//std::cout << "y = " << y<< std::endl;
				}
				if (x >= K1) {
					x = (rand() + rand()) % K1;
					//std::cout << "x = " << x << std::endl;
				}
				if (y == -1) {
					y = (x+rand())%K1;
				}
				if (x >=K1 || y >= K1) std::cout << x << " " << y << std::endl;
				double r = rand() % (END - START + 1) + x%END;
				
				file_in << x << " " << y << " " << r << "\n";
			}
		}
	}
};
