#pragma once
#include <iostream>
#include <vector>

//const double INFINIT = 100000000;
template<class T>
class CMP {
public:
	virtual bool operator()(T elem1, T elem2) {
		return elem1 < elem2;
	}
	virtual ~CMP() {}
};

template<class T>
class CompareForMinHeap : public CMP<T> {
public:
	bool operator()(T elem1, T elem2) {
		return elem1 < elem2;
	}
};
template<class T>
class CompareForMaxHeap :public CMP<T> {
public:
	bool operator()(T elem1, T elem2) {
		return elem1 > elem2;
	}
};
template<class T1, class T2>
class BinHeap2 {
	std::vector<std::pair<T1, T2>> data; // приоритет, значение
	std::vector<int> index; //нумерация с 0
	size_t size;
	CMP<T1>* comp;
public:
	void Print() {
		for (int i = 0; i < size; i++) {
			std::cout << "{" << data[i].first << "," << data[i].second << "}\n";
		}
		std::cout << "\n";
	}
	int leftChild(int index) {
		if (2 * index + 1 < size) {
			return 2 * index + 1;
		}
		else return -1;
	}
	int rightChild(int index) {
		if (2 * index + 2 < size) {
			return 2 * index + 2;
		}
		else return -1;
	}
	int minChild(int index) {
		int rC = rightChild(index);
		int lC = leftChild(index);
		if (rC != -1 && lC != -1) {
			if (comp->operator()(data[lC].first, data[rC].first)) {
				return lC;
			}
			else return rC;
		}
		else {
			if (rC == -1) {
				return lC;
			}
			else return rC;
		}
	}
	int parent(int index) {
		if ((index != 0) && (index < size)) {
			return (index - 1) / 2;
		}
		else return -1;
	}
	void swap(int index1, int index2) {
		if (index1 != index2) {
			std::pair<T1, T2> tmp = data[index1];
			int tmpIndex = index[data[index1].second];
			index[data[index1].second] = index[data[index2].second];
			index[data[index2].second] = tmpIndex;

			data[index1] = data[index2];
			data[index2] = tmp;
		}
	}
	size_t getSize() {
		return size;
	}
	void diving(int index) {
		int j1 = index;
		int j2 = minChild(index);
		while ((j2 != -1) && (!comp->operator()(data[j1].first, data[j2].first))) { //(data[j1].first > data[j2].first)
			swap(j1, j2);
			j1 = j2;
			j2 = minChild(j1);
		}
	}
	void emersion(int index) {
		int j1 = index;
		int j2 = parent(index);
		while ((j2 != -1) && (comp->operator()(data[j1].first, data[j2].first))) {
			swap(j1, j2);
			j1 = j2;
			j2 = parent(j1);
		}
	}
	std::pair<T1, T2> getMin() {
		if (size > 0) {
			return data[0];
		}
		else return { INFINIT,INFINIT };
	}
	void push(std::pair<T1, T2> pair) {
		data.push_back(pair);
		index[pair.second] = size;
		size++;
		emersion(size - 1);
	}
	void pop() {
		swap(0, size - 1);
		index[data[size - 1].second] = INFINIT;
		data.pop_back();
		size--;
		diving(0);
	}
	BinHeap2(const std::vector<std::pair<T1, T2>>& data_, CMP<T1>* comp_, int N) {
		index = std::vector<int>(N, INFINIT);
		size = data_.size();
		comp = comp_;
		data = data_;
		for (int i = size - 1; i >= 0; i--) {
			index[data[i].second] = i;
		}
		makeHeap();
	}
	void  makeHeap() {
		for (int i = size - 1; i >= 0; i--) {
			diving(i);
		}
	}
	BinHeap2& operator=(const BinHeap2& B) {
		data = B.data;
		index = B.index;
		size = B.size;
		comp = B.comp;
		return *this;
	}
	void printIndex() {
		for (int i = 0; i < index.size(); i++) {
			std::cout << i << "-" << index[i] << " ";
		}
		std::cout << std::endl;
	}
	int getIndex(int i) {
		return index[i];
	}
	bool decreaseKey(int i, double newValue) {
		if (index[i] != INFINIT) { // вынести проверку наверх и добавить функцию проверки, есть ли индекс в очереди 
			data[index[i]].first = newValue;
			emersion(index[i]);
			return 1;
		}
		return 0;
	}
	bool deleteKey(int i) {
		if (index[i] != INFINIT) {
			int tmp = index[i];
			swap(index[i], size - 1);
			data.pop_back();
			index[i] = INFINIT;
			size--;
			diving(tmp);
			emersion(tmp);
			return 1;
		}
		return 0;
	}
	/*bool deleteKey(int i) {
		if (index[i] != INFINIT) {
			int tmp = index[i];
			swap(index[i], size - 1);
			data.pop_back();
			index[i] = INFINIT;
			size--;
			while ((tmp != -1) || (tmp >= 0)) { 
				diving(tmp);
				tmp = parent(tmp);
			}
			return 1;
		}
		return 0;
	}*/
};
template<class T1, class T2>
class priority_queue2 {
	BinHeap2<T1, T2> data;
	size_t size;
	CMP<T1>* comp;
public:
	bool decreaseKey(int i, double newValue) { //заодно показывает, есть ли вершина в очереди или нет
		return data.decreaseKey(i, newValue);
	}
	bool deleteKey(int i) { //заодно показывает, есть ли вершина в очереди или нет
		size--;
		return data.deleteKey(i);
	}
	size_t Size() {
		return size;
	}
	priority_queue2(std::vector<std::pair<T1, T2>> data_= std::vector<std::pair<T1, T2>>(0), CMP<T1>* comp_ = new CompareForMinHeap<double>(), int N = 0) :data(BinHeap2<T1, T2>(data_, comp_, N)) {
		comp = comp_;
		size = data_.size();
	}
	priority_queue2(priority_queue2&& Q) {
		data = Q.data;
		size = Q.size;
		comp = Q.comp;
	}
	void push(std::pair<T1, T2> pair) {
		data.push(pair);
		size++;
	}
	std::pair<T1, T2> top() {
		return data.getMin();
	}
	void pop() {
		data.pop();
		size--;
	}
	void printQueue() {
		data.Print();
	}
	priority_queue2& operator=(const priority_queue2& Q) {
		data = Q.data;
		size = Q.size;
		comp = Q.comp;
		return *this;
	}
	priority_queue2& operator=(priority_queue2&& Q) {
		data = Q.data;
		size = Q.size;
		comp = Q.comp;
		return *this;
	}
};

