#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <vector>
#include <cmath>
#include <time.h>

#define N 400000

using namespace std;
struct TreeNode {
	int number;
	TreeNode* left;
	TreeNode* right;
	
	TreeNode(int n) {
		number = n;
		left = right = nullptr;
	}

	static void destroy(TreeNode* node) {
		if (node) {
			destroy(node->left);
			destroy(node->right);
			delete node;
		}
	}
};

class Ordenamiento {
public:
	static void bubbleSort(int a[], int n);
	static void cocktailSort(int a[], int n);
	static void insertionSort(int a[], int n);
	static void bucketSort(int a[], int n);
	static void countingSort(vector <int>& a);
	static void countingSort(int a[], int n, int exp);
	static void mergeSort(int a[], int begin, int end, int n);
	static void merge(int a[], int begin, int middle, int end, int n);
	static void binaryTreeSort(int a[], int n);
	static void radixSort(int a[], int n);
	static void shellSort(int a[], int n);
	static void selectionSort(int a[], int n);
	static void heapSort(int a[], int n);
	static void quickSort(int a[], int lo, int hi);
	static int getMax(int a[], int n);
private:
	static void insert(TreeNode*& node, int number);
	static void inOrder(TreeNode* node, vector<int>& a);
	static void heapify(int a[], int n);
	static void siftDown(int a[], int start, int end);
	static int partition(int a[], int lo, int hi);
};

void Ordenamiento::bubbleSort(int a[], int n) {
	bool swapped = true;
	int j = 0;

	while (swapped) {
		swapped = false;
		j++;
		for (int i = 0; i < n - j; ++i) {
			if (a[i] > a[i + 1]) {
				swap(a[i], a[i + 1]);
				swapped = true;
			}
		}
	}
}

void Ordenamiento::cocktailSort(int a[], int n) {
	bool swapped = true;

	while (swapped) {
		swapped = false;
		
		for (int i = 0; i < n - 1; ++i) {
			if (a[i] > a[i + 1]) {
				std::swap(a[i], a[i + 1]);
				swapped = true;
			}
		}

		if (!swapped)
			break;
		swapped = false;

		for (int i = n - 1; i > 0; --i) {
			if (a[i - 1] > a[i]) {
				swap(a[i], a[i - 1]);
				swapped = true;
			}
		}
	}
}

void Ordenamiento::insertionSort(int a[], int n) {
	for (int i = 1; i < n; ++i) {
		for (int j = i; j > 0 && a[j - 1] > a[j]; --j)
			swap(a[j], a[j - 1]);
	}
}

void Ordenamiento::bucketSort(int* a, int n)
{
    int minValue = a[0];
    int maxValue = a[0];
    
    for (int i = 1; i < n; i++)
    {
        if (a[i] > maxValue)
            maxValue = a[i];
        if (a[i] < minValue)
            minValue = a[i];
    }
    
    int bucketLength = maxValue - minValue + 1;
    vector<int>* bucket = new vector<int>[bucketLength];
    
    for (int i = 0; i < bucketLength; i++)
    {
        bucket[i] = vector<int>();
    }
    
    for (int i = 0; i < n; i++)
    {
        bucket[a[i] - minValue].push_back(a[i]);
    }
    
    int k = 0;
    for (int i = 0; i < bucketLength; i++)
    {
        int bucketSize = bucket[i].size();
        
        if (bucketSize > 0)
        {
            for (int j = 0; j < bucketSize; j++)
            {
                a[k] = bucket[i][j];
                k++;
            }
        }
    }
}

void Ordenamiento::countingSort(vector <int>& arr)
{
    int max = *max_element(arr.begin(), arr.end());
    int min = *min_element(arr.begin(), arr.end());
    int range = max - min + 1;
    
    vector<int> count(range), output(arr.size());
    for(int i = 0; i < arr.size(); i++)
        count[arr[i]-min]++;
    
    for(int i = 1; i < count.size(); i++)
        count[i] += count[i-1];
    
    for(int i = arr.size()-1; i >= 0; i--)
    {
        output[ count[arr[i]-min] -1 ] = arr[i];
        count[arr[i]-min]--;
    }
    
    for(int i=0; i < arr.size(); i++)
        arr[i] = output[i];
}
void Ordenamiento::countingSort(int a[], int n, int exp)
{
    int output[n]; // output array
    int i, count[10] = {0};

    for (i = 0; i < n; i++)
        count[ (a[i]/exp)%10 ]++;

    for (i = 1; i < 10; i++)
        count[i] += count[i - 1];
    
    for (i = n - 1; i >= 0; i--)
    {
        output[count[ (a[i]/exp)%10 ] - 1] = a[i];
        count[ (a[i]/exp)%10 ]--;
    }
    

    for (i = 0; i < n; i++)
        a[i] = output[i];
}

void Ordenamiento::mergeSort(int a[], int begin, int end, int n) {
	if (end - begin < 2)
		return;

	int middle = (begin + end) / 2;
	mergeSort(a, begin, middle, n);
	mergeSort(a, middle, end, n);
	merge(a, begin, middle, end, n);
}

void Ordenamiento::merge(int a[], int begin, int middle, int end, int n) {
	int i = begin, j = middle;
	int aux[n];

	for (int k = begin; k < end; ++k) {
		if (i < middle && (j >= end || a[i] <= a[j]))
			aux[k] = a[i++];
		else
			aux[k] = a[j++];
	}

	for (int k = begin; k < end; ++k)
		a[k] = aux[k];
}

void Ordenamiento::binaryTreeSort(int a[], int n) {
	TreeNode* root = nullptr;
	for (int i = 0; i < n; ++i)
		insert(root, a[i]);

	vector<int> aux;
	inOrder(root, aux);

	TreeNode::destroy(root);

	for (int i = 0; i < aux.size(); ++i)
		a[i] = aux[i]; 
}

void Ordenamiento::insert(TreeNode*& node, int number) {
	if (!node)
		node = new TreeNode(number);
	else if (number < node->number)
		insert(node->left, number);
	else
		insert(node->right, number);
}

void Ordenamiento::inOrder(TreeNode* node, vector<int>& a) {
	if (node) {
		inOrder(node->left, a);
		a.push_back(node->number);
		inOrder(node->right, a);
	}
}
int Ordenamiento::getMax(int a[], int n)
{
    int mx = a[0];
    for (int i = 1; i < n; i++)
        if (a[i] > mx)
            mx = a[i];
    return mx;
}

void Ordenamiento::radixSort(int a[], int n)
{
    int m = getMax(a, n);
    for (int exp = 1; m/exp > 0; exp *= 10)
        countingSort(a, n, exp);
}

void Ordenamiento::shellSort(int a[], int n) {
	int gaps[8] = {701, 301, 132, 57, 23, 10, 4, 1};
	
	for (int i = 0; i < 8; ++i) {
		int gap = gaps[i];
		for (int j = gap; j < n; ++j) {
			int temp = a[j], k;
			for (k = j; k >= gap && a[k - gap] > temp; k -= gap)
				a[k] = a[k - gap];
			a[k] = temp;
		}
	} 
}

void Ordenamiento::selectionSort(int a[], int n) {
	for (int i = 0; i < n-1; ++i) {
		int min = i;
		for (int j = i+1; j < n; ++j) {
			if (a[j] < a[min])
				min = j;
		}
		
		if (min != i)
			swap(a[i], a[min]);
	}
}

void Ordenamiento::heapSort(int a[], int n) {
	heapify(a, n);
	int end = n - 1;
	while (end > 0) {
		swap(a[end--], a[0]);
		siftDown(a, 0, end);
	}
}

void Ordenamiento::heapify(int a[], int n) {
	int start = floor((n-2) / 2);
	while (start >= 0) {
		siftDown(a, start, n-1);
		--start;
	}
}

void Ordenamiento::siftDown(int a[], int start, int end) {
	int root = start;
	while (root * 2 + 1 <= end) {
		int child = root * 2 + 1;
		int sw = root;
		if (a[sw] < a[child])
			sw = child;
		if (child+1 <= end && a[sw] < a[child+1])
			sw = child + 1;
		if (sw == root)
			return;
		else {
			swap(a[root], a[sw]);
			root = sw;
		}
	}
}

void Ordenamiento::quickSort(int a[], int lo, int hi) {
	if (lo < hi) {
		int p = partition(a, lo, hi);
		quickSort(a, lo, p);
		quickSort(a, p+1, hi);
	}
}

int Ordenamiento::partition(int a[], int lo, int hi) {
	int pivot = a[lo];
	int i = lo - 1;
	int j = hi + 1;

	while (true) {
		do {
			--j;
		} while (a[j] > pivot);

		do {
			++i;
		} while (a[i] < pivot);
		
		if (i < j)
			swap(a[i], a[j]);
		else
			return j;
	}
}


int main() {
	int a[N], b[N];
	int i;
	vector<int> arr;
	for (int i = 0; i < N; ++i)
	{
		a[i] = rand() % N;
		arr.push_back(a[i]);
		b[i] = a[i];
	}
	printf("Con un arreglo de %d", N);


	clock_t t;


	//Quick
	t = clock();
	Ordenamiento::quickSort(a,0, N-1);
	t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Quick Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Heap
	t = clock();
	Ordenamiento::heapSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Heap Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}


	//Selection
	t = clock();
	Ordenamiento::selectionSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Selection Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Shell
	t = clock();
	Ordenamiento::shellSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Shell Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Radix
	t = clock();
	Ordenamiento::radixSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Radix Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Binarytree
	t = clock();
	Ordenamiento::binaryTreeSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Binary Tree Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}


	//Merge
	t = clock();
	Ordenamiento::mergeSort(a,0,N-1,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Merge Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Counting
	t = clock();
	Ordenamiento::countingSort(arr);
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Counting Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Bucket
	t = clock();
	Ordenamiento::bucketSort(a,N-1);//
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Bucket Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Insertion
	t = clock();
	Ordenamiento::insertionSort(a,N-1);
	t = clock() - t;
	time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Inserton Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Cocktail
	t = clock();
	Ordenamiento::cocktailSort(a,N-1);
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Cocktail Sort: %f ms \n", time_taken*1000);

	for (i = 0; i < N; ++i)
	{
		a[i] = b[i];
	}

	//Bubble
	t = clock();
	Ordenamiento::bubbleSort(a,N-1);
	t = clock() - t;
    time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("Bubble Sort: %f ms \n", time_taken*1000);

	return 0;
}
