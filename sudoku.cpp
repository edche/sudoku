#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

#define N 9
#define BT 0
#define BTFC 1
#define BTFCH 2
#define UNASSIGNED 0
#define BOX 3

#define EASY 0
#define MEDIUM 1
#define HARD 2
#define EVIL 3

typedef vector<pair<int, int>>::iterator candIter;

int nodesExpanded = 0;

void printGrid(int grid[N][N]) {
	for (int row = 0; row < N; ++row) {
		for (int col = 0; col < N; ++col) {
			cout << grid[row][col] << " ";	
		}
		cout << endl;
	}
}

vector<pair<int, int>> getUnassigned(int grid[N][N]) {
	vector<pair<int, int>> unassigned;
	for (int row = 0; row < N; ++row) {
		for (int col = 0; col < N; ++col) {
			if (grid[row][col] == UNASSIGNED) {
				unassigned.push_back(make_pair(row, col));
			}
		}
	}
	random_shuffle(unassigned.begin(), unassigned.end());
	return unassigned;
}	

bool checkValid(int grid[N][N], int row, int col, int num) {
	// Check rows and columns for duplicate
	for (int i = 0; i < N; ++i) {
		if (grid[row][i] == num || grid[i][col] == num) {
			return false;
		}
	}

	// Check box
	int start_row, start_col;
	start_row = (row/BOX)*BOX;
	start_col = (col/BOX)*BOX;
	for (int i = 0; i < BOX; ++i) {
		for (int j =0; j < BOX; ++j) {
			if (grid[start_row + i][start_col + j] == num) {
				return false;
			}
		}
	}
	return true;
}

bool solve(int grid[N][N], candIter start, candIter end) {
	if (start == end) {
		return true;
	}
	int num, row, col;
	row = start->first;
	col = start->second;

	vector<int> candidates;
	for (int i = 0; i < N; ++i) {
		candidates.push_back(i+1);
	}
	random_shuffle(candidates.begin(), candidates.end());

	for (int i = 0; i < candidates.size(); ++i) {
		num = candidates[i];
		if (checkValid(grid, row, col, num)) {
			grid[row][col] = num;
			nodesExpanded++;
			if (solve(grid, start + 1, end)) {
				return true;
			}

			// Failed somewhere
			grid[row][col] = UNASSIGNED;
		}
	}
	return false;

}
int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Please type in 0, 1, 2, or 3 for difficulty level" << endl;
		exit(0);
	}
	int difficulty = atoi(argv[1]);

	clock_t start;
	start = clock();

	if (difficulty == EASY) {
		int grid[N][N] = {{0, 6, 1, 0, 0, 0, 0, 5, 2},
				  {8, 0, 0, 0, 0, 0, 0, 0, 1},
				  {7, 0, 0, 5, 0, 0, 4, 0, 0},
				  {9, 0, 3, 6, 0, 2, 0, 4, 7},
				  {0, 0, 6, 7, 0, 1, 5, 0, 0},
				  {5, 7, 0, 9, 0, 3, 2, 0, 6},
				  {0, 0, 4, 0, 0, 9, 0, 0, 5},
				  {1, 0, 0, 0, 0, 0, 0, 0, 8},
				  {6, 2, 0, 0, 0, 0, 9, 3, 0}};

		//printGrid(grid);
		vector<pair<int, int>> unassigned = getUnassigned(grid);

		if (solve(grid, unassigned.begin(), unassigned.end())) {
			printGrid(grid);
		} else {
			cout << "No Valid Solution" << endl;
		}
		
		/*
		for (vector<pair<int, int>>::iterator it = unassigned.begin(); it != unassigned.end(); ++it) {
			cout << ' ' << it->first << ' ' << it->second << endl;
		}
		*/
	}
	cout << "Time: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;	
	cout << "Nodes expanded: " << nodesExpanded << endl;
}
