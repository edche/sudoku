#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>
#include <unistd.h>

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

bool forwardCheck(int grid[N][N], int numPossible[N][N]) {
	for (int row = 0; row < N; ++row){
		for (int col = 0; col < N; ++col) {
			if (grid[row][col] == UNASSIGNED && numPossible[row][col] == 0) {
				return false;
			}
		}
	}
	return true;
}
bool solveBT(int grid[N][N], candIter start, candIter end) {
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

	vector<int>::iterator cIter;
	for (cIter = candidates.begin(); cIter != candidates.end(); ++cIter) {
		num = *cIter;
		if (checkValid(grid, row, col, num)) {
			grid[row][col] = num;
			nodesExpanded++;
			if (solveBT(grid, start + 1, end)) {
				return true;
			}

			// Failed somewhere
			grid[row][col] = UNASSIGNED;
		}
	}
	return false;
}

bool solveBTFC(int grid[N][N], bool possibleValues[N][N][N], int numPossible[N][N], candIter start, candIter end) {
	if (start == end) {
		return true;
	}
	int num, row, col;
	row = start->first;
	col = start->second;
	
	vector<int> candidates;
	for (int i = 0; i < N; ++i) {
		if (possibleValues[row][col][i])
			candidates.push_back(i+1);
	}
	random_shuffle(candidates.begin(), candidates.end());

	vector<int>::iterator cIter;
	for (cIter = candidates.begin(); cIter != candidates.end(); ++cIter) {
		num = *cIter;
		grid[row][col] = num;
		nodesExpanded++;

		// Update possible values table					
		// Update Rows and Cols
		for (int ii = 0; ii < N; ++ii) {
			if (possibleValues[row][ii][num-1]) {
				numPossible[row][ii] -= 1;
				possibleValues[row][ii][num-1] = false;
			}

			if (possibleValues[ii][col][num-1]) {
				numPossible[ii][col] -= 1;	
				possibleValues[ii][col][num-1] = false;
			}
		}

		// Update Box
		int startRow = (row/BOX)*BOX;
		int startCol = (col/BOX)*BOX;

		for (int ii = 0; ii < BOX; ++ii) {
			for (int jj = 0; jj < BOX; ++jj) {
				if (possibleValues[startRow + ii][startCol + jj][num-1]) {
					numPossible[startRow + ii][startCol + jj] -= 1;
				}
				possibleValues[startRow + ii][startCol + jj][num-1] = false;
			}
		}
		if (forwardCheck(grid, numPossible)) {
			if (solveBTFC(grid, possibleValues, numPossible, start + 1, end)) {
				return true;
			}
		}
		// Failed somewhere
		grid[row][col] = UNASSIGNED;
		
		// Add back num to the possible values
		for (int ii = 0; ii < N; ++ii) {
			if (checkValid(grid, row, ii, num)) { 
				if (!possibleValues[row][ii][num-1]) {
					numPossible[row][ii] += 1;
				}
				possibleValues[row][ii][num-1] = true;			
			}
			if (checkValid(grid, ii, col, num)) {
				if (!possibleValues[ii][col][num-1]) {
					numPossible[ii][col] += 1;
				}
				possibleValues[ii][col][num-1] = true;
				}
		}

		for (int ii = 0; ii < BOX; ++ii) {
			for (int jj = 0; jj < BOX; ++jj) {
				if (checkValid(grid, startRow + ii, startCol + jj, num)) {
				        if (!possibleValues[startRow + ii][startCol + jj][num-1]) {
						numPossible[startRow + ii][startCol + jj] += 1;
					}
					possibleValues[startRow + ii][startCol + jj][num-1] = true;
				}
			}
		}
	}
	return false;
}

bool solveBTFCH(int grid[N][N], bool possibleValues[N][N][N], int numPossible[N][N], vector<pair<int, int>> unassigned) {

	if (unassigned.empty()) {
		return true;
	}
	// Find most constrained variable

	int num, row, col;
	vector<pair<int, int>>::iterator uaIter;
	vector<pair<int, int>>::iterator iterIdx;
	pair<int, int> mcv;
	int maxConstr = 999;
	for (uaIter = unassigned.begin(); uaIter != unassigned.end(); ++uaIter) {
		row = uaIter->first;
		col = uaIter->second;
		if (numPossible[row][col] < maxConstr) {
			maxConstr = numPossible[row][col];
			iterIdx = uaIter;			
		}
	}

	mcv = *iterIdx;
	row = mcv.first;
	col = mcv.second;

	vector<int> candidates;
	for (int i = 0; i < N; ++i) {
		if (possibleValues[row][col][i]) {
			candidates.push_back(i+1);
		}
	}
	
	vector<int>::iterator cIter;
	
	// Now perform BTFC 
	for (cIter = candidates.begin(); cIter != candidates.end(); ++cIter) {
		num = *cIter;

		grid[row][col] = num;
		nodesExpanded++;

		// Update possible values table					
		// Update Rows and Cols
		for (int ii = 0; ii < N; ++ii) {
			if (possibleValues[row][ii][num-1]) {
				numPossible[row][ii] -= 1;
				possibleValues[row][ii][num-1] = false;
			}

			if (possibleValues[ii][col][num-1]) {
				numPossible[ii][col] -= 1;	
				possibleValues[ii][col][num-1] = false;
			}
		}

		// Update Box
		int startRow = (row/BOX)*BOX;
		int startCol = (col/BOX)*BOX;

		for (int ii = 0; ii < BOX; ++ii) {
			for (int jj = 0; jj < BOX; ++jj) {
				if (possibleValues[startRow + ii][startCol + jj][num-1]) {
					numPossible[startRow + ii][startCol + jj] -= 1;
				}
				possibleValues[startRow + ii][startCol + jj][num-1] = false;
			}
		}
		unassigned.erase(iterIdx);
		if (forwardCheck(grid, numPossible)) {
			if (solveBTFCH(grid, possibleValues, numPossible, unassigned)) {
				return true;
			}
		}
		// Failed somewhere
		grid[row][col] = UNASSIGNED;		
		unassigned.push_back(mcv);
		iterIdx = unassigned.end()-1;

		// Add back num to the possible values
		for (int ii = 0; ii < N; ++ii) {
			if (checkValid(grid, row, ii, num)) { 
				if (!possibleValues[row][ii][num-1]) {
					numPossible[row][ii] += 1;
				}
				possibleValues[row][ii][num-1] = true;			
			}
			if (checkValid(grid, ii, col, num)) {
				if (!possibleValues[ii][col][num-1]) {
					numPossible[ii][col] += 1;
				}
				possibleValues[ii][col][num-1] = true;
			}
		}

		for (int ii = 0; ii < BOX; ++ii) {
			for (int jj = 0; jj < BOX; ++jj) {
				if (checkValid(grid, startRow + ii, startCol + jj, num)) {
				        if (!possibleValues[startRow + ii][startCol + jj][num-1]) {
						numPossible[startRow + ii][startCol + jj] += 1;
					}
					possibleValues[startRow + ii][startCol + jj][num-1] = true;
				}
			}
		}
	}
	return false;

}	

int main(int argc, char** argv) {
	if (argc < 4) {
		cout << "Please type in 0, 1, 2, or 3 for difficulty level (0 Easy, 1 Medium, 2 Hard, 3 Evil)" << endl;
		cout << "followed by the number of trials you wish to perform" << endl;
		cout << "followed by the 0, 1, 2 to identify which method you wish to use." << endl;
		cout << "0 is Bactracking. 1 is BT with FC. 2 is BT with FC and Heuristics." << endl;
		exit(0);
	}
	
	int difficulty = atoi(argv[1]);
	int numTrials = atoi(argv[2]);
	int method = atoi(argv[3]);

	cout << "Difficulty = " << difficulty << " numTrials = " << numTrials << " method = " << method << endl;

	vector<int> nodeCount;
	vector<double> timeCount;
	
	srand(time(0));
	for (int trial = 0; trial < numTrials; ++trial) {
		nodesExpanded = 0;
		clock_t start;
		int grid[N][N];
		if (difficulty == EASY) {
			int filledGrid[N][N] = 	  {{0, 6, 1, 0, 0, 0, 0, 5, 2},
						  {8, 0, 0, 0, 0, 0, 0, 0, 1},
						  {7, 0, 0, 5, 0, 0, 4, 0, 0},
						  {9, 0, 3, 6, 0, 2, 0, 4, 7},
						  {0, 0, 6, 7, 0, 1, 5, 0, 0},
						  {5, 7, 0, 9, 0, 3, 2, 0, 6},
						  {0, 0, 4, 0, 0, 9, 0, 0, 5},
						  {1, 0, 0, 0, 0, 0, 0, 0, 8},
						  {6, 2, 0, 0, 0, 0, 9, 3, 0}};

			
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					grid[i][j] = filledGrid[i][j];
				}
			}
					
		} else if (difficulty == MEDIUM) {
			int filledGrid[N][N] = {{5, 0, 0, 6, 1, 0, 0, 0, 0},
					  {0, 2, 0, 4, 5, 7, 8, 0, 0},
					  {1, 0, 0, 0, 0, 0, 5, 0, 3},
					  {0, 0, 0, 0, 2, 1, 0, 0, 0},
					  {4, 0, 0, 0, 0, 0, 0, 0, 6},
					  {0, 0, 0, 3, 6, 0, 0, 0, 0},
					  {9, 0, 3, 0, 0, 0, 0, 0, 2},
					  {0, 0, 6, 7, 3, 9, 0, 8, 0},
					  {0, 0, 0, 0, 8, 6, 0, 0, 5}};
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					grid[i][j] = filledGrid[i][j];
				}
			}		


		} else if (difficulty == HARD) {
			int filledGrid[N][N] = {{0, 4, 0, 0, 2, 5, 9, 0, 0},
					  {0, 0, 0, 0, 3, 9, 0, 4, 0},
					  {0, 0, 0, 0, 0, 0, 0, 6, 1},
					  {0, 1, 7, 0, 0, 0, 0, 0, 0},
					  {6, 0, 0, 7, 5, 4, 0, 0, 9},
					  {0, 0, 0, 0, 0, 0, 7, 3, 0},
					  {4, 2, 0, 0, 0, 0, 0, 0, 0},
					  {0, 9, 0, 5, 4, 0, 0, 0, 0},
					  {0, 0, 8, 9, 6, 0, 0, 5, 0}};
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					grid[i][j] = filledGrid[i][j];
				}
			}	
		} else if (difficulty == EVIL) {
			int filledGrid[N][N] = {{0, 6, 0, 8, 2, 0, 0, 0, 0},
					  {0, 0, 2, 0, 0, 0, 8, 0, 1},
					  {0, 0, 0, 7, 0, 0, 0, 5, 0},
					  {4, 0, 0, 5, 0, 0, 0, 0, 6},
					  {0, 9, 0, 6, 0, 7, 0, 3, 0},
					  {2, 0, 0, 0, 0, 1, 0, 0, 7},
					  {0, 2, 0, 0, 0, 9, 0, 0, 0},
					  {8, 0, 4, 0, 0, 0, 7, 0, 0},
					  {0, 0, 0, 0, 4, 8, 0, 2, 0}};
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) {
					grid[i][j] = filledGrid[i][j];
				}
			}	
		}
		vector<pair<int, int>> unassigned = getUnassigned(grid);	
		start = clock();
		if (method == BT) {
			if (solveBT(grid, unassigned.begin(), unassigned.end())) {
				printGrid(grid);
			} else {
				cout << "No Valid Solution" << endl;
			}

		} else {
			// Construct Possible Values Table
			bool possibleValues[N][N][N] = {{false}};
			int numPossible[N][N] = {{0}};
			for (int row = 0; row < N; ++row) {
				for (int col = 0; col < N; ++col) {
					for (int num = 1; num <= N; ++num) {
						if (checkValid(grid, row, col, num)) {
							possibleValues[row][col][num-1] = true;
							numPossible[row][col] += 1;
						}
					}
				}
			}


			if (method == BTFC) {
				if (solveBTFC(grid, possibleValues, numPossible, unassigned.begin(), unassigned.end())) {
					printGrid(grid);
				} else {
					cout << "No Valid Solution" << endl;
				}

			} else if (method == BTFCH) {
				if (solveBTFCH(grid, possibleValues, numPossible, unassigned)) {
					printGrid(grid);
				} else {
					cout << "No Valid Solution" << endl;
				}

			}
		}

		cout << "Finished Trial #" << trial << endl;
		double runtime = (clock() - start) / (double)(CLOCKS_PER_SEC); 
		cout << "Time: " << runtime << " s" << endl;	
		cout << "Nodes expanded: " << nodesExpanded << endl;

		nodeCount.push_back(nodesExpanded);
		timeCount.push_back(runtime);		
	}
	double meanNode = accumulate(nodeCount.begin(), nodeCount.end(), 0.0) /nodeCount.size();
	double sd = 0;

	double meanTime = accumulate(timeCount.begin(), timeCount.end(), 0.0) / timeCount.size();

	for (int i = 0; i < numTrials; ++i) {
		sd += pow(nodeCount[i] - meanNode,2); 
	}

	double sdtime = 0;

	for (int i = 0; i < numTrials; ++i) {
		sdtime += pow(timeCount[i]- meanTime,2);		
	}
	sdtime = sqrt(sdtime/numTrials);
	sd = sqrt(sd/numTrials);
	cout << "********************************************" << endl;
	cout << "Finished! Average runtime " << meanTime << " s Std Dev: " << sdtime << endl;
	cout << "Average Nodes Expanded " << meanNode << " Std Dev: " << sd<< endl;
}
