#include <bits/stdc++.h>

using namespace std;

bool isMagic(const vector<vector<int>>& s)
{
	// Verify that elements are not repeated
	vector<int> aux = {s[0][0], s[0][1], s[0][2], s[1][0], s[1][1], s[1][2], s[2][0], s[2][1], s[2][2]};

	for (int i = 0; i < 8; i++)
		for (int j = i + 1; j < 9; j++)
			if (aux[i] == aux[j]) return false;

	int sum0 = s[0][0] + s[0][1] + s[0][2];

	int sum = s[1][0] + s[1][1] + s[1][2]; if (sum != sum0) return false;
	sum = s[2][0] + s[2][1] + s[2][2]; if (sum != sum0) return false;

	sum = s[0][0] + s[1][0] + s[2][0]; if (sum != sum0) return false;
	sum = s[0][1] + s[1][1] + s[2][1]; if (sum != sum0) return false;
	sum = s[0][2] + s[1][2] + s[2][2]; if (sum != sum0) return false;

	sum = s[0][0] + s[1][1] + s[2][2]; if (sum != sum0) return false;
	sum = s[2][0] + s[1][1] + s[0][2]; if (sum != sum0) return false;

	return true;
}

unsigned int costBetweenSquares(const vector<vector<int>>& s1,
	       	                const vector<vector<int>>& s2)
{
	unsigned int cost{};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			cost += abs(s1[i][j] - s2[i][j]);
	return cost;
}

unsigned int formingMagicSquare(vector<vector<int>> s)
{
	unsigned int minCost{1000};
	for (int a1 = 1; a1 <= 9; a1++) {
	for (int a2 = 1; a2 <= 9; a2++) {
	for (int a3 = 1; a3 <= 9; a3++) {
	for (int a4 = 1; a4 <= 9; a4++) {
	for (int a5 = 1; a5 <= 9; a5++) {
	for (int a6 = 1; a6 <= 9; a6++) {
	for (int a7 = 1; a7 <= 9; a7++) {
	for (int a8 = 1; a8 <= 9; a8++) {
	for (int a9 = 1; a9 <= 9; a9++) {
		vector<vector<int>> trialSquare = {
			{a1, a2, a3},
			{a4, a5, a6},
			{a7, a8, a9}};
		bool isMagicF = isMagic(trialSquare);
		if (isMagicF) {
			unsigned int cost = costBetweenSquares(trialSquare, s);
			if (cost < minCost) {
				minCost = cost;
			}
		}
	}}}}}}}}}

	return minCost;
}

int main()
{
	vector<vector<int>> square1 = {{5, 3, 4}, {1, 5, 8}, {6, 4, 2}};
	//vector<vector<int>> square2 = {{4, 9, 2}, {3, 5, 7}, {8, 1, 5}};
	//vector<vector<int>> square3 = {{4, 8, 2}, {4, 5, 7}, {6, 1, 6}};

	cout << formingMagicSquare(square1) << endl;
	return 0;
}
