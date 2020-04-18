#include <bits/stdc++.h>

using namespace std;

bool isMagic(const int s[9])
{
	const int sum0 = s[0] + s[1] + s[2];

	if (s[3] + s[4] + s[5] != sum0) return false;
	if (s[6] + s[7] + s[8] != sum0) return false;

	if (s[0] + s[3] + s[6] != sum0) return false;
	if (s[1] + s[4] + s[7] != sum0) return false;
	if (s[2] + s[5] + s[8] != sum0) return false;

	if (s[0] + s[4] + s[8] != sum0) return false;
	if (s[6] + s[4] + s[2] != sum0) return false;

	return true;
}

unsigned int costBetweenSquares(const int s1[9], const int s2[9])
{
	unsigned int cost{};
	for (int i = 0; i < 9; i++)
		cost += abs(s1[i] - s2[i]);
	return cost;
}

unsigned int formingMagicSquare(const int s[9])
{
	unsigned int minCost{1000};
	for (int a1 = 1; a1 <= 9; a1++) {
	for (int a2 = 1; a2 <= 9; a2++) {
		if (a2 != a1)
	for (int a3 = 1; a3 <= 9; a3++) {
		if (a3 != a2 && a3 != a1)
	for (int a4 = 1; a4 <= 9; a4++) {
		if (a4 != a3 && a4 != a2 && a4 != a1)
	for (int a5 = 1; a5 <= 9; a5++) {
		if (a5 != a4 && a5 != a3 && a5 != a2 && a5 != a1)
	for (int a6 = 1; a6 <= 9; a6++) {
		if (a6 != a5 && a6 != a4 && a6 != a3 && a6 != a2 && a6 != a1)
	for (int a7 = 1; a7 <= 9; a7++) {
		if (a7 != a6 && a7 != a5 && a7 != a4 && a7 != a3 && a7 != a2 && a7 != a1)
	for (int a8 = 1; a8 <= 9; a8++) {
		if (a8 != a7 && a8 != a6 && a8 != a5 && a8 != a4 && a8 != a3 && a8 != a2 && a8 != a1)
	for (int a9 = 1; a9 <= 9; a9++) {
		if (a9 != a8 && a9 != a7 && a9 != a6 && a9 != a5 && a9 != a4 && a9 != a3 && a9 != a2 && a9 != a1) {
		const int trialSquare[9] = {a1, a2, a3, a4, a5, a6, a7, a8, a9};
		bool isMagicF = isMagic(trialSquare);
		if (isMagicF) {
			unsigned int cost = costBetweenSquares(trialSquare, s);
			if (cost < minCost) {
				minCost = cost;
			}
		}
		}
	}}}}}}}}}

	return minCost;
}

void printSquare(const int s[9])
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
			cout << s[i * 3 + j] << " "; 
	        cout << endl;
	}
	cout << endl;
}


int main()
{

	const int square1[9] = {5, 3, 4, 1, 5, 8, 6, 4, 2};
	//const int square2[9] = {4, 9, 2, 3, 5, 7, 8, 1, 5};
	//const int square3[9] = {4, 8, 2, 4, 5, 7, 6, 1, 6};

	cout << formingMagicSquare(square1) << endl;
	return 0;
}
