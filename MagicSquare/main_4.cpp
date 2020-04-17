#include <bits/stdc++.h>

using namespace std;

bool isMagic(const int s[9])
{
	// Verify that elements are not repeated
	for (int i = 0; i < 8; i++)
		for (int j = i + 1; j < 9; j++)
			if (s[i] == s[j]) return false;

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
	for (int a3 = 1; a3 <= 9; a3++) {
	for (int a4 = 1; a4 <= 9; a4++) {
	for (int a5 = 1; a5 <= 9; a5++) {
	for (int a6 = 1; a6 <= 9; a6++) {
	for (int a7 = 1; a7 <= 9; a7++) {
	for (int a8 = 1; a8 <= 9; a8++) {
	for (int a9 = 1; a9 <= 9; a9++) {
		const int trialSquare[9] = {
			a1, a2, a3,
			a4, a5, a6,
			a7, a8, a9};
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
	const int square1[9] = {5, 3, 4, 1, 5, 8, 6, 4, 2};
	//const int square2[9] = {4, 9, 2, 3, 5, 7, 8, 1, 5};
	//const int square3[9] = {4, 8, 2, 4, 5, 7, 6, 1, 6};

	cout << formingMagicSquare(square1) << endl;
	return 0;
}
