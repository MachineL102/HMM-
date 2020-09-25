#pragma once

#include <iostream>
#include <vector>

typedef int real;

using namespace std;

class HMM
{
	int N;	// 状态种类
	int M;	// 观测种类
	int Tmax;
	vector<vector<int>> A;		// N * N
	vector<vector<int>> B;		// N * M
	vector<int> P0;				// N , 每个状态的初始概率
	vector<vector<int>> alpha;	// Tmax * N		第t时刻状态为j的前向变量
	vector<vector<int>> beta;	// Tmax * N		第t时刻状态为j的后向变量
	
	vector<vector<int>>* pA;
	vector<vector<int>>* pB;
	vector<int>* pP0;
	vector<vector<int>>* pbeta;
	vector<vector<int>>* palpha;

public:
	HMM();
	HMM(vector<vector<int>>& TransProb, vector<vector<int>>& ObsProb, vector<int>& InitProb);
	~HMM();
	real GetProbabilityF(vector<int>& ObsIdxs);
	real GetProbabilityB(vector<int>& ObsIdxs);
	vector<int> ViterbiStateIdxs(vector<int>& ObsIdxs);
	vector<int> PosteriorDecodingIdxs(vector<int>& ObsIdxs);
};
