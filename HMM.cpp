#include "HMM.h"

HMM::HMM()
{
	N = 0;
	M = 0;
	Tmax = 0;
}

HMM::HMM(vector<vector<int>>& TransProb, vector<vector<int>>& ObsProb, vector<int>& InitProb)
{
	N = TransProb.size();
	M = ObsProb[0].size();
	/*
	  vector<vector<int>>  A(N,N);
	  vector<vector<int>>  B(N,M);
	  FArr1D  P0(N);
	*/
	pA = new vector<vector<int>>(TransProb);
	pB = new vector<vector<int>>(ObsProb);
	pP0 = new vector<int>(InitProb);

	A = *pA;
	B = *pB;
	P0 = *pP0;
}

HMM:: ~HMM()
{
	delete pA;
	delete pB;
	delete pP0;
	delete pbeta;
	delete palpha;
}

// Forward algorithm
real HMM::GetProbabilityF(vector<int>& ObsIdxs)
{
	// 给定长度为Tmax第观测数据，输出该次观测的概率P(O|u)
	real sum = 0.0;
	Tmax = ObsIdxs.size();

	palpha = new vector<vector<int>>();
	palpha->reserve(Tmax);  // 省去了初始化为0第构造函数调用
	alpha = *palpha;  // 使用数组下标访问
	// Initialization
	int t = 0;
	for (int i = 0; i < N; i++)
		alpha[0].push_back(P0[i]*B[i][ObsIdxs[t]]);

	// Recursion
	for (t = 0; t < Tmax-1; t++)
	{
		for (int j = 0; j < N; j++)
		{
			sum = 0.0;
			for (int i = 0; i < N; i++) {
				sum += alpha[t][i] * A[i][j];  // 在t时刻的状态i上求和
			}
			alpha[t+1].push_back(sum * B[j][ObsIdxs[t+1]]);
		};
	}
	// Termination
	t = Tmax-1;
	sum = 0.0;
	for (int i = 0; i < N; i++)
		sum += alpha[t][i];

	return sum;
}

// Backward algorithm
real HMM::GetProbabilityB(vector<int>& ObsIdxs)
{
	real sum = 0.0;
	Tmax = ObsIdxs.size();
	pbeta = new vector<vector<int>>();
	palpha->reserve(Tmax);  // 省去了初始化为0第构造函数调用
	beta = *pbeta;

	// Initialization
	int t = Tmax-1;
	for (int i = 0; i < N; i++)
		beta[t][i] = 1.0;

	// Recursion
	for (t = Tmax - 2; t >= 0; t--)
	{
		for (int i = 1; i <= N; i++)
		{
			sum = 0.0;
			for (int j = 1; j <= N; j++) {
				sum += A[i][j] * B[j][ObsIdxs[t + 1]] * beta[t + 1][j];
			}
			beta[t][i] = sum;
		};
	}

	// Termination
	// B[i][ObsIdxs[t]] : 状态i发射第t种观测值的概率
	t = 0;
	sum = 0.0;
	for (int i = 1; i <= N; i++)
		sum += P0[i] * B[i][ObsIdxs[t]] * beta[t][i];

	return sum;
}

// Viterbi algorithm: 求出给定观测值条件下的最优序列 P(Q|O)
vector<int> HMM::ViterbiStateIdxs(vector<int>& ObsIdxs)
{
	real p, q;
	int t, k, i;
	vector<int> StateIdxs(Tmax);
	vector<vector<int>> delta;
	delta.reserve(Tmax);
	vector<vector<int>> psi;
	psi.reserve(Tmax);

	// Initialization
	t = 0;
	for (int i = 0; i < N; i++)
	{
		delta[t].push_back(P0[i] * B[i][ObsIdxs[t]]);
		psi[t].push_back(0); // Outside limits - not used
	};

	// Recursion
	for (t = 1; t < Tmax; t++)
	{
		for (int j = 0; j < N; j++)
		{
			p = 0.0;
			k = 0; // Outside limits, must be redefined below
			for (int i = 1; i <= N; i++)
			{
				q = delta[t-1][i] * A[i][j];
				if (q > p)
				{
					p = q;
					k = i;
				};
			};
			delta[t].push_back(p * B[j][ObsIdxs[t]]);
			psi[t].push_back(k);
		};
	};

	// Termination
	t = Tmax-1;
	p = 0.0;
	k = 0; // Outside limits, must be redefined below
	for (int i = 0; i < N; i++)
	{
		q = delta[t][i];
		if (q > p)
		{
			p = q;
			k = i;
		};
	};
	StateIdxs[t] = k;  // q* in Rabiner's paper

	// Path (state sequence) backtracking
	for (t = Tmax - 2; t >= 0; t--)
	{
		StateIdxs[t] = psi[t+1][StateIdxs[t+1]];  // t+1时刻使下一个状态为k的最有可能的状态
	};

	return StateIdxs;
}
