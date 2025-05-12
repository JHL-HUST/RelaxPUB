#pragma once
#include "Utility.h"
#include "Timer.h"
using namespace std;
long long nodeCnt, edgeCnt;
long long subMax, subSum, subCnt;
class MSearcher {
public:
    ui n;

    char *matrix;
    ui *commNeis;

    ui *neiInG;
    ui *neiInP;

    ui K;
    ui *bestSol;
    ui bestSz;

    ui *neighbors;
    ui *nonneighbors;

    std::vector<ui> addList;
	ui addListSz=0;
    std::vector<std::pair<ui, ui> > removed;
    ui removeSz=0;

    ui *PC; // union of P and C, where P is at the front
    ui *PC_rid; // reverse ID for PC
    std::queue<ui> Qv;
    std::queue<std::pair<ui,ui> > Qe;
    ui *lvlID;

	ui P_end;
	ui C_end;
	ui lvl;

	std::vector<std::pair<ui,ui> > vp;

	std::vector<std::vector<ui>> disnode;
	std::vector<ui> colored_nodes;
	std::vector<ui> record;
	ui *nbnn;
	//ui *has_include;
	ui *nbnns;
	ui *nbnns2;
	ui *ranked;
	bool *colored;
	bool *conflicted;
	float *profits;

	bool enable=true;
	bool enableInput=true;

public:
    MSearcher() {
        n = 0;
        matrix = NULL;
        commNeis = NULL;

        neiInG = neiInP = NULL;
        
        bestSol = NULL;
        K = bestSz = 0;

        neighbors = nonneighbors = NULL;

        PC = PC_rid = NULL;
        lvlID = NULL;
        removeSz = 0;

		nbnn = NULL;
		//has_include = NULL;
		nbnns = NULL;
		nbnns2 = NULL;
		ranked = NULL;
		colored = NULL;
		conflicted = NULL;
		profits = NULL;
    }

    ~MSearcher() {
        if(matrix != NULL) {
            delete[] matrix;
            matrix = NULL;
        }
        if(commNeis != NULL) {
        	delete[] commNeis;
        	commNeis = NULL;
        }
        if(neiInG != NULL) {
            delete[] neiInG;
            neiInG = NULL;
        }
        if(neiInP != NULL) {
            delete[] neiInP;
            neiInP = NULL;
        }
        if(bestSol != NULL) {
            delete[] bestSol;
            bestSol = NULL;
        }
        if(PC != NULL) {
            delete[] PC;
            PC = NULL;
        }
        if(PC_rid != NULL) {
            delete[] PC_rid;
            PC_rid = NULL;
        }
        if(neighbors != NULL) {
            delete[] neighbors;
            neighbors = NULL;
        }
        if(nonneighbors != NULL) {
            delete[] nonneighbors;
            nonneighbors = NULL;
        }
        if(lvlID != NULL) {
        	delete[] lvlID;
        	lvlID = NULL;
        }
		if(nbnn != NULL) {
        	delete[] nbnn;
        	nbnn = NULL;
        }
		if(nbnns != NULL) {
        	delete[] nbnns;
        	nbnns = NULL;
        }
		if(nbnns2 != NULL) {
        	delete[] nbnns2;
        	nbnns2 = NULL;
        }
		if(ranked != NULL) {
        	delete[] ranked;
        	ranked = NULL;
        }
		if(colored != NULL) {
        	delete[] colored;
        	colored = NULL;
        }
		if(conflicted != NULL) {
        	delete[] conflicted;
        	conflicted = NULL;
        }
		if(profits != NULL) {
        	delete[] profits;
        	profits = NULL;
        }
    }

    void load(ui _n, const std::vector<std::pair<int,int> > &vp) {
		n = _n;
        matrix = new char[n*n];
        commNeis = new ui[n*n];
        
        neiInG = new ui[n];
        neiInP = new ui[n];
        bestSol = new ui[n];
        PC = new ui[n];
        PC_rid = new ui[n];
        neighbors = new ui[n];
        nonneighbors = new ui[n];
        lvlID = new ui[n];
        memset(commNeis, 0, sizeof(ui)*n*n);
        memset(matrix, 0, sizeof(char)*n*n);
        for(ui i = 0; i < n; i++) neiInG[i] = 0;
        for(ui i = 0;i < vp.size();i ++) {
            assert(vp[i].first >= 0&&vp[i].first < n&&vp[i].second >= 0&&vp[i].second < n);
        	ui a = vp[i].first, b = vp[i].second;
            neiInG[a] ++;
            neiInG[b] ++;
            matrix[a*n + b] = matrix[b*n + a] = 1;
        }

		nbnn = new ui[n];
		//has_include = NULL;
		nbnns = new ui[n];
		nbnns2 = new ui[n+1];
		ranked = new ui[n];
		colored = new bool[n];
		memset(colored, false, sizeof(bool) * n);
		conflicted = new bool[n+1];
		profits = new float[n];
    }

    void run(ui K_, std::vector<ui> &res) {
        K = K_;
        bestSz = res.size();
		P_end=C_end=lvl=0;
        init();
        if(PC_rid[0] < C_end){
			if(!moveC2P_R(0)){
				if(C_end>=bestSz){
					enable=enableInput;
					subMax=max(subMax,(long long)(C_end-bestSz));
					subCnt++;
					subSum+=(long long)(C_end-bestSz);
					lvl++,reduce(0, 0),lvl--;
					enable=true;
				}
			}
		}
        if(bestSz > res.size()) {
            res.clear();
            for(int i = 0;i < bestSz;i ++) res.push_back(bestSol[i]);
        }
    }

    void runWHOLE(ui K_, std::vector<ui> &res) {
        K = K_;
        bestSz = res.size();
		P_end=C_end=lvl=0;
        init();
		enable=enableInput;
		subMax=max(subMax,(long long)(C_end-bestSz));
		subCnt++;
		subSum+=(long long)(C_end-bestSz);
		lvl++,reduce(0,0),lvl--;
		enable=true;
        if(bestSz > res.size()) {
            res.clear();
            for(int i = 0;i < bestSz;i ++) res.push_back(bestSol[i]);
        }
    }

private:
    void init() {
    	// the following computes a degeneracy ordering
		while(!Qv.empty()) Qv.pop();
    	while(!Qe.empty()) Qe.pop();
		ui *seq = neighbors;
    	ui *core = nonneighbors;
    	ui *vis = PC;
    	memset(vis, 0, sizeof(ui)*n);
    	ui max_core = 0, UB = 0, idx = n;
    	for(ui i = 0;i < n;i ++) {
			ui u, min_neiInG = n;
			for(ui j = 0;j < n;j ++) if(!vis[j]&&neiInG[j] < min_neiInG) {
				u = j;
				min_neiInG = neiInG[j];
			}
			if(min_neiInG > max_core) max_core = min_neiInG;
			core[u] = max_core;
			seq[i] = u;
			vis[u] = 1;

			ui t_UB = core[u] + K;
			if(n - i < t_UB) t_UB = n - i;
			if(t_UB > UB) UB = t_UB;

			if(idx == n&&min_neiInG + K >= n - i) idx = i;

			for(ui j = 0;j < n;j ++) if(!vis[j]&&matrix[u*n + j]) -- neiInG[j];
		}

    	if(n - idx > bestSz) {
    		bestSz = n - idx;
    		for(ui i = idx;i < n;i ++) bestSol[i-idx] = seq[i];
    		printf("Degen find a solution of size %u\n", bestSz);
    	}

        memset(neiInP, 0, sizeof(int)*n);
        C_end = 0;
        for(ui i = 0;i < n;i ++) PC_rid[i] = n;
        for(ui i = 0;i < n;i ++) if(core[i] + K > bestSz) {
            PC[C_end] = i; PC_rid[i] = C_end;
            ++ C_end;
        }

        if(PC_rid[0] == n) {
        	C_end = 0;
        	return ;
        }

        for(ui i = 0;i < C_end;i ++) {
        	ui u = PC[i];
        	neiInG[u] = 0;
        	for(ui j = 0;j < C_end;j ++) if(matrix[u*n + PC[j]]) ++ neiInG[u];
        }

        for(ui i = 0;i < C_end;i ++) {
        	ui neighbors_n = 0;
        	char *t_matrix = matrix + PC[i]*n;
        	for(ui j = 0;j < C_end;j ++) if(t_matrix[PC[j]]) neighbors[neighbors_n ++] = PC[j];
        	for(ui j = 0;j < neighbors_n;j ++) for(ui k = j+1;k < neighbors_n;k ++) {
        		++ commNeis[neighbors[j]*n + neighbors[k]];
        		++ commNeis[neighbors[k]*n + neighbors[j]];
        	}
        }

        memset(lvlID, 0, sizeof(ui)*n);
        for(ui i = 0;i < C_end;i ++) lvlID[PC[i]] = n;

        for(ui i = 0;i < C_end;i ++) for(ui j = i+1;j < C_end;j ++) {
        	if(matrix[PC[i]*n + PC[j]]&&upper_bound_based_prune(0, PC[i], PC[j])) {
        		Qe.push(std::make_pair(PC[i], PC[j]));
        	}
        }

        removeSz = 0;
        if(CTCP()) C_end = 0;
    }

	bool checkNei(int u) {
		int notSat = 0;
		char *t_matrix = matrix + u*n;
		for(int i = P_end; i < C_end; ++i) {
			int ele = PC[i];
			if(u != ele && !t_matrix[ele] && K+neiInG[ele] < C_end) {
				if(notSat) return false;
				notSat++;
			}
		}
		return true;
	}

	bool checkTriangle(int u) {
		std::vector<int> NotNei;
		char *t_matrix = matrix + u*n;
		for(int i = P_end; i < C_end; ++i) {
			int ele = PC[i];
			if(u != ele && !t_matrix[ele]) NotNei.emplace_back(ele);
		}
		return !matrix[NotNei[0]*n+NotNei[1]];
	}

	void branch(ui &begIdx, ui &endIdx){
		nodeCnt++, edgeCnt++; endIdx=begIdx;
		ui minnei=0x3f3f3f3f; ui pivot;
		char *t_matrix = matrix + 0*n;
		for(ui i = P_end;i < C_end;i ++) {
			ui v = PC[i];
			if(!t_matrix[v] && COND(v)){ //HOP2 first
				if(addListSz == endIdx) {
					addList.push_back(v);
					++ endIdx;
					++ addListSz;
				}
				else addList[endIdx++] = v;
				return;
			}
			if (neiInG[v] < minnei)
			{
				minnei = neiInG[v];
				pivot = v;
			}
		}
		t_matrix = matrix + pivot*n;
		for(ui i = P_end;i < C_end;i ++) if(!t_matrix[PC[i]]) {
			if(addListSz == endIdx) {
				addList.push_back(PC[i]);
				++ endIdx;
				++ addListSz;
			}
			else addList[endIdx++] = PC[i];
		}
		std::sort(addList.data()+begIdx,addList.data()+endIdx,[&](int a,int b){return neiInG[a]>neiInG[b];});
	}

    void reduce(ui begIdx, ui endIdx) {
        if(C_end <= bestSz){
			return;
		}
        if(P_end == C_end) {
            bestSz = P_end;
            for(ui i = 0;i < bestSz;i ++) bestSol[i] = PC[i];
			nodeCnt++;
            return ;
        }
        bool isKPlex = true;
        for(ui i = 0;i < C_end;i ++) if(neiInG[PC[i]] + K < C_end) {
            isKPlex = false;
            break;
        }
        if(isKPlex) {
            bestSz = C_end;
            for(ui i = 0;i < bestSz;i ++) bestSol[i] = PC[i];
			nodeCnt++;
        	return ;
        }
        else if(C_end == bestSz + 1){
			return ;
		}
		branching_num++;
		//reduction for high neiInG vertices
        ui pivot = n;
		for(ui i = P_end;i < C_end;i ++) {
			ui v = PC[i];
			const int NotNei=C_end-neiInG[v];
			if(NotNei <= 2) {
				pivot = v;
				break;
			}
			else if(neiInP[v] == P_end) {
				if(NotNei == 3 && checkTriangle(v)) {
					pivot = v;
					break;
				}else if(NotNei <= K+1 && checkNei(v)){
					pivot = v;
					break;
				}
			}
		}

		//reduction for the whole P
		if(pivot == n){
			//int cursor=bound();
			int cursor=SeesawUB();
			if(cursor>=C_end){
				return;
			}
			if(cursor+1==C_end) pivot=PC[cursor];
		}
		
		//reduction for high neiInG vertices
		if(pivot < n){
       		ui oldRemoveSz = removeSz;
			if(!moveC2P_R(pivot))lvl++, reduce(begIdx, endIdx), lvl--;
       		unR(oldRemoveSz);
			moveP2C();
			return;
		}

		//no more un-reduced lazy branches
		if(begIdx >= endIdx || PC_rid[addList[endIdx-1]] >= C_end || PC_rid[addList[endIdx-1]] < P_end)
			branch(begIdx,endIdx); //branch in lazy way for better reduction of generated branches

		pivot = addList[-- endIdx];

        ui oldBestSz = bestSz;
        ui oldRemoveSz = removeSz;
        // the first branch includes u into S
		
        if(!moveC2P_R(pivot)) edgeCnt++, lvl++, reduce(begIdx, endIdx), lvl--;
        unR(oldRemoveSz);

        while(!Qe.empty()) Qe.pop();
        // the second branch exclude u from S

		if(removed.size() == removeSz) {
			removed.push_back(std::make_pair(-1,pivot));
			++ removeSz;
		}
		else removed[removeSz ++] = std::make_pair(-1,pivot);
		bool pruned = moveP2X_R();
        if(!pruned && bestSz > oldBestSz) pruned = updateR();
        if(!pruned && !CTCP()) lvl++, reduce(endIdx, endIdx), lvl--;
        unR(oldRemoveSz);
    }

    bool moveC2P_R(ui u) {
		swap_pos(P_end, PC_rid[u]); ++ P_end;
        char *t_matrix = matrix + u*n;

        ui neighbors_n = 0, nonneighbors_n = 0;
        for(ui i = 0;i < C_end;i ++) if(i != P_end-1) {
            if(t_matrix[PC[i]]) neighbors[neighbors_n++] = PC[i];
            else nonneighbors[nonneighbors_n++] = PC[i];
        }

        for(ui i = 0;i < neighbors_n;i ++) ++ neiInP[neighbors[i]];

        if(neiInP[u] + K == P_end) { // only neighbors of u in R can be candidates --- RR2
        	ui i = 0;
        	while(i < nonneighbors_n&&PC_rid[nonneighbors[i]] < P_end) ++ i;
            for(;i < nonneighbors_n;i ++) { // remove non-neighbors from R
            	assert(lvlID[nonneighbors[i]] > lvl);
            	lvlID[nonneighbors[i]] = lvl;
                Qv.push(nonneighbors[i]);
            }
        }
        else { // only non-neighbors of u may change their allowance --- RR1
        	ui i = 0;
        	while(i < nonneighbors_n&&PC_rid[nonneighbors[i]] < P_end) ++ i;
            for(;i < nonneighbors_n;i ++) if(P_end - neiInP[nonneighbors[i]] >= K) {
            	assert(lvlID[nonneighbors[i]] > lvl);
            	lvlID[nonneighbors[i]] = lvl;
                Qv.push(nonneighbors[i]);
            }
        }

        // RR2
        for(ui i = 0;i < nonneighbors_n&&PC_rid[nonneighbors[i]] < P_end;i ++) if(neiInP[nonneighbors[i]] + K == P_end) {
            char *tt_matrix = matrix + nonneighbors[i]*n;
            for(ui j = P_end;j < C_end;j ++) if(lvlID[PC[j]] > lvl&&!tt_matrix[PC[j]]) {
            	lvlID[PC[j]] = lvl;
                Qv.push(PC[j]);
            }
        }
if(enable){
        // RR4
        for(ui i = 0;i < nonneighbors_n;i ++) {
            int v = nonneighbors[i];
            assert(!t_matrix[v]);
            if(PC_rid[v] < P_end||lvlID[v] == lvl||t_matrix[v]) continue;
            if(upper_bound_based_prune(P_end, u, v)) {
            	lvlID[v] = lvl;
                Qv.push(v);
            }
        }

        // update cn(.,.)
        for(ui i = 0;i < neighbors_n;i ++) { // process common neighbors of u
            for(ui j = i+1;j < neighbors_n;j ++) {
#ifndef NDEBUG
            	if(!commNeis[neighbors[i]*n + neighbors[j]]) {
            		printf("cn[neighbors[i]*n + neighbors[j]]: %u %u\n", commNeis[neighbors[i]*n + neighbors[j]], commNeis[neighbors[j]*n + neighbors[i]]);
            	}
#endif
                assert(commNeis[neighbors[i]*n + neighbors[j]]);
                -- commNeis[neighbors[i]*n + neighbors[j]];
                -- commNeis[neighbors[j]*n + neighbors[i]];
            }
        }
        int new_n = 0;
        for(ui i = 0;i < nonneighbors_n;i ++) if(lvlID[nonneighbors[i]] > lvl) nonneighbors[new_n ++] = nonneighbors[i];
        nonneighbors_n = new_n;
        for(ui i = 1;i < nonneighbors_n;i ++) { // process common non-neighbors of u
            ui w = nonneighbors[i];
            for(ui j = 0;j < i;j ++) {
            	ui v = nonneighbors[j];
            	if(!upper_bound_based_prune(P_end, v, w)) continue;
            	if(PC_rid[w] < P_end) return true; // v, w \in S --- UB2
            	else if(PC_rid[v] >= P_end) { // v, w, \in R --- RR5
            		if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
            	}
            	else { // RR4
            		assert(lvlID[w] > lvl);
            		lvlID[w] = lvl;
            		Qv.push(w);
            		break;
            	}
            }
        }
}
        return CTCP();
    }

    bool CTCP() {
		if(!enable){
			while(!Qe.empty())Qe.pop();
		}
        while(!Qv.empty()||!Qe.empty()) {
        	while(!Qv.empty()) {
				ui u = Qv.front(); Qv.pop(); // remove u
				assert(PC[PC_rid[u]] == u);
				assert(PC_rid[u] >= P_end&&PC_rid[u] < C_end);
				-- C_end;
				swap_pos(PC_rid[u], C_end);
				bool terminate = false;
				ui neighbors_n = 0;
				char *t_matrix = matrix + u*n;
				for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
					ui w = PC[i];
					neighbors[neighbors_n++] = w;
					-- neiInG[w];
					if(neiInG[w] + K <= bestSz) {
						if(i < P_end) terminate = true; // UB1
						else if(lvlID[w] > lvl) { // RR3
							lvlID[w] = lvl;
							Qv.push(w);
						}
					}
				}
				// UB1
				if(terminate) {
					for(ui i = 0;i < neighbors_n;i ++) ++ neiInG[neighbors[i]];
					lvlID[u] = n;
					++ C_end;
					return true;
				}
				else{
					if(removed.size() == removeSz) {
						removed.push_back(std::make_pair(-1,u));
						++ removeSz;
					}
					else removed[removeSz ++] = std::make_pair(-1,u);
				}

if(enable){
				for(ui i = 1;i < neighbors_n;i ++) {
					ui w = neighbors[i];
					for(ui j = 0;j < i;j ++) {
						ui v = neighbors[j];
						assert(commNeis[v*n+w]);
#ifndef NDEBUG
						ui common_neighbors = 0;
        	for(ui k = P_end;k <= C_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors;
        	assert(commNeis[v*n + w] == common_neighbors);
        	assert(commNeis[w*n + v] == common_neighbors);
#endif
						-- commNeis[v*n + w];
						-- commNeis[w*n + v];
#ifndef NDEBUG
						common_neighbors = 0;
        	for(ui k = P_end;k < C_end;k ++) if(matrix[PC[k]*n + v]&&matrix[PC[k]*n + w]) ++ common_neighbors;
        	assert(commNeis[v*n + w] == common_neighbors);
        	assert(commNeis[w*n + v] == common_neighbors);
#endif

						if(!upper_bound_based_prune(P_end, v, w)) continue;

						if(PC_rid[w] < P_end) terminate = true; // v, w \in S --- UB2
						else if(PC_rid[v] >= P_end) { // v, w, \in R --- RR5
							if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
						}
						else if(lvlID[w] > lvl) { // RR4
							lvlID[w] = lvl;
							Qv.push(w);
						}
					}
				}
				if(terminate) {
					return true;
				}
}
        	}
        	
			if(Qe.empty()) break;

if(enable){
        	ui v = Qe.front().first, w =  Qe.front().second; Qe.pop();
        	if(lvlID[v] <= lvl||lvlID[w] <= lvl||!matrix[v*n + w]) continue;
        	assert(PC_rid[v] >= P_end&&PC_rid[v] < C_end&&PC_rid[w] >= P_end&&PC_rid[w] < C_end);

        	if(neiInG[v] + K <= bestSz + 1) {
        		lvlID[v] = lvl;
        		Qv.push(v);
        	}
        	if(neiInG[w] + K <= bestSz + 1) {
        		lvlID[w] = lvl;
        		Qv.push(w);
        	}
        	if(!Qv.empty()) continue;
        	assert(matrix[v*n + w]);
        	matrix[v*n + w] = matrix[w*n + v] = 0;
        	-- neiInG[v]; -- neiInG[w];

        	if(removed.size() == removeSz) {
        		removed.push_back(std::make_pair(v,w));
        		++ removeSz;
        	}
        	else removed[removeSz ++] = std::make_pair(v,w);

        	char *t_matrix = matrix + v*n;
        	for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
        		-- commNeis[w*n + PC[i]];
        		-- commNeis[PC[i]*n + w];
        		if(!upper_bound_based_prune(P_end, w, PC[i])) continue;
        		if(i < P_end) {
        			if(lvlID[w] > lvl) {
        				lvlID[w] = lvl;
        				Qv.push(w);
        			}
        		}
        		else if(matrix[w*n + PC[i]]) Qe.push(std::make_pair(w, PC[i]));
        	}
        	t_matrix = matrix + w*n;
        	for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
        		-- commNeis[v*n + PC[i]];
        		-- commNeis[PC[i]*n + v];
        		if(!upper_bound_based_prune(P_end, v, PC[i])) continue;
        		if(i < P_end) {
        			if(lvlID[v] > lvl) {
        				lvlID[v] = lvl;
        				Qv.push(v);
        			}
        		}
        		else if(matrix[v*n + PC[i]]) Qe.push(std::make_pair(v, PC[i]));
        	}
}
        }
        return false;
    }

    void unR(ui oldRemoveSz) {
        while(!Qv.empty()) {
            ui u = Qv.front(); Qv.pop();
            assert(lvlID[u] == lvl&&PC_rid[u] < C_end);
            lvlID[u] = n;
        }
		while(removeSz>oldRemoveSz){
        	ui v = removed[--removeSz].first, w = removed[removeSz].second;
			if(v==-1){
				ui u = w; // remove u
				assert(PC[PC_rid[u]] == u && PC_rid[u] == C_end);
				C_end++; lvlID[u] = n;
				ui neighbors_n = 0;
				char *t_matrix = matrix + u*n;
				for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
					ui w = PC[i];
					neighbors[neighbors_n++] = w;
					++ neiInG[w];
				}
				if(enable){
					for(ui i = 1;i < neighbors_n;i ++) {
						ui w = neighbors[i];
						for(ui j = 0;j < i;j ++) {
							ui v = neighbors[j];
							++ commNeis[v*n + w];
							++ commNeis[w*n + v];
						}
					}
				}
			}
			else{
				assert(PC_rid[v] >= P_end&&PC_rid[v] < C_end&&PC_rid[w] >= P_end&&PC_rid[w] < C_end);
				assert(lvlID[v] == n && lvlID[w] == n && !matrix[v*n + w]);
				matrix[v*n + w] = matrix[w*n + v] = 1;
				++ neiInG[v]; ++ neiInG[w];

				char *t_matrix = matrix + v*n;
				for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
					++ commNeis[w*n + PC[i]];
					++ commNeis[PC[i]*n + w];
				}
				t_matrix = matrix + w*n;
				for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) {
					++ commNeis[v*n + PC[i]];
					++ commNeis[PC[i]*n + v];
				}
			}
        }
    }

    void moveP2C() {
    	assert(P_end);
        ui u = PC[-- P_end];
        ui neighbors_n = 0;
        char *t_matrix = matrix + u*n;
        for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) neighbors[neighbors_n ++] = PC[i];
        for(ui i = 0;i < neighbors_n;i ++) -- neiInP[neighbors[i]];

if(enable){
        for(ui i = 0;i < neighbors_n;i ++) {
        	ui v = neighbors[i];
        	for(ui j = i+1;j < neighbors_n;j ++) {
        		ui w = neighbors[j];
        		++ commNeis[v*n + w];
        		++ commNeis[w*n + v];
        	}
        }
}
    }

    bool moveP2X_R() {
    	assert(P_end);
		ui u = PC[P_end-1];
		-- P_end; -- C_end;
		swap_pos(P_end, C_end);
		lvlID[u] = lvl;

		bool ret = false;
        ui neighbors_n = 0;
        char *t_matrix = matrix + u*n;
        for(ui i = 0;i < C_end;i ++) if(t_matrix[PC[i]]) neighbors[neighbors_n ++] = PC[i];
        for(ui i = 0;i < neighbors_n;i ++) {
        	-- neiInP[neighbors[i]];
        	-- neiInG[neighbors[i]];
        	if(neiInG[neighbors[i]] + K <= bestSz) {
        		if(PC_rid[neighbors[i]] < P_end) ret =  true;
        		else {
        			assert(lvlID[neighbors[i]] > lvl);
        			lvlID[neighbors[i]] = lvl;
        			Qv.push(neighbors[i]);
        		}
        	}
        }
        if(ret) return true;

if(enable){
        for(ui i = 1;i < neighbors_n;i ++) if(lvlID[neighbors[i]] > lvl) {
        	ui w = neighbors[i];
        	for(ui j = 0;j < i;j ++) {
        		ui v = neighbors[j];
        		if(!upper_bound_based_prune(P_end, v, w)) continue;

        		if(PC_rid[w] < P_end) return true; // v, w \in S
				else if(PC_rid[v] >= P_end) { // v, w, \in R
					if(matrix[v*n + w]) Qe.push(std::make_pair(v,w));
				}
				else {
					assert(lvlID[w] > lvl);
					lvlID[w] = lvl;
					Qv.push(w);
					break;
				}
        	}
        }
}
		return false;
	}

    bool updateR() {
    	for(ui i = 0;i < P_end;i ++) if(neiInG[PC[i]] + K <= bestSz) return true;

if(enable){
    	for(ui i = 0;i < P_end;i ++) for(ui j = i+1;j < P_end;j ++) {
    		if(upper_bound_based_prune(P_end, PC[i], PC[j])) return true;
    	}
}

    	for(ui i = P_end;i < C_end;i ++) if(lvlID[PC[i]] > lvl){
    		if(P_end - neiInP[PC[i]] >= K||neiInG[PC[i]] + K <= bestSz) {
    			assert(lvlID[PC[i]] > lvl);
    			lvlID[PC[i]] = lvl;
    			Qv.push(PC[i]);
    			continue;
    		}
    		char *t_matrix = matrix + PC[i]*n;
    		for(ui j = 0;j < P_end;j ++) {
    			if((P_end - neiInP[PC[j]] == K&&!t_matrix[PC[j]])||(enable&&upper_bound_based_prune(P_end, PC[i], PC[j])))
    			{
    				assert(lvlID[PC[i]] > lvl);
    				lvlID[PC[i]] = lvl;
    				Qv.push(PC[i]);
    				break;
    			}
    		}
    	}

if(enable){
    	for(ui i = P_end;i < C_end;i ++) if(lvlID[PC[i]] > lvl) {
    		for(ui j = i+1;j < C_end;j ++) if(lvlID[PC[i]] < lvl&&matrix[PC[i]*n + PC[j]]) {
    			if(upper_bound_based_prune(P_end, PC[i], PC[j])) Qe.push(std::make_pair(PC[i], PC[j]));
    		}
    	}
}
        return false;
    }

    int bound() {
    	vp.clear();
    	for(ui i = 0;i < P_end;i ++) vp.push_back(std::make_pair(K-(P_end-neiInP[PC[i]]), PC[i]));
		// for(ui i = 0;i < P_end;i ++) vp.push_back(std::make_pair(-(neiInG[PC[i]]-neiInP[PC[i]]), PC[i]));
    	sort(vp.begin(), vp.end());
    	ui UB = P_end, cursor = P_end;
    	for(ui i = 0;i < (ui)vp.size(); i++) {
    		ui u = vp[i].second;
    		if(vp[i].first == 0) continue;
    		ui count = 0;
    		char *t_matrix = matrix + u*n;
    		for(ui j = cursor;j < C_end;j ++) if(!t_matrix[PC[j]]) {
    			if(j != cursor + count) swap_pos(j, cursor+count);
    			++ count;
    		}
    		ui t_ub = count;
    		if(vp[i].first < t_ub) t_ub = vp[i].first;
    		if(UB + t_ub <= bestSz) {
    			UB += t_ub;
    			cursor += count;
    		}
    		else {
    			return cursor + (bestSz - UB);
    		}
    	}
		cursor+=(bestSz-UB);
		if(cursor>C_end)cursor=C_end;
    	return cursor;
		//[0, cursor - 1] contain at most bestSz k-plex
    }

	void qsort_profits(ui s[], int l, int r)
	{
		if (l < r)
		{
			int i = l, j = r, x = s[l];
			while (i < j)
			{
				while (i < j && profits[s[j]] >= profits[x]) // ���������ҵ�һ��С��x����
					j--;
				if (i < j)
					s[i++] = s[j];
				while (i < j && profits[s[i]] < profits[x]) // ���������ҵ�һ�����ڵ���x����
					i++;
				if (i < j)
					s[j--] = s[i];
			}
			s[i] = x;
			qsort_profits(s, l, i - 1); // �ݹ����
			qsort_profits(s, i + 1, r);
		}
	}

		void qsort_value(ui s[], int l, int r)
	{
		if (l < r)
		{
			int i = l, j = r, x = s[l];
			while (i < j)
			{
				while (i < j && s[j] < x) // ���������ҵ�һ��С��x����
					j--;
				if (i < j)
					s[i++] = s[j];
				while (i < j && s[i] >= x) // ���������ҵ�һ�����ڵ���x����
					i++;
				if (i < j)
					s[j--] = s[i];
			}
			s[i] = x;
			qsort_value(s, l, i - 1); // �ݹ����
			qsort_value(s, i + 1, r);
		}
	}
	
	ui con_num2(ui i, ui color_num)
	{
		ui cnum = 0;
		char *t_matrix = matrix + PC[i] * n;
		for (ui j = 0; j < color_num; j++)
		{
			if (t_matrix[PC[colored_nodes[j]]])
				cnum++;
		}
		return cnum;
	}

	ui con_num(ui i, ui color_num, bool *conflicted, int color_ub)
	{
		record.clear();
		char *t_matrix = matrix + PC[i] * n;
		//ui c_num = 1;
		for (ui j = 0; j < color_num; j++)
		{
			if (nbnns2[j] > color_ub)
				continue;
			if (t_matrix[PC[colored_nodes[j]]] && !conflicted[j])
			{
				//record[c_num - 1] = j;
				record.push_back(j);
				//c_num++;
			}
		}
		return record.size() + 1;
	}

	bool can_color(ui i, ui color_num)
	{
		char *t_matrix = matrix + PC[i] * n;
		for (ui j = 0; j < color_num; j++)
		{
			if (t_matrix[PC[colored_nodes[j]]])
			{
				return false;
			}
		}
		return true;
	}

	void try_color(ui &color_num, ui &color_ub, ui lb, ui P_end, ui C_end, ui *has_include) //Ⱦɫʱbound���ܳ���lb����ζ��ֻ����Ⱦɫ������
	{
		ui large_num = 0;
		colored_nodes.clear();
		//for (ui u = C_end - 1; u >= P_end; u--)
		for (ui u = P_end; u < C_end; u++)
		{
			ui i = ranked[u - P_end];
			if (has_include[i])
				continue;
			ui invest = nbnn[i];
			if (invest > lb)
			{
				large_num++;
				if (large_num > lb)
					continue;
			}
			//	continue;
			if (can_color(i, color_num))
			{
				colored[i] = true;
				nbnns[color_num] = invest;
				//colored_nodes[color_num++] = i;
				colored_nodes.push_back(i);
				color_num++;
			}
		}
		ui color_num_tmp = color_num;
		//color_num = colored_nodes.size();
		if (color_num > 0){
			memset(nbnns2, 0, sizeof(int) * (C_end - P_end + 1));
			memcpy(nbnns2, nbnns, sizeof(int) * color_num);
			qsort_value(nbnns, 0, color_num - 1);
			for (ui i = 0; i < color_num; i++)
			{
				if (nbnns[i] > color_ub)
					color_ub++;
				else
					break;
			}
			
			int can_add = color_ub;
			memset(conflicted, false, sizeof(bool) * (C_end - P_end + 1));
			for (ui i = 0; i < color_num; i++)
			{
				if (nbnns[i] > color_ub)
				{
					//conflicted[i] = true;
					can_add--;
				}
				else
					break;
			}

			if (can_add > 0){
				for (ui u = P_end; u < C_end; u++)
				{
					ui i = ranked[u - P_end];
					if (colored[i] || has_include[i] || nbnn[i] < color_ub)
						continue;

					ui add_num = con_num(i, color_num, conflicted, color_ub);
					if (add_num <= can_add)
					{
						can_add -= add_num;
						for (ui j = 0; j < add_num - 1; j++)
							conflicted[record[j]] = true;
						//insert[insert_num++] = i;
						colored_nodes.push_back(i);
						conflicted[color_num] = true;
						color_num++;
						if (can_add == 0)
							break;
					}
				}
			}

			for (ui u = P_end; u < C_end; u++)
			{
				ui i = ranked[u - P_end];
				if (colored[i] || has_include[i] || nbnn[i] >= color_ub)
					continue;
				ui cnum = con_num2(i, color_num);
				if (cnum <= color_ub - nbnn[i])
				{
					colored_nodes.push_back(i);
					color_num++;
				}
			}
			
			for (ui i = 0; i < color_num_tmp; i++)
				colored[colored_nodes[i]] = false;
		}
	}

	int SeesawUB() {
    	// vp.clear();
    	// for(ui i = 0;i < P_end;i ++) vp.push_back(std::make_pair(K-(P_end-neiInP[PC[i]]), PC[i]));
		// // for(ui i = 0;i < P_end;i ++) vp.push_back(std::make_pair(-(neiInG[PC[i]]-neiInP[PC[i]]), PC[i]));
    	// sort(vp.begin(), vp.end());
    	// ui UB = P_end, cursor = P_end;

		disnode.clear();
		ui *has_include = neighbors;
		memset(has_include, 0, sizeof(ui) * C_end);
		ui out_num = C_end - P_end;
		ui lb = bestSz - P_end;
		ui r_num = out_num;
		for (ui i = 0; i < P_end; i++)
		{
			disnode.push_back(std::vector<ui>());
			char *t_matrix = matrix + PC[i] * n;
			for (ui j = P_end; j < C_end; j++)
				if (!t_matrix[PC[j]])
					disnode[i].push_back(j);
		}
		for (ui i = P_end; i < C_end; i++)
		{
			ranked[i - P_end] = i;
			nbnn[i] = K - (P_end - neiInP[PC[i]]);
			profits[i] = (float) (out_num - neiInG[PC[i]] + neiInP[PC[i]] - 1) / (nbnn[i] - 1 + 0.00000001);
		}
		qsort_profits(ranked, 0, out_num - 1);
		float last_pub_profit = out_num;

		while (lb > 0)
		{
			ui color_num = 0, color_ub = 0;
			try_color(color_num, color_ub, lb, P_end, C_end, has_include);
			float color_profit_rate = (float) (color_num) / color_ub;
			//float color_profit_rate = -1;
			ui saved_invest = 0, saved_count = 0;
			int investor_idx = -1;
			float max_profit_rate = -1;
			if (last_pub_profit > color_profit_rate && 1 < 0)
			{
				for (ui i = 0; i < P_end; i++)
				{
					ui ub = K - (P_end - neiInP[PC[i]]);
					if (ub == 0 || disnode[i].empty())
						continue;
					ui count = 0;
					for (ui j : disnode[i])
						if (has_include[j] == 0)
							count++;
					ui invest = count > ub ? ub : count;
					if (invest == 0)
					{
						disnode[i].clear();
						continue;
					}
					if (invest > lb)
						continue;
					float profit_rate = (float) count / invest;
					if (profit_rate > max_profit_rate || (profit_rate == max_profit_rate && invest > saved_invest))
					{
						saved_invest = invest;
						saved_count = count;
						investor_idx = i;
						max_profit_rate = profit_rate;
					}
				}
				last_pub_profit = max_profit_rate;
			}
			
			if (color_profit_rate > max_profit_rate || (color_profit_rate == max_profit_rate && color_ub > saved_invest))
			{
				color_time++;
				for (ui i = 0; i < color_num; i++)
				{
					has_include[colored_nodes[i]] = 1;
				}
				lb -= color_ub;
				r_num -= (color_num);
				if (r_num <= lb)
					return C_end;
			}
			else
			{
				if (investor_idx < 0)
					break;
				pub_time++;
				for (ui j : disnode[investor_idx])
					has_include[j] = 1;
				disnode[investor_idx].clear();
				lb -= saved_invest;
				r_num -= saved_count;
				if (r_num <= lb)
					return C_end;
			}		
		}

		ui brunching_idx = C_end;
		r_num -= lb;
		for (ui i = 0; i < out_num; i++)
		{
			ui u = ranked[out_num - 1 - i];
			//assert(u >= P_end && u < C_end);
			// wish put u to tail
			if (has_include[u] == 0)
			{
				brunching_idx--;
				has_include[u] = 1;
				if (u != brunching_idx)
				{
					swap_pos(u, brunching_idx);
					i--;
				}
				if (brunching_idx == C_end - r_num)
					break;
			}
		}
		return brunching_idx;
    	// for(ui i = 0;i < (ui)vp.size(); i++) {
    	// 	ui u = vp[i].second;
    	// 	if(vp[i].first == 0) continue;
    	// 	ui count = 0;
    	// 	char *t_matrix = matrix + u*n;
    	// 	for(ui j = cursor;j < C_end;j ++) if(!t_matrix[PC[j]]) {
    	// 		if(j != cursor + count) swap_pos(j, cursor+count);
    	// 		++ count;
    	// 	}
    	// 	ui t_ub = count;
    	// 	if(vp[i].first < t_ub) t_ub = vp[i].first;
    	// 	if(UB + t_ub <= bestSz) {
    	// 		UB += t_ub;
    	// 		cursor += count;
    	// 	}
    	// 	else {
    	// 		return cursor + (bestSz - UB);
    	// 	}
    	// }
		// cursor+=(bestSz-UB);
		// if(cursor>C_end)cursor=C_end;
    	// return cursor;
		
		//[0, cursor - 1] contain at most bestSz k-plex
    }

    bool upper_bound_based_prune(ui P_end, ui u, ui v) {
    	// ui ub = P_end + 2*K - (P_end - neiInP[u]) - (P_end - neiInP[v]) + cn[u*n + v];
    	ui ub = 2*K + neiInP[u] - P_end + neiInP[v] + commNeis[u*n + v];
    	if(PC_rid[u] >= P_end) {
    		-- ub; // P_end ++
    		if(matrix[u*n+v]) ++ ub; // neiInP[v] ++
    	}
    	if(PC_rid[v] >= P_end) {
    		-- ub;
    		if(matrix[v*n+u]) ++ ub;
    	}
    	return ub <= bestSz;
    }

    void swap_pos(ui i, ui j) {
        std::swap(PC[i], PC[j]);
        PC_rid[PC[i]] = i;
        PC_rid[PC[j]] = j;
    }

	bool COND(ui v){
		return (neiInG[v]+K<=bestSz+1)||(neiInP[v]+K==P_end+1)||(neiInP[0]+K==P_end+1);
	}
};
