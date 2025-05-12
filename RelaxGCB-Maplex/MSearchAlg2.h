#pragma once
#include <algorithm>
#include <assert.h>
#include "MGraph.h"
#include "LinearHeap.h"
using namespace std;

class FastExactSearch2{
public:
	const MCsrGraph &cg;
    MBitGraph *bg;
	MBitGraph *rbg;
    int K;
    int lb;
	int maxsec;
	int *core;
	ui *seq;
	int* triangles;

	ui* stkvtx;	//Ӧ����Candidate�еĵ�
	ui* stkvtx_tmp;
	// int* stkdegree;
	int* stkcon2p;
	int* stkbound;
	int* stkpcolor;

	// int** degree;
	int** con2p; 
	int** bound;
	int** pcolors;
	
	int *plexpos;
	int *candpos;
	int depth;

	int *nbnn;
	int *nbnns;
	int *nbnns2;
	int *ranked;
	float *profits;
	// bool *has_include;
	bool *has_include2;
	bool *colored;
	bool *conflicted;
	int *maxcolors;
	int* stkcbound;
	int **cbound;
	//vector<ui> insert;
	vector<ui> colored_nodes;
	vector<ui> record;
	// vector<ui> removed;

	clock_t startclk;
	clock_t endclk;
	uint64_t nnodes;
	
	int ub;
    ui *sol;
    int szsol;
	int isoptimal;

	typedef struct {
		MBitSet64 bitcolor;
		ui *vertices;
		int sz;
	}ColorType;
	ColorType *colors;
	int maxcolor;
#define PLSTART(d)	(plexpos[d])
#define PLEND(d)	(candpos[d])
#define CASTART(d) (candpos[d])
#define CAEND(d)	(plexpos[d+1])	
#define SZ_PLEX(d) (PLEND(d) - PLSTART(d))	//number of vertices in P
#define SZ_CAND(d)	(CAEND(d) - CASTART(d)) //number of vertices in C
#define SZ_VTX(d) (CAEND(d) - PLSTART(d)) // number of verices in P union C

    FastExactSearch2(const MCsrGraph &csrg, int valueK, int maxsec, int lb=0):
	cg(csrg){
    
    printf("%d %d\n",csrg.nbedges, csrg.nbvtx);
    //for (int ii = 0; ii < csrg.nbedges; ii++){
    //  if (csrg.edges[ii] == 0)
    //    printf("!!!!!!!!\n");
    //}
    /*
    for (int ii = 0; ii < csrg.nbvtx; ii++)
      printf("%d ",csrg.pstart[ii]);
    printf("\n\n");
    for (int ii = 0; ii < csrg.nbedges; ii++)
      printf("%d ",csrg.edges[ii]);
    printf("\n");
    */
		bg = new MBitGraph(csrg);
		rbg = new MBitGraph(csrg);
		rbg->reverseGraph();
		this->K = valueK;
		this->lb = lb;
		this->maxsec = maxsec;
		this->nnodes = 0UL;
		this->sol =  new ui[bg->nbvtx];
		this->szsol = lb;		
		cout << lb <<endl;

		//Initialize bound
		core = new int[bg->nbvtx];
		seq = new ui[bg->nbvtx];

		int* vpos = new int[bg->nbvtx];
		//triangles = new int[cg.nbedges];
		int maxcore = coreDecomposition(cg, core, seq, vpos);
		//countingTriangles(cg, vpos, triangles);
		delete[] vpos;
		
		ub =  maxcore+K + 1;
		nnodes = 0;
		isoptimal = 1;

		stkvtx = new ui[csrg.nbvtx * ub];
		stkvtx_tmp = new ui[csrg.nbvtx * ub];
		stkcon2p = new int[csrg.nbvtx * ub];
		stkbound = new int[csrg.nbvtx * ub];
		stkpcolor = new int[(csrg.nbvtx + 1) * ub];
		// stkdegree = new int[csrg.nbvtx * ub];
		stkcbound = new int[csrg.nbvtx * ub];
		
		con2p = new int*[ub];
		bound = new int*[ub];
		pcolors = new int*[ub];
		// degree = new int*[ub];
		cbound = new int*[ub];
		for (int d = 0; d < ub; d++) {
			con2p[d] = &(stkcon2p[csrg.nbvtx * d]);
			bound[d] = &(stkbound[csrg.nbvtx * d]);
			pcolors[d] = &(stkpcolor[(csrg.nbvtx + 1) * d]);
			// degree[d] = &(stkdegree[csrg.nbvtx * d]);
			cbound[d] = &(stkcbound[csrg.nbvtx * d]);
		}

		plexpos = new int[ub + 1];
		candpos = new int[ub + 1];
		//Init color
		colors = new ColorType[bg->nbvtx+1];
		for (int i = 0; i < bg->nbvtx + 1; i++) {
			colors[i].bitcolor = MBitSet64(bg->nbvtx);
			colors[i].vertices = new ui[bg->nbvtx];
			colors[i].sz = 0;			
		}		
		pcolors[0][0] = 0;
		
		nbnn = new int[csrg.nbvtx];
		nbnns = new int[csrg.nbvtx+1];
		nbnns2 = new int[csrg.nbvtx+1];
		ranked = new int[csrg.nbvtx];
		profits = new float[csrg.nbvtx];
		// has_include = new bool[csrg.nbvtx];
		has_include2 = new bool[csrg.nbvtx];
		colored = new bool[csrg.nbvtx];
		memset(colored, false, sizeof(bool) * csrg.nbvtx);
		conflicted = new bool[csrg.nbvtx+1];
		maxcolors = new int[ub];

		// nbnn = new int[csrg.nbvtx * ub];
		// nbnns = new int[csrg.nbvtx * ub];
		// ranked = new int[csrg.nbvtx * ub];
		// profits = new float[csrg.nbvtx * ub];
		// has_include = new int[csrg.nbvtx * ub];
		// colored = new bool[csrg.nbvtx * ub];
		// memset(colored, false, sizeof(bool) * csrg.nbvtx * ub);
		// conflicted = new bool[csrg.nbvtx * ub];

		depth = 0;
		PLSTART(depth) = 0;
		PLEND(depth) = 0;
		CASTART(depth) = 0;
		int pos = CASTART(depth);
		for (int i = 0; i < bg->nbvtx; i++) {
			ui u = seq[bg->nbvtx - i - 1];		
			stkvtx[pos++] = u;
			bound[depth][u] = core[u] + K;
			//potent[depth][u] = core[u];
			con2p[depth][u] = 0;	
			// degree[depth][u] = cg.degree(u);
		}
		CAEND(depth) = pos;
		PLSTART(depth+1) = CAEND(depth);
		

		//candSortByDescPortential();
		#ifdef DEBUG
			showStack();
		#endif
	}

	#define OP_BIT_THRESH(n) (n >> 6)
	inline int intsecColor(ui u, int c) {
		return colors[c].bitcolor.intersect(*(bg->rows[u]));
	}
	
	void bitColorSort() {
		int j = 0; 
		maxcolor = 1;
		//int threshold = szSolution - szPlex;
		colors[1].sz = 0;
		colors[1].bitcolor.clear();
		colors[2].sz = 0;
		colors[2].bitcolor.clear();
		for (int i = CASTART(depth); i < CAEND(depth); i++) {
			int u = stkvtx[i];
			int c = 1;
			int sum = 0;
			while (1) {
				int con = intsecColor(u, c);
				if (con == 0) {
					break;
				}else
					c++;
			}
			colors[c].bitcolor.set(u);
			colors[c].vertices[colors[c].sz] = u;
			colors[c].sz++;		
			//Open the next color to ensure correctness
			if (c > maxcolor) {
				maxcolor = c;
				colors[maxcolor + 1].bitcolor.clear();
				colors[maxcolor + 1].sz = 0;
			}
		}
		int pos = CASTART(depth);
		pcolors[depth][0] = pos;
		int accbound = 0;
		for (int c = 1; c <= maxcolor; c++) {
			for (int i = 0; i < colors[c].sz; i++) {
				ui u = colors[c].vertices[i];
				stkvtx[pos++] = u;
				bound[depth][u] = accbound + min(i + 1, K);		//bound[depth][u]ָ�ľ��ǵ�u������bound		
			}
			accbound += min(colors[c].sz, K);
			pcolors[depth][c] = pos;
		}
		assert(pos == CAEND(depth));
		pcolors[depth][maxcolor + 1] = cg.nbvtx + 1; //maximum color
	}

	void qsort_profits(int s[], int l, int r)
	{
		if (l < r)
		{
			int i = l, j = r, x = s[l];
			while (i < j)
			{
				while (i < j && profits[s[j]] >= profits[x]) // ?????������????????????x????
					j--;
				if (i < j)
					s[i++] = s[j];
				while (i < j && profits[s[i]] < profits[x]) // ??����?��???????????��??????x????
					i++;
				if (i < j)
					s[j--] = s[i];
			}
			s[i] = x;
			qsort_profits(s, l, i - 1); // ???��?��??
			qsort_profits(s, i + 1, r);
		}
	}
	
	void qsort_value(int s[], int l, int r)
	{
		if (l < r)
		{
			int i = l, j = r, x = s[l];
			while (i < j)
			{
				while (i < j && s[j] < x) // ??????????????��??x????
					j--;
				if (i < j)
					s[i++] = s[j];
				while (i < j && s[i] >= x) // ?????????????????????x????
					i++;
				if (i < j)
					s[j--] = s[i];
			}
			s[i] = x;
			qsort_value(s, l, i - 1); // ??????
			qsort_value(s, i + 1, r);
		}
	}

	ui con_num(ui v, ui color_num, bool *conflicted, int c, int color_ub) {
		record.clear();
		//char *t_matrix = matrix + SR[i] * n;
		//ui c_num = 1;
		for (ui j = 0; j < color_num; j++)
		{
			if (nbnns2[j] > color_ub)
				continue;
			ui u = colors[c].vertices[j];
			if (bg->rows[v]->test(u) && !conflicted[j])
			{
				//record[c_num - 1] = j;
				record.push_back(j);
				//c_num++;
			}
		}
		return record.size() + 1;
	}

	ui con_num2(ui v, ui color_num, int c) {
		//record.clear();
		//char *t_matrix = matrix + SR[i] * n;
		ui cnum = 0;
		for (ui j = 0; j < color_num; j++)
		{
			ui u = colors[c].vertices[j];
			if (bg->rows[v]->test(u))
				cnum++;
		}
		return cnum;
	}

	void try_color(ui &color_ub, ui lb, int c) {
		ui large_num = 0;
		colored_nodes.clear();
		colors[c].sz = 0;
		colors[c].bitcolor.clear();
		
		for (int u = CAEND(depth) - 1; u >= CASTART(depth); u--) {
			//ui i = ranked[u - CASTART(depth)];
			ui i = ranked[CAEND(depth) - 1 - u];
			//ui i = CAEND(depth) - 1 - u;
			ui v = stkvtx[i + CASTART(depth)];
			if (has_include2[i])
				continue;
			ui invest = nbnn[i];
			if (invest > lb) {
				large_num++;
				if (large_num > lb)
					continue;
			}
			if (intsecColor(v, c) == 0) {
				colored[i] = true;
				nbnns[colors[c].sz] = invest;
				colors[c].bitcolor.set(v);
				colors[c].vertices[colors[c].sz] = v;
				colors[c].sz++;	
				colored_nodes.push_back(i);	
				// pcolors[depth][v] = c;
			}
		}
		ui color_num_tmp = colors[c].sz;
		if (colors[c].sz > 0) {
			memset(nbnns2, 0, sizeof(int) * (CAEND(depth) - CASTART(depth) + 1));
			memcpy(nbnns2, nbnns, sizeof(int) * colors[c].sz);
			qsort_value(nbnns, 0, colors[c].sz - 1);
			for (ui j = 0; j < colors[c].sz; j++)
			{
				if (nbnns[j] > color_ub)
					color_ub++;
				else
					break;
			}
			
			memset(conflicted, false, sizeof(bool) * (CAEND(depth) - CASTART(depth) + 1));
			int can_add = color_ub;
			for (ui j = 0; j < colors[c].sz; j++)
			{
				if (nbnns[j] > color_ub){
					//conflicted[j] = true;
					can_add--;
				}
				else
					break;
			}

			if (can_add > 0){
				for (int u = CASTART(depth); u < CAEND(depth); u++)
				{
					ui i = ranked[u - CASTART(depth)];
					ui v = stkvtx[i + CASTART(depth)];
					if (colored[i] || has_include2[i] || nbnn[i] < color_ub)
						continue;

					ui add_num = con_num(v, colors[c].sz, conflicted, c, color_ub);
					if (add_num <= can_add)
					{
						can_add -= add_num;
						for (ui j = 0; j < add_num - 1; j++)
							conflicted[record[j]] = true;
						//insert[insert_num++] = i;
						//insert.push_back(i);
						colors[c].bitcolor.set(v);
						colors[c].vertices[colors[c].sz] = v;
						conflicted[colors[c].sz] = true;
						colors[c].sz++;	
						colored_nodes.push_back(i);	
						// pcolors[depth][v] = c;
						if (can_add == 0)
							break;
					}
				}
			}

			for (int u = CASTART(depth); u < CAEND(depth); u++)
			{
				ui i = ranked[u - CASTART(depth)];
				ui v = stkvtx[i + CASTART(depth)];
				if (colored[i] || has_include2[i] || nbnn[i] >= color_ub)
					continue;
				//ui cnum = con_num2(v, colors[c].sz, c);
				ui cnum = intsecColor(v, c);
				if (cnum <= color_ub - nbnn[i])
				{
					colors[c].bitcolor.set(v);
					colors[c].vertices[colors[c].sz] = v;
					colors[c].sz++;	
					colored_nodes.push_back(i);	
					// pcolors[depth][v] = c;
				}
			}

			for (ui i = 0; i < color_num_tmp; i++)
				colored[colored_nodes[i]] = false;
		}
		// if (colors[c].sz == color_ub) {
		// 	color_ub = 0;
		// 	colors[c].sz = 0;
		// }
	}
	/*
	void bitColorSort2() {
		ui lb = szsol - SZ_PLEX(depth);
		ui out_num = CAEND(depth) - CASTART(depth);
		ui r_num = out_num;
		memset(has_include, false, sizeof(bool) * out_num);
		memset(has_include2, false, sizeof(bool) * out_num);
		for (int i = CASTART(depth); i < CAEND(depth); i++) {
			ui v = stkvtx[i];
			ranked[i - CASTART(depth)] = i - CASTART(depth);
			nbnn[i - CASTART(depth)] = K - (SZ_PLEX(depth) - con2p[depth][v]);
			profits[i - CASTART(depth)] = (float) (out_num - degree[depth][v] + con2p[depth][v]) / nbnn[i - CASTART(depth)];
		}
		qsort_profits(ranked, 0, out_num - 1);
		
		int c = 0;
		while (lb > 0) {
			ui color_ub = 0, insert_num = 0;
			c++;
			try_color(color_ub, insert_num, lb, c);
			for (ui i = 0; i < colors[c].sz; i++) {
				has_include2[colored_nodes[i]] = true;
				//has_include2[colors[c].vertices[i]] = true;
			}
			if (color_ub == colors[c].sz)
				break;
			for (ui i = 0; i < colors[c].sz; i++) {
				has_include[colored_nodes[i]] = true;
				//has_include[colors[c].vertices[i]] = true;
			}
			lb -= color_ub;
			r_num -= colors[c].sz;
			if (r_num <= lb) {
				bound[depth][stkvtx[CAEND(depth) - 1]] = 0;
				return;
			}
		}

		int include_color = lb == 0 ? c : c - 1;
		ui r_num2 = lb == 0 ? r_num : r_num - colors[c].sz;

		// while (r_num2 > 0) {
		// 	ui color_ub = 0, insert_num = 0;
		// 	c++;
		// 	try_color(color_ub, insert_num, out_num, c);
		// 	for (ui i = 0; i < colors[c].sz; i++) {
		// 		has_include2[colored_nodes[i]] = true;
		// 		//has_include2[colors[c].vertices[i]] = true;
		// 	}
		// 	r_num2 -= colors[c].sz;
		// }

		//order colored nodes

		int pos2 = CAEND(depth);
		r_num -= lb;
		for (ui i = 0; i < out_num; i++) {
			ui u = ranked[out_num - 1 - i];
			//ui v = stkvtx[u + CASTART(depth)];
			if (has_include[u] == false) {
				pos2--;
				if (u + CASTART(depth) < pos2)
					swap(stkvtx[u + CASTART(depth)], stkvtx[pos2]);
				bound[depth][stkvtx[pos2]] = r_num + szsol - SZ_PLEX(depth) - (CAEND(depth) - 1 - pos2);
				if (pos2 == CAEND(depth) - r_num)
					break;
			}
		}

		int pos = CASTART(depth);
		pcolors[depth][0] = pos;
		for (int i = 1; i <= include_color; i++) {
			for (ui j = 0; j < colors[i].sz; j++) {
				ui u = colors[i].vertices[j];
				stkvtx[pos++] = u;
				bound[depth][u] = 0;
			}
			pcolors[depth][i] = pos;
		}

		cout << pos << " " << pos2 << "  ";

		for (int i = pos + 1; i < CAEND(depth); i++) {
			pcolors[depth][++include_color] = i;
		}
		pcolors[depth][++include_color] = cg.nbvtx + 1;

		// int pos = CASTART(depth);
		// pcolors[depth][0] = pos;
		// int accbound = 0;
		// for (int c = 1; c <= maxcolor; c++) {
		// 	for (int i = 0; i < colors[c].sz; i++) {
		// 		ui u = colors[c].vertices[i];
		// 		stkvtx[pos++] = u;
		// 		bound[depth][u] = accbound + min(i + 1, K);		//bound[depth][u]ָ�ľ��ǵ�u������bound		
		// 	}
		// 	accbound += min(colors[c].sz, K);
		// 	pcolors[depth][c] = pos;
		// }
		// assert(pos == CAEND(depth));
		// pcolors[depth][maxcolor + 1] = cg.nbvtx + 1; //maximum color

		//bound last should be r_num - lb + szsol - SZ_PLEX(depth)

		// int pos = CAEND(depth);
		// r_num -= lb;
		// for (ui i = 0; i < out_num; i++) {
		// 	ui u = ranked[out_num - 1 - i];
		// 	//ui v = stkvtx[u + CASTART(depth)];
		// 	if (has_include[u] == false) {
		// 		pos--;
		// 		if (u + CASTART(depth) < pos)
		// 			swap(stkvtx[u + CASTART(depth)], stkvtx[pos]);
		// 		bound[depth][stkvtx[pos]] = r_num + szsol - SZ_PLEX(depth) - (CAEND(depth) - 1 - pos);
		// 		if (pos == CAEND(depth) - r_num)
		// 			break;
		// 	}
		// }
		// for (ui i = CASTART(depth); i < pos; i++)
		// 	bound[depth][stkvtx[i]] = 0;
	}
	*/

	void bitColorSort3() {
		ui lb = szsol - SZ_PLEX(depth);
		ui out_num = CAEND(depth) - CASTART(depth);
		ui r_num = out_num;
		//memset(has_include, false, sizeof(bool) * out_num);
		memset(has_include2, false, sizeof(bool) * out_num);
		for (int i = CASTART(depth); i < CAEND(depth); i++) {
			ui v = stkvtx[i];
			ranked[i - CASTART(depth)] = i - CASTART(depth);
			nbnn[i - CASTART(depth)] = K - (SZ_PLEX(depth) - con2p[depth][v]);
			// profits[i - CASTART(depth)] = (float) (out_num - degree[depth][v] + con2p[depth][v]) / (nbnn[i - CASTART(depth)] - 1 + 0.00000001);
			profits[i - CASTART(depth)] = (float) (cg.nbvtx - cg.degree(v) + con2p[depth][v]) / (nbnn[i - CASTART(depth)] - 1 + 0.00000001);
			// profits[i - CASTART(depth)] = (float) (out_num) / (nbnn[i - CASTART(depth)] - 1 + 0.00000001);
		}
		qsort_profits(ranked, 0, out_num - 1);
		
		int c = 0;
		ui color_ub;
		int pos = CASTART(depth);
		pcolors[depth][0] = pos;
		while (lb > 0) {
			color_ub = 0;
			c++;
			try_color(color_ub, lb, c);
			cbound[depth][c] = color_ub;
			for (ui i = 0; i < colors[c].sz; i++)
				has_include2[colored_nodes[i]] = true;
			// for (ui i = 0; i < insert_num; i++)
			// 	has_include2[insert[i]] = true;
			// if (color_ub == (colors[c].sz + insert_num))
			// 	break;
			// for (ui i = 0; i < colors[c].sz; i++)
			// 	has_include[colored_nodes[i]] = true;
			// for (ui i = 0; i < insert_num; i++)
			// 	has_include[insert[i]] = true;

			for (int j = 0; j < colors[c].sz; j++) {
				ui u = colors[c].vertices[j];
				stkvtx_tmp[pos++] = u;
				bound[depth][u] = 0;		//bound[depth][u]ָ�ľ��ǵ�u������bound		
			}
			pcolors[depth][c] = pos;

			lb -= color_ub;
			r_num -= (colors[c].sz);
			if (r_num <= lb) {
				// for (int i = CASTART(depth); i < CAEND(depth); i++)
				// 	bound[depth][stkvtx[i]] = 0;
				bound[depth][stkvtx[CAEND(depth) - 1]] = 0;
				return;
			}
		}
		// memcpy(has_include,has_include2,sizeof(bool) * out_num);
		// int include_color = lb == 0 ? c : c - 1;
		ui r_num2 = lb == 0 ? r_num : r_num - colors[c].sz;

		int accbound = szsol - SZ_PLEX(depth);
		while (r_num2 > 0) {
			color_ub = 0;
			c++;
			try_color(color_ub, out_num, c);
			cbound[depth][c] = color_ub;
			for (ui i = 0; i < colors[c].sz; i++)
				has_include2[colored_nodes[i]] = true;
			// for (ui i = 0; i < insert_num; i++)
			// 	has_include2[insert[i]] = true;
			r_num2 -= (colors[c].sz);

			for (int i = 0; i < colors[c].sz; i++) {
				ui u = colors[c].vertices[i];
				stkvtx_tmp[pos++] = u;
				bound[depth][u] = i+1 < color_ub ? accbound + i+1 : accbound + color_ub;  //bound[depth][u]ָ�ľ��ǵ�u������bound		
			}
			accbound += color_ub;
			pcolors[depth][c] = pos;
		}
		maxcolors[depth] = c;

		assert(pos == CAEND(depth));
		pcolors[depth][c + 1] = cg.nbvtx + 1; //maximum color
	
		// if (r_num2 != 0)
		// 	cout << "r_num2: " << r_num2 << endl;

		for (int i = CASTART(depth); i < CAEND(depth); i++)
			stkvtx[i] = stkvtx_tmp[i];
		//cout << "A: ";
		// int num = 0;
		// memset(has_include2, false, sizeof(bool) * out_num);
		// for (ui i = 0; i < out_num; i++) {
		// 	ui u = ranked[out_num - 1 - i];
		// 	//ui v = stkvtx[u + CASTART(depth)];
		// 	if (has_include[u] == false) {
		// 		has_include2[u] = true;
		// 		num++;
		// 		//cout << stkvtx[u + CASTART(depth)] << " ";
		// 		if (num == r_num)
		// 			break;
		// 	}
		// }
		//cout << endl;

		//cout << "B: ";

		//bound last should be r_num - lb + szsol - SZ_PLEX(depth)
		
		// for (int i = CAEND(depth) - 1; i >= CASTART(depth); i--) {
		// 	// ui v = stkvtx[i];
		// 	// ranked[i - CASTART(depth)] = i - CASTART(depth);
		// 	ui u = stkvtx[ranked[i - CASTART(depth)] + CASTART(depth)];
		// 	if (has_include[ranked[i - CASTART(depth)]] == false) {
		// 		pos--;
		// 		stkvtx[pos] = u;
		// 		bound[depth][u] = 
		// 	}
		// }

		
		

		// for (ui i = 0; i < out_num; i++) {
		// 	//ui u = ranked[out_num - 1 - i];
		// 	//ui v = stkvtx[u + CASTART(depth)];
		// 	ui u = out_num - 1 - i;
		// 	if (has_include[u] == false) {
		// 		pos--;
		// 		if (u + CASTART(depth) < pos){
		// 			swap(stkvtx[u + CASTART(depth)], stkvtx[pos]);
		// 			// i--;
		// 		}
					
		// 		bound[depth][stkvtx[pos]] = r_num + szsol - SZ_PLEX(depth) - (CAEND(depth) - 1 - pos);
		// 		if (pos == CAEND(depth) - r_num)
		// 			break;
		// 	}
		// }

		// // for (ui i = pos; i < CAEND(depth); i++)
		// // 	cout << stkvtx[i] << " ";
		// // cout << endl;

		// for (ui i = CASTART(depth); i < pos; i++)
		// 	bound[depth][stkvtx[i]] = 0;
	}
	

	void showcolor(int vpos){
		printf("vertex: %d\n", stkvtx[vpos]);
		printf("Pos: \n" );
		for (int i = CASTART(depth); i < vpos; i++){
			printf("%d\t",i);
		}
		printf("cand [%d]: \n", vpos - CASTART(depth));
		for (int i = CASTART(depth); i < vpos; i++){
			printf("%d\t",stkvtx[i]);
		}
		printf("\n");
		printf("connection: \n");
		for (int i = CASTART(depth);i < vpos; i++){
			if (bg->rows[stkvtx[vpos]]->test(stkvtx[i]))
				printf("1\t");
			else printf("0\t");
		}
		printf("\n");
		printf("colors start: %d ", CASTART(depth));
		for (int c = 1; pcolors[depth][c] != cg.nbvtx + 1; c++){
			printf("C%d: %d | ", c, pcolors[depth][c]);
		}
		printf("\n");

	}
	int lookahead(int vpos){
		int sum = 0;
		int c = 1;
		int p = CASTART(depth);
		int intesec = 0;		
		while (pcolors[depth][c] <= vpos){			
			if (bg->rows[stkvtx[vpos]]->test(stkvtx[p])){
				intesec++;
			}
			if (intesec == K){				
				sum += K;
				p = pcolors[depth][c++]; // jump to next color
				intesec = 0;
			}else if(p == pcolors[depth][c] - 1){
				sum += intesec;
				++c, ++p;
				intesec = 0;
			}else{
				p++;
			}
			
		}
		return sum;
	}
	
	int lookahead2(int vpos){
		int sum = 0;
		int c = 1;
		int p = CASTART(depth);
		int intesec = 0;		
		// while (pcolors[depth][c] <= vpos){	
		while (p < vpos) {		
			if (bg->rows[stkvtx[vpos]]->test(stkvtx[p])){
				intesec++;
			}
			if (intesec == cbound[depth][c]){				
				sum += cbound[depth][c];
				p = pcolors[depth][c++]; // jump to next color
				intesec = 0;
			}else if(p == pcolors[depth][c] - 1){
				sum += intesec;
				++c, ++p;
				intesec = 0;
			}else{
				p++;
			}
			
		}
		sum += intesec;
		return sum;
	}

	int lookahead3(int vpos){
		int sum = 0;
		int intesec_num[maxcolors[depth]];
		for (int i = 1; i <= maxcolors[depth]; i++)
			intesec_num[i] = 0;
		for (int i = CASTART(depth); i < vpos; i++) {
			if (bg->rows[stkvtx[vpos]]->test(stkvtx[i])) {
				int c = pcolors[depth][stkvtx[i]];
				intesec_num[c]++;
			}
		}
		for (int i = 1; i <= maxcolors[depth]; i++)
			sum += min(intesec_num[i], cbound[depth][i]);
		
		return sum;
	}

	int interrupt() {
		if ((double)(clock() - startclk) / CLOCKS_PER_SEC > maxsec) {
			return 1;
		}
		return 0;
	}

	void refineCandidate(){
		if (SZ_PLEX(depth) == 0) return;
		MBitSet64 mask(cg.nbvtx);
		for (int i = PLSTART(depth); i < CAEND(depth); i++){
			mask.set(stkvtx[i]);
		}
		ui lastu = stkvtx[PLEND(depth)-1];
		mask &= *(bg->rows[lastu]);
		int pos = CASTART(depth);
		
		for (int i = CASTART(depth); i < CAEND(depth); i++){
			ui v = stkvtx[i];
			int com = bg->rows[v]->intersect(mask);
			if (mask.test(v) && com + 2* K > szsol){
				stkvtx[pos++] = v; // keep v
			}else if (!mask.test(v) && com + 2*K -2 > szsol){
				stkvtx[pos++] = v; //keep v
			}else{
				if (mask.test(v)) mask.set(v);	
				// removed.push_back(v);		
			}
		}
		CAEND(depth) = pos;
	}

	/** 
	 * move the last vertex in stkcand to plex
	 * increase the stack level
	 * update stkvtx, stkcon2p, stkdeg.
	 * */
	void selectLast() {
		//We only consider the last candidate vertex under current config	
		ui u = stkvtx[CAEND(depth) - 1];
		MBitSet64 gmark(cg.nbvtx);		
		//update the best known		
		set<ui> satu;
		int pos = PLSTART(depth+1);		
	#ifdef DEBUG
		printf("--------------Take %d and deepen ---------------\n", u);
	#endif
		//Update the saturated vertices in the growing Plex
		//stkcon2p[depth+1] = stkcon2p[depth]; 
		for (ui i = PLSTART(depth); i < PLEND(depth); i++){
			ui v = stkvtx[i];
			int con = con2p[depth][v];
			if (bg->rows[v]->test(u)){
				con++;
			}			
			if (con == depth + 1 - K){
				satu.insert(v); // satu vertex
			}
			con2p[depth+1][v] = con;
			stkvtx[pos++] = v;			
			gmark.set(v);
		}
		if (con2p[depth][u] == depth + 1 - K)
			satu.insert(u);
		stkvtx[pos++] = u;		
		//gmark.set(u);
		con2p[depth+1][u] = con2p[depth][u];
		
		PLSTART(depth + 1) = CAEND(depth);
		PLEND(depth + 1) =  pos;
		assert(depth + 1 == SZ_PLEX(depth+1));
		assert(CASTART(depth + 1) == PLEND(depth+1));

		//check every candidate vertex except the last one.
		for (int i = CASTART(depth); i < CAEND(depth) - 1; i++){
			ui v = stkvtx[i];
			int kept = 1;
			int con = con2p[depth][v];
			//check if v is kept in the new solution			
			for (auto x: satu){
				//x is not connected to v, remove
				if (!bg->rows[x]->test(v)){
					kept = 0;
					break;
				}
			}			
			if (kept){				
				if (bg->rows[v]->test(u)){
					con++;
				}
				if (con + K < depth + 2){
					kept = 0;
				}
			}

			if (kept){
				//push into stack without changing the order
				con2p[depth+1][v] = con;
				stkvtx[pos++] = v;				
			}
			// else
			// 	removed.push_back(v);
		}
		CAEND(depth+1) = pos;
	}
	

#ifdef DEBUG
	void showStack() {
		printf("depth: %d\n", depth);
		printf("PLEX SZ %d : ", SZ_PLEX(depth));
		for (int i = PLSTART(depth); i < PLEND(depth); i++){
			int v = stkvtx[i];
			printf("%d(toP:%d, deg%d)\n  ", v, stkcon2p[i], stkdeg[i][v]);
		}
		printf("\n CAND SZ %d: ", SZ_CAND(depth));
		for (int i = CASTART(depth); i < CAEND(depth); i++){
			ui v = stkvtx[i];
			printf("%d(toP: %d, deg%d)\n ", v, stkcon2p[i], stkdeg[depth][stkvtx[i]]);
		}
		printf("\n");
	}
	#else
	void showStack() {}
	#endif
	// void checkDegree(){
	// 	int *deg = new int[cg.nbvtx];
	// 	bool exit1 = false;
	// 	for (int i = PLSTART(depth); i < CAEND(depth); i++){
	// 		ui u = stkvtx[i];
	// 		deg[u] = 0;
	// 		for (int j = PLSTART(depth); j < CAEND(depth); j++){
	// 			ui v = stkvtx[j];
	// 			if (bg->rows[u]->test(v)){
	// 				deg[u]++;
	// 			}
	// 		}
	// 		if (deg[u] != degree[depth][u]){
	// 			exit1 = true;
	// 		}
	// 		assert(deg[u] == degree[depth][u]);
	// 	}	
	// 	if (exit1){
	// 		cout << "!!!!!!!!!\n";
	// 		exit(1);
	// 	}
	// 	delete[] deg;
	// }

	void expand() {
		//int stop = 0;
		nnodes++;
		assert(depth == SZ_PLEX(depth));
		//checkDegree();
		if (SZ_PLEX(depth) > szsol){
			copy(stkvtx + PLSTART(depth), stkvtx + PLEND(depth), sol);
			szsol = SZ_PLEX(depth);
		}		
	#ifdef DEBUG
			showStack();
	#endif		
		if (SZ_VTX(depth) <= szsol)
			return;
		refineCandidate();
		if (SZ_VTX(depth) <= szsol)
			return;
		
		// if (depth > 0){
		// 	for (int i = PLSTART(depth); i < CAEND(depth); i++){
		// 		ui v = stkvtx[i];
		// 		degree[depth][v] = degree[depth-1][v];
		// 		for (ui u : removed){
		// 			if (bg->rows[v]->test(u))
		// 				degree[depth][v]--;
		// 		}
		// 	}
		// 	removed.clear();
		// }

		//checkDegree();
		bitColorSort3();
		while (SZ_CAND(depth) > 0) {
			if (interrupt()) { //stop due to interruption
				isoptimal = 0;
				break;
			}
			ui last = stkvtx[CAEND(depth) - 1];
			//ui lastv = stkvtx[CAEND(depth) - 1];
			if (SZ_PLEX(depth) + bound[depth][last] <= szsol || SZ_PLEX(depth) + SZ_CAND(depth) <= szsol) {
				return;
			}
			/*  color bound*/			
			int tol = K - SZ_PLEX(depth) + con2p[depth][last];
			int potential = lookahead2(CAEND(depth) - 1);
			if (SZ_PLEX(depth) + potential + tol <= szsol ) {
				// for (int i = PLSTART(depth); i < CAEND(depth); i++){
				// 	ui v = stkvtx[i];
				// 	if (bg->rows[v]->test(last))
				// 		degree[depth][v]--;
				// }
				//assert(potent[depth][last] == lookahead(last));			
				--CAEND(depth);
				continue;
			}
			//update Plex
			selectLast();	
			++depth;	
				
			expand();
			// removed.clear();
			//drop the last vertex
			--depth;	
			--CAEND(depth);
			// for (int i = PLSTART(depth); i < CAEND(depth); i++){
			// 	ui v = stkvtx[i];
			// 	if (bg->rows[v]->test(last))
			// 		degree[depth][v]--;
			// }
		}
	}

    void search(){
		startclk = clock();
		//start to search the optimal solution.
		expand(); 
		endclk = clock();
	}

    double getRunningTime(){
		return Utility::elapse_seconds(startclk, endclk);
	}
	uint64_t getRunningNodes(){
		return nnodes;
	}
    void getBestSolution(ui *vtx, ui &sz){
		sz = szsol;
		memcpy(vtx, sol, sizeof(ui) * szsol);
	}
	
};





class FastHeuSearch{
public:
	MCsrGraph &g;
	int K;
	int *core;
	int issort;
	ui *seq;
	int *pos;
	
	ui* sol;
	int szsol;
	
	ui* curP;
	int szP;
	int* degP;
	const int maxdepth = 1;

	FastHeuSearch(MCsrGraph &csrg, int valueK):	g(csrg), K(valueK){
		core = new int[g.nbvtx];
		seq = new ui[g.nbvtx];
		pos = new int[g.nbvtx];
		szsol = 0;
		sol = new ui[g.nbvtx];
		issort= 0;
	}

	~FastHeuSearch(){
		delete[] core;
		delete[] seq;
		delete[] pos;
		delete[] sol;
 	}

	void Add2P(ui u){
		for (ui j = g.pstart[u]; j < g.pstart[u+1]; j++){
			ui nei = g.edges[j];
			degP[nei]++;
		}
	}
	void DelFrP(ui u){
		for (ui j = g.pstart[u]; j < g.pstart[u+1]; j++){
			ui nei = g.edges[j];
			degP[nei]--;
		}
	}
	int isFeasible(ui u){
		if (degP[u] + K < szP + 1)	return 0;
		int isf = 1;
		for (int i = 0; i < szP; i++){
			int v = curP[i];
			if (degP[v] == szP - K){ //satu
				ui *it = find(g.edges+ g.pstart[u], g.edges+g.pstart[u+1], v);
				if (it == g.edges + g.pstart[u+1]){
					isf = 0;
					break;
				}
			}
		}
		return isf;
	}
	void expand(ui* cand, ui szc){
		while(szc > 0){
			if (szP + szc <= szsol || szP + core[cand[szc-1]] + K <= szsol){ 
				break;
			}
			ui u = cand[--szc];//pop u
			curP[szP++] = u;
			Add2P(u);
			ui tsz = szc;
			szc = 0;
			for (int j = 0; j < tsz; j++){
				if (isFeasible(cand[j])){
					cand[szc++] = cand[j];
				}
			}
			if (szP > szsol){
				szsol = szP;
				memcpy(sol, curP, szsol * sizeof(ui));
			}
		}
	}

	void heuPolish(){
		if (!issort){
			coreDecomposition(g, core, seq, pos);
			issort = 1;
		} 
		curP = new ui[g.nbvtx];
		szP = 0;
		degP = new int[g.nbvtx];
		
		
		ui *cand = new ui[g.nbvtx];
		ui szc = 0;
		
		for (int i = 0; i < g.nbvtx; i++){
			ui u = seq[g.nbvtx - i -1];
			if (core[u] + K <= szsol) break;
			curP[0] = u;
			szP = 1;
			szc = g.nbvtx - i - 1;	
			memcpy(cand, seq, sizeof(ui) * szc);		
			memset(degP, 0, sizeof(int) * g.nbvtx);
			for (ui j = g.pstart[u]; j < g.pstart[u+1]; j++)
				degP[g.edges[j]] = 1; 
			expand(cand, szc);			

			szP = 0;
		}
		delete[] cand;
		printf("The maximum solution after polish %d\n", szsol);
	}
	void coreHeurSolution(){    
		assert(core != nullptr && seq != nullptr && pos != nullptr && sol != nullptr);

		int *deg = new int[g.nbvtx];
		szsol = 0;
		MBitSet64 *issol = new MBitSet64(g.nbvtx);
		int maxcore = -1;
		int hit = 0;

		for (ui i = 0; i < g.nbvtx; i++) {		
			deg[i] = g.degree(i);
			seq[i] = i;
			core[i] = -1;
		}
		ListLinearHeap *linear_heap = new ListLinearHeap(g.nbvtx, g.nbvtx);
		linear_heap->init(g.nbvtx, g.nbvtx - 1, seq, (ui*)deg);
		for (ui i = 0; i < g.nbvtx; i++) {
			ui u, key;
			linear_heap->pop_min(u, key);
			if ((int)key > maxcore) 
				maxcore = key;		
			seq[i] = u;
			core[u] = maxcore;
			pos[u] = i;
			if (key >= g.nbvtx-i -K && !hit){
				ui sz = i+1;			
				linear_heap->get_ids_of_larger_keys(seq, sz, key);							
				for (ui j = i+1; j < g.nbvtx; j++){//vertices after i is not ordered
					//core[seq[j]] = max_core;				
					sol[szsol++] = seq[j];				
					issol->set(seq[j]);
				}
				hit = 1;
				//break;
			}
		
			for (ui j = g.pstart[u]; j < g.pstart[u + 1]; j++){
				if (core[g.edges[j]] == -1)
					linear_heap->decrement(g.edges[j]);
			}
		}	
		issort = 1;
		//extend the heuristic solution.
		//compute degrees
		memset(deg, 0, sizeof(int) * g.nbvtx);	
		//MBitSet64 *issatu = new MBitSet64(og_nbvtx);
		int nsatu = 0;
		for (int i = 0; i < szsol; i++){
			ui u = sol[i];
			//printf("%d ", u);
			for (ui j = g.pstart[u]; j < g.pstart[u+1]; j++){
				ui nei = g.edges[j];
				if (issol->test(nei)){
					deg[u]++;
				}
			}
			if (deg[u] == szsol - K){
				nsatu ++;
			}
		}
		//extending the remaining vertices
		for (int i = g.nbvtx - szsol - 1; i >= 0; i--){
			ui u = seq[i];
			if (core[u] + K <= szsol) break;
			int cntnei = 0, cntsatu = 0;
			for (int j = g.pstart[u]; j < g.pstart[u+1]; j++){
				ui nei = g.edges[j];
				if (issol->test(nei)){
					cntnei++;
					if (deg[nei] == szsol - K){ //nei is satu
						cntsatu++;
					}
				}			
			}
			//extend u
			if (cntnei >= szsol - K + 1 && cntsatu == nsatu){ 
				printf("extend %d\n",u);
				assert(deg[u] == 0);
				deg[u] = cntnei;
				sol[szsol++] = u;			
				issol->set(u);
				//update number of saturated vertices.			
				for (ui j = g.pstart[u]; j < g.pstart[u+1]; j++){
					ui nei = g.edges[j];
					if (issol->test(nei)){
						deg[nei]++;
					}
				}
				//recount satu vertices
				nsatu = 0;
				for (int j = 0; j < szsol; j++){				
					if (deg[sol[j]] == szsol - K){
						nsatu++;
					}
				}
			}
		}
		
		delete linear_heap;
		delete issol;
		delete[] deg;
		printf("The maximum heuristc solution %d\n", szsol);

	}	
};



