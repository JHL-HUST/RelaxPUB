
/*
 *
 * g++ -O3 DiseMKP.cpp -o disemkp -std=c++11 -DGOP
 * usage: ./disemkp instance -x k
 *
 */

#include "MBitSet.h"
#include "MGraph.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <limits.h>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>
#include <assert.h>
#include <bitset>
#include <iostream>

#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define MAX_NODE 80000000
#define max_expand_depth 100000
#define pop(stack) stack[--stack ## _fill_pointer]
#define push(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer
//#define is_neibor(i,j) matrice[i][j]

#define CUR_KPX_SIZE Clique_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]
//#define MIN(a,b) a<=b?a:b
//#define BIT_MAP_SIZE 4097


using namespace std;

static int TIME_OUT, CUT_OFF=0;
static double BEST_SOL_TIME;

static int FORMAT = 1, NB_NODE, NB_NODE_O, ADDED_NODE, NB_EDGE, NB_EDGE_O,
MAX_KPX_SIZE, MAX_ISET_SIZE, INIT_KPX_SIZE, HEUR_KPX_SIZE,
NB_BACK_CLIQUE, MATRIX_ROW_WIDTH, MAX_VERTEX_NO, K_CORE_G = 0, Reasoning_Point=0,MAX_MATRIX_SIZE=2000;


static int Max_Degree = 0, UPPER_BOUND=0;
static int Node_Degree[MAX_NODE];
static int *CNN,*Node_Count;


std::bitset<MAX_NODE> Node_State;
std::bitset<MAX_NODE> Node_State2;
std::bitset<MAX_NODE> Node_State3;

static char instance[1024]={'\0'};

static int **Node_Neibors;

static int Candidate_Stack_fill_pointer = 0;
static int Candidate_Stack[MAX_NODE * 2];
static int NBNN_Stack_fill_pointer = 0;
static char NBNN_Stack[MAX_NODE*2];
static int Clique_Stack_fill_pointer,Energy_Stack_fill_pointer;
static int  *Clique_Stack, *MaxCLQ_Stack, *Energy_Stack,*Index_Stack;
static int Cursor_Stack[max_expand_depth];
static int Cursor_Stack_fill_pointer = 0;
static char *Removed;
static char *Touched;

static unsigned char * Adj_Matrix;


static int Rollback_Point;
static int Branching_Point;

static int *Old_Name;
static int *Second_Name;
static int NB_CANDIDATE = 0, FIRST_INDEX;
static int START_MAXSAT_THD = 15;

static int Extra_Node_Stack_fill_pointer = 0;
static int *Extra_Node_Stack;

static int Last_Idx = 0;
static int cut_ver = 0, total_cut_ver = 0;
static int cut_inc = 0, total_cut_inc = 0;
static int cut_iset = 0, total_cut_iset = 0;
static int cut_satz = 0, total_cut_satz = 0;
static long long Branches_Nodes[6];
static int LAST_IN;
static float Dynamic_Radio = 0.70;
static int REBUILD_MATRIX = FALSE;
static int CUR_MAX_NODE;
static int Branches[1200];
static int * Init_Adj_List;
static int BLOCK_COUNT = 0;
static int *BLOCK_LIST[100];
static double READ_TIME, INIT_TIME, SEARCH_TIME;
static long long N0_0 = 0, N0_1 = 0, N1_0 = 0, N1_1 = 0, L1 = 0;
static double D0 = 0, D1 = 0, D2 = 0, Dt = 0;

typedef struct {
	MBitSet64 bitcolor;
	int *vertices;
	int sz;
}ColorType;
ColorType *colors;

MCsrGraph kernalG;
MBitGraph *bg;
int *target_normal, *target_abnormal;
int max_node_index;

int *ranked;
int *nbnn;
float *profits;
bool *colored;
bool *conflicted;
int *nbnns;
int *nbnns2;
//vector<int> insert;
vector<int> colored_nodes;
vector<int> record;
//vector<int> removed;

int *nondegree_temp;

// FILE* pp = fopen("result.xls","a+");

struct Paramerters{
	int KX;
	int ORDER;
	int FORMAT;
	int START_MAXSAT_THD;
	int CUT_OFF;
	int TIME_OUT;
}PARA;

static double get_utime() {
	struct rusage utime;
	getrusage(RUSAGE_SELF, &utime);
	return (double) (utime.ru_utime.tv_sec
					 + (double) utime.ru_utime.tv_usec / 1000000);
}

static int int_cmp_desc(const void * a, const void * b) {
	return *((int *) b) - *((int *) a);
}

static int int_cmp_asc(const void * a, const void * b) {
	return *((int *) a) - *((int *) b);
}
static int is_adjacent(int node1, int node2) {
	int neibor, *neibors;
	int end;
	if(node1>node2){
		neibors = Node_Neibors[node1];
		end=node1;
	}else{
		neibors = Node_Neibors[node2];
		node2=node1;
		end=node2;
	}
	for (neibor = *neibors; neibor <node2 && neibor!=end ; neibor = *(neibors+=2));
	return neibor==node2;
}

static int nbE=0;

static void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge,int offset) {
	int i, block_size = 40960000, free_size = 0;
	Init_Adj_List = (int *) malloc((2 * nb_edge + nb_node) * sizeof(int));
	if (Init_Adj_List == NULL ) {
		for (i = 1; i <= NB_NODE; i++) {
			if (Node_Degree[i - offset] + 1 > free_size) {
				Node_Neibors[i] = (int *) malloc(block_size * sizeof(int));
				BLOCK_LIST[BLOCK_COUNT++] = Node_Neibors[i];
				free_size = block_size - (Node_Degree[i - offset] + 1);
			} else {
				Node_Neibors[i] = Node_Neibors[i - 1]
				+ Node_Degree[i - 1 - offset] + 1;
				free_size = free_size - (Node_Degree[i - offset] + 1);
			}
		}
	} else {
		BLOCK_COUNT = 1;
		BLOCK_LIST[BLOCK_COUNT - 1] = Init_Adj_List;
		Node_Neibors[1] = Init_Adj_List;
		for (i = 2; i <= NB_NODE; i++) {
			Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1 - offset]
			+ 1;
		}
	}
}
static int read_graph_node_node(char *input_file, int format) {
	int j, l_node, r_node, nb_edge = 0, max_node = NONE, offset = 0;
	int node = 1;
	char line[128], word[10];
	FILE* fp_in = fopen(input_file, "r");
	
	if (fp_in == NULL ) {
		printf("R can not open file %s\n", input_file);
		return FALSE;
	}
	
	if (format == 1)
		printf("R reading file <e n1 n2> ...\n");
	else
		printf("R reading file <n1 n2> ...\n");
	
	memset(Node_Degree, 0, MAX_NODE*sizeof(int));
	
	while (fgets(line, 127, fp_in) != NULL ) {
		if ((format == 1 && line[0] == 'e')
			|| (format == 2 && line[0] != '%')) {
			if (format == 1)
			  sscanf(line, "%s%d%d", word, &l_node, &r_node);
			else
			  sscanf(line, "%d%d", &l_node, &r_node);
			 
			if (l_node >= 0 && r_node >= 0 && l_node != r_node) {
				
				nb_edge++;
				
				if (l_node > max_node)
					max_node = l_node;
				if (r_node > max_node)
					max_node = r_node;
				
				if (offset ==0 && (l_node == 0 || r_node == 0)){
					offset = 1;
				}
				
				if (max_node+offset>=MAX_NODE) {
					printf("! The graph goes beyond the maximum size (%d) can be processed.\n",MAX_NODE);
					printf("! Please modify the definition of the variable MAX_NODE to fit the size.\n");
					exit(1);
				}
				
				Node_Degree[l_node]++;
				Node_Degree[r_node]++;
				
			}
		}
	}
	NB_NODE = max_node;
	NB_NODE = NB_NODE + offset;
	
	printf("R the graph size is %d\n", NB_NODE);
	
	Node_Neibors = (int **) malloc((NB_NODE + 1) * sizeof(int *));
	allcoate_memory_for_adjacency_list(NB_NODE, nb_edge, offset);
	memset(Node_Degree, 0, (NB_NODE + 1) * sizeof(int));
	
	nb_edge = 0;
	fseek(fp_in, 0L, SEEK_SET);
	while (fgets(line, 127, fp_in) != NULL ) {
		if ((format == 1 && line[0] == 'e')
			|| (format == 2 && line[0] != '%')) {
			if (format == 1)
				sscanf(line, "%s%d%d", word, &l_node, &r_node);
			else
				sscanf(line, "%d%d", &l_node, &r_node);
			if (l_node >= 0 && r_node >= 0 && l_node != r_node) {
				if(offset){
					l_node +=offset;
					r_node +=offset;
				}
				for (j = 0; j < Node_Degree[l_node]; j++) {
					if (Node_Neibors[l_node][j] == r_node)
						break;
				}
				if (j == Node_Degree[l_node]) {
					Node_Neibors[l_node][Node_Degree[l_node]] = r_node;
					Node_Neibors[r_node][Node_Degree[r_node]] = l_node;
					Node_Degree[l_node]++;
					Node_Degree[r_node]++;
					nb_edge++;
				}
			}
		}
	}
	NB_EDGE = nb_edge;
	Max_Degree = 0;
	N0_0 = NB_NODE;
	for (node = 1; node <= NB_NODE; node++) {
		Node_Neibors[node][Node_Degree[node]] = NONE;		
		if (Node_Degree[node] > Max_Degree)
			Max_Degree = Node_Degree[node];
	}
	UPPER_BOUND=Max_Degree+PARA.KX;
	return TRUE;
}

static int build_simple_graph_instance(char *input_file) {
	printf("# reading instance ...\n");
	const char * fileStyle="clq";
	if(strrchr(input_file, '.')!=NULL)
		fileStyle = strrchr(input_file, '.') + 1;
	
	if (strcmp(fileStyle, "clq") == 0) {
		read_graph_node_node(input_file, 1);
	} else if (strcmp(fileStyle, "edges") == 0) {
		read_graph_node_node(input_file, 2);
	} else if (strcmp(fileStyle, "mtx") == 0) {
		read_graph_node_node(input_file, 2);
	} else if (FORMAT == 1) {
		read_graph_node_node(input_file, 1);
	} else if (FORMAT == 2) {
		read_graph_node_node(input_file, 2);
	} else {
		read_graph_node_node(input_file, 1);
	}
	printf("R Instance Information: #node=%d #edge=%d density=%9.8f\n", NB_NODE,
		   NB_EDGE, ((float) NB_EDGE * 2 / NB_NODE / (NB_NODE - 1)));
	NB_NODE_O = NB_NODE;
	NB_EDGE_O = NB_EDGE;
	D0 = ((float) NB_EDGE * 2 / NB_NODE / (NB_NODE - 1));
	
	READ_TIME = get_utime();
	printf("R the reading time is %4.2lf \n", READ_TIME);
	return TRUE;
}

static int sort_by_degeneracy_ordering() {
	int *degree_counter, *where;
	int p1, i, node = NONE, neibor, *neibors, t, j, h, k;
	int cur_degree = 0;
	INIT_KPX_SIZE = 0;
	printf("I computing an initial k-plex...\n");
	
	where = Candidate_Stack + NB_NODE + 1;
	degree_counter=(int *)calloc(Max_Degree + 1,sizeof(int));
	
	
	for (node = 1; node <= NB_NODE; node++) {
	  assert(Node_Degree[node]<=Max_Degree);
	  degree_counter[Node_Degree[node]]++;
	}
	j = 0;
	for (i = 0; i <= Max_Degree; i++) {
		k = degree_counter[i];
		degree_counter[i] = j;
		j += k;
	}
	
	for (node = 1; node <= NB_NODE; node++) {
		Candidate_Stack[t = degree_counter[Node_Degree[node]]++] = node;
		where[node] = t;
	}
	
	for (i = Max_Degree; i > 0; i--) {
		degree_counter[i] = degree_counter[i - 1];
	}
	degree_counter[0] = 0;
	
	Candidate_Stack[NB_NODE] = DELIMITER;
	ptr(Candidate_Stack) = NB_NODE + 1;
	
	p1 = 0;
	cur_degree = Node_Degree[Candidate_Stack[p1]];
	while (p1 < NB_NODE) {
		node = Candidate_Stack[p1];
		
		if (cur_degree > K_CORE_G)
			K_CORE_G = cur_degree;
		
		if (p1 < NB_NODE - 1
			&& Node_Degree[node] == Node_Degree[Candidate_Stack[p1 + 1]]) {
			degree_counter[Node_Degree[node]] = p1 + 1;
		}
		if (Node_Degree[node] > MAX_VERTEX_NO)
			MAX_VERTEX_NO = Node_Degree[node];
		
		if (Node_Degree[node] >= NB_NODE - p1 - PARA.KX) { //when KX=1, it is a clique.
			BEST_SOL_TIME=get_utime();
			MAX_KPX_SIZE=HEUR_KPX_SIZE=INIT_KPX_SIZE = NB_NODE - p1;
			printf("I the upper bound of k-plex %d ...\nI the initial %d-plex  %d...\n", UPPER_BOUND, PARA.KX,INIT_KPX_SIZE);
			MaxCLQ_Stack = (int *) malloc((UPPER_BOUND+1) * sizeof(int));
			Clique_Stack = (int *) malloc((UPPER_BOUND+1) * sizeof(int));
			Energy_Stack = (int *) malloc((UPPER_BOUND+1) * sizeof(int));
			Index_Stack  = (int *) malloc((UPPER_BOUND+1) * sizeof(int));
			assert(MaxCLQ_Stack);assert(Clique_Stack);assert(Energy_Stack);assert(Index_Stack);	 
			memcpy(MaxCLQ_Stack, Candidate_Stack + p1,INIT_KPX_SIZE * sizeof(int));
			break;
		}		
		neibors = Node_Neibors[node];
		for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
			if (where[neibor] > p1) {
				t = where[neibor];
				h = degree_counter[Node_Degree[neibor]];
				
				k = Candidate_Stack[h];
				
				Candidate_Stack[h] = neibor;
				where[neibor] = h;
				
				Candidate_Stack[t] = k;
				where[k] = t;
				
				degree_counter[Node_Degree[neibor]]++;
				
				Node_Degree[neibor]--;
				if (Node_Degree[neibor]
					!= Node_Degree[Candidate_Stack[h - 1]]) {
					degree_counter[Node_Degree[neibor]] = h;
				}
			}
		}
		p1++;
	}
	free(degree_counter);
	
	if (UPPER_BOUND == INIT_KPX_SIZE) {
		MAX_KPX_SIZE = INIT_KPX_SIZE;
		printf("I find the maximum %d-plex in initial phase!\n",PARA.KX);
		return TRUE;
	} else {
		return FALSE;
	}
}




static long long BRANCHING_COUNT = 0;

static void reduce_instance_with_unsupport_property();
static int rebuild_instance(int pre);
static void	 check_solution();

static void store_maximum_clique(int node, int silent) {
	if (Reasoning_Point == 0)
		push(node, Clique_Stack);
	else
		push(Second_Name[node], Clique_Stack);
	
	MAX_KPX_SIZE = ptr(Clique_Stack);
	BEST_SOL_TIME=get_utime();
	memcpy(MaxCLQ_Stack, Clique_Stack, MAX_KPX_SIZE * sizeof(int));
	
	ptr(NBNN_Stack) = ptr(Candidate_Stack) = NB_NODE + 1;
	ptr(Cursor_Stack) = 1;
	ptr(Clique_Stack) = 0;
	ptr(Energy_Stack) =0;
	Rollback_Point = 0;
	Reasoning_Point=0;
	//	Vertex_UB[CURSOR]=MAX_KPX_SIZE;
	if (silent==0)
		printf("C %4d |%7d |%14lld| %8.2lf |       ", MAX_KPX_SIZE, CURSOR,BRANCHING_COUNT,BEST_SOL_TIME);
	total_cut_ver += cut_ver;
	cut_ver = 0;
	total_cut_inc += cut_inc;
	cut_inc = 0;
	total_cut_iset += cut_iset;
	cut_iset = 0;
	total_cut_satz += cut_satz;
	cut_satz = 0;
	Last_Idx = CURSOR;

	
	check_solution();
#ifndef NOIP
	int cur_node=Candidate_Stack[CURSOR];
	  
	reduce_instance_with_unsupport_property();
	if(rebuild_instance(0)){
	  CURSOR=0;
	}else{
	  CURSOR=ptr(Candidate_Stack)-1;
	  for(int i=ptr(Candidate_Stack)-2, node=Candidate_Stack[i];i>0;node=Candidate_Stack[--i]){
	    if(node<=cur_node){							     
	      CURSOR=i;
	      if(node==cur_node)
		break;
	    }else{
	      break;
	    }
	  }
	}
#endif
}

static long long cut_level=0;
static long long cut_count=0;

static float BMTHD=0.3;

static inline void reset_state(){
	for(int i=CURSOR+1, cn=Candidate_Stack[i];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		Node_State.reset(cn);
		Node_State2.reset(cn);
	}
}


static inline void reset_state1(){
	for(int i=0;i<ptr(Clique_Stack);i++){
		if(Energy_Stack[i]<0){
			Energy_Stack[i]=-Energy_Stack[i]-1;
		}
	}
	for(int i=CURSOR+1, cn=Candidate_Stack[i];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		Node_State.reset(cn);
		Node_State2.reset(cn);
	}
}
static inline void reset_state2(){
	for(int i=0;i<ptr(Clique_Stack);i++){
		if(Energy_Stack[i]<0){
			Energy_Stack[i]=-Energy_Stack[i]-1;
		}
	}
}


static int produce_subgraph0 () {
	int i = CURSOR,j=0, neibor, max = 0, *neibors,*neibors2;	
	int start=ptr(Candidate_Stack);	
	int bnode=Candidate_Stack[CURSOR];	
	assert(ptr(Candidate_Stack) == ptr(NBNN_Stack));
	
	for(int cn=Candidate_Stack[i=CURSOR+1];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		Node_State.reset(cn);
		Node_State2.set(cn);
	}
	//cout << Candidate_Stack[i=CURSOR+1] <<" "<< i - CURSOR - 1 << " ";
	//mark neighbors of current branching vertex,and set state	
	neibors = Node_Neibors[bnode];
	for (neibor = *neibors; neibor != bnode; neibor = *(neibors+=2)) {
	  assert(neibor>0);
		if(Node_State2[neibor])
			Node_State.set(neibor);
	}
	//removed.clear();
	NB_CANDIDATE=0;
	for(int cn=Candidate_Stack[i=CURSOR+1];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		assert(cn!=bnode);
		if(Node_State[cn]){
			push(cn, Candidate_Stack);
			push(NBNN_Stack[i],NBNN_Stack);
			NB_CANDIDATE++;
		}else if((NBNN_Stack[CURSOR]<PARA.KX-1) &&  ( NBNN_Stack[i]<PARA.KX-1)){
			push(cn, Candidate_Stack);
			push(NBNN_Stack[i]+1,NBNN_Stack);
			NB_CANDIDATE++;
		}
		//else
		//	removed.push_back(cn);
	}

	push(DELIMITER, Candidate_Stack);
	push(DELIMITER, NBNN_Stack);
	// check_candidates3();
	
	if(NB_CANDIDATE+CUR_KPX_SIZE<=MAX_KPX_SIZE){
		reset_state();
		return 0;
	}
	
	for(int cn=Candidate_Stack[i=CURSOR+1];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		Node_State.reset(cn);
		Node_State2.reset(cn);
	}
	
	for(int cn=Candidate_Stack[j=start];cn!=DELIMITER;cn=Candidate_Stack[++j]){
		Node_State2.set(cn);
	}
	
	//-----------------------------//
	
	int count=0;
	for(i=0;i<ptr(Clique_Stack)-1;i++){
		if(is_adjacent(Clique_Stack[i],bnode)==FALSE){
			assert(Energy_Stack[i]<PARA.KX-1);
			Energy_Stack[i]=-(Energy_Stack[i]+1);
			count++;
		}
	}
	assert(count==NBNN_Stack[CURSOR]);
	
	assert(ptr(Candidate_Stack) == ptr(NBNN_Stack));
	
	for(i=0;i<ptr(Clique_Stack)-1;i++){
		assert(Energy_Stack[i]<PARA.KX);
		if((Energy_Stack[i]<0 && (Energy_Stack[i] == 1-PARA.KX))){
			int nn=Clique_Stack[i];
			
			neibors = Node_Neibors[nn];
			for (neibor = *neibors; neibor != nn; neibor = *(neibors+=2)) {
				if(Node_State2[neibor])
					Node_State.set(neibor);
			}
			
			for(int cn=Candidate_Stack[j=start];cn!=DELIMITER;cn=Candidate_Stack[++j]){
				assert(nn!=cn || nn+cn!=0);
				if(cn>0 && Node_State[cn]==0){
					Candidate_Stack[j]=-cn;
					NB_CANDIDATE--;
				}
			}
			
			neibors = Node_Neibors[nn];
			for (neibor = *neibors; neibor != nn; neibor = *(neibors+=2)) {
				if(Node_State2[neibor])
					Node_State.reset(neibor);
			}
			
			if(NB_CANDIDATE+CUR_KPX_SIZE<=MAX_KPX_SIZE){
				for(int cn=Candidate_Stack[j=start];cn!=DELIMITER;cn=Candidate_Stack[++j]){
					if(cn>0)
						Node_State2.reset(cn);
					else{
						//removed.push_back(-cn);
						Node_State2.reset(-cn);
					}
				}
				for(int i=0;i<ptr(Clique_Stack);i++){
					if(Energy_Stack[i]<0){
						Energy_Stack[i]=-Energy_Stack[i]-1;
					}
				}				
				return 0;
			}
		}
	}	
	int _count=0;
	i=j=start;
	for(int cn=Candidate_Stack[i];cn!=DELIMITER;cn=Candidate_Stack[++i]){
		if(cn>0){
#if !defined(COMM)  && !defined(GOP) && !defined(SEC)
			Node_State2.reset(cn);
#endif
			Candidate_Stack[j]=cn;
			NBNN_Stack[j]=NBNN_Stack[i];
			_count++;
			j++;
		}else{
			//removed.push_back(-cn);
			Node_State2.reset(-cn);
		}
	}
	assert(_count==NB_CANDIDATE);
	Candidate_Stack[j]=DELIMITER;
	NBNN_Stack[j]=DELIMITER;
	ptr(Candidate_Stack)=j+1;
	ptr(NBNN_Stack)=j+1;
	
	if(NB_CANDIDATE+CUR_KPX_SIZE>MAX_KPX_SIZE){
		
#if !defined(COMM)  && !defined(GOP) && !defined(SEC)
		for(i=0;i<ptr(Clique_Stack);i++){
			if(Energy_Stack[i]<0){
				Energy_Stack[i]=-Energy_Stack[i];
			}
		}
#endif	
		return ptr(Candidate_Stack)-1-(MAX_KPX_SIZE-CUR_KPX_SIZE);
	}else{
		reset_state1();
		return 0;
	}
}

static void qsort_profits(int s[], int l, int r)
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

static void qsort_value(int s[], int l, int r){
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

static int cut_by_iteration_partition_SEC();

static int con_num(int tnode, int color_num, bool *conflicted, int color_ub) {
		record.clear();
		//char *t_matrix = matrix + SR[i] * n;
		//ui c_num = 1;
		for (int j = 0; j < color_num; j++)
		{
			if (nbnns2[j] > color_ub)
				continue;
			int u = colors[1].vertices[j];
			if (bg->rows[target_normal[tnode]-1]->test(u) && !conflicted[j])
			{
				//record[c_num - 1] = j;
				record.push_back(j);
				//c_num++;
			}
		}
		return record.size() + 1;
	}

static int con_num2(int tnode, int color_num) {
		//char *t_matrix = matrix + SR[i] * n;
		int cnum = 0;
		for (int j = 0; j < color_num; j++)
		{
			int u = colors[1].vertices[j];
			if (bg->rows[target_normal[tnode]-1]->test(u))
				cnum++;
		}
		return cnum;
	}

inline int intsecColor(unsigned int u, int c) {
	return colors[c].bitcolor.intersect(*(bg->rows[u]));
}

static void try_color(int &color_ub, int lb, int start, int end){
	int large_num = 0;
	colored_nodes.clear();

	for (int u = start; u <= end; u++){
	//for (int u = end; u >= start; u--){
		int i = ranked[u - start];
		int tnode = Candidate_Stack[i + start];
		if (Node_State[tnode] == 0)
			continue;
		int invest = nbnn[i];
		if (invest > lb){
			large_num++;
			if (large_num > lb)
				continue;
		}
		//colors[1].bitcolor.intersect(*(bg->rows[target_normal[tnode]-1]));
		//if (intsecColor(target_normal[tnode]-1, 1) == 0){
		//int con = intsecColor(target_normal[tnode]-1, 1);
		//if (colors[1].bitcolor.intersect(*(bg->rows[target_normal[tnode]-1])) == 0){
		//if (can_color(i, color_num)){
		
		if (intsecColor(target_normal[tnode]-1, 1) == 0){
			colored[i] = true;
			nbnns[colors[1].sz] = invest;
			colors[1].bitcolor.set(target_normal[tnode]-1);
			colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
			colors[1].sz++;	
			colored_nodes.push_back(i);	
		}
	}

	int color_num_tmp = colors[1].sz;

	if (colors[1].sz > 0) {
		memset(nbnns2, 0, sizeof(int) * (end - start + 1));
		memcpy(nbnns2, nbnns, sizeof(int) * colors[1].sz);
		qsort_value(nbnns, 0, colors[1].sz - 1);
		for (int j = 0; j < colors[1].sz; j++)
		{
			if (nbnns[j] > color_ub)
				color_ub++;
			else
				break;
		}
		
		memset(conflicted, false, sizeof(bool) * (end - start + 1));
		int can_add = color_ub;
		for (int j = 0; j < colors[1].sz; j++)
		{
			if (nbnns[j] > color_ub){
				//conflicted[j] = true;
				can_add--;
			}
			else
				break;
		}

		if (can_add > 0){
			for (int u = start; u <= end; u++)
			{
				int i = ranked[u - start];
				int tnode = Candidate_Stack[i + start];
				if (colored[i] || Node_State[tnode] == 0 || nbnn[i] < color_ub)
					continue;
				int add_num = con_num(tnode, colors[1].sz, conflicted, color_ub);
				if (add_num <= can_add)
				{
					can_add -= add_num;
					for (int j = 0; j < add_num - 1; j++)
						conflicted[record[j]] = true;
					//insert[insert_num++] = i;
					//insert.push_back(i);
					colors[1].bitcolor.set(target_normal[tnode]-1);
					colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
					conflicted[colors[1].sz] = true;
					colors[1].sz++;	
					colored_nodes.push_back(i);	
					if (can_add == 0)
						break;
				}
			}
		}

		for (int u = start; u <= end; u++)
		{
			int i = ranked[u - start];
			int tnode = Candidate_Stack[i + start];
			if (colored[i] || Node_State[tnode] == 0 || nbnn[i] >= color_ub)
				continue;
			int cnum = con_num2(tnode, colors[1].sz);
			if (cnum <= color_ub - nbnn[i])
			{
				colors[1].bitcolor.set(target_normal[tnode]-1);
				colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
				colors[1].sz++;	
				colored_nodes.push_back(i);	
			}
		}

		// if (color_ub == nbnns[0])
		// {
		// 	//insert.clear();
		// 	int can_add = color_ub;
		// 	memset(conflicted, false, sizeof(bool) * (end - start + 1));

		// 	for (int u = start; u <= end; u++)
		// 	//for (int u = end; u >= start; u--)
		// 	{
		// 		int i = ranked[u - start];
		// 		int tnode = Candidate_Stack[i + start];
		// 		if (colored[i] || Node_State[tnode] == 0 || nbnn[i] < color_ub)
		// 			continue;
		// 		int add_num = con_num(tnode, colors[1].sz, conflicted);
		// 		if (add_num <= can_add)
		// 		{
		// 			can_add -= add_num;
		// 			for (int j = 0; j < add_num - 1; j++)
		// 				conflicted[record[j]] = true;
		// 			//insert[insert_num++] = i;
		// 			//insert.push_back(i);
		// 			colors[1].bitcolor.set(target_normal[tnode]-1);
		// 			colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
		// 			conflicted[colors[1].sz] = true;
		// 			colors[1].sz++;	
		// 			colored_nodes.push_back(i);	
		// 			if (can_add == 0)
		// 				break;
		// 		}
		// 	}

		// 	for (int u = start; u <= end; u++)
		// 	//for (int u = end; u >= start; u--)
		// 	{
		// 		int i = ranked[u - start];
		// 		int tnode = Candidate_Stack[i + start];
		// 		if (colored[i] || Node_State[tnode] == 0 || nbnn[i] >= color_ub)
		// 			continue;
		// 		int cnum = con_num2(tnode, colors[1].sz);
		// 		if (cnum <= color_ub - nbnn[i])
		// 		{
		// 			colors[1].bitcolor.set(target_normal[tnode]-1);
		// 			colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
		// 			colors[1].sz++;	
		// 			colored_nodes.push_back(i);	
		// 		}
		// 	}
		// 	//insert_num = insert.size();

		// 	// for (int u = start; u <= end; u++)
		// 	// {
		// 	// 	int i = ranked[u - start];
		// 	// 	int tnode = Candidate_Stack[i + start];
		// 	// 	if (colored[i] || Node_State[tnode] == 0 || nbnn[i] >= color_ub)
		// 	// 		continue;
		// 	// 	int cnum = con_num2(tnode, colors[1].sz);
		// 	// 	if (cnum <= color_ub - nbnn[i])
		// 	// 	{
		// 	// 		colors[1].bitcolor.set(target_normal[tnode]-1);
		// 	// 		colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
		// 	// 		colors[1].sz++;	
		// 	// 		colored_nodes.push_back(i);	
		// 	// 	}
		// 	// }
		// }
		// else
		// {
		// 	for (int u = start; u <= end; u++)
		// 	//for (int u = end; u >= start; u--)
		// 	{
		// 		int i = ranked[u - start];
		// 		int tnode = Candidate_Stack[i + start];
		// 		if (colored[i] || Node_State[tnode] == 0 || nbnn[i] >= color_ub)
		// 			continue;
		// 		int cnum = con_num2(tnode, colors[1].sz);
		// 		if (cnum <= color_ub - nbnn[i])
		// 		{
		// 			colors[1].bitcolor.set(target_normal[tnode]-1);
		// 			colors[1].vertices[colors[1].sz] = target_normal[tnode]-1;
		// 			colors[1].sz++;	
		// 			colored_nodes.push_back(i);	
		// 		}
		// 	}
		// }

		for (ui i = 0; i < color_num_tmp; i++)
			colored[colored_nodes[i]] = false;
	}
}

static int cut_by_iteration_partition_dise3(){
  //step 0: Mark the candidate vertices
  	int lb=MAX_KPX_SIZE - CUR_KPX_SIZE;
	// if(CUR_KPX_SIZE<5 || lb<PARA.KX)
    // 	return cut_by_iteration_partition_SEC();
  	int end=ptr(Candidate_Stack)-2,start=end;
  	for(int pnode=Candidate_Stack[start];pnode!=DELIMITER;pnode=Candidate_Stack[--start]){
		//cout << pnode << " ";
    	Node_State.set(pnode);
		Node_State3.set(pnode);
  	}
	//cout << endl;
  	start++;
	
	int out_num = end + 1 - start;
	int r_num = out_num;
		
  	//step 1: multiset-partition w.r.t the solution
	int idx=ptr(Candidate_Stack);
	for(int i=0;i<ptr(Clique_Stack);i++){
		Index_Stack[i]=idx;
		int ub1=PARA.KX-1 - Energy_Stack[i];
		if(ub1==0)continue;		
		int cnode=Clique_Stack[i];
		int *neibors=Node_Neibors[cnode];
		//clear marks of neighbors of cnode
		for(int neibor = *neibors; neibor!=cnode; neibor = *(neibors+=2)){
			Node_State.reset(neibor);
		}
		//all marked nodes are non-neighbors
		for(int j=end,pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
			if(Node_State[pnode]){
				Candidate_Stack[idx++]=pnode;
			}else{
				Node_State.set(pnode);
			}
		}
		Candidate_Stack[idx++]=NONE;
	}	
	// maximum profit:try to use lb cost to buy as most as possible vertices 

	memset(nondegree_temp, 0, sizeof(int) * (NB_NODE+1));
	for(int i = start; i <= end; i++){
		int pnode1 = Candidate_Stack[i];
		int *neibors=Node_Neibors[pnode1];
		for(int neibor = *neibors; neibor!=pnode1; neibor = *(neibors+=2)){
			Node_State3.reset(neibor);
		}
		int ndeg = 0;
		//all marked nodes are non-neighbors
		for (int j = start; j <= end; j++){
			int pnode2 = Candidate_Stack[j];
			if (Node_State3[pnode2]){
				// nondegree_temp[target_normal[pnode1]]++;
				// nondegree_temp[target_normal[pnode2]]++;
				ndeg++;
			}
			else
				Node_State3.set(pnode2);
		}
		nondegree_temp[target_normal[pnode1]] = ndeg;
	}

	for (int i = start; i <= end; i++){
		ranked[i - start] = i - start;
		nbnn[i - start] = PARA.KX - NBNN_Stack[i];
		//profits[i - start] = (float) (nondegree_temp[target_normal[Candidate_Stack[i]]]) / nbnn[i - start];
		profits[i - start] = (float) (nondegree_temp[target_normal[Candidate_Stack[i]]]) / (nbnn[i - start] - 1 + 0.00000001);
	}
	qsort_profits(ranked, 0, out_num - 1);
	// for (int i = start; i <= end; i++)
	// 	cout << profits[ranked[i - start]] << " ";
	// cout << endl;

	//step 2: abstract the pub and color partition with the most profit
	float last_pub_profit = out_num;
	while (lb > 0){
		int color_ub = 0;
		colors[1].sz = 0;
		colors[1].bitcolor.clear();
		try_color(color_ub, lb, start, end);
		float color_profit_rate = (float) (colors[1].sz) / color_ub;
		//float color_rrate = (float) (r_num - colors[1].sz - insert_num) / (lb - color_ub);
		// if (color_profit_rate <= 1)
		// 	color_profit_rate = -10000;
		//color_profit_rate = -10000;
		//cout << color_ub << " " << colors[1].sz + insert_num << "   ";

		int saved_invest,investor_idx=-1,saved_count;
		float max_profit_rate = -1;		
		//float saved_rrate = -1;
		if (last_pub_profit > color_profit_rate){
			for(int i=0;i<ptr(Clique_Stack);i++){
				int ub1=PARA.KX-1 - Energy_Stack[i];   
				if(ub1==0) continue;
				idx=Index_Stack[i];
				//if(idx<0 || ub1 >= lb) continue;
				if(idx<0) continue;
				int cnode=Clique_Stack[i],count=0;
				for(int pnode=Candidate_Stack[idx++];pnode!=NONE;pnode=Candidate_Stack[idx++]){
					if(Node_State[pnode])count++;
				}				
				//save the maximum profile partition
				int invest=count>ub1? ub1:count;
				if(invest==0){
					Index_Stack[i]=-1;
					continue;
				}
				if(invest>lb)
					continue;      
				// assert(invest>0 && invest<=lb);      
				float profit_rate = (float)(count) / (invest);
				//float rrate = (float)(r_num - count) / (lb - invest);
				//if(profit_rate>max_profit_rate||(profit_rate==max_profit_rate && rrate<saved_rrate)){
				if(profit_rate>max_profit_rate||(profit_rate==max_profit_rate && invest>saved_invest)){
					saved_invest=invest;
					investor_idx=i;
					max_profit_rate=profit_rate;  
					saved_count = count;
					//saved_rrate = rrate;
				}	   
			}
			last_pub_profit = max_profit_rate;
		}
		//if (color_profit_rate > max_profit_rate || (color_profit_rate == max_profit_rate && color_rrate < saved_rrate))
		if (color_profit_rate > max_profit_rate || (color_profit_rate == max_profit_rate && color_ub > saved_invest))
		//if (color_profit_rate > max_profit_rate && color_rrate < saved_rrate)
		{
			for (int i = 0; i < colors[1].sz; i++){
				Node_State.reset(Candidate_Stack[colored_nodes[i] + start]);
			}
			// for (ui i = 0; i < insert_num; i++){
			// 	Node_State.reset(Candidate_Stack[insert[i] + start]);
			// }
			lb -= color_ub;
			r_num -= (colors[1].sz);
			if (r_num <= lb)
				return start;
		}
		else{
			if (investor_idx < 0)
				break;
			for(int j=Index_Stack[investor_idx],pnode=Candidate_Stack[j];pnode!=NONE;pnode=Candidate_Stack[++j]){
				Node_State.reset(pnode);
			}
			Index_Stack[investor_idx]=-1;
			lb -= saved_invest;
			r_num -= saved_count;
			if (r_num <= lb)
				return start;
		}	
	}

	// cout << "A: ";
	// int num = 0;
	// for (int i = 0; i < out_num; i++) {
	// 	int u = ranked[out_num - 1 - i];
	// 	//int u = i;
	// 	if (Node_State[Candidate_Stack[u + start]]){
	// 		num++;
	// 		cout << Candidate_Stack[u + start] << " ";
	// 		if (num == r_num){
	// 			cout << endl;
	// 			break;
	// 		}
	// 	}
	// }

	int branching_idx = start - 1;
	r_num -= lb;
	for (int i = 0; i < out_num; i++) {
		int u = ranked[out_num - 1 - i];
		//int u = i;
		if (Node_State[Candidate_Stack[u + start]]){
			//Node_State3.reset(Candidate_Stack[u + start]);
			Node_State.reset(Candidate_Stack[u + start]);
			branching_idx++;
			if (u + start != branching_idx){
				swap(Candidate_Stack[u + start], Candidate_Stack[branching_idx]);
				swap(NBNN_Stack[u + start], NBNN_Stack[branching_idx]);
				i--;
			}
			if (branching_idx == start - 1 + r_num){
				// for (int j = i + 1; j < out_num; j++) {
				// 	int v = ranked[out_num - 1 - j];
				// 	if (Node_State[Candidate_Stack[v + start]])
				// 		Node_State.reset(Candidate_Stack[v + start]);
				// }
				break;
			}
		}
	}

	// cout << "B: ";
	// for (int i = 0; i < r_num; i++) {
	// 	int u = i + start;
	// 	cout << Candidate_Stack[u] << " ";
	// }
	// cout << endl;

	return branching_idx + 1;

	branching_idx = start - 1;
	for (int i = start; i <= end; i++) {
		int pnode = Candidate_Stack[i];
		if (Node_State3[pnode] == 0) {
			branching_idx++;
			swap(Candidate_Stack[i], Candidate_Stack[branching_idx]);
			swap(NBNN_Stack[i], NBNN_Stack[branching_idx]);
			if (branching_idx == start - 1 + r_num)
				break;
		}
	}

	return branching_idx + 1;

	// for (ui i = 0; i < out_num; i++) {
	// 	ui u = ranked[out_num - 1 - i];
	// 	//ui v = stkvtx[u + CASTART(depth)];
	// 	if (has_include[u] == false) {
	// 		pos--;
	// 		if (u + CASTART(depth) < pos)
	// 			swap(stkvtx[u + CASTART(depth)], stkvtx[pos]);
	// 		if (pos == CAEND(depth) - r_num)
	// 			break;
	// 	}
	// }

	// ui brunching_idx = R_end;
	// r_num -= lb;
	// for (ui i = 0; i < out_num; i++){
	// 	ui u = ranked[out_num - 1 - i];
	// 	// wish put u to tail
	// 	if (has_include[u] == 0){
	// 		brunching_idx--;
	// 		if (u < brunching_idx)
	// 			swap_pos(u, brunching_idx);
	// 		if (brunching_idx == R_end - r_num)
	// 			break;
	// 	}
	// }
	// return brunching_idx;

	int j=end,p=end,temp;
	for(int pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
		if(Node_State[pnode]==0){
		if(j<p){
			Candidate_Stack[j]=Candidate_Stack[p];
			Candidate_Stack[p]=pnode;
			temp=NBNN_Stack[j];
			NBNN_Stack[j]=NBNN_Stack[p];
			NBNN_Stack[p]=temp;
		}
		p--;
		}
		else{
			Node_State.reset(pnode);
		}
	}
	p++; 

	//if (lb != 0)
	if (branching_idx + 1 != p - lb)
		cout << start << " " << end << " " << branching_idx + 1 << " " << p - lb << " " << lb << endl;

	return (p-lb < start)? start:p-lb;
}



static int cut_by_iteration_partition_dise(){
  //step 0: Mark the candidate vertices
  int lb=MAX_KPX_SIZE - CUR_KPX_SIZE; assert(lb>=0);
  if(CUR_KPX_SIZE<5 || lb<PARA.KX)
    return cut_by_iteration_partition_SEC();
  int end=ptr(Candidate_Stack)-2,start=end;
  for(int pnode=Candidate_Stack[start];pnode!=DELIMITER;pnode=Candidate_Stack[--start]){
    Node_State.set(pnode);
  }
  start++;
  //step 1: multiset-partition w.r.t the solution
  int idx=ptr(Candidate_Stack);
  for(int i=0;i<ptr(Clique_Stack);i++){
    Index_Stack[i]=idx;
    int ub1=PARA.KX-1 - Energy_Stack[i];
    assert(ub1>=0);
    if(ub1==0)continue;		
    int cnode=Clique_Stack[i];
    int *neibors=Node_Neibors[cnode];
    //clear marks of neighbors of cnode
    for(int neibor = *neibors; neibor!=cnode; neibor = *(neibors+=2)){
      Node_State.reset(neibor);
    }
    //all marked nodes are non-neighbors
    for(int j=end,pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
      if(Node_State[pnode]){
	Candidate_Stack[idx++]=pnode;
      }else{
	Node_State.set(pnode);
      }
    }
    Candidate_Stack[idx++]=NONE;
  }	
  // maximum profit:try to use lb cost to buy as most as possible vertices 
	
  //step 2: abstract the partition with the most profit
  while(lb>0){
    int saved_idx,saved_invest,investor_idx=-1;
    float max_profit_rate=-1;		
    for(int i=0;i<ptr(Clique_Stack);i++){
      int ub1=PARA.KX-1 - Energy_Stack[i];      
      if(ub1==0) continue;
      idx=Index_Stack[i];
      if(idx<0) continue;
	   
      int cnode=Clique_Stack[i],count=0;
      for(int pnode=Candidate_Stack[idx++];pnode!=NONE;pnode=Candidate_Stack[idx++]){
	if(Node_State[pnode])count++;
      }				
      //save the maximum profile partition
      int invest=count>ub1? ub1:count;
      if(invest==0){
	Index_Stack[i]=-1;
	continue;
      }
      if(invest>lb)
	continue;      
      assert(invest>0 && invest<=lb);      
      float profit_rate=(count)/((float)invest);
      if(profit_rate>max_profit_rate||(profit_rate==max_profit_rate && invest>saved_invest)){
	saved_invest=invest;
	investor_idx=i;
	max_profit_rate=profit_rate;  
      }	   
    }
    if(investor_idx<0)break;
   for(int j=Index_Stack[investor_idx],pnode=Candidate_Stack[j];pnode!=NONE;pnode=Candidate_Stack[++j]){
     Node_State.reset(pnode);
   }
  Index_Stack[investor_idx]=-1;
      lb=lb-saved_invest;

  }
  int j=end,p=end,temp;
  for(int pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
    if(Node_State[pnode]==0){
      if(j<p){
	Candidate_Stack[j]=Candidate_Stack[p];
	Candidate_Stack[p]=pnode;
	temp=NBNN_Stack[j];
	NBNN_Stack[j]=NBNN_Stack[p];
	NBNN_Stack[p]=temp;
      }
      p--;
    }else{
      Node_State.reset(pnode);
    }
  }
  p++; 
  return (p-lb < start)? start:p-lb;
}

static int cut_by_iteration_partition_SEC(){
  int ub=0,k=ptr(Candidate_Stack)-1;
  int lb=MAX_KPX_SIZE - CUR_KPX_SIZE;
  int min_k;
  assert(Candidate_Stack[k]==DELIMITER);
  int total=0;
  for(int j=k-1,pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
    Node_State.set(pnode);
    total++;
  }
  int cut=0;
#ifdef INVE
  for(int i=ptr(Clique_Stack)-1;i>=0;i--){
#else
    for(int i=0;i<ptr(Clique_Stack);i++){
#endif
      int cnode=Clique_Stack[i];
      int ub1=PARA.KX-1 - Energy_Stack[i];
      assert(ub1>=0);
      if(ub1==0)continue;
      int j=k-1,p=k,count=0;

      int neibor,*neibors = Node_Neibors[cnode];
      for (neibor = *neibors; neibor != cnode; neibor = *(neibors+=2)) {
	Node_State.reset(neibor);
      }
      for(int pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
	if(Node_State[pnode]){
	  if(j<p-1){
	    int temp=Candidate_Stack[p-1];
	    Candidate_Stack[p-1]=pnode;
	    Candidate_Stack[j]=temp;
	    temp=NBNN_Stack[p-1];
	    NBNN_Stack[p-1]=NBNN_Stack[j];
	    NBNN_Stack[j]=temp;
	  }
	  p--;
	  count++;
	}else{
	  Node_State.set(pnode);
	}
      }
	
      assert(j<p);
      assert(Candidate_Stack[j]==DELIMITER);
			
      min_k=j;
      if(count<ub1)
	ub1=count;
			
      if(ub+ub1<=lb){
	cut+=count;
	ub=ub+ub1;
	k=p;
	if(k==min_k+1||ub==lb)
	  break;
      }else{
	for(int j=ptr(Candidate_Stack)-2,pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
	  Node_State.reset(pnode);
	}
	return k-(lb-ub) > min_k ? k-(lb-ub): min_k+1;
      }
    }
    for(int j=ptr(Candidate_Stack)-2,pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
      Node_State.reset(pnode);
    }
    if(k==min_k+1){
      return k;
    }else{
      assert(k>min_k+1);
      int x= k-(lb-ub) > min_k ? k-(lb-ub): min_k+1;
      return x;
    }
  }

static int cut_by_iteration_partition(){
	int ub=0,k=ptr(Candidate_Stack)-1;
	int lb=MAX_KPX_SIZE - CUR_KPX_SIZE;
	int min_k;
	assert(Candidate_Stack[k]==DELIMITER);
	
#ifdef INVE
	for(int i=ptr(Clique_Stack)-1;i>=0;i--){
#else
           for(int i=0;i<ptr(Clique_Stack);i++){
#endif
			int cnode=Clique_Stack[i];
			int ub1=PARA.KX-1 - Energy_Stack[i];
			assert(ub1>=0);
			if(ub1==0)continue;
			int j=k-1,p=k,count=0;
			for(int pnode=Candidate_Stack[j];pnode!=DELIMITER;pnode=Candidate_Stack[--j]){
				if(is_adjacent(cnode,pnode)==FALSE){
					if(j<p-1){
						int temp=Candidate_Stack[p-1];
						Candidate_Stack[p-1]=pnode;
						Candidate_Stack[j]=temp;
						temp=NBNN_Stack[p-1];
						NBNN_Stack[p-1]=NBNN_Stack[j];
						NBNN_Stack[j]=temp;
					}
					p--;
					count++;
				}
			}
			assert(j<p);
			assert(Candidate_Stack[j]==DELIMITER);
			
			min_k=j;
			
			if(count<ub1)
				ub1=count;
			
			if(ub+ub1<=lb){
				ub=ub+ub1;
				k=p;
				if(k==min_k+1)
					break;
			}else{
				return k-(lb-ub) > min_k ? k-(lb-ub): min_k+1;
			}
		}
		if(k==min_k+1){
			return k;
		}
		else{
			assert(k>min_k+1);
			int x= k-(lb-ub) > min_k ? k-(lb-ub): min_k+1;
			return x;
		}
	}

	static int cut_by_common_neibor(int start){
		int i,j,neibor, max = 0, *neibors;
		int bnode=Candidate_Stack[CURSOR];
		
		if(CUR_KPX_SIZE>2){
			for(int cn=Candidate_Stack[i=start];cn!=DELIMITER;cn=Candidate_Stack[++i]){
				Node_State2.reset(cn);
				Node_State.reset(cn);
			}
			for(i=0;i<ptr(Clique_Stack);i++){
				if(Energy_Stack[i]<0){
					Energy_Stack[i]=-Energy_Stack[i];
				}
			}
			return ptr(Candidate_Stack)-1-(MAX_KPX_SIZE-(CUR_KPX_SIZE));
		}
		
		neibors = Node_Neibors[bnode];
		for (neibor = *neibors; neibor != bnode; neibor = *(neibors+=2)) {
			if(Node_State2[neibor])
				Node_State.set(neibor);
		}
		
		for(int cn=Candidate_Stack[i=start];cn!=DELIMITER;cn=Candidate_Stack[++i]){
			CNN[cn]=0;
			
			neibors = Node_Neibors[cn];
			for (neibor = *neibors; neibor != cn; neibor = *(neibors+=2)) {
				if(Node_State[neibor])
					CNN[cn]++;
			}
			
			int flag= Node_State[cn] ? 0:1;
			int min_NN=MAX_KPX_SIZE - (CUR_KPX_SIZE+1)+1 - CNN[cn];
			
			if(NBNN_Stack[CURSOR] + NBNN_Stack[i]+min_NN+flag > 2*PARA.KX-2){
				Node_State.reset(cn);
				Candidate_Stack[i]=-cn;
				NB_CANDIDATE--;
			}
		}
		
		for(int cn=Candidate_Stack[i=start];cn!=DELIMITER;cn=Candidate_Stack[++i]){
			if(cn>0){
				Node_State2.reset(cn);
				Node_State.reset(cn);
			}else{
				Node_State2.reset(-cn);
				Node_State.reset(-cn);
			}
		}
		
		if(NB_CANDIDATE+CUR_KPX_SIZE>MAX_KPX_SIZE){
			int j=start,count=0;
			float cutted=0.0;
			for(int cn=Candidate_Stack[i=start];cn!=DELIMITER;cn=Candidate_Stack[++i]){
				if(cn>0){
					Candidate_Stack[j]=cn;
					NBNN_Stack[j]=NBNN_Stack[i];
					count++;
					j++;
				}else{
					cutted++;
					//removed.push_back(-cn);
				}
			}
			Candidate_Stack[j]=DELIMITER;
			ptr(Candidate_Stack)=j+1;
			NBNN_Stack[j]=DELIMITER;
			ptr(NBNN_Stack)=j+1;
			assert(count==NB_CANDIDATE);
			
			for(i=0;i<ptr(Clique_Stack);i++){
				if(Energy_Stack[i]<0){
					Energy_Stack[i]=-Energy_Stack[i];
				}
			}
			return ptr(Candidate_Stack)-1-(MAX_KPX_SIZE-(CUR_KPX_SIZE));
		}else{
			for(int i=0;i<ptr(Clique_Stack);i++){
				if(Energy_Stack[i]<0){
					Energy_Stack[i]=-Energy_Stack[i]-1;
				}
			}
			return 0;
		}
	}
	
	
	static void init_for_search() {
		int i, node;
		cut_ver = 0;
		cut_inc = 0;
		cut_iset = 0;
		cut_satz = 0;
		total_cut_ver = 0;
		total_cut_inc = 0;
		total_cut_iset = 0;
		total_cut_satz = 0;
		
		Branches_Nodes[0] = 0;
		Branches_Nodes[1] = 0;
		Branches_Nodes[2] = 0;
		Branches_Nodes[3] = 0;
		Branches_Nodes[4] = 0;
		Branches_Nodes[5] = 0;
		
		Last_Idx = NB_NODE;
		NB_BACK_CLIQUE = 0;
		MAX_KPX_SIZE = 0;
		ptr(Clique_Stack) = 0;
		ptr(Energy_Stack) =0;
		ptr(Cursor_Stack) = 0;
		Rollback_Point = 0;
		Reasoning_Point=0;
		
		
		// push(NB_NODE - INIT_KPX_SIZE - 1, Cursor_Stack);
		push(NB_NODE - INIT_KPX_SIZE, Cursor_Stack);
		MAX_KPX_SIZE = INIT_KPX_SIZE;
		
		//	assert(Candidate_Stack[ptr(Candidate_Stack)- 2]==1);
		assert(ptr(NBNN_Stack)==ptr(Candidate_Stack));
		
		for (i = 0; i < ptr(Candidate_Stack) - 1; i++) {
			node = Candidate_Stack[i];
			//	Vertex_UB[i] = Node_Degree[node] + PARA.KX;
			NBNN_Stack[i]=0;
		      
		}
		
		NBNN_Stack[NB_NODE]=DELIMITER;
		ptr(NBNN_Stack)=ptr(Candidate_Stack);
	}
	
	static void allocate_memory() {
		Second_Name = (int *) malloc((MAX_VERTEX_NO + 1) * sizeof(int));
		
		CNN=(int *)malloc((MAX_VERTEX_NO+1)*sizeof(int));
		Node_Count=(int *)malloc((MAX_VERTEX_NO+1)*sizeof(int));
	}
	

static void bnb_search(int cutoff, int silent) {
		//removed.clear();
		int bnode;
		init_for_search();
		BRANCHING_COUNT = 0;
		if (silent==0) {
			printf("C  -------------------------------------------------------------------------------\n");
			printf("C  Size|   Index|   NB_Branches|   Time(s)| R-dep|   Nodes|    Edges|    Density|\n");
		}
		while (CURSOR> 0) {
			//if(CUT_OFF>0 && get_utime()> CUT_OFF){
			if(get_utime()> 2000){
				TIME_OUT=TRUE;
				break;
			}
			bnode=Candidate_Stack[--CURSOR];
#ifdef ISET
			if(Reasoning_Point && bnode>0)
				continue;
#endif
			
			if(bnode==DELIMITER) {
				ptr(NBNN_Stack)=ptr(Candidate_Stack)=CURSOR+1;
				ptr(Cursor_Stack)--;
				ptr(Clique_Stack)--;
				ptr(Energy_Stack)--;
				int pop_node=Clique_Stack[ptr(Clique_Stack)];

				for(int i=0;i<ptr(Clique_Stack);i++){
					if(is_adjacent(Clique_Stack[i],pop_node)==FALSE){
						Energy_Stack[i]--;
					}
				}
				if(pop_node==Reasoning_Point){
					Reasoning_Point=0;
				}
			} else {
				if(MAX_KPX_SIZE==CUR_KPX_SIZE) {
					store_maximum_clique(bnode,silent);
				
				}else {
					BRANCHING_COUNT++;
					Rollback_Point=ptr(Candidate_Stack);
					push(bnode, Clique_Stack);					
					push(NBNN_Stack[CURSOR],Energy_Stack);
					
#ifdef COMM
					if((Branching_Point=produce_subgraph0())==0
					   || (Branching_Point=cut_by_common_neibor(Rollback_Point))==0)
#elif  defined(PART)
						if((Branching_Point=produce_subgraph0())==0
						   || (Branching_Point=cut_by_iteration_partition())==0
						   )
#elif defined(PART11)
						 if((Branching_Point=produce_subgraph0())==0
						   || (Branching_Point=cut_by_iteration_partition_SEC())==0
						   )						 
#elif defined(SEC)
						   if((Branching_Point=produce_subgraph0())==0
						      || (Branching_Point=cut_by_common_neibor(Rollback_Point))==0
						   || (Branching_Point=cut_by_iteration_partition_SEC())==0
						   )
#elif defined(GOP)
						   if((Branching_Point=produce_subgraph0())==0
						      || (Branching_Point=cut_by_common_neibor(Rollback_Point))==0
						   || (Branching_Point=cut_by_iteration_partition_dise3())==0
						   )
#elif defined(ALL)
							if((Branching_Point=produce_subgraph0())==0
							   || (Branching_Point=cut_by_common_neibor(Rollback_Point))==0
							   || (Branching_Point=cut_by_iteration_partition())==0
							   )
#else
								if((Branching_Point=produce_subgraph0())==0)
#endif
									{
									ptr(NBNN_Stack)=ptr(Candidate_Stack)=Rollback_Point;
									ptr(Clique_Stack)--;
									ptr(Energy_Stack)--;
									continue;
									}
					
					push(Branching_Point,Cursor_Stack);
					//check_solution2();
				}
			}
		}
		
		SEARCH_TIME = get_utime();
		if (silent==0) {
			/*printf(
			 "C  -----------------------------------------------------------------------------\n");
			 printf("C %4d |%7d |%8d %10d %10d %10d|%14lld %8.2lf\n", MAX_KPX_SIZE, CURSOR,cut_ver,cut_inc, cut_iset, cut_satz,BRANCHING_COUNT,SEARCH_TIME);
			 total_cut_ver += cut_ver;
			 total_cut_inc += cut_inc;
			 total_cut_iset += cut_iset;
			 total_cut_satz += cut_satz;*/
			printf(
				   "C --------------------------------------------------------------------------------\n");
			printf("C %4d |%7d |%14lld| %8.2lf\n", MAX_KPX_SIZE, CURSOR,BRANCHING_COUNT,SEARCH_TIME);
		}
	}
	
	static int * Adj_List;
	
#define New_Name Node_Degree
	
	
	static void free_block() {
		int i = 0;
		for (i = 0; i < BLOCK_COUNT; i++)
			free(BLOCK_LIST[i]);
	}


	static void compute_numbers_of_common_neighbors(){
	  int i,j,p1,p2,p3,node, *neibors1,*neibors2,*ptr;
	  printf("computing the numbers of common neighbors....\n");
	 
	  //list all triangles <p1>p2>p3>, compute common neibors for each edge
	  for(p1= Candidate_Stack[i=0];p1!=DELIMITER; p1= Candidate_Stack[++i]){
	    for(neibors1 = Node_Neibors[p1],p2 = *neibors1; p2<p1; p2 = *(neibors1+=2)) {
	      Node_State.set(p2);
	    }			 
	    for(neibors1 = Node_Neibors[p1],p2 = *neibors1; p2<p1; p2 = *(neibors1+=2)) {
	      assert(Node_State[p2]);
	      //assert(Node_Degree[p2]+PARA.KX > MAX_KPX_SIZE);			   
	      for (neibors2=Node_Neibors[p2],p3= *neibors2; p3<p2; p3 = *(neibors2+=2)) {
		if(Node_State[p3]){
		  Node_State2.set(p3);
		  *(neibors1+1)=*(neibors1+1)+1;
		  *(neibors2+1)=*(neibors2+1)+1;
		}
	      }			    
	      for (neibors2=Node_Neibors[p1],p3= *neibors2;p3< p2; p3 = *(neibors2+=2)) {
		if(Node_State2[p3]){
		  *(neibors2+1)=*(neibors2+1)+1;
		  Node_State2.reset(p3);
		}
	      }
	    }			  
	    for(neibors1 = Node_Neibors[p1],p2 = *neibors1; p2<p1; p2 = *(neibors1+=2)) {
	      Node_State.reset(p2);
	    }
	  }		   
	}

	static inline void update_cnn(int p1,int p2){
	  int pi,pj;
	  int *neibors1,*neibors2;
	   if(Touched[p1]==0){
	     Candidate_Stack[ptr(Candidate_Stack)-1]=p1;
	     push(DELIMITER,Candidate_Stack);
	     Touched[p1]=1;
	   }
	   if(Touched[p2]==0){
	     Candidate_Stack[ptr(Candidate_Stack)-1]=p2;
	     push(DELIMITER,Candidate_Stack);
	     Touched[p2]=1;
	   }
	  for(neibors1 = Node_Neibors[p1],pi = *neibors1; pi!=p1; pi = *(neibors1+=2)) {
	    if(pi>0 && !Removed[pi])
	      Node_State.set(pi);
	  }
	  for(neibors1 = Node_Neibors[p2],pi = *neibors1; pi!=p2; pi = *(neibors1+=2)) {
	    if(pi>0 && !Removed[pi] && Node_State[pi])
	      Node_State2.set(pi);
	  }
	  for(int i=0,px=p1;i<2;px=p2,i++){
	    for(neibors1 = Node_Neibors[px],pi = *neibors1; pi!=px; pi = *(neibors1+=2)) {
	      if(pi>0 && !Removed[pi] && Node_State2[pi]){
		if(Touched[pi]==0){
		  Candidate_Stack[ptr(Candidate_Stack)-1]=pi;
		  push(DELIMITER,Candidate_Stack);
		  Touched[pi]=1;
		}		  
		if(pi<px){
		  *(neibors1+1)=*(neibors1+1)-1;
		}else{
		  for(neibors2 = Node_Neibors[pi],pj = *neibors2; pj<pi; pj = *(neibors2+=2)){
		    if(pj==px){
		      *(neibors2+1)=*(neibors2+1)-1;
		      break;
		    }
		  }
		}
	      }
	    }
	  }
	  for(neibors1 = Node_Neibors[p1],pi = *neibors1; pi!=p1; pi = *(neibors1+=2)) {
	    if(pi>0 && !Removed[pi]){
	      Node_State.reset(pi);
	      Node_State2.reset(pi);
	    }
	  }
	}

	static void reduce_instance_with_unsupport_property(){
	  int level=0;
	  int i,j,p1,p2,p3,node, *neibors1,*neibors2,*ptr;		
	  Extra_Node_Stack=&Candidate_Stack[NB_NODE+2];	
	  do{
	    printf("\b\b\b\b\b%3d |",++level);fflush(stdout);
	    for(p1= Candidate_Stack[i=0],j=0;p1!=DELIMITER; p1= Candidate_Stack[++i]){
	      if(!Removed[p1]){
		Candidate_Stack[j++]=p1;
		Touched[p1]=1;
	      }
	    }
	    ptr(Candidate_Stack)=j+1;
	    Candidate_Stack[j]=DELIMITER;
	    
	    for(p1= Candidate_Stack[i=0];p1!=DELIMITER; p1= Candidate_Stack[++i]){
	      Touched[p1]=0;
	      for (neibors1=Node_Neibors[p1],p2=*neibors1;(p2<p1)&&(p1+p2>0);p2=*(neibors1+=2)){
		if(p2>0 && !Removed[p2] && *(neibors1+1)+2*PARA.KX<=MAX_KPX_SIZE){	
		  Node_Degree[p1]--;
		  *neibors1=-p2;		 
		  for (neibors2=Node_Neibors[p2],p3 = *neibors2; p3 != p2; p3 = *(neibors2+=2)) {
		    if(p3==p1){
		      *neibors2=-p1;
		      Node_Degree[p2]--;
		      break;
		    }
		  }
		  update_cnn(p1,p2);
		}
	      }
	    }
	    Candidate_Stack[j]=DELIMITER;
			  
	    ptr(Extra_Node_Stack)=0;
	    for(p1= Candidate_Stack[i=0];p1!=DELIMITER; p1= Candidate_Stack[++i]){
	      if(Node_Degree[p1]+PARA.KX<=MAX_KPX_SIZE){			   
		push(p1,Extra_Node_Stack);
	      }
	    }
	    push(DELIMITER,Extra_Node_Stack);
			
	    //removed vertices recursively
	    for(p1=Extra_Node_Stack[i=0];p1!=DELIMITER;p1=Extra_Node_Stack[++i]){
	      assert(p1>0);
	      assert(Node_Degree[p1]<=MAX_KPX_SIZE-PARA.KX);
	      assert(Removed[p1]==0);
	      Removed[p1]=1;			
	      for (neibors1=Node_Neibors[p1],p2 = *neibors1; p2 != p1; p2 = *(neibors1+=2)) {
		if(p2>0 && !Removed[p2] && Node_Degree[p2]+PARA.KX>MAX_KPX_SIZE){
		  Node_Degree[p2]--;
		  if(Node_Degree[p2]+PARA.KX<=MAX_KPX_SIZE){
		    Extra_Node_Stack[ptr(Extra_Node_Stack)-1]=p2;
		    push(DELIMITER,Extra_Node_Stack);
		  }
		}
	      }				
	    }
	    //update cnn after removing vertices
	    for(p1=Extra_Node_Stack[i=0];p1!=DELIMITER;p1=Extra_Node_Stack[++i]){
	      for (neibors1=Node_Neibors[p1],p2 = *neibors1; p2 != p1; p2 = *(neibors1+=2)) {
		if(p2>0 && !Removed[p2])
		  Node_State.set(p2);
	      }
	      for (neibors1=Node_Neibors[p1],p2 = *neibors1; p2 != p1; p2 = *(neibors1+=2)) {
		if(p2>0 && !Removed[p2]){
		  for (neibors2=Node_Neibors[p2],p3 = *neibors2;p3 <p2 && p3+p2>0; p3 = *(neibors2+=2)) {
		    if(p3>0 && Node_State[p3]){
		      *(neibors2+1)-=1;
		    }
		  }
		}			     
	      }
	      for (neibors1=Node_Neibors[p1],p2 = *neibors1; p2 != p1; p2 = *(neibors1+=2)) {
		if(p2>0 && !Removed[p2])
		  Node_State.reset(p2);
	      }
	    }
	  }while(ptr(Extra_Node_Stack)>1);
	}
	

	static int rebuild_instance(int pre){
		int i=0,j=0,nb_edge=0;
		int  node, *neibors, *neibors2, *addr;
		MAX_VERTEX_NO=0;
		int ret=FALSE;
		for(int p1= Candidate_Stack[i=0];p1!=DELIMITER; p1= Candidate_Stack[++i]){
			if(!Removed[p1]){
				Candidate_Stack[j++] = p1;
				Node_State.set(p1);
				Node_Degree[p1]=0;				
			}
		}
		NB_NODE = j;
		
		if(NB_NODE<=MAX_KPX_SIZE){
		  ret=TRUE;
		}
		
		NBNN_Stack[j]=Candidate_Stack[j] = DELIMITER;
		ptr(NBNN_Stack)=ptr(Candidate_Stack) = j + 1;
		
		NB_EDGE = 0;
		neibors2=Adj_List;
		for(int p1= Candidate_Stack[i=0];p1!=DELIMITER; p1= Candidate_Stack[++i]){
			neibors = Node_Neibors[p1];
			for (node = *neibors; node != p1; node = *(neibors+=2)) {			
				if(node>0 && Node_State[node]) {
					Node_Degree[p1]++;
					*neibors2++ = node;
					*neibors2++=*(neibors+1);
					NB_EDGE++;
				}
			}
			(*neibors2) = p1;
			neibors2++;
		}
		Node_Neibors[Candidate_Stack[0]]=Adj_List;
		for(int p0=Candidate_Stack[0],p1= Candidate_Stack[i=1];p1!=DELIMITER;p0=p1,p1= Candidate_Stack[++i]){
			Node_Neibors[p1]=Node_Neibors[p0]+ 2*Node_Degree[p0]+1;
			Node_State.reset(p1);
		}
		
		NB_EDGE=NB_EDGE/2;
		
		MAX_VERTEX_NO=Candidate_Stack[0];
		if(pre)
		     printf("after reducing #node=%d #edges=%d #density=%.8f\n", j,
		       NB_EDGE, NB_NODE<=1? 0:((float)2* NB_EDGE  / NB_NODE / (NB_NODE-1)));
		else
		printf("%8d|%9d| %.8f|\n", j,
		       NB_EDGE, NB_NODE<=1? 0:((float)2* NB_EDGE  / NB_NODE / (NB_NODE-1)));
		
		max_node_index = j - 1;

		return ret;
	}
	
	
	static int reduce_instance_with_degree() {
		
		int flag=1;
		int i=0,nb=0, node, *neibors, *neibors2, *addr;
		for (i = 0,nb=0; i < NB_NODE; i++) {
			node = Candidate_Stack[i];
			if (flag && Node_Degree[node]+PARA.KX <= INIT_KPX_SIZE) {
				Node_State[node] = PASSIVE;
			} else {
				flag=0;
				Candidate_Stack[nb++] = node;
				Node_State[node] = ACTIVE;
			}
		}
		NB_NODE = nb;
		
		if(NB_NODE<=INIT_KPX_SIZE){
			printf("I find the optimal solution in level-1 reduction.\n");
			return TRUE;
		}
		
		N0_1 = NB_NODE;
		NBNN_Stack[nb]=Candidate_Stack[nb] = DELIMITER;
		ptr(NBNN_Stack)=ptr(Candidate_Stack) = nb + 1;
		
		Old_Name = (int *)malloc((NB_NODE + 1) * sizeof(int));
		Removed=(char *)malloc((NB_NODE+1)*sizeof(char));
		Touched=(char *)malloc((NB_NODE+1)*sizeof(char));
		for (i = 0; i < NB_NODE; i++) {
			Old_Name[NB_NODE - i] = Candidate_Stack[i];
			New_Name[Candidate_Stack[i]] = NB_NODE - i;
			Candidate_Stack[i] = NB_NODE - i;
		}
		
		NB_EDGE = 0;
		for (i = NB_NODE; i > 0; i--) {
		  neibors = Node_Neibors[Old_Name[i]];
			neibors2 = neibors;
			nb = 0;
			for (node = *neibors; node != NONE; node = *(++neibors)) {
				if (Node_State[node]==ACTIVE) {
					(*neibors2) = New_Name[node];
					neibors2++;
					nb++;
				}
			}
			(*neibors2) = NONE;
			NB_EDGE += nb;
			qsort(Node_Neibors[Old_Name[i]], nb, sizeof(int), int_cmp_asc);
			assert(nb+PARA.KX > INIT_KPX_SIZE);
		}
		NB_EDGE=NB_EDGE/2;
		
		Adj_List = (int *) malloc((4*NB_EDGE + NB_NODE) * sizeof(int));
		addr = Adj_List;
		
		if (Adj_List == NULL ) {
			printf("can't allocate enough memory for Adj_List!\n");
			exit(1);
		}
		for (MAX_VERTEX_NO = 0,i = NB_NODE; i > 0; i--) {
		  Removed[i]=0;
		  Touched[i]=1;
		  Node_Degree[i] = 0;
			
		  Node_State.reset(i);
		  Node_State2.reset(i);
		  neibors = Node_Neibors[Old_Name[i]];
		  for (node = *neibors; node != NONE; node = *(++neibors)) {
		    assert(i!=node);
		    *(addr++) = node;
		    *(addr++) = 0;
		    Node_Degree[i]++;
		  }
		  *(addr++) = i;
		  if (Node_Degree[i] > MAX_VERTEX_NO)
		    MAX_VERTEX_NO = Node_Degree[i];
		  assert(Node_Degree[i]+PARA.KX> INIT_KPX_SIZE);
		}
		free_block();
		Node_Neibors[NB_NODE] = Adj_List;
		for (i = NB_NODE - 1; i > 0; i--) {
			Node_Neibors[i] = Node_Neibors[i + 1] + 2*Node_Degree[i + 1] + 1;
		}
		
		D1 = ((float) NB_EDGE*2  / NB_NODE / (NB_NODE - 1));
		printf("I the L1-Reduced graph #node %d #edge %d #density %9.8f\n", NB_NODE,
			   NB_EDGE, ((float) NB_EDGE*2  / NB_NODE / (NB_NODE - 1)));
		//	printf("I the largest subgraph is %d\n", MAX_VERTEX_NO);
		return FALSE;
	}
	
	static int initialize() {
	  int r = sort_by_degeneracy_ordering();
	  if (r == FALSE) {
	    r=reduce_instance_with_degree();
	    if(r==FALSE){
	      compute_numbers_of_common_neighbors();
	      reduce_instance_with_unsupport_property();
	      r=rebuild_instance(1);
	    }
			
	  }
	  INIT_TIME = get_utime() - READ_TIME;
	  if(r==FALSE)
	    compute_numbers_of_common_neighbors();
	  if(r==TRUE){
	    SEARCH_TIME=READ_TIME+INIT_TIME;
	    printf("Solution: ");
	    for(int i=0;i<MAX_KPX_SIZE;i++){
	      printf("%d ", MaxCLQ_Stack[i]);
	    }
	    printf("\n");
	  }
	  printf("I the initial time is %4.2lf \n", INIT_TIME);
	  return r;
	}
	
	static char * getInstanceName(char *s) {
		if (strrchr(s, '/') == NULL )
			return s;
		else
			return strrchr(s, '/') + 1;
	}

void print_compile_options(){
  printf("compiled at %s,%s\n",__TIME__,__DATE__);

#ifdef COMM
	  printf("compiled with --COMM");
#endif
		
#ifdef  GOP
	  printf("compiled with --GOP");
#endif

#ifdef  SEC
	  printf("compiled with --SEC");
	      
#endif
	  printf("\n");
}
	
	void parse_parameters(int argc,char* argv[]){
		PARA.ORDER =-1;
		PARA.FORMAT=1;
		PARA.START_MAXSAT_THD=15;
		PARA.CUT_OFF=0;
		PARA.TIME_OUT=0;
		for (int i = 2; i < argc; i++) {
			if (strcmp(argv[i],"-x")==0){
				sscanf(argv[++i], "%d", &PARA.KX);
			}else if (strcmp(argv[i], "-o") == 0){
				sscanf(argv[++i], "%d", &PARA.ORDER);
			}else if (strcmp(argv[i], "-f") == 0){
				sscanf(argv[++i], "%d", &PARA.FORMAT);
			} else if (strcmp(argv[i], "-i") == 0) {
				sscanf(argv[++i], "%d", &PARA.START_MAXSAT_THD);
			}else if (strcmp(argv[i], "-c") == 0) {
				sscanf(argv[++i], "%d", &PARA.CUT_OFF);
			}else if(strcmp(argv[i],"-m")==0){
				sscanf(argv[++i], "%f", &BMTHD);
			}else if(strcmp(argv[i],"-ms")==0){
				sscanf(argv[++i],"%d",&MAX_MATRIX_SIZE);
			}
		}
		strcpy(instance,argv[1]);
		printf("# Solving %d-plex in %s\n\n",PARA.KX,instance);
	}

	static void check_solution(){
		int count=0,count1=0,node1;
		if(MAX_KPX_SIZE>INIT_KPX_SIZE){
			for(int i=0; i<MAX_KPX_SIZE ; i++){
				count=0;count1=0;
				node1=MaxCLQ_Stack[i];
				for(int j=0;j<MAX_KPX_SIZE;j++){
					if(MaxCLQ_Stack[j]==node1)
						count1++;
					else if(is_adjacent(node1,MaxCLQ_Stack[j])==FALSE)
						count++;
				}
				if(count>PARA.KX-1)
					std::cout<<"count "<<count<<" "<<PARA.KX-1<<endl;
				assert(count1==1);
				assert(count<=PARA.KX-1);
				//assert(count==Energy_Stack[i]);
			}
		}
	}
	
	static void print_solution(){
		printf("Solution: ");
		for(int i=0;i<MAX_KPX_SIZE;i++){
			if(MAX_KPX_SIZE>INIT_KPX_SIZE){
				printf("%d ",Old_Name[MaxCLQ_Stack[i]]);
			}else{
				printf("%d ",MaxCLQ_Stack[i]);
			}
		}
		printf("\n");
	}

	static void read_color_graph(){
		colors = new ColorType[2];
		colors[1].bitcolor = MBitSet64(max_node_index + 1);
		colors[1].vertices = new int[NB_NODE];
		colors[1].sz = 0;	
		target_normal = new int[Candidate_Stack[0]+1];
      	target_abnormal = new int[NB_NODE+1];
		int normal_node = 1;
		//cout << max_node_index << endl;
		for(int i = max_node_index; i >= 0; i--){
			int temp_n = Candidate_Stack[i];
			target_normal[temp_n] = normal_node;
			target_abnormal[normal_node] = temp_n;
			normal_node++;
			//cout << normal_node << " " << temp_n << "    ";
		}
		kernalG.nbvtx = NB_NODE; kernalG.nbedges = 2*NB_EDGE;
		kernalG.pstart = new int[NB_NODE+1];
		kernalG.edges = new unsigned int[2*NB_EDGE+2];
		kernalG.pstart[0] = 0;
		for (int i = 0; i < NB_NODE; i++)
			kernalG.pstart[i+1] = kernalG.pstart[i] + Node_Degree[target_abnormal[i+1]];
		int edge_num = 0;
		//cout<<edge_num<<" "<<NB_NODE<<endl;
		for (int i = 1; i <= NB_NODE; i++){
			int neibor, *neibors = Node_Neibors[target_abnormal[i]];
			for (neibor = *neibors; neibor!=target_abnormal[i]; neibor = *(neibors+=2)){
				// if (i == NB_NODE)
				// 	cout << neibor << " ";
				kernalG.edges[edge_num++] = target_normal[neibor]-1;
			}
		}
		cout<<edge_num<<" "<<NB_NODE<<endl;
		bg = new MBitGraph(kernalG);

		ranked = new int[NB_NODE+1];
		nbnn = new int[NB_NODE+1];
		profits = new float[NB_NODE+1];
		colored = new bool[NB_NODE+1];
		memset(colored, false, sizeof(bool) * (NB_NODE+1));
		conflicted = new bool[NB_NODE+1];
		nondegree_temp = new int[NB_NODE+1];
		nbnns = new int[NB_NODE+1];
		nbnns2 = new int[NB_NODE+1];
	} 

	static void free_memory(){
		delete[] target_normal;
		delete[] target_abnormal;
		delete[] colors;
		//delete[] bg; 

		delete[] ranked;
		delete[] nbnn;
		delete[] profits;
		delete[] colored;
		delete[] conflicted;
		delete[] nondegree_temp;
		delete[] nbnns;
		delete[] nbnns2;
	}

	int main(int argc, char *argv[]) {
	  print_compile_options();
	  parse_parameters(argc,argv);
	  if (build_simple_graph_instance(argv[1])) {
	    if (initialize() == FALSE) {
			if (NB_NODE >= 500000){
				// fprintf(pp,"%s\t%d\t0\tOut of Memory!\n",argv[1],PARA.KX);
				// fclose(pp);
				return 0;
			}
	      allocate_memory();
		  read_color_graph();
	      bnb_search(0, 0);
	      print_solution();
		  free_memory();
	    }
	  }
	//   if (TIME_OUT != TRUE)
	//   	fprintf(pp,"%s\t%d\t1\t%d\t%lld\t%0.3f\n",argv[1],PARA.KX,MAX_KPX_SIZE,BRANCHING_COUNT,SEARCH_TIME);
	//   else
	// 	fprintf(pp,"%s\t%d\t0\t%d\t%lld\t%0.3f\n",argv[1],PARA.KX,MAX_KPX_SIZE,BRANCHING_COUNT,SEARCH_TIME);
	// 	fclose(pp);
	    printf(
		 ">>%s @ %s |V| %d |E| %d K %d MaxKPX %d InitKPX %d Tree %lld Read_Time %4.2lf Init_Time %4.2lf Search_Time %4.2lf Total %4.2lf \\\\\n",
		 argv[1],getInstanceName(argv[1]),NB_NODE_O, NB_EDGE_O, PARA.KX,
		 MAX_KPX_SIZE, HEUR_KPX_SIZE, BRANCHING_COUNT, READ_TIME, INIT_TIME,
		 SEARCH_TIME - READ_TIME - INIT_TIME, SEARCH_TIME);	  
	    		
	  return 0;
	}
	
