#include "Graph.h"
#include "Utility.h"
#include "Timer.h"
#include <signal.h>

using namespace std;
#define LEN_LIMIT (1<<10)
char filename[LEN_LIMIT];

// FILE* fp;

static void SIGINT_exit(int signum) {
	printf("Time Out!\n");
	//fp = fopen("result.txt","a+");
	// fprintf(fp,"0\tTimeOut!\n");
	// fclose(fp);
	exit(1);
}

int main(int argc, char *argv[]) {
	Timer TT;
	// fp = fopen("result.txt","a+");
	// fprintf(fp,"%s\t%d\t",argv[1],atoi(argv[2]));
	signal(SIGTERM, SIGINT_exit);
	signal(SIGINT, SIGINT_exit);
	branching_num = 0;
	color_time = 0;
	pub_time = 0;

	printf("\n-----------------------------------------------------------------------------------------\n");
	
	if (argc == 3) {
		strncpy(filename, argv[1], LEN_LIMIT);
		int k = atoi(argv[2]);
		Graph *graph = new Graph(filename, k);
		graph->read();
		graph->search();
		graph->write();
		// delete graph; // there are some bugs in releasing memory
	}
	else printf("[usage]: exe file k\n");
	printf("-----------------------------------------------------------------------------------------\n\n");

	printf("Branching Number: %lld; Color times: %lld; Pub times: %lld\n",branching_num,color_time,pub_time);
	// fprintf(fp,"1\t%d\t%lld\t%lld\t%lld\t%s\n",Best_Size,branching_num,color_time,pub_time,Utility::integer_to_string(TT.elapsed()).c_str());
	// fclose(fp);
}