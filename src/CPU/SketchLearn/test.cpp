#include "SmallSketch.h"
#include "MultiLevelSketch.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <unordered_map>

#include <iostream>

using namespace std;

#define START_FILE_NO 1
#define END_FILE_NO 1

struct FIVE_TUPLE{	char key[13];	};
typedef vector<FIVE_TUPLE> TRACE;
TRACE traces[END_FILE_NO - START_FILE_NO + 1];

void ReadInTraces(const char *trace_prefix)
{
	for(int datafileCnt = START_FILE_NO; datafileCnt <= END_FILE_NO; ++datafileCnt)
	{
		char datafileName[100];
		sprintf(datafileName, "%s%d.dat", trace_prefix, datafileCnt - 1);
		FILE *fin = fopen(datafileName, "rb");

		FIVE_TUPLE tmp_five_tuple;
		traces[datafileCnt - 1].clear();
		while(fread(&tmp_five_tuple, 1, 13, fin) == 13)
		{
			traces[datafileCnt - 1].push_back(tmp_five_tuple);
		}
		fclose(fin);

		printf("Successfully read in %s, %ld packets\n", datafileName, traces[datafileCnt - 1].size());
	}
	printf("\n");
}

int main()
{
	ReadInTraces("../../../data/");

	#define KEY_LEN 4
	#define SK_D 3
	#define THETA 0.5
	MultiLevelSketch<KEY_LEN, SK_D> *sketchlearn = NULL;

	for(int datafileCnt = START_FILE_NO; datafileCnt <= END_FILE_NO; ++datafileCnt){
		unordered_map<string, int> Real_Freq;
		sketchlearn = new MultiLevelSketch<KEY_LEN, SK_D>(10 * 1024); // 600 * 1024 / 12

		int packet_cnt = (int)traces[datafileCnt - 1].size();
		for(int i = 0; i < packet_cnt; ++i)
		{
			sketchlearn->insert((uint8_t*)(traces[datafileCnt - 1][i].key));
			string str((const char*)(traces[datafileCnt - 1][i].key), 4);
			Real_Freq[str]++;
		}

		// sketchlearn -> model_inference(THETA);
	}

}	
