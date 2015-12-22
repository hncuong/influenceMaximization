#include "simpath.hpp"
#include <unistd.h>
using namespace std;


time_t startTime;
float getTime(){
	time_t curTime;
	time(&curTime);

	float min = ((float)(curTime - startTime));
	return min;
}

float getCurrentMemoryUsage() {

	string pid = intToStr(unsigned(getpid()));
	string outfile = "temp/tmp_" + pid + ".txt";
	string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
	system(command.c_str());

	string mem_str;
	ifstream ifs(outfile.c_str());
	std::getline(ifs, mem_str);
	ifs.close();

	mem_str = mem_str.substr(0, mem_str.size()-1);
	float mem = (float)strToInt(mem_str);

	return mem/1024;
	
	return 0;
}

int main(int argc, char * argv[])
{

	SimPath *s = new SimPath();
	double _pruneVal = 0.001; 
	int _topL = 4;
	int _k = 5;
	if (argc == 2){
		_k = strToInt(argv[1]);		
	}
	string file = "netHEP.txt";
	s->setParameter(_pruneVal,_topL,_k);
	s->readDataSets(file);
	
	time(&startTime);
	s->findVertexCover();
	s->simPathCore();
    	cout << "Total time taken : " << getTime() << " seconds" << endl;
	cout << "All memory released, current usage : " << getCurrentMemoryUsage() << " M" << endl;
	return 0;
}
