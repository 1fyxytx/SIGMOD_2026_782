#include "LCHQuery.h"

using namespace std;


int main(int argc, char* argv[]){
	init_workers();
	SetLCR_Construction(  "/mnt/data/zyy/LCHIndex/uk2002/index_",
						"/mnt/data/zyy/part/uk2002/p10",
						"/home/hnu/Disk0/zyy_dataset/",
				1049, 1200, 20); //  距离查询的输入
	

	worker_finalize();

	return 0;
}

