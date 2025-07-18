#include "Dish.h"
// #include "Compare.h"

using namespace std;

int main(int argc, char* argv[]){
	init_workers();
	SetLCR_Construction(  "/mnt/data/zyy/Full/indochina.graph",
				          "/mnt/data/zyy/part/indochina/p10",
			              "/mnt/data/zyy/LCHIndex/it2004/index_",
			        1049, 123980, 32); //  Distri4hop 索引构建 的输入

	worker_finalize();
	

	return 0;
}
