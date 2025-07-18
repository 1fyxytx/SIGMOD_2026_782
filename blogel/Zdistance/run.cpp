#include "Dish.h"
// #include "Compare.h"

using namespace std;

int main(int argc, char* argv[]){
	init_workers();
	SetLCR_Construction(  "/blogel/test.graph",
				          "/blogel/p2",
			              "/blogel/index_",
			        1049, 123980, 32); //  Distri4hop 索引构建 的输入

	worker_finalize();
	

	return 0;
}
