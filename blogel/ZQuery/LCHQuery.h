#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <queue>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <ext/hash_map>
#include <string.h>
#include <map>
#include <numeric>
#include <utility>
#include <time.h>
#include <deque>
#include <set>
#include <algorithm>
#include <bitset>

using namespace std;


struct STLCRVertexValue {
	vector<pair<unsigned, unsigned> > indL, indE; // 两种index
	int src, dst;

    // void set_Edge(unsigned val, int flg){

	// 	if (flg == 1)
	// 		indE.push_back(val);
	// 	else
	// 		indL.push_back(val);
    // }

	void Initial(){
		src = 1000, dst = -1000;
	}

	void empty(){
		vector<pair<unsigned, unsigned> >().swap(indL);
		vector<pair<unsigned, unsigned> >().swap(indE);
	}

    friend ibinstream & operator<<(ibinstream & m, const STLCRVertexValue & v){
		m<<v.indL;
		m<<v.indE;
		m<<v.src;
		m<<v.dst;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, STLCRVertexValue & v){
		m>>v.indL;
    	m>>v.indE;
		m>>v.src;
		m>>v.dst;

		return m;
    }
};


struct MsgInf {
	int dis;

	MsgInf(){dis = -1;}
    
	MsgInf(int d){dis = d;}

    friend ibinstream& operator<<(ibinstream& m, const MsgInf& idm)
    {
        m << idm.dis;

        return m;
    }

    friend obinstream& operator>>(obinstream& m, MsgInf& idm)
    {
        m >> idm.dis;

        return m;
    }
};


class LCRVertex : public BVertex<VertexID, STLCRVertexValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	void Inner_Compute(MessageContainer& messages){
		if (step_num() == 1){ // 添加元素，剩下的在block里面做
			if (id == src) 
				value().src = 0;
			if (id == dst) 
				value().dst = 0;
		}else{ // 第二个超步，将message中的顶点添加到send中就完成
			for (MessageIter it = messages.begin(); it != messages.end(); it++) {
				if ((*it).dis > 0 and value().src > (*it).dis)
					value().src = (*it).dis;
				
				if ((*it).dis < 0 and value().dst < (*it).dis)
					value().dst = (*it).dis;
			}
			
			if (value().src < 1000 and value().dst > -1000){
				if (khop > value().src-value().dst)
					khop = value().src-value().dst;
			}
		}

		vote_to_halt();
    }
};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {
public:
	//
	virtual void compute(MessageContainer &messages, VertexContainer &vertexes){
		// == 先执行 内部的传递 ==

	}
};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	vector<vector<vector<unsigned> > > Lab1, Lab2;
	vector<int> D1, D2;
	vector<unsigned> *label, *Clabel, *outlabel;
	vector<unsigned> v2p, p2v, v2p_Global, p2v_Global;
	vector<pair<int, int> > QueryPairs;
	unsigned MAXDIS, MAXMOV, MASK;
	unsigned MAXDIS1, MAXMOV1, MASK1;
	vector<int> ActiveV;

	vector<pair<unsigned, int> > Psrc, Pdst, Csrc, Cdst;
	long long InSize = 0, BoundSize = 0;


	int BoundNum = 0, candvalue = 0;
	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}


	virtual void blockInit(VertexContainer& vertexList, BlockContainer& blockList){
		
		vector<int> srcL, dstL;
		int start = 0, tgt = 1000;
		for (int i=start; ; ++i){
			
			if (v2copy[p2v[i]] == 0) continue;

			if (_my_rank % 2 == 1){
				srcL.push_back(p2v[i]);
			}else{
				dstL.push_back(p2v[i]);
			}

			if (srcL.size() >= tgt or dstL.size() >= tgt) break;
		}

		vector<vector<int> > SrcL(_num_workers), DstL(_num_workers);
		for (int i=0; i<_num_workers; ++i){
			SrcL[i] = srcL; DstL[i] = dstL;
		}

		all_to_all_cat(SrcL, DstL);
		worker_barrier();

		srcL.clear(), dstL.clear();

		for (int i=0; i<_num_workers; ++i){
			srcL.insert(srcL.end(), SrcL[i].begin(), SrcL[i].end());
			dstL.insert(dstL.end(), DstL[i].begin(), DstL[i].end());
		}

		vector<vector<int> >().swap(SrcL);
		vector<vector<int> >().swap(DstL);

		sort(srcL.begin(), srcL.end());
		sort(dstL.begin(), dstL.end());

		for (int i=0; i<srcL.size(); ++i){
			if (srcL[i] != dstL[i] and v2part[srcL[i]] != v2part[dstL[i]]){
				// if (v2copy[srcL[i]] == 1 and v2copy[dstL[i]] == 1){
					pair<int, int> elem(srcL[i], dstL[i]);
					QueryPairs.push_back(elem);
				// }
			}
		}

	}


	void active_LCR_vcompute(){ }

	void all_LCR_vcompute(){ }


	void all_bcompute_LCR(){ }


	void InLRead(string indexName){
		
		int p2v_cnt = 0;

		// ==  ActiveV[_my_rank] 代表本分区顶点的数量  ==
		label = new vector<unsigned>[ActiveV[_my_rank]];
		
		ifstream fin((indexName).c_str(), ios::binary);

		fin.read(reinterpret_cast<char*>(&MAXMOV), sizeof(unsigned));
		fin.read(reinterpret_cast<char*>(&MASK), sizeof(unsigned));

		// ===== 读取 p2v =====
		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		p2v.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(p2v.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		v2p.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(v2p.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容
		

		// === v2copy ===
		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		v2copy.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(v2copy.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容


		// ===== 读取label =====
		fin.read(reinterpret_cast<char*>(&BoundNum), sizeof(unsigned));

		for (int i=BoundNum; i<ActiveV[_my_rank]; i++){
			unsigned labcnt = 0;
			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
			label[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(label[i].data()), labcnt * sizeof(unsigned)); // 读内层内容

			InSize += labcnt;
		}

		fin.close();
	}


	void BdLRead(string indexName){
		int p2v_cnt = 0, AllNode = 0;

		ifstream fin((indexName).c_str(), ios::binary);

		fin.read(reinterpret_cast<char*>(&MAXMOV1), sizeof(unsigned));
		fin.read(reinterpret_cast<char*>(&MASK1), sizeof(unsigned));

		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		p2v_Global.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(p2v_Global.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		fin.read(reinterpret_cast<char*>(&p2v_cnt), sizeof(unsigned));
		v2p_Global.resize(p2v_cnt);
		fin.read(reinterpret_cast<char*>(v2p_Global.data()), p2v_cnt * sizeof(unsigned)); // 读内层内容

		fin.read(reinterpret_cast<char*>(&AllNode), sizeof(unsigned));

		Clabel = new vector<unsigned>[AllNode];

		long long ccc = 0;
		for (int i=0; i<AllNode; ++i){			
			unsigned labcnt = 0;
			fin.read(reinterpret_cast<char*>(&labcnt), sizeof(unsigned));
			Clabel[i].resize(labcnt);
			fin.read(reinterpret_cast<char*>(Clabel[i].data()), labcnt * sizeof(unsigned)); // 读内层内容
		}

		fin.close();

		long long maxlabs = 0;
		int ccnt = 0;

		for (int i=0; i<AllNode; ++i){
			int orig_id = p2v_Global[i];

			if (v2part[orig_id] != _my_rank) continue;

			int self_id = v2p[orig_id];

			label[self_id].insert(label[self_id].end(), Clabel[i].begin(), Clabel[i].end());

			BoundSize += Clabel[i].size();

			vector<unsigned>().swap(Clabel[i]);
		}

		delete [] Clabel;
	}


	void initial(int src, int dst){
		
		if (v2part[src] == _my_rank){
			unsigned vsrc = v2p[src];

			if (vsrc < BoundNum){
				if (v2copy[src] == 1) Csrc.push_back(make_pair(src, 0));
				else                  Psrc.push_back(make_pair(src, 0));
			}else{
				for (auto& elem:label[vsrc]){
					unsigned vid = elem>>MAXMOV, dd = elem&MASK;
					if (v2copy[p2v[vid]] == 1) 
						Csrc.push_back(make_pair(p2v[vid], dd));
					else
						Psrc.push_back(make_pair(p2v[vid], dd));
				}
			}
		}

		
		if (v2part[dst] == _my_rank){
			unsigned vdst = v2p[dst];

			if (vdst < BoundNum){
				if (v2copy[dst] == 1) Cdst.push_back(make_pair(dst, 0));
				else                  Pdst.push_back(make_pair(dst, 0));
			}else{
				for (auto& elem:label[vdst]){
					unsigned vid = elem>>MAXMOV, dd = elem&MASK;
					if (v2copy[p2v[vid]] == 1) 
						Cdst.push_back(make_pair(p2v[vid], dd));
					else
						Pdst.push_back(make_pair(p2v[vid], dd));
				}
			}
		}

	}


	int get_min_sum(const vector<pair<unsigned int, int> >& A, const vector<pair<unsigned int, int> >& B) {
		size_t i = 0, j = 0;
		int min_sum = 1000;
		while (i < A.size() && j < B.size()) {
			if (A[i].first < B[j].first) {
				++i;
			} else if (A[i].first > B[j].first) {
				++j;
			} else {
				int sum = A[i].second + B[j].second;
				if (sum < min_sum) min_sum = sum;
				++i;
				++j;
			}
		}
		return min_sum;
	}


	void Refresh(vector<pair<unsigned, int> >& A) {
		if (A.empty()) return;
		sort(A.begin(), A.end()); // 按 first 升序，second 升序
		size_t write = 0;
		for (size_t i = 0; i < A.size(); ) {
			A[write++] = A[i];
			unsigned cur_first = A[i].first;
			// 跳过所有 first 相同的
			size_t j = i + 1;
			while (j < A.size() && A[j].first == cur_first) ++j;
			i = j;
		}
		A.resize(write);
	}


	void run_LCR(const WorkerParams& params){
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		src = params.src, dst = params.dst;
		khop = 10000;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string ss;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		ActiveV.resize(_num_workers);

		while(getline(infile, ss)){
			char* part = new char[strlen(ss.c_str())+1];
			strcpy(part, ss.c_str());
			v2part.push_back( atoi(part) );
			ActiveV[atoi(part)] += 1;
			v2degree.push_back(0);
			delete part;
			nnn += 1;
		}

		totalV = nnn;
		Lab1.resize(_num_workers), Lab2.resize(_num_workers);
        
		string indexName = params.input_path+"In_"+ to_string(_my_rank);
		InLRead(indexName);

		indexName = params.input_path+"Bd";
		BdLRead(indexName);

        worker_barrier();

		long long TotalIn = all_sum_LL(InSize), Innum = all_sum_LL(ActiveV[_my_rank]);		
		long long TotalBound = all_sum_LL(BoundSize), Bnum = all_sum_LL(BoundNum);

		if (_my_rank == 0){
			cout<<"Average InSize: "<<TotalIn/Innum<<endl;
			cout<<"Average BoundSize: "<<TotalBound/Bnum<<endl;
		}


		blockInit(vertexes, blocks); // === 随机生成查询任务队列 ===
// =================================================================
		float TotalTime = 0;
		long long TotalMessage = 0;
		int cnt1 = 0;
	

		long long global_msg_num = 0;
		long long msg_num = 0;

		for (int ii=0; ii<QueryPairs.size(); ++ii){

			vector<vector<pair<unsigned, int> > > MsgSrc(_num_workers), MsgDst(_num_workers);
			
			src = QueryPairs[ii].first, dst = QueryPairs[ii].second;
			
			float s_1 = clock();

			// == 初始化 ==
			initial(src, dst);

			int dup = 100000, ddd = -1;

			for (int k=0; k<2; ++k){
								
				if (k == 0){ // Psrc and Pdst

					for (pair<unsigned, int>& elem: Psrc){
						unsigned vid = elem.first; // 原本的id，vid至少是边界点
						int dis = elem.second; // 确认是源于src还是dst

						for (unsigned& lab: label[v2p[vid]]){
							unsigned tgt = lab>>MAXMOV1, dd = lab&MASK1;
							// == tgt基于p2v_global回归到原始id ==
							unsigned orig_id = p2v_Global[tgt], vpart = v2part[orig_id];
							int ddis = dis+dd; // 确认距离的正负

							MsgSrc[vpart].push_back(make_pair(orig_id, ddis));
							
							if (vpart != _my_rank) msg_num += 1;
						}	
					}

					for (pair<unsigned, int>& elem: Pdst){
						unsigned vid = elem.first; // 原本的id，vid至少是边界点
						int dis = elem.second; // 确认是源于src还是dst

						for (unsigned& lab: label[v2p[vid]]){
							unsigned tgt = lab>>MAXMOV1, dd = lab&MASK1;
							// == tgt基于p2v_global回归到原始id ==
							unsigned orig_id = p2v_Global[tgt], vpart = v2part[orig_id];
							int ddis = dis+dd; // 确认距离的正负

							MsgDst[vpart].push_back(make_pair(orig_id, ddis));

							if (vpart != _my_rank) msg_num += 1;
						}	
					}
				
				}else{  // Csrc and Cdst
					
					for (pair<unsigned, int>& elem: Csrc){
						unsigned vid = elem.first; // 原本的id，vid至少是边界点
						int dis = elem.second; // 确认是源于src还是dst

						if (dis+2 > dup) continue; // 下限约束

						for (unsigned& lab: label[v2p[vid]]){
							unsigned tgt = lab>>MAXMOV1, dd = lab&MASK1;
							// == tgt基于p2v_global回归到原始id ==
							unsigned orig_id = p2v_Global[tgt], vpart = v2part[orig_id];
							int ddis = dis+dd; // 确认距离的正负

							if (ddis + 1 > dup) break; // 后续的label也不会满足距离约束

							Psrc.push_back(make_pair(orig_id, ddis));
						}	
					}

					for (pair<unsigned, int>& elem: Cdst){
						unsigned vid = elem.first; // 原本的id，vid至少是边界点
						int dis = elem.second; // 确认是源于src还是dst

						if (dis+2 > dup) continue; // 下限约束

						for (unsigned& lab: label[v2p[vid]]){
							unsigned tgt = lab>>MAXMOV1, dd = lab&MASK1;
							// == tgt基于p2v_global回归到原始id ==
							unsigned orig_id = p2v_Global[tgt], vpart = v2part[orig_id];
							int ddis = dis+dd; // 确认距离的正负

							if (ddis + 1 > dup) break; // 后续的label也不会满足距离约束

							Pdst.push_back(make_pair(orig_id, ddis));
						}	
					}

					Refresh(Psrc), Refresh(Pdst);

					for (pair<unsigned int, int>& elem: Psrc){
						int vpart = v2part[elem.first];
						MsgSrc[vpart].push_back(elem);
						if (vpart != _my_rank) msg_num += 1;
					}

					for (pair<unsigned int, int>& elem: Pdst){
						int vpart = v2part[elem.first];
						MsgDst[vpart].push_back(elem);
						if (vpart != _my_rank) msg_num += 1;
					}
				}
				
				all_to_all_cat(MsgSrc, MsgDst);
				worker_barrier();

				for (int i=0; i<_num_workers; ++i){
					Csrc.insert(Csrc.end(), MsgSrc[i].begin(), MsgSrc[i].end());
					MsgSrc[i].clear();
					Cdst.insert(Cdst.end(), MsgDst[i].begin(), MsgDst[i].end());
					MsgDst[i].clear();
				}

				Psrc.clear(); Pdst.clear();
				Refresh(Csrc);
				Refresh(Cdst);
				
				if (k == 0){
					ddd = get_min_sum(Csrc, Cdst);
					dup = all_min(ddd);
				}
			}
			
			worker_barrier();
			float s_2 = clock();


			TotalTime += (s_2-s_1);

			global_msg_num = all_sum_LL(msg_num);

			if (_my_rank == 0 and ii % 500 == 0){
				cout<<" LCH d:"<<ii<<"  Average: "<<float(TotalTime)/(ii*CLOCKS_PER_SEC)<<" s , Comm: "<<float(global_msg_num)/(ii*1024)<<endl;
			}


			Psrc.clear(), Pdst.clear();
			Csrc.clear(), Cdst.clear();
			vector<vector<pair<unsigned, int> > >().swap(MsgSrc);
			vector<vector<pair<unsigned, int> > >().swap(MsgDst);
		}
		
	}

};


void SetLCR_Construction(string in_path, string partition_path, string out_path, int src, int dst, int khop){
	WorkerParams param;
	param.input_path=in_path;
	param.partition_path = partition_path;
	param.output_path=out_path;
	param.src = src;
	param.dst = dst;
	param.khop = khop;
	param.force_write=true;
	param.native_dispatcher=false;
	LCRBlockWorker worker;
	worker.set_compute_mode(LCRBlockWorker::VB_COMP);

	worker.run_LCR(param);
};

