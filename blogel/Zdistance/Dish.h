#include "../utils/communication.h"
#include "../blogel/BVertex.h"
#include "../blogel/Block.h"
#include "../blogel/BWorker.h"
#include "../blogel/BGlobal.h"
#include "../blogel/BType.h"
#include <queue>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <unordered_map>
#include <string.h>
#include <map>
#include <numeric>
#include <utility>
#include <time.h>
#include <deque>
#include <set>
#include <algorithm>
#include <bitset>
#include <vector>
#include <omp.h>

// 核心点：并行化树分解，如何确认顶点之间的互不影响性
// 结合分区和树分解的性质，进一步加速分解过程
#define MAXINT ((unsigned) 4294967295)

using namespace std;

struct EnumValue {
	vector<int> Local, Edge; 

    void set_Edge(int adj_id, int id1){

		if (v2part[id1] != v2part[adj_id])
			Edge.push_back(adj_id);
		else
			Local.push_back(adj_id);
    }

	void empty(){
		vector<int>().swap(Local);
		vector<int>().swap(Edge);
	}

    friend ibinstream & operator<<(ibinstream & m, const EnumValue & v){
		m<<v.Local;
		m<<v.Edge;

    	return m;
    }


    friend obinstream & operator>>(obinstream & m, EnumValue & v){
		m>>v.Local;
    	m>>v.Edge;
    	
		return m;
    }
};


struct MsgInf {
	vector<treeinf> Neig;

	MsgInf(){}

    friend ibinstream& operator<<(ibinstream& m, const MsgInf& idm)
    {
        return m;
    }

    friend obinstream& operator>>(obinstream& m, MsgInf& idm)
    {
        return m;
    }
};


class LCRVertex : public BVertex<VertexID, EnumValue, MsgInf>{
public:
	virtual void compute(MessageContainer& messages){}

	// === 构建 dist ==

	


	// =======边界索引构建策略=========、
	void Bound_Compute(MessageContainer& messages){ // 按道理讲，这个过程不需要接收其他分区顶点的信息

		vote_to_halt();
	}

};


class LCRBlock : public Block<char, LCRVertex, MsgInf> {
public:

	void Boundary_compute(VertexContainer &vertexes){ // 边界图计算,不接收信息

		vote_to_halt();
	}



	virtual void compute(MessageContainer &messages, VertexContainer &vertexes){}

};


class LCRBlockWorker : public BWorker<LCRBlock>{
public:
	unsigned MAXDIS, MAXMOV, MASK;
	unsigned MAXDIS1, MAXMOV1, MASK1;
	vector<unsigned> *label, *root, *flags;
	vector<unsigned> *Clabel, *Cflags;


	vector<TriElem> *treelab;
	vector<vector<unsigned> > con, conB, conI;
	vector<unsigned> InnerCopy, InnerBound;
	vector<int> idList;
	bool *usd_bp;
	bool* is_indep;
	BPLabel *label_bp;
    int BoundNum = 0, topk = 0, LocalActive = 0;
	long long BoundEdge = 0;
	long long cntt = 0;
	// ==================================
	int n_org, n, n_core;
	long long m, m_core;
	
	vector<vector<int>> nbr, cost;
	vector<int> ord, rank;
	vector<int> ActiveV;
	vector<TreeNode> tree;
	double t;
	

	virtual LCRVertex *toVertex(char *line){
		LCRVertex *v;
		return v;
	}


	void ParameterInitial(){
		unsigned actV = *max_element(ActiveV.begin(), ActiveV.end());
		MAXDIS = 2; MAXMOV = 1;
		while( MAXINT / (actV * 2) >= MAXDIS ) {
			MAXDIS *= 2;
			++MAXMOV;
		}
		MASK = MAXDIS - 1;

		v2p.resize(totalV, -1);
		p2v.resize(totalV, -1);
		v2p_Global.resize(totalV, -1);
		p2v_Global.resize(totalV, -1);	

		InnerCopy.resize(totalV, -1);
		InnerBound.resize(totalV, -1);
		
	}

	
    virtual void blockInit(VertexContainer &vertexes, BlockContainer &blocks){

		ParameterInitial();

		SelfCopy = 0;

		for (int i=0; i<v2degree.size(); ++i){
			int flg = v2degree[i] > 0? 1 : 0;
			int valB = v2copy[i]*000000000+ flg * 10000000 + v2degree_R[i]; // 优先copy，其次bound，
			int valC = v2copy[i]*900000000 + flg * 10000000 + v2degree_R[i];

			if (v2part[i] == _my_rank){
				sortList.push_back(make_pair(valB, -i));
				if (v2copy[i] == 1) SelfCopy += 1;
			}

			if (v2degree[i] > 0)
				sortTotal.push_back(make_pair(valC, -i)); // pos就是id					
		}
		
		sort(sortList.rbegin(), sortList.rend());
		
		for (int i=0; i<sortList.size(); ++i){
			v2p[-sortList[i].second] = i;
			p2v[i] = -sortList[i].second; // new 2 old

			if (v2copy[-sortList[i].second] == 1) InnerCopy[i] = 1;

			if (v2degree[-sortList[i].second] > 0) InnerBound[i] = 1;
		}
		
		sort(sortTotal.rbegin(), sortTotal.rend());
		
		for (int i=0; i<sortTotal.size(); ++i){
			v2p_Global[-sortTotal[i].second] = i;
			p2v_Global[i] = -sortTotal[i].second; // new 2 old
		}

		// ======================================
        con.resize(sortList.size());
		conB.resize(sortTotal.size());

		MAXDIS1 = 2; MAXMOV1 = 1;
		while( MAXINT / (BdV * 2) >= MAXDIS1 ) {
			MAXDIS1 *= 2;
			++MAXMOV1;
		}
		MASK1 = MAXDIS1 - 1;


		for (int i=0; i<sortList.size(); ++i){
            int ovid = p2v[i];
            vector<int>& local = vertexes[vert2place[ovid]]->value().Local;
            
			for( int p = 0; p < local.size(); ++p ) {
                int j = v2p[local[p]];

                con[j].push_back(i); // 内部存储的
		    }
			
			// ===================================
			vector<int>& edge = vertexes[vert2place[ovid]]->value().Edge;
			if (edge.size() > 0){
				int newPos = v2p_Global[ovid]; // 在边界图上的定位
				for( int p = 0; p < edge.size(); ++p ) {
					int newVid = v2p_Global[edge[p]]; // 在边界图上的定位

					unsigned elem = newVid << MAXMOV1 | 1;

					conB[newPos].push_back(elem); // 存的是原始id信息
				}
			}

			vertexes[vert2place[ovid]]->value().empty();
        }

		vector<pair<int, int>>().swap(sortList);
		vector<LCRVertex *>().swap(vertexes);
		vector<pair<int, int>>().swap(sortTotal);	

		long long edges = 0;
		for (int i=0; i<con.size(); ++i){
			edges += con[i].size();
		}

	}


	void AddVertex(char *line, int va){
		// 三个元素  V_A, V_B, Label
		int vb, pa, pb;
		// pa = v2part[va];
		LCRVertex* v = new LCRVertex;
		v->id = va; v->bid = 0;
		load_vertex(v);
		vert2place[va] = vertexes.size()-1;
		char* s1 = strtok(line," ");

		while(s1){
			totalEdges += 1;
			vb = atoi(s1) - 1;
			v->value().set_Edge(vb, va);			
			s1=strtok(NULL," ");
		}
		
		v2degree[va] = v->value().Edge.size();
		v2degree_R[va] = v->value().Edge.size() + v->value().Local.size();

		if (v->value().Edge.size() > 0) BoundNum += 1;
	}


	void VertexCut_CutEdge(int batch){ // == 每个分区各自选择 topk 的copy点 ==
		vector<pair<int, int> > RankV; // 贡献值 id
		vector<vector<int> > VList;

		v2copy.resize(totalV, 0);
		for (int i=0; i<vertexes.size(); ++i){
			if (vertexes[i]->value().Edge.size() > 0){
				idList.push_back(vertexes[i]->id);
			}
		}
				
		while(1){
			for (int i=0; i<idList.size(); ++i){
				int vid = idList[i];

				if (v2copy[vid] == 1) continue;

				vector<int>& edges = vertexes[vert2place[vid]]->value().Edge;
				int val = 0;
				
				if (v2degree[vid] == 1){
					if (v2degree[edges[0]] == 1)
						if (v2degree_R[vid] > v2degree_R[edges[0]] or 
						    (v2degree_R[vid] == v2degree_R[edges[0]] and vid < edges[0])) 
							val += 1;
				}else{
					for (int j=0; j<edges.size(); ++j)
						if (v2copy[edges[j]] == 0)
							val += 1;
				}

				if (val > 0) // val = 0 表示没有切边了，自身也不是边界点了
					RankV.push_back(make_pair(val, -vid));
			}

			int cnn11 = all_sum(RankV.size());
			if (cnn11 == 0)
				break;

			sort(RankV.rbegin(), RankV.rend());

			// ==== 批量选择顶点来成为复写点 ====
			if (batch > RankV.size())  batch = RankV.size();

			VList.resize(_num_workers);
			
			for (int i=0; i<batch; ++i)
				for (int j=0; j<_num_workers; ++j)
					VList[j].push_back(-RankV[i].second);
				
			all_to_all(VList);
			worker_barrier();

			for (int j=0; j<_num_workers; ++j)
				for (int kk=0; kk<VList[j].size(); ++kk)
					v2copy[VList[j][kk]] = 1;
	
			vector<vector<int > >().swap(VList);
			vector<pair<int, int> >().swap(RankV);
		}
		
		CopyNum = accumulate(v2copy.begin(), v2copy.end(), 0);

		if (_my_rank == 0){
			cout<<"B: "<<BdV<<"  C: "<<CopyNum<<endl;
		}
	}


	bool can_update(int v, int dis, char *nowdis) {
		for( int i = 0; i < (int) label[v].size(); ++i ) {
			int w = label[v][i]>>MAXMOV, d = label[v][i]&MASK;
			if( nowdis[w] >= 0 && nowdis[w] + d <= dis ) return false;
		}
		return true;
	}

	bool can_update_C(int v, int dis, char *nowdis) {
		for( int i = 0; i < (int) Clabel[v].size(); ++i ) {
			int w = Clabel[v][i]>>MAXMOV1, d = Clabel[v][i]&MASK1;
			if( nowdis[w] >= 0 && nowdis[w] + d <= dis ) return false;
		}
		return true;
	}
	
	int prune_by_root(int u, int v, int d) {
		// 依托 root 的labels 进行剪枝
		int mind = 1000000;
		for (int i=0; i<topk; ++i){
			int dd = root[u][i]+root[v][i];
			if (mind > dd)
				mind = dd;
		}

		if (d >= mind) 
			return 0;
		else
			return 1;
    }
	
	void BoundGraph(){
		
		vector<vector<unsigned> >().swap(con);

		for (int i=0; i<BoundNum; ++i){
			// obnd_id 是在边界图里面的 位置
			int oid = p2v[i], 
			    ob_id = v2p_Global[oid];
			
			vector<unsigned>& lbs = label[i];
			vector<unsigned>& fls = flags[i];

			for (int j=0; j<lbs.size(); ++j){
				
				if (fls[j] != 0) continue; // 不影响距离

				unsigned vid = lbs[j] >> MAXMOV,
						 dis = lbs[j] & MASK;

				int ovid = p2v[vid], ob_vid = v2p_Global[ovid];
				
				if (ob_id == ob_vid or dis == 0) continue;

				// ============= 
				unsigned elem1 = ob_vid << MAXMOV1 | dis,
					     elem2 =  ob_id << MAXMOV1 | dis;
				
				conB[ob_id].push_back(elem1);
				conB[ob_vid].push_back(elem2);
			}

			vector<unsigned>().swap(lbs);
			vector<unsigned>().swap(fls);

		}

		
		// delete[] label;
		// delete[] flags;


		vector<vector<vector<unsigned> > > CB(_num_workers);
		
		long long Edgess, aa = 0;
		for (int i=0; i<conB.size(); ++i){
			aa += conB[i].size();
		}
		
		Edgess = all_sum_LL(aa);
		worker_barrier();

		if (_my_rank == 0)
		cout<<"V: "<<CopyNum<<"   "<<BdV<<"   E: "<<Edgess/2<<endl;

		// == 收集 CopyNum 顶点的标签信息  ==
		int start = 0, batch = 100000;
		if (start+batch > BdV)
			batch = BdV-start;

		while (1){
			for (int i=start; i<start+batch; ++i){
				
				if (conB[i].size() == 0) continue;

				conB[i].push_back(i); 
				for (int j=0; j<_num_workers; ++j)
					if (j != _my_rank)
						CB[j].push_back(conB[i]);
				conB[i].pop_back();
			}
			
			all_to_all(CB);
			worker_barrier();

			for(int i=0; i<_num_workers; ++i){
				
				for (int j=0; j<CB[i].size(); ++j){
					
					vector<unsigned>& inf = CB[i][j];
					int vid = inf[inf.size()-1];
					inf.pop_back();

					conB[vid].insert(conB[vid].end(), inf.begin(), inf.end());
					
					vector<unsigned>().swap(inf);
				}
			}
			
			for (int i=0; i<_num_workers; ++i)
				vector<vector<unsigned> >().swap(CB[i]);

			vector<vector<vector<unsigned> > >(CB).swap(CB);

			start += batch;

			int cntt = all_sum(start);

			// if (_my_rank == 0)
			// 	cout<<"Cnt: "<<cntt<<"  "<< _num_workers*BdV<<endl;

			if (cntt == _num_workers*BdV) break;

			if (start+batch > BdV)
				batch = BdV-start;
		}


		for (int i=0; i<conB.size(); ++i){
			
			if (conB[i].size() == 0) continue;

			vector<unsigned>& edges = conB[i];
			sort(edges.begin(), edges.end());
			int p = 1;
			for( int j = 1; j < (int) edges.size(); ++j ){
				unsigned id1 = edges[j-1] >> MAXMOV1,
						 id2 = edges[j] >> MAXMOV1;
				if( id1 != id2 ) 
					edges[p++] = edges[j];
			}
			edges.resize(p);
		}
	}


	unsigned Query(int u, int v){

		if( u == v ) return 0;
		unsigned lu = (unsigned)label[u].size(), lv = (unsigned)label[v].size(), dis = MAXD;
		for( int i = 0, j = 0; i < lu && j < lv; ++i ) {	
			for( ; j < lv && label[v][j]>>MAXMOV < label[u][i]>>MAXMOV; ++j ) ;
			if( j < lv && label[v][j]>>MAXMOV == label[u][i]>>MAXMOV) 
				dis = min(dis, (label[u][i]&MASK) + (label[v][j]&MASK));
		}
		return dis;
	}


	void Part2hop(){

		totalV = con.size();

		omp_set_num_threads(threads);
		double t = omp_get_wtime();

		int **pos = new int*[totalV];
		for( int i = 0; i < totalV; ++i ) 
			pos[i] = new int[MAXDIS];

		label = new vector<unsigned>[totalV]; // unsigned 整合了 id+dis
		flags = new vector<unsigned>[totalV]; 

		for(int i = 0; i < totalV; ++i){
			if (i < BoundNum){ // LocalActive
				label[i].push_back((((unsigned)i)<<MAXMOV) | 0);
				flags[i].push_back(0); // 0表示这个index是必要的
				pos[i][0] = 1;
			}else{
				pos[i][0] = 0;
			}

			for (int j=0; j<con[i].size(); ++j){
				if (con[i][j] < i && con[i][j] < BoundNum){
					unsigned elem = con[i][j]<<MAXMOV | 1;
					label[i].push_back(elem);
					flags[i].push_back(0);
				}
			}
			
			pos[i][1] = (int) label[i].size();
		}

		int dis = 2;
		for( long long cnt = 1; cnt && dis <= MAXDIS; ++dis ){
			
			cnt = 0;
			vector<unsigned> *label_new = new vector<unsigned>[totalV];
			vector<unsigned> *flg_new = new vector<unsigned>[totalV];
			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0;
				unsigned char *used = new unsigned char[totalV/8+1];
				memset( used, 0, sizeof(unsigned char) * (totalV/8+1) );
				vector<int> cand, candFg, flgs(totalV, 0);
				
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
	
				for( int u = pid; u < totalV; u += np ){

					cand.clear();
					candFg.clear();
					
					for( int i = 0; i < con[u].size(); ++i ){
						
						int w = con[u][i];

						for( int j = pos[w][dis-2]; j < pos[w][dis-1]; ++j ) {
							
							int v = label[w][j] >> MAXMOV;
							
							if( v >= u ) break;

							int ff = flags[w][j];
							
							if (ff < 2){
								if (InnerCopy[w] == 1) // 似乎没有了意义
									ff = 2;
								else if(InnerBound[w] == 1 and InnerCopy[w] != 1){
										ff = 1;
								}
							}
						
							if (ff != 0){
								if (flgs[v] == 0) candFg.push_back(v);
								if (flgs[v] < ff) flgs[v] = ff;
							}
	
							if( !(used[v/8]&(1<<(v%8))) ) {
								used[v/8] |= (1<<(v%8)), cand.push_back(v);
							}
					
						}
					}

					int n_cand = 0;
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i] >> MAXMOV] = label[u][i] & MASK;
					
					for( int i = 0; i < (int) cand.size(); ++i ) {
						used[cand[i]/8] = 0;

						if( can_update(cand[i], dis, nowdis) ) 
							cand[n_cand++] = cand[i]; 
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());
					for( int i = 0; i < (int) cand.size(); ++i ) {
						label_new[u].push_back((((unsigned)cand[i])<<MAXMOV) | (unsigned) dis), ++local_cnt;
						flg_new[u].push_back(flgs[cand[i]]);
					}
					for( int i = 0; i < (int) label[u].size(); ++i ) 
						nowdis[label[u][i]>>MAXMOV] = -1;

					for (int vid: candFg) flgs[vid] = 0;
				}
				#pragma omp critical
				{
					cnt += local_cnt;
				}
				delete[] used; delete[] nowdis;
			}

			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){
					label[u].insert(label[u].end(), label_new[u].begin(), label_new[u].end());
					flags[u].insert(flags[u].end(), flg_new[u].begin(), flg_new[u].end());

					vector<unsigned>(label[u]).swap(label[u]);
					vector<unsigned>(flags[u]).swap(flags[u]);

					vector<unsigned>().swap(label_new[u]);
					vector<unsigned>().swap(flg_new[u]);

					pos[u][dis] = (int) label[u].size();
				}
			}

			delete[] label_new,flg_new;
		}

		delete[] pos;

	}

	

	void Core2hop(){
		totalV = conB.size();

		omp_set_num_threads(threads);

		vector<int>* pos = new vector<int>[totalV];

		Clabel = new vector<unsigned>[totalV]; // unsigned 整合了 id+dis
		Cflags = new vector<unsigned>[totalV]; 

		for(int i = 0; i < totalV; ++i){
			if (i < CopyNum){
				Clabel[i].push_back((((unsigned)i)<<MAXMOV1) | 0);
				Cflags[i].push_back(0); // 0表示这个index是必要的
				pos[i].push_back(1);
			}else{
				pos[i].push_back(0);
			}
		}

		int dis = 1;
		for( long long cnt = 1; cnt && dis <= MAXDIS; ++dis ){
			cnt = 0;
			vector<unsigned> *label_new = new vector<unsigned>[totalV];
			vector<unsigned> *flg_new = new vector<unsigned>[totalV];

			#pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				long long local_cnt = 0;
				unsigned char *used = new unsigned char[totalV/8+1];
				memset( used, 0, sizeof(unsigned char) * (totalV/8+1) );
				
				char *nowdis = new char[totalV];
				memset( nowdis, -1, sizeof(char) * totalV);
				vector<int> cand, candFg, flgs(totalV, 0);

				for( int u = pid; u < totalV; u += np ){

					cand.clear();
					candFg.clear();

					for( int i = 0; i < conB[u].size(); ++i ){
						int w = conB[u][i] >> MAXMOV1, 
						    d = conB[u][i] &  MASK1;

						if (d == 0 or d > dis) continue;

						for( int j = d==dis? 0:pos[w][dis-d-1]; j < pos[w][dis-d]; ++j ) {

							int v = Clabel[w][j] >> MAXMOV1;
							
							if( v >= u ) break;

							int ff = Cflags[w][j];

							if (dis >= 2 and w < CopyNum) ff = 3;

							if (ff != 0){
								if (flgs[v] == 0) candFg.push_back(v);
								if (flgs[v] < ff) flgs[v] = ff;
							}

							if( !(used[v/8]&(1<<(v%8))) ) 
								used[v/8] |= (1<<(v%8)), cand.push_back(v);
						}
					}

					int n_cand = 0;
					for( int i = 0; i < (int) Clabel[u].size(); ++i ) 
						nowdis[Clabel[u][i] >> MAXMOV1] = Clabel[u][i] & MASK1;
					
					for( int i = 0; i < (int) cand.size(); ++i ) {
						used[cand[i]/8] = 0;
						// if( !prune_by_bp(u, cand[i], dis) )
						if( can_update_C(cand[i], dis, nowdis) ) 
							cand[n_cand++] = cand[i]; 
					}

					cand.resize(n_cand);
					sort(cand.begin(), cand.end());
					for( int i = 0; i < (int) cand.size(); ++i ) {
						label_new[u].push_back((((unsigned)cand[i])<<MAXMOV1) | (unsigned) dis), ++local_cnt;
						flg_new[u].push_back(flgs[cand[i]]);
					}
					for( int i = 0; i < (int) Clabel[u].size(); ++i ) 
						nowdis[Clabel[u][i]>>MAXMOV1] = -1;

					for (int vid: candFg) flgs[vid] = 0;
				}
				#pragma omp critical
				{
					cnt += local_cnt;
				}
				delete[] used; delete[] nowdis;
			}

			#	pragma omp parallel
			{
				int pid = omp_get_thread_num(), np = omp_get_num_threads();
				for( int u = pid; u < totalV; u += np ){
					Clabel[u].insert(Clabel[u].end(), label_new[u].begin(), label_new[u].end());
					Cflags[u].insert(Cflags[u].end(), flg_new[u].begin(), flg_new[u].end());

					vector<unsigned>().swap(label_new[u]);
					vector<unsigned>().swap(flg_new[u]);

					pos[u].push_back((int) Clabel[u].size());
				}
			}

			// if (_my_rank == 0)
			cout<<"Distance: "<<dis<<"   Cnt: "<<cnt<<endl;

			delete[] label_new;
		}

		delete[] pos;
	}


	void DHI_Store(string new_filename){
		// FILE *fout = fopen( (new_filename).c_str(), "wb" );
		
		ofstream fout((new_filename).c_str(), ios::binary);

		// == 先把记录label的两个参数存储  ==
		fout.write(reinterpret_cast<const char*>(&MAXMOV), sizeof(unsigned));
		fout.write(reinterpret_cast<const char*>(&MASK), sizeof(unsigned)); 
		
		// == 存储 p2v ==
		int p2v_cnt = p2v.size();
		fout.write(reinterpret_cast<const char*>(&p2v_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(p2v.data()), sizeof(unsigned)*p2v_cnt);

		int v2p_cnt = v2p.size();
		fout.write(reinterpret_cast<const char*>(&v2p_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(v2p.data()), sizeof(unsigned)*v2p_cnt);

		int copy_cnt = v2copy.size();
		fout.write(reinterpret_cast<const char*>(&copy_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(v2copy.data()), sizeof(unsigned)*copy_cnt);


		// == 先写label[i]的size ==
		fout.write(reinterpret_cast<const char*>(&BoundNum), sizeof(unsigned)); // 内部起始位置

		for (int i=BoundNum; i<con.size(); ++i){
			int len = (int) label[i].size(), c = 0;

			for (int j=0; j<len; ++j){
				if (flags[i][j] == 0){
					label[i][c] = label[i][j];
					c += 1;
				}
			}

			label[i].resize(c);

			fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
			fout.write(reinterpret_cast<const char*>(label[i].data()), c*sizeof(unsigned));
		}

		fout.close();
	}


	void DHBound_Store(string new_filename){
		
		ofstream fout((new_filename).c_str(), ios::binary);

		// == 先把记录label的两个参数存储  ==
		fout.write(reinterpret_cast<const char*>(&MAXMOV), sizeof(unsigned));
		fout.write(reinterpret_cast<const char*>(&MASK), sizeof(unsigned)); 
		
		// == 存储 p2v_Global ==
		int p2v_cnt = p2v_Global.size();
		fout.write(reinterpret_cast<const char*>(&p2v_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(p2v_Global.data()), sizeof(unsigned)*p2v_cnt);

		int v2p_cnt = v2p_Global.size();
		fout.write(reinterpret_cast<const char*>(&v2p_cnt), sizeof(unsigned)); 
		fout.write(reinterpret_cast<const char*>(v2p_Global.data()), sizeof(unsigned)*v2p_cnt);

		// == 先写label[i]的size ==
		int AllNode = conB.size();
		fout.write(reinterpret_cast<const char*>(&AllNode), sizeof(unsigned)); // 内部起始位置

		for (int i=0; i<AllNode; ++i){
			int len = (int) Clabel[i].size(), c;

			if (i < CopyNum){
				fout.write(reinterpret_cast<const char*>(&len), sizeof(unsigned)); // label[i]的size
				fout.write(reinterpret_cast<const char*>(Clabel[i].data()), len*sizeof(unsigned));
			}else{
				c = 0;
				for (int j=0; j<len; ++j){
					if (Cflags[i][j] == 0){
						Clabel[i][c] = Clabel[i][j];
						c += 1;
					}
				}
				Clabel[i].resize(c);

				fout.write(reinterpret_cast<const char*>(&c), sizeof(unsigned)); // label[i]的size
				fout.write(reinterpret_cast<const char*>(Clabel[i].data()), c*sizeof(unsigned));
			}
		}

		fout.close();
		
	}


	void PartitionRewrite(string new_filename, int cc){
		vector< vector< pair<int, int> > > Parts(_num_workers);
		for (int i=0; i<vertexes.size(); ++i){
			int vid = vertexes[i]->id;
			
			if (v2degree[vid] == v2degree_R[vid] and v2degree[vid] <= cc){
				
				// int flg = 0, p1 = v2part[vertexes[i]->value().Edge[0]];
				vector<int> neigs(_num_workers);

				for (int ii=0; ii<vertexes[i]->value().Edge.size(); ++ii){
					int p2 = v2part[vertexes[i]->value().Edge[ii]];
					neigs[p2] += 1;
				}

				auto max_it = max_element(neigs.begin(), neigs.end());
				int p1 = distance(neigs.begin(), max_it);
				
				pair<int, int> aa(vid, p1);
				for (int ii=0; ii<_num_workers; ++ii){
					Parts[ii].push_back(aa);
				}
			}

		}
		all_to_all(Parts);

		for (int i=0; i<_num_workers; ++i){
			for (pair<int, int> aa:Parts[i]){
				v2part[aa.first] = aa.second;
			}
		}

		if (_my_rank == 0){
			const char *file = new_filename.c_str();
			fstream outfileX;
			outfileX.open(file, ios::out);

			for (int i=0; i<v2part.size(); ++i){
				string Inf1 = to_String(v2part[i]); // con[i][0] is the number of in neighbors 
				outfileX<<Inf1<<endl;
			}
		}

	}
	

	void check(int cc){
		int cntt = 0, cc1 = 0;
		for (int i=0; i<vertexes.size(); ++i){
			int vid = vertexes[i]->id;
			if (v2degree[vid] == v2degree_R[vid] and v2degree[vid] <= cc){
				cc1 += 1;
				int flg = 0, p1 = v2part[vertexes[i]->value().Edge[0]];
				
				// if (vertexes[i]->value().Edge.size() != cc) cout<<"??"<<endl;
				for (int ii=1; ii<vertexes[i]->value().Edge.size(); ++ii){
					int p2 = v2part[vertexes[i]->value().Edge[ii]];
					if (p1 != p2){
						flg = 1;
						break;
					}
				}

				if (flg == 0){
					cntt += 1;
				}
			}
		}

		cout<<_my_rank<<"    "<<cntt<<"  "<<cc1<<endl;
	}


	void related_compute(int deg){
		vector<int> Val(con.size(), 0);
		int cnnt = 0;
		for (int i=0; i<CopyNum; ++i){
			if (con[i].size()<deg){
				cnnt += 1;
				for (int neig: con[i]){
					if (v2copy[neig] != 1){
						Val[neig] = 1;
					}
				}
			}
		}

		int vval = accumulate(Val.begin(), Val.end(), 0);
		cout<<_my_rank<<"  "<<cnnt<<"   "<<float(vval)/cnnt<<endl;
	}

	
	void run_LCR(const WorkerParams& params){
		
		ifstream infile, infile1, infile2; // 先把分区文件读取
		
		threads = params.khop;

		const char *part_path = params.partition_path.c_str(); 
		infile.open(part_path);
		string s;
		if(!infile.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int nnn = 0;
		ActiveV.resize(_num_workers);

		while(getline(infile, s)){
			char* part = new char[strlen(s.c_str())+1];
			strcpy(part, s.c_str());
			v2part.push_back( atoi(part) ); // 
			v2degree.push_back(0); 
			delete part;
			nnn += 1;
		}
		v2degree_R.resize(v2degree.size());
		
		const char *filepath = params.input_path.c_str(); //"data//Amazon_New.txt";

		// ======================
		infile1.open(filepath);
		if(!infile1.is_open()){
			cout<<"No such file!"<<endl;
			exit(-1);
		}

		int xxx = 0;
		while(getline(infile1, s)){
			if (xxx > 0 and v2part[xxx-1] == _my_rank){
				char* strc = new char[strlen(s.c_str())+1];
				strcpy(strc, s.c_str());
				AddVertex(strc, xxx-1);
				delete strc;
			}
			xxx += 1;
		}

		
		vector<vector<int> > V2D(_num_workers), V2DR(_num_workers);
		for(int i=0; i<_num_workers; ++i){
			if (i == _my_rank) 
				continue;
			V2D[i] = v2degree; // edge.size()
			V2DR[i] = v2degree_R; // edge.size() + local.size()
		}

		all_to_all_cat(V2D, V2DR);
		worker_barrier();

		for (int i=0; i<_num_workers; ++i){
			for (int j=0; j<V2D[i].size(); ++j){
				v2degree[j]   += V2D[i][j];
				v2degree_R[j] += V2DR[i][j]; 
			}
		}


		BdV = all_sum(BoundNum);
		totalV = all_sum(vertexes.size());
		ALLVertices = totalV;
		long long totalE = all_sum_LL(totalEdges);
		long long cutE   = accumulate(v2degree.begin(), v2degree.end(), 0);
		

		VertexCut_CutEdge(200);
		
		ActiveV[_my_rank] = vertexes.size();

		blockInit(vertexes, blocks);


		float t = omp_get_wtime();

		Part2hop();
		
		float t1 = omp_get_wtime();

		cout<<"rank: "<<_my_rank<<"   time: "<<t1-t<<" s   "<<endl;
		
		// string new_filename = params.output_path + "In_"+ to_String(_my_rank);
		// DHI_Store(new_filename);		

		BoundGraph();

		worker_barrier();
		
		if (_my_rank == 0){

			float t2 = omp_get_wtime();
			
			Core2hop();
			
			float t3 = omp_get_wtime();

			cout<<"Core time: "<<t3-t2<<endl;

			for (int i=0; i<conB.size(); ++i){
				if (i < CopyNum) 
					aaa += Clabel[i].size();
				else{
					for (int j=0; j<Clabel[i].size(); ++j)
						if (Cflags[i][j] == 0)
							aaa += 1;
				}
			}
			cout<<"Total  size: "<<float(aaa*4)/(1024*1024)<<" MB"<<endl;

			// new_filename = params.output_path + "Bd";
			// DHBound_Store(new_filename);

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

