#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<map>
#include<set>
#include<cstdlib>

using namespace std;

struct Node{
	int id; //结点编号
	multimap<double,Node*> *inNodes; //结点入结点
	multimap<double,Node*> *outNodes; //结点出结点
	int inDegree; //结点入度
	int outDegree; //结点出度
	bool onCurrentPath; //是否出现在当前路径中——用于Optimizing Simpath
	double pp; //记录从种子结点到该结点的路径概率——用于计算影响传播值
	double inf; //记录从种子结点到该结点的影响——用于计算影响传播值
};

/*将字符串转化为int*/
inline unsigned int strToInt(string s) {
	unsigned int i;
	istringstream myStream(s);

	if (myStream>>i) {
		return i;
	} else {
		cout << "String " << s << " is not a number." << endl;
		return atoi(s.c_str());
	}
	return i;
}

/*将字符串转化为double*/
inline double strToDouble(string s) {
    return atof(s.c_str());
    double i;
    istringstream myStream(s);

    if (myStream>>i) {
        return i;
    } else {
        cout << "String " << s << " is not a float." << endl;
		return atof(s.c_str());
    }
    return i;
}

class SimPath{
private:
	map<int,Node*> *nodes; //存放所有结点信息（包括出度、入度等）
	int edges; //保存边的总数
	set<int> nodeIds; //存放所有顶点id
	set<int> vcNodeIds; //存放最小顶点覆盖的顶点id
	set<int> nvcNodeIds; //存放不在最小顶点覆盖中的顶点id
	set<int> seedNodeIds; //存放种子结点id即最终结果
	map<int,double> nvcInf; //存放不在VC中的结点的影响值
	map<int,double> infSUp; //用于simPathSpreadNormal 存放inf{V-x}(S) x属于U
	map<int,double> infSDown;  //存放inf{V-S}(x)
	double pruneVal; //路径剪枝阈值
	int topL; //设置候选结点数目
	int k; //设置需要选出的种子节点集合
public:
	SimPath();//构造函数
	void setParameter(double _pruneVal,int _topL, int _k);//设置参数
	void readDataSets(string file); //读取文件构造结点集合
	void findVertexCover();//第一轮迭代前寻找最小顶点覆盖
	void simPathCore(); //算法核心，依次寻找k个种子结点保存在seedNodeIds中
	double simPathSpreadFirst(int id); //第一轮迭代计算结点覆盖集中的点
	double simPathSpreadNormal(set<int> S,set<int> U); //其他轮迭代计算影响值Inf(S)同时计算Inf{V-x}(S)，x属于U
	double backtrackNormal(int id,set<int> S,set<int> U); //计算Inf{V-S+u-x}(u),x属于U
	double backtrackSimple(int id,set<int> S); //计算Inf{V-S}(x)
};

/*构造函数*/
SimPath::SimPath(){
	nodes = new map<int,Node*>();
	pruneVal = 0.001;
	topL = 4;
}

/*设置参数*/
void SimPath::setParameter(double _pruneVal,int _topL, int _k){
	pruneVal = _pruneVal;
	topL = _topL;
	k = _k;
}

/*读取文件构造结点集合*/
void SimPath::readDataSets(string file){
    cout << "Reading \"" << file << "\"data file... " <<endl;
	ifstream dataFile (file.c_str(), ios::in);

    if (dataFile.is_open()) {
		int id1,id2,b; //读取结点编号1，结点编号2以及影响因子b
		string split = " \t";	

        while (!dataFile.eof() )	{
            string line;
			getline (dataFile,line);
			if (line.empty()) {continue;} //读入空行跳过
			edges++;
			if (edges == 0) {continue;} //读入第一行跳过
			string::size_type pos = line.find_first_of(split);
			int	prevpos = 0;

			// 获取第一个结点编号
			int id1 = strToInt(line.substr(prevpos, pos-prevpos));

			// 获取第二个结点编号
			prevpos = line.find_first_not_of(split, pos);
			pos = line.find_first_of(split, prevpos);
			int id2 = strToInt(line.substr(prevpos, pos-prevpos));

			// 获取边权重影响因子
			double b = 0;
			prevpos = line.find_first_not_of(split, pos);
			pos = line.find_first_of(split, prevpos);
			if (pos == string::npos) 
				b = strToDouble(line.substr(prevpos));
			else
				b = strToDouble(line.substr(prevpos, pos-prevpos));

			// 两结点相等或b=0均为无效边
			if (id1==id2||b == 0) {
                edges--; 
                continue;
			}

			// 根据id获取相应Node
            Node *node1 = NULL;
            Node *node2 = NULL;
			
            pair<set<int>::iterator,bool> r1, r2; //记录插入结果从而判断该结点是否出现过
			r1 = nodeIds.insert(id1);
            r2 = nodeIds.insert(id2);
			
            if (r1.second == true) // id1插入成功，说明结点未出现过
			{   
                node1 = new Node();
				node1->id = id1;
				node1->inNodes = new multimap<double, Node *>();
				node1->outNodes = new multimap<double, Node *>();
                nodes->insert(make_pair(id1, node1));
            }else{ // id1插入失败，说明结点出现过
                map<int, Node *>::iterator it1 = nodes->find(id1);
                if (it1!= nodes->end()) 
                    node1 = it1->second;
                else 
                    cout << "Find Node" << id1 <<"Error!"<<endl;
            } 

            if (r2.second == true) // id2插入成功，说明结点未出现过
			{   
                node2 = new Node();
				node2->id = id2;
				node2->inNodes = new multimap<double, Node *>();
				node2->outNodes = new multimap<double, Node *>();
                nodes->insert(make_pair(id2, node2));
            }else{ // id2插入失败，说明结点出现过
                map<int, Node *>::iterator it2 = nodes->find(id2);
                if (it2!= nodes->end()) 
                    node2 = it2->second;
                else 
                     cout << "Find Node" << id2 <<"Error!"<<endl;
            } 
			
			//更新两个结点信息
			node1->outNodes->insert(make_pair(b, node2));
			node2->inNodes->insert(make_pair(b, node1));
        } 
       dataFile.close();

	} else {
        cout << "Cannot open file! File Name is \"" << file <<"\""<<endl;
		exit(1);
    } 

    //保存各结点入度出度数用以计算顶点覆盖
	map<int,Node*>::iterator it = nodes->begin();
	while( it != nodes->end() )
	{
		Node *node = it->second;
		node->inDegree = node->inNodes->size();
		node->outDegree = node->outNodes->size();
		it++;
	}

	//输出提示信息
    cout << "NodeIds size: " << nodeIds.size()<< endl;
	cout << "Nodes size: " << nodes->size() << endl;
    cout << "Edges num: " << edges << endl;
	cout << "Read Data Set Complete." << endl;
}

/*使用最大入度贪心规则找到结点覆盖*/
void SimPath::findVertexCover(){
	cout << "Calculating Vertex Cover..." << endl;
	
	//按序根据入度大小建立结点队列（map默认从小到大）
    multimap<int, Node*> degreeQueue;
	map<int, Node*>::iterator it = nodes->begin();
	while( it != nodes->end() )
	{
        Node *node = it->second;  
		degreeQueue.insert(make_pair(node->inDegree, node));
		it++;
    } 
    
	//按照最大入度贪心原则，从大到小顺序读取并依次加入VC or NonVC
    multimap<int, Node *>::reverse_iterator rit = degreeQueue.rbegin(); 
	while (vcNodeIds.size() + nvcNodeIds.size() < nodeIds.size()) {
        Node *node = rit->second;
        ++rit;
		int id = node->id;

        //若当前node还未被加入NonVC集合中则加入VC集合（贪心原则）
		if (nvcNodeIds.find(id) == nvcNodeIds.end()) {
			vcNodeIds.insert(id);
        } 
        //更新与之相邻的结点信息，即将相邻边去除
		multimap<double,Node*> *in = node->inNodes;
		multimap<double,Node*>::iterator itIn = in->begin();
		while( itIn != in->end() )
		{
			Node *inNode = itIn->second;
			inNode->outDegree--;           
			if ( inNode->outDegree == 0 && inNode->inDegree == 0) {
				if ( vcNodeIds.find(inNode->id) == vcNodeIds.end()) {
					nvcNodeIds.insert(inNode->id);
                } 
            } 
			itIn++;
        }
		multimap<double,Node*> *out = node->outNodes;
		multimap<double,Node*>::iterator itOut = out->begin();
		while( itOut != out->end() )
		{
			Node *outNode = itOut->second;
			outNode->inDegree--;          
			if ( outNode->inDegree == 0 && outNode->outDegree == 0) {
				if ( vcNodeIds.find(outNode->id) == vcNodeIds.end()) {
					nvcNodeIds.insert(outNode->id);
                } 
            } 
			itOut++;
        }
	}

	//输出提示信息
	cout << "VC size: " << vcNodeIds.size() << endl;
    cout << "NonVC size: " << nvcNodeIds.size() << endl;
}

/*算法核心，依次寻找k个种子结点保存在seedNodeIds中*/
void SimPath::simPathCore(){
	cout << "Start Calculating..." << endl;

    multimap<double, int> celfQueue; //根据marginal gain按序保存结点id
	double totalInf = 0;//保存总影响值
	int count = 0;//记录每个种子结点需要调用Spread次数

    seedNodeIds.clear();

    //初始化所有不在结点覆盖集合中的点的影响值为0
	//map<int,double> nvcInf;
	for(set<int>::iterator it = nvcNodeIds.begin();it!=nvcNodeIds.end();it++){
		nvcInf.insert(make_pair(*it,1.0));
	}
    
    // 计算所有在结点覆盖集合中的点的影响值，同时根据Theorem2更新不在结点覆盖集合中的点的影响值
	for(set<int>::iterator it = vcNodeIds.begin();it!=vcNodeIds.end();it++){
		int id = *it;
        double inf = simPathSpreadFirst(id);
        celfQueue.insert(make_pair(inf, id));
		count++;
	}

    // 计算所有不在结点覆盖集合中的点的影响值
	for (map<int, double>::iterator it = nvcInf.begin(); it != nvcInf.end(); it++) {
        int id = it->first;
        double inf = it->second; //根据Theorem 2
        celfQueue.insert(make_pair(inf, id));
		count++;
    }
    
	//CELF队列中MG最大的即为第一个Seed Node
	multimap<double,int>::iterator it = celfQueue.end();
	it--;
	totalInf = it->first;
	seedNodeIds.insert(it->second);
	celfQueue.erase(it);
	cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<*seedNodeIds.begin()<<", mg = "<<totalInf<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;

    // 首种子结点选择完成，按照Look Ahead Optimization原则进行接下来种子结点选取
    count = 0;
    set<int> examinedNodeIds; //基于当前S重新计算Margin Gain的结点
	set<int> topLNodeIds; // CELF队列中MG最大的L个结点
    int realTopL;  

    //随着S更新反复计算Margin Gain选择种子结点
    while (seedNodeIds.size() < k) {
		topLNodeIds.clear();
		if (celfQueue.size() < topL) {
			realTopL = celfQueue.size()-1;
		}else{
			realTopL = topL;
		}
		multimap<double, int>::iterator itCelf = celfQueue.end();
		itCelf--;

		int compareId;
		double compareMg;
		bool findExamined = false;
		bool topOneExamined = false;
        for (int i = 0; i < realTopL; i++) {
            int xId = itCelf->second;
			// 当前结点未被重新计算过
            if (examinedNodeIds.find(xId) == examinedNodeIds.end()) {
                topLNodeIds.insert(xId);
				examinedNodeIds.insert(xId);//即将被重新计算！！！
                count++;
            } else {
                // 在topL中找到未被重新计算的结点
                findExamined = true;
                compareId = xId;
                compareMg = itCelf->first;
                if (i == 0) {
                    topOneExamined = true; 
                }
                break;
            } 
            itCelf--; 
        }

        // 候选结点集合第一个就被重新计算，则直接加入种子集合中,更新总影响值，更新CELF队列
        if (topOneExamined == true) {
            totalInf += compareMg;
			seedNodeIds.insert(compareId);
			cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<compareId<<", mg = "<<compareMg<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;
            count = 0;
			examinedNodeIds.clear(); //种子集合变化，所有结点均需重新计算
            celfQueue.erase(itCelf);
            continue;
        }
		// 若Top L中的所有结点都未重新计算，需要与下一个L的最大结点进行比较
        if (findExamined == false) {
            itCelf--;
            compareId = itCelf->second;
            compareMg = itCelf->first;
        }

        // 将所有TopL候选结点从CELF队列中删除
		for (int i = 0; i < topLNodeIds.size(); i++) {
             multimap<double, int>::iterator iter = celfQueue.end();
             iter--;
             celfQueue.erase(iter);
        }
        
		
		int maxId = 0;
		double maxMg = 0;

		//计算 inf{V-x}(S)
		infSUp.clear();
		simPathSpreadNormal(seedNodeIds, topLNodeIds); //更新infSUp
        
		//计算 inf{V-S}(x)并将Margin Gain计算结果保存其中
		infSDown.clear();
		for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
            int xId = *it;
			double xInf = backtrackSimple(xId, seedNodeIds); //计算inf{V-S}(x)
            // compute node x's marginal gain;
            double xMg = infSUp.find(xId)->second + xInf - totalInf;
            infSDown.insert(make_pair(xId, xMg)); // pair(ID, newMG)
			if(xMg > maxMg)
			{
				maxId = xId;
				maxMg = xMg;
			}
        }
        if (maxMg >= compareMg) {  // TopL中存在经过重新计算且最大MG的结点，加入种子集合中,更新总影响值，更新CELF队列
            totalInf += maxMg;
            seedNodeIds.insert(maxId);
			cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<maxId<<", mg = "<<maxMg<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;
            count = 0;
            examinedNodeIds.clear();

            //将TopL中其他所有结点重新插入CELF队列
            for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
                int xId = *it;
				if ( xId != maxId) {
                    double xMg = infSDown.find(xId)->second;
                    celfQueue.insert(make_pair(xMg, xId));
                } 
            } 
            continue;
        } else {   // TopL中所有结点都不如下一个L的最大结点的MG大，则将所有结点重新插入CELF队列
            for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
                int xId = *it;
                double xMg = infSDown.find(xId)->second;
                celfQueue.insert(make_pair(xMg, xId));
            }

        }
	}
}

/*计算编号为id的结点的影响值，同时根据Theorem2更新不在结点覆盖集合中的点的影响值*/
double SimPath::simPathSpreadFirst(int id){

	Node *uNode = nodes->find(id)->second; // 根据id找到结点u
    uNode->inf = 1;
    uNode->pp = 1;
	uNode->onCurrentPath = true;
	int uId = uNode->id;
	int uOutNum = uNode->outNodes->size();

	//构造U（U = {x|x不在结点覆盖集合中且是u的入结点）
	set<int> UIds;
	map<int, Node *> UNodes;
	map<int, double> bMap; // 记录u到其出结点的影响因子值
	multimap<double, Node*> *in = uNode->inNodes; 
	for (multimap<double, Node *>::iterator it = in->begin(); it != in->end(); it++) {
		Node *v = it->second;
		int vId = v->id;
		if (nvcNodeIds.find(vId) != nvcNodeIds.end()) {
			UIds.insert(vId);
			bMap.insert(make_pair(vId, it->first));
			v->inf = 1; //记录u在V-v上的影响值
			v->pp = 1;  //记录u到v的路径概率
			v->onCurrentPath = false;
			UNodes.insert(make_pair(vId, v));
        } 
	} 

    double inf = 1; // 记录u结点的影响值
    double pp = 1; // 记录路径概率 
	vector<Node *> Q; //记录当前路径结点
    set<int> idQ; //记录当前路径结点id
    map<int, set<int> > D; //记录结点被访问过的邻居

    Q.push_back(uNode);
    idQ.insert(uId);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // 进行Forward即寻找路径
            Node *xNode = Q.back(); // 获取当前路径最后一个结点
            int xId = xNode->id; 
            if (xId == lastNodeId) { //当前结点路径搜索完毕跳出while循环进行Backtrace
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // 获取u到当前结点x的路径概率

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //当前结点x没有出结点则跳过
            }

			//遍历当前结点x的所有出结点以寻找路径
            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y就是u时跳过
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y 在路径上已经存在，跳过防止形成环
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y 已经访问过，跳过防止重复计算
                    continue;
                }

                double ppNext = pp * it->first; // 找到新传播路径，更新路径概率

                //判断是否达到阈值，若低于阈值则停止搜索路径
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); //将y加入x的已访问列表中
                    continue;
                } 

				//若高于阈值则沿当前结点继续向下找
                inf += ppNext;  //更新u的总影响值
                yNode->pp = ppNext; // 保存u至结点y的路径概率值
				yNode->onCurrentPath = true; //更改y在当前路径中
                Q.push_back(yNode); // 将y结点加至当前路径结点Q中
                idQ.insert(yId); // 将y结点编号加至当前路径结点idQ中
                D[xId].insert(yId); // 标记x结点的邻居y已被访问过

                // 更新影响值Inf{V-v}(u)
                for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
					if (it->second->onCurrentPath == false) {
                        it->second->inf += ppNext;
                    }
                }
                break; //跳出当前for循环，将y当作最后一个结点寻找路径
            } 
		}

        //当前路径最后一个结点已搜索完毕，通过Backtrack继续寻找新路径。
        Node *lastNode = Q.back();
		int lastId = lastNode->id;

        if (Q.size() == 1) {   
            if(lastId != uId) {
                cout << "The only remaining node in Q is: " << lastId << ", but not: " << uId << endl;
                exit(1);
            } 
            if(D[lastId].size() < uOutNum) {
                continue;
            } 
        }

        //将最后一个结点从当前路径删除
		map<int, Node*>::iterator it = UNodes.find(lastId);
        if (it != UNodes.end()) {
			it->second->onCurrentPath = false;
        }
        idQ.erase(lastId);
        D.erase(lastId);
        Q.pop_back();
    }


    // 更新U的影响值（U = {x|x不在结点覆盖集合中且是u的入结点）
    for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
        int xId = it->first;
        Node *xNode = it->second;
        nvcInf[xId] += bMap[xId] * xNode->inf; //根据Theorem2
    } 

    return inf;
}
/*计算影响值Inf(S)，同时计算Inf{V-x}(S)，x属于U*/
double SimPath::simPathSpreadNormal(set<int> S, set<int> U)
{
	double inf = 0;
	infSUp.clear();
	for (set<int>::iterator it = U.begin(); it != U.end(); ++it){
		infSUp.insert(make_pair(*it, 0.0));
    } 
        
    for (set<int>::iterator it = S.begin(); it != S.end(); ++it) {
		int u = *it;
		inf +=  backtrackNormal(u,S,U);
	} 
	return inf;
}

/*计算Inf{V-S+u}(u)，同时更新Inf{V-S+u-x}(u)，x属于U*/
double SimPath::backtrackNormal(int id, set<int> S,set<int> U){
	Node *uNode = nodes->find(id)->second; // 根据id找到结点u
	int uId = uNode->id;
	int uOutNum = uNode->outNodes->size();
    uNode->inf = 1;
    uNode->pp = 1;
	uNode->onCurrentPath = true;

	//构造U
	set<int> UIds;
	map<int, Node *> UNodes;
	for (set<int>::iterator it = U.begin(); it != U.end(); it++) {
		Node *v = nodes->find(*it)->second;
		int vId = v->id;
		v->inf = 1; //记录u在V-v上的影响值
		v->pp = 1;
		v->onCurrentPath = false;
		UNodes.insert(make_pair(vId, v));
		UIds.insert(vId);
	} 

    double inf = 1; // cov(u), initially 1 (counting itself first)
    double pp = 1; // path prob. 
	vector<Node *> Q; //记录当前路径结点
    set<int> idQ; //记录当前路径结点id
    map<int, set<int> > D; //记录结点被访问过的邻居

    Q.push_back(uNode);
    idQ.insert(id);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // FORWARD starts here!!!
            Node *xNode = Q.back(); // 获取当前路径最后一个结点
            int xId = xNode->id; 
    
            if (xId == lastNodeId) {
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // get the current path prob. till this node

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //当前结点x没有邻居
            }

            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y就是u时
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y 在路径上已经存在，跳过防止形成环
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y 已经计算过，跳过防止重复计算
                    continue;
				} else if (S.find(yId)!=S.end()){  //y不在W（V-S）中
					continue;
				}

                double ppNext = pp * it->first; // pp = pp * b(x,y)
                // 判断是否达到阈值
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); // y is explored
                    continue;
                } 

				//若高于阈值未被剪掉则说明找到一条新路径
                inf += ppNext;  //更新u的总影响值
                yNode->pp = ppNext; // 更新u至结点y的路径概率值
				yNode->onCurrentPath = true; //此时y在当前路径中 -->可改

                Q.push_back(yNode); // 将y结点加至当前路径结点Q中
                idQ.insert(yId); // 将y结点id加至当前路径结点idQ中
                D[xId].insert(yId); // 记录x结点的邻居y已被访问过

                // 更新影响值（V-v）
                for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
					if (it->second->onCurrentPath == false) {
                        it->second->inf += ppNext;
                    }
                }
                break;
            } 

		}

         //当前路径最后一个结点已搜索完毕，通过Backtrack继续寻找新路径。
        Node *lastNode = Q.back();
		int lastId = lastNode->id;

        if (Q.size() == 1) {   
            if(lastId != uId) {
                cout << "The only remaining node in Q is: " << lastId << ", but not: " << uId << endl;
                exit(1);
            } 
            if(D[lastId].size() < uOutNum) {
                continue;
            } 
        }

        //将最后一个结点从当前路径删除
		map<int, Node*>::iterator it = UNodes.find(lastId);
        if (it != UNodes.end()) {
			it->second->onCurrentPath = false;
        }
        idQ.erase(lastId);
        D.erase(lastId);
        Q.pop_back();

    }


    // update partial coverage for U 
    for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
        int xId = it->first;
        Node *xNode = it->second;
		infSUp[xId] += xNode->inf; 
    } 
    return inf;
}

/*根据编号id计算Inf{V-S}(x)的值*/
double SimPath::backtrackSimple(int id, set<int> S){
	Node *uNode = nodes->find(id)->second; // 根据id找到结点u
	int uId = uNode->id;
    uNode->inf = 1;
    uNode->pp = 1;

    double inf = 1; // cov(u), initially 1 (counting itself first)
    double pp = 1; // path prob. 
	vector<Node *> Q; //记录当前路径结点
    set<int> idQ; //记录当前路径结点id
    map<int, set<int> > D; //记录结点被访问过的邻居

    Q.push_back(uNode);
    idQ.insert(id);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // FORWARD starts here!!!
            Node *xNode = Q.back(); // 获取当前路径最后一个结点
            int xId = xNode->id; 
    
            if (xId == lastNodeId) {
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // get the current path prob. till this node

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //当前结点x没有邻居
            }

            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y就是u时
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y 在路径上已经存在，跳过防止形成环
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y 已经计算过，跳过防止重复计算
                    continue;
				} else if (S.find(yId)!=S.end()){  //y不在W（V-S）中
					D[xId].insert(yId);
					continue;
				}

                double ppNext = pp * it->first; // pp = pp * b(x,y)
                // 判断是否达到阈值
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); // y is explored
                    continue;
                } 

				//若高于阈值未被剪掉则说明找到一条新路径
                inf += ppNext;  //更新u的总影响值
                yNode->pp = ppNext; // 更新u至结点y的路径概率值
				yNode->onCurrentPath = true; //此时y在当前路径中 -->可改

                Q.push_back(yNode); // 将y结点加至当前路径结点Q中
                idQ.insert(yId); // 将y结点id加至当前路径结点idQ中
                D[xId].insert(yId); // 记录x结点的邻居y已被访问过

                break;
            } 

		}

        //进行Backtrace，即寻找下一个可探索路径的点
        Node *lastNode = Q.back(); 
		int lastId = lastNode->id;
        idQ.erase(lastId);
        D.erase(lastId);
        Q.pop_back();
    }

    return inf;
}

int main()
{
	SimPath *s = new SimPath();
	double _pruneVal = 0.001; 
	int _topL = 4;
	int _k = 5;
	string file = "hep_LT2.inf";
	s->setParameter(_pruneVal,_topL,_k);
	s->readDataSets(file);
	s->findVertexCover();
	s->simPathCore();
	return 0;
}