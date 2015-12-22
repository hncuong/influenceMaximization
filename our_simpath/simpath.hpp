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
	int id; //½áµã±àºÅ
	multimap<double,Node*> *inNodes; //½áµãÈë½áµã
	multimap<double,Node*> *outNodes; //½áµã³ö½áµã
	int inDegree; //½áµãÈë¶È
	int outDegree; //½áµã³ö¶È
	bool onCurrentPath; //ÊÇ·ñ³öÏÖÔÚµ±Ç°Â·¾¶ÖÐ¡ª¡ªÓÃÓÚOptimizing Simpath
	double pp; //¼ÇÂ¼´ÓÖÖ×Ó½áµãµ½¸Ã½áµãµÄÂ·¾¶¸ÅÂÊ¡ª¡ªÓÃÓÚ¼ÆËãÓ°Ïì´«²¥Öµ
	double inf; //¼ÇÂ¼´ÓÖÖ×Ó½áµãµ½¸Ã½áµãµÄÓ°Ïì¡ª¡ªÓÃÓÚ¼ÆËãÓ°Ïì´«²¥Öµ
};

inline string intToStr(int i) {
	stringstream ss;
	ss << i;
	return ss.str();
}

/*½«×Ö·û´®×ª»¯Îªint*/
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

/*½«×Ö·û´®×ª»¯Îªdouble*/
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
	map<int,Node*> *nodes; //´æ·ÅËùÓÐ½áµãÐÅÏ¢£¨°üÀ¨³ö¶È¡¢Èë¶ÈµÈ£©
	int edges; //±£´æ±ßµÄ×ÜÊý
	set<int> nodeIds; //´æ·ÅËùÓÐ¶¥µãid
	set<int> vcNodeIds; //´æ·Å×îÐ¡¶¥µã¸²¸ÇµÄ¶¥µãid
	set<int> nvcNodeIds; //´æ·Å²»ÔÚ×îÐ¡¶¥µã¸²¸ÇÖÐµÄ¶¥µãid
	set<int> seedNodeIds; //´æ·ÅÖÖ×Ó½áµãid¼´×îÖÕ½á¹û
	map<int,double> nvcInf; //´æ·Å²»ÔÚVCÖÐµÄ½áµãµÄÓ°ÏìÖµ
	map<int,double> infSUp; //ÓÃÓÚsimPathSpreadNormal ´æ·Åinf{V-x}(S) xÊôÓÚU
	map<int,double> infSDown;  //´æ·Åinf{V-S}(x)
	double pruneVal; //Â·¾¶¼ôÖ¦ãÐÖµ
	int topL; //ÉèÖÃºòÑ¡½áµãÊýÄ¿
	int k; //ÉèÖÃÐèÒªÑ¡³öµÄÖÖ×Ó½Úµã¼¯ºÏ
public:
	SimPath();//¹¹Ôìº¯Êý
	void setParameter(double _pruneVal,int _topL, int _k);//ÉèÖÃ²ÎÊý
	void readDataSets(string file); //¶ÁÈ¡ÎÄ¼þ¹¹Ôì½áµã¼¯ºÏ
	void findVertexCover();//µÚÒ»ÂÖµü´úÇ°Ñ°ÕÒ×îÐ¡¶¥µã¸²¸Ç
	void simPathCore(); //Ëã·¨ºËÐÄ£¬ÒÀ´ÎÑ°ÕÒk¸öÖÖ×Ó½áµã±£´æÔÚseedNodeIdsÖÐ
	double simPathSpreadFirst(int id); //µÚÒ»ÂÖµü´ú¼ÆËã½áµã¸²¸Ç¼¯ÖÐµÄµã
	double simPathSpreadNormal(set<int> S,set<int> U); //ÆäËûÂÖµü´ú¼ÆËãÓ°ÏìÖµInf(S)Í¬Ê±¼ÆËãInf{V-x}(S)£¬xÊôÓÚU
	double backtrackNormal(int id,set<int> S,set<int> U); //¼ÆËãInf{V-S+u-x}(u),xÊôÓÚU
	double backtrackSimple(int id,set<int> S); //¼ÆËãInf{V-S}(x)
};

/*¹¹Ôìº¯Êý*/
SimPath::SimPath(){
	nodes = new map<int,Node*>();
	pruneVal = 0.001;
	topL = 4;
}

/*ÉèÖÃ²ÎÊý*/
void SimPath::setParameter(double _pruneVal,int _topL, int _k){
	pruneVal = _pruneVal;
	topL = _topL;
	k = _k;
}

/*¶ÁÈ¡ÎÄ¼þ¹¹Ôì½áµã¼¯ºÏ*/
void SimPath::readDataSets(string file){
    cout << "Reading \"" << file << "\"data file... " <<endl;
	ifstream dataFile (file.c_str(), ios::in);

    if (dataFile.is_open()) {
		int id1,id2,b; //¶ÁÈ¡½áµã±àºÅ1£¬½áµã±àºÅ2ÒÔ¼°Ó°ÏìÒò×Ób
		string split = " \t";	

        while (!dataFile.eof() )	{
            string line;
			getline (dataFile,line);
			if (line.empty()) {continue;} //¶ÁÈë¿ÕÐÐÌø¹ý
			edges++;
			if (edges == 0) {continue;} //¶ÁÈëµÚÒ»ÐÐÌø¹ý
			string::size_type pos = line.find_first_of(split);
			int	prevpos = 0;

			// »ñÈ¡µÚÒ»¸ö½áµã±àºÅ
			int id1 = strToInt(line.substr(prevpos, pos-prevpos));

			// »ñÈ¡µÚ¶þ¸ö½áµã±àºÅ
			prevpos = line.find_first_not_of(split, pos);
			pos = line.find_first_of(split, prevpos);
			int id2 = strToInt(line.substr(prevpos, pos-prevpos));

			// »ñÈ¡±ßÈ¨ÖØÓ°ÏìÒò×Ó
			double b = 0;
			prevpos = line.find_first_not_of(split, pos);
			pos = line.find_first_of(split, prevpos);
			if (pos == string::npos) 
				b = strToDouble(line.substr(prevpos));
			else
				b = strToDouble(line.substr(prevpos, pos-prevpos));

			// Á½½áµãÏàµÈ»òb=0¾ùÎªÎÞÐ§±ß
			if (id1==id2||b == 0) {
                edges--; 
                continue;
			}

			// ¸ù¾Ýid»ñÈ¡ÏàÓ¦Node
            Node *node1 = NULL;
            Node *node2 = NULL;
			
            pair<set<int>::iterator,bool> r1, r2; //¼ÇÂ¼²åÈë½á¹û´Ó¶øÅÐ¶Ï¸Ã½áµãÊÇ·ñ³öÏÖ¹ý
			r1 = nodeIds.insert(id1);
            r2 = nodeIds.insert(id2);
			
            if (r1.second == true) // id1²åÈë³É¹¦£¬ËµÃ÷½áµãÎ´³öÏÖ¹ý
			{   
                node1 = new Node();
				node1->id = id1;
				node1->inNodes = new multimap<double, Node *>();
				node1->outNodes = new multimap<double, Node *>();
                nodes->insert(make_pair(id1, node1));
            }else{ // id1²åÈëÊ§°Ü£¬ËµÃ÷½áµã³öÏÖ¹ý
                map<int, Node *>::iterator it1 = nodes->find(id1);
                if (it1!= nodes->end()) 
                    node1 = it1->second;
                else 
                    cout << "Find Node" << id1 <<"Error!"<<endl;
            } 

            if (r2.second == true) // id2²åÈë³É¹¦£¬ËµÃ÷½áµãÎ´³öÏÖ¹ý
			{   
                node2 = new Node();
				node2->id = id2;
				node2->inNodes = new multimap<double, Node *>();
				node2->outNodes = new multimap<double, Node *>();
                nodes->insert(make_pair(id2, node2));
            }else{ // id2²åÈëÊ§°Ü£¬ËµÃ÷½áµã³öÏÖ¹ý
                map<int, Node *>::iterator it2 = nodes->find(id2);
                if (it2!= nodes->end()) 
                    node2 = it2->second;
                else 
                     cout << "Find Node" << id2 <<"Error!"<<endl;
            } 
			
			//¸üÐÂÁ½¸ö½áµãÐÅÏ¢
			node1->outNodes->insert(make_pair(b, node2));
			node2->inNodes->insert(make_pair(b, node1));
        } 
       dataFile.close();

	} else {
        cout << "Cannot open file! File Name is \"" << file <<"\""<<endl;
		exit(1);
    } 

    //±£´æ¸÷½áµãÈë¶È³ö¶ÈÊýÓÃÒÔ¼ÆËã¶¥µã¸²¸Ç
	map<int,Node*>::iterator it = nodes->begin();
	while( it != nodes->end() )
	{
		Node *node = it->second;
		node->inDegree = node->inNodes->size();
		node->outDegree = node->outNodes->size();
		it++;
	}

	//Êä³öÌáÊ¾ÐÅÏ¢
    cout << "NodeIds size: " << nodeIds.size()<< endl;
	cout << "Nodes size: " << nodes->size() << endl;
    cout << "Edges num: " << edges << endl;
	cout << "Read Data Set Complete." << endl;
}

/*Ê¹ÓÃ×î´óÈë¶ÈÌ°ÐÄ¹æÔòÕÒµ½½áµã¸²¸Ç*/
void SimPath::findVertexCover(){
	cout << "Calculating Vertex Cover..." << endl;
	
	//°´Ðò¸ù¾ÝÈë¶È´óÐ¡½¨Á¢½áµã¶ÓÁÐ£¨mapÄ¬ÈÏ´ÓÐ¡µ½´ó£©
    multimap<int, Node*> degreeQueue;
	map<int, Node*>::iterator it = nodes->begin();
	while( it != nodes->end() )
	{
        Node *node = it->second;  
		degreeQueue.insert(make_pair(node->inDegree, node));
		it++;
    } 
    
	//°´ÕÕ×î´óÈë¶ÈÌ°ÐÄÔ­Ôò£¬´Ó´óµ½Ð¡Ë³Ðò¶ÁÈ¡²¢ÒÀ´Î¼ÓÈëVC or NonVC
    multimap<int, Node *>::reverse_iterator rit = degreeQueue.rbegin(); 
	while (vcNodeIds.size() + nvcNodeIds.size() < nodeIds.size()) {
        Node *node = rit->second;
        ++rit;
		int id = node->id;

        //Èôµ±Ç°node»¹Î´±»¼ÓÈëNonVC¼¯ºÏÖÐÔò¼ÓÈëVC¼¯ºÏ£¨Ì°ÐÄÔ­Ôò£©
		if (nvcNodeIds.find(id) == nvcNodeIds.end()) {
			vcNodeIds.insert(id);
        } 
        //¸üÐÂÓëÖ®ÏàÁÚµÄ½áµãÐÅÏ¢£¬¼´½«ÏàÁÚ±ßÈ¥³ý
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

	//Êä³öÌáÊ¾ÐÅÏ¢
	cout << "VC size: " << vcNodeIds.size() << endl;
    cout << "NonVC size: " << nvcNodeIds.size() << endl;
}

/*Ëã·¨ºËÐÄ£¬ÒÀ´ÎÑ°ÕÒk¸öÖÖ×Ó½áµã±£´æÔÚseedNodeIdsÖÐ*/
void SimPath::simPathCore(){
	cout << "Start Calculating..." << endl;

    multimap<double, int> celfQueue; //¸ù¾Ýmarginal gain°´Ðò±£´æ½áµãid
	double totalInf = 0;//±£´æ×ÜÓ°ÏìÖµ
	int count = 0;//¼ÇÂ¼Ã¿¸öÖÖ×Ó½áµãÐèÒªµ÷ÓÃSpread´ÎÊý

    seedNodeIds.clear();

    //³õÊ¼»¯ËùÓÐ²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐµÄµãµÄÓ°ÏìÖµÎª0
	//map<int,double> nvcInf;
	for(set<int>::iterator it = nvcNodeIds.begin();it!=nvcNodeIds.end();it++){
		nvcInf.insert(make_pair(*it,1.0));
	}
    
    // ¼ÆËãËùÓÐÔÚ½áµã¸²¸Ç¼¯ºÏÖÐµÄµãµÄÓ°ÏìÖµ£¬Í¬Ê±¸ù¾ÝTheorem2¸üÐÂ²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐµÄµãµÄÓ°ÏìÖµ
	for(set<int>::iterator it = vcNodeIds.begin();it!=vcNodeIds.end();it++){
		int id = *it;
        double inf = simPathSpreadFirst(id);
        celfQueue.insert(make_pair(inf, id));
		count++;
	}

    // ¼ÆËãËùÓÐ²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐµÄµãµÄÓ°ÏìÖµ
	for (map<int, double>::iterator it = nvcInf.begin(); it != nvcInf.end(); it++) {
        int id = it->first;
        double inf = it->second; //¸ù¾ÝTheorem 2
        celfQueue.insert(make_pair(inf, id));
		count++;
    }
    
	//CELF¶ÓÁÐÖÐMG×î´óµÄ¼´ÎªµÚÒ»¸öSeed Node
	multimap<double,int>::iterator it = celfQueue.end();
	it--;
	totalInf = it->first;
	seedNodeIds.insert(it->second);
	celfQueue.erase(it);
	cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<*seedNodeIds.begin()<<", mg = "<<totalInf<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;

    // Ê×ÖÖ×Ó½áµãÑ¡ÔñÍê³É£¬°´ÕÕLook Ahead OptimizationÔ­Ôò½øÐÐ½ÓÏÂÀ´ÖÖ×Ó½áµãÑ¡È¡
    count = 0;
    set<int> examinedNodeIds; //»ùÓÚµ±Ç°SÖØÐÂ¼ÆËãMargin GainµÄ½áµã
	set<int> topLNodeIds; // CELF¶ÓÁÐÖÐMG×î´óµÄL¸ö½áµã
    int realTopL;  

    //Ëæ×ÅS¸üÐÂ·´¸´¼ÆËãMargin GainÑ¡ÔñÖÖ×Ó½áµã
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
			// µ±Ç°½áµãÎ´±»ÖØÐÂ¼ÆËã¹ý
            if (examinedNodeIds.find(xId) == examinedNodeIds.end()) {
                topLNodeIds.insert(xId);
				examinedNodeIds.insert(xId);//¼´½«±»ÖØÐÂ¼ÆËã£¡£¡£¡
                count++;
            } else {
                // ÔÚtopLÖÐÕÒµ½Î´±»ÖØÐÂ¼ÆËãµÄ½áµã
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

        // ºòÑ¡½áµã¼¯ºÏµÚÒ»¸ö¾Í±»ÖØÐÂ¼ÆËã£¬ÔòÖ±½Ó¼ÓÈëÖÖ×Ó¼¯ºÏÖÐ,¸üÐÂ×ÜÓ°ÏìÖµ£¬¸üÐÂCELF¶ÓÁÐ
        if (topOneExamined == true) {
            totalInf += compareMg;
			seedNodeIds.insert(compareId);
			cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<compareId<<", mg = "<<compareMg<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;
            count = 0;
			examinedNodeIds.clear(); //ÖÖ×Ó¼¯ºÏ±ä»¯£¬ËùÓÐ½áµã¾ùÐèÖØÐÂ¼ÆËã
            celfQueue.erase(itCelf);
            continue;
        }
		// ÈôTop LÖÐµÄËùÓÐ½áµã¶¼Î´ÖØÐÂ¼ÆËã£¬ÐèÒªÓëÏÂÒ»¸öLµÄ×î´ó½áµã½øÐÐ±È½Ï
        if (findExamined == false) {
            itCelf--;
            compareId = itCelf->second;
            compareMg = itCelf->first;
        }

        // ½«ËùÓÐTopLºòÑ¡½áµã´ÓCELF¶ÓÁÐÖÐÉ¾³ý
		for (int i = 0; i < topLNodeIds.size(); i++) {
             multimap<double, int>::iterator iter = celfQueue.end();
             iter--;
             celfQueue.erase(iter);
        }
        
		
		int maxId = 0;
		double maxMg = 0;

		//¼ÆËã inf{V-x}(S)
		infSUp.clear();
		simPathSpreadNormal(seedNodeIds, topLNodeIds); //¸üÐÂinfSUp
        
		//¼ÆËã inf{V-S}(x)²¢½«Margin Gain¼ÆËã½á¹û±£´æÆäÖÐ
		infSDown.clear();
		for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
            int xId = *it;
			double xInf = backtrackSimple(xId, seedNodeIds); //¼ÆËãinf{V-S}(x)
            // compute node x's marginal gain;
            double xMg = infSUp.find(xId)->second + xInf - totalInf;
            infSDown.insert(make_pair(xId, xMg)); // pair(ID, newMG)
			if(xMg > maxMg)
			{
				maxId = xId;
				maxMg = xMg;
			}
        }
        if (maxMg >= compareMg) {  // TopLÖÐ´æÔÚ¾­¹ýÖØÐÂ¼ÆËãÇÒ×î´óMGµÄ½áµã£¬¼ÓÈëÖÖ×Ó¼¯ºÏÖÐ,¸üÐÂ×ÜÓ°ÏìÖµ£¬¸üÐÂCELF¶ÓÁÐ
            totalInf += maxMg;
            seedNodeIds.insert(maxId);
			cout<<"No "<<seedNodeIds.size()<<" seed: id = "<<maxId<<", mg = "<<maxMg<<", total= "<<totalInf<<", calculate margin gain num ="<<count<<endl;
            count = 0;
            examinedNodeIds.clear();

            //½«TopLÖÐÆäËûËùÓÐ½áµãÖØÐÂ²åÈëCELF¶ÓÁÐ
            for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
                int xId = *it;
				if ( xId != maxId) {
                    double xMg = infSDown.find(xId)->second;
                    celfQueue.insert(make_pair(xMg, xId));
                } 
            } 
            continue;
        } else {   // TopLÖÐËùÓÐ½áµã¶¼²»ÈçÏÂÒ»¸öLµÄ×î´ó½áµãµÄMG´ó£¬Ôò½«ËùÓÐ½áµãÖØÐÂ²åÈëCELF¶ÓÁÐ
            for (set<int>::iterator it = topLNodeIds.begin(); it != topLNodeIds.end(); it++) {
                int xId = *it;
                double xMg = infSDown.find(xId)->second;
                celfQueue.insert(make_pair(xMg, xId));
            }

        }
	}
}

/*¼ÆËã±àºÅÎªidµÄ½áµãµÄÓ°ÏìÖµ£¬Í¬Ê±¸ù¾ÝTheorem2¸üÐÂ²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐµÄµãµÄÓ°ÏìÖµ*/
double SimPath::simPathSpreadFirst(int id){

	Node *uNode = nodes->find(id)->second; // ¸ù¾ÝidÕÒµ½½áµãu
    uNode->inf = 1;
    uNode->pp = 1;
	uNode->onCurrentPath = true;
	int uId = uNode->id;
	int uOutNum = uNode->outNodes->size();

	//¹¹ÔìU£¨U = {x|x²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐÇÒÊÇuµÄÈë½áµã£©
	set<int> UIds;
	map<int, Node *> UNodes;
	map<int, double> bMap; // ¼ÇÂ¼uµ½Æä³ö½áµãµÄÓ°ÏìÒò×ÓÖµ
	multimap<double, Node*> *in = uNode->inNodes; 
	for (multimap<double, Node *>::iterator it = in->begin(); it != in->end(); it++) {
		Node *v = it->second;
		int vId = v->id;
		if (nvcNodeIds.find(vId) != nvcNodeIds.end()) {
			UIds.insert(vId);
			bMap.insert(make_pair(vId, it->first));
			v->inf = 1; //¼ÇÂ¼uÔÚV-vÉÏµÄÓ°ÏìÖµ
			v->pp = 1;  //¼ÇÂ¼uµ½vµÄÂ·¾¶¸ÅÂÊ
			v->onCurrentPath = false;
			UNodes.insert(make_pair(vId, v));
        } 
	} 

    double inf = 1; // ¼ÇÂ¼u½áµãµÄÓ°ÏìÖµ
    double pp = 1; // ¼ÇÂ¼Â·¾¶¸ÅÂÊ 
	vector<Node *> Q; //¼ÇÂ¼µ±Ç°Â·¾¶½áµã
    set<int> idQ; //¼ÇÂ¼µ±Ç°Â·¾¶½áµãid
    map<int, set<int> > D; //¼ÇÂ¼½áµã±»·ÃÎÊ¹ýµÄÁÚ¾Ó

    Q.push_back(uNode);
    idQ.insert(uId);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // ½øÐÐForward¼´Ñ°ÕÒÂ·¾¶
            Node *xNode = Q.back(); // »ñÈ¡µ±Ç°Â·¾¶×îºóÒ»¸ö½áµã
            int xId = xNode->id; 
            if (xId == lastNodeId) { //µ±Ç°½áµãÂ·¾¶ËÑË÷Íê±ÏÌø³öwhileÑ­»·½øÐÐBacktrace
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // »ñÈ¡uµ½µ±Ç°½áµãxµÄÂ·¾¶¸ÅÂÊ

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //µ±Ç°½áµãxÃ»ÓÐ³ö½áµãÔòÌø¹ý
            }

			//±éÀúµ±Ç°½áµãxµÄËùÓÐ³ö½áµãÒÔÑ°ÕÒÂ·¾¶
            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y¾ÍÊÇuÊ±Ìø¹ý
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y ÔÚÂ·¾¶ÉÏÒÑ¾­´æÔÚ£¬Ìø¹ý·ÀÖ¹ÐÎ³É»·
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y ÒÑ¾­·ÃÎÊ¹ý£¬Ìø¹ý·ÀÖ¹ÖØ¸´¼ÆËã
                    continue;
                }

                double ppNext = pp * it->first; // ÕÒµ½ÐÂ´«²¥Â·¾¶£¬¸üÐÂÂ·¾¶¸ÅÂÊ

                //ÅÐ¶ÏÊÇ·ñ´ïµ½ãÐÖµ£¬ÈôµÍÓÚãÐÖµÔòÍ£Ö¹ËÑË÷Â·¾¶
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); //½«y¼ÓÈëxµÄÒÑ·ÃÎÊÁÐ±íÖÐ
                    continue;
                } 

				//Èô¸ßÓÚãÐÖµÔòÑØµ±Ç°½áµã¼ÌÐøÏòÏÂÕÒ
                inf += ppNext;  //¸üÐÂuµÄ×ÜÓ°ÏìÖµ
                yNode->pp = ppNext; // ±£´æuÖÁ½áµãyµÄÂ·¾¶¸ÅÂÊÖµ
				yNode->onCurrentPath = true; //¸ü¸ÄyÔÚµ±Ç°Â·¾¶ÖÐ
                Q.push_back(yNode); // ½«y½áµã¼ÓÖÁµ±Ç°Â·¾¶½áµãQÖÐ
                idQ.insert(yId); // ½«y½áµã±àºÅ¼ÓÖÁµ±Ç°Â·¾¶½áµãidQÖÐ
                D[xId].insert(yId); // ±ê¼Çx½áµãµÄÁÚ¾ÓyÒÑ±»·ÃÎÊ¹ý

                // ¸üÐÂÓ°ÏìÖµInf{V-v}(u)
                for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
					if (it->second->onCurrentPath == false) {
                        it->second->inf += ppNext;
                    }
                }
                break; //Ìø³öµ±Ç°forÑ­»·£¬½«yµ±×÷×îºóÒ»¸ö½áµãÑ°ÕÒÂ·¾¶
            } 
		}

        //µ±Ç°Â·¾¶×îºóÒ»¸ö½áµãÒÑËÑË÷Íê±Ï£¬Í¨¹ýBacktrack¼ÌÐøÑ°ÕÒÐÂÂ·¾¶¡£
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

        //½«×îºóÒ»¸ö½áµã´Óµ±Ç°Â·¾¶É¾³ý
		map<int, Node*>::iterator it = UNodes.find(lastId);
        if (it != UNodes.end()) {
			it->second->onCurrentPath = false;
        }
        idQ.erase(lastId);
        D.erase(lastId);
        Q.pop_back();
    }


    // ¸üÐÂUµÄÓ°ÏìÖµ£¨U = {x|x²»ÔÚ½áµã¸²¸Ç¼¯ºÏÖÐÇÒÊÇuµÄÈë½áµã£©
    for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
        int xId = it->first;
        Node *xNode = it->second;
        nvcInf[xId] += bMap[xId] * xNode->inf; //¸ù¾ÝTheorem2
    } 

    return inf;
}
/*¼ÆËãÓ°ÏìÖµInf(S)£¬Í¬Ê±¼ÆËãInf{V-x}(S)£¬xÊôÓÚU*/
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

/*¼ÆËãInf{V-S+u}(u)£¬Í¬Ê±¸üÐÂInf{V-S+u-x}(u)£¬xÊôÓÚU*/
double SimPath::backtrackNormal(int id, set<int> S,set<int> U){
	Node *uNode = nodes->find(id)->second; // ¸ù¾ÝidÕÒµ½½áµãu
	int uId = uNode->id;
	int uOutNum = uNode->outNodes->size();
    uNode->inf = 1;
    uNode->pp = 1;
	uNode->onCurrentPath = true;

	//¹¹ÔìU
	set<int> UIds;
	map<int, Node *> UNodes;
	for (set<int>::iterator it = U.begin(); it != U.end(); it++) {
		Node *v = nodes->find(*it)->second;
		int vId = v->id;
		v->inf = 1; //¼ÇÂ¼uÔÚV-vÉÏµÄÓ°ÏìÖµ
		v->pp = 1;
		v->onCurrentPath = false;
		UNodes.insert(make_pair(vId, v));
		UIds.insert(vId);
	} 

    double inf = 1; // cov(u), initially 1 (counting itself first)
    double pp = 1; // path prob. 
	vector<Node *> Q; //¼ÇÂ¼µ±Ç°Â·¾¶½áµã
    set<int> idQ; //¼ÇÂ¼µ±Ç°Â·¾¶½áµãid
    map<int, set<int> > D; //¼ÇÂ¼½áµã±»·ÃÎÊ¹ýµÄÁÚ¾Ó

    Q.push_back(uNode);
    idQ.insert(id);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // FORWARD starts here!!!
            Node *xNode = Q.back(); // »ñÈ¡µ±Ç°Â·¾¶×îºóÒ»¸ö½áµã
            int xId = xNode->id; 
    
            if (xId == lastNodeId) {
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // get the current path prob. till this node

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //µ±Ç°½áµãxÃ»ÓÐÁÚ¾Ó
            }

            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y¾ÍÊÇuÊ±
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y ÔÚÂ·¾¶ÉÏÒÑ¾­´æÔÚ£¬Ìø¹ý·ÀÖ¹ÐÎ³É»·
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y ÒÑ¾­¼ÆËã¹ý£¬Ìø¹ý·ÀÖ¹ÖØ¸´¼ÆËã
                    continue;
				} else if (S.find(yId)!=S.end()){  //y²»ÔÚW£¨V-S£©ÖÐ
					continue;
				}

                double ppNext = pp * it->first; // pp = pp * b(x,y)
                // ÅÐ¶ÏÊÇ·ñ´ïµ½ãÐÖµ
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); // y is explored
                    continue;
                } 

				//Èô¸ßÓÚãÐÖµÎ´±»¼ôµôÔòËµÃ÷ÕÒµ½Ò»ÌõÐÂÂ·¾¶
                inf += ppNext;  //¸üÐÂuµÄ×ÜÓ°ÏìÖµ
                yNode->pp = ppNext; // ¸üÐÂuÖÁ½áµãyµÄÂ·¾¶¸ÅÂÊÖµ
				yNode->onCurrentPath = true; //´ËÊ±yÔÚµ±Ç°Â·¾¶ÖÐ -->¿É¸Ä

                Q.push_back(yNode); // ½«y½áµã¼ÓÖÁµ±Ç°Â·¾¶½áµãQÖÐ
                idQ.insert(yId); // ½«y½áµãid¼ÓÖÁµ±Ç°Â·¾¶½áµãidQÖÐ
                D[xId].insert(yId); // ¼ÇÂ¼x½áµãµÄÁÚ¾ÓyÒÑ±»·ÃÎÊ¹ý

                // ¸üÐÂÓ°ÏìÖµ£¨V-v£©
                for (map<int, Node *>::iterator it = UNodes.begin(); it != UNodes.end(); ++it) {
					if (it->second->onCurrentPath == false) {
                        it->second->inf += ppNext;
                    }
                }
                break;
            } 

		}

         //µ±Ç°Â·¾¶×îºóÒ»¸ö½áµãÒÑËÑË÷Íê±Ï£¬Í¨¹ýBacktrack¼ÌÐøÑ°ÕÒÐÂÂ·¾¶¡£
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

        //½«×îºóÒ»¸ö½áµã´Óµ±Ç°Â·¾¶É¾³ý
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

/*¸ù¾Ý±àºÅid¼ÆËãInf{V-S}(x)µÄÖµ*/
double SimPath::backtrackSimple(int id, set<int> S){
	Node *uNode = nodes->find(id)->second; // ¸ù¾ÝidÕÒµ½½áµãu
	int uId = uNode->id;
    uNode->inf = 1;
    uNode->pp = 1;

    double inf = 1; // cov(u), initially 1 (counting itself first)
    double pp = 1; // path prob. 
	vector<Node *> Q; //¼ÇÂ¼µ±Ç°Â·¾¶½áµã
    set<int> idQ; //¼ÇÂ¼µ±Ç°Â·¾¶½áµãid
    map<int, set<int> > D; //¼ÇÂ¼½áµã±»·ÃÎÊ¹ýµÄÁÚ¾Ó

    Q.push_back(uNode);
    idQ.insert(id);
   
    while (Q.empty() == false) {
        int lastNodeId = 0;
        while (true) {
            // FORWARD starts here!!!
            Node *xNode = Q.back(); // »ñÈ¡µ±Ç°Â·¾¶×îºóÒ»¸ö½áµã
            int xId = xNode->id; 
    
            if (xId == lastNodeId) {
                break;
            }
            lastNodeId = xId;

            pp = xNode->pp; // get the current path prob. till this node

			multimap<double,Node *> *outFromX = xNode->outNodes;
            if (outFromX == NULL || outFromX->empty() == true) {
                continue; //µ±Ç°½áµãxÃ»ÓÐÁÚ¾Ó
            }

            for (multimap<double,Node*>::iterator it = outFromX->begin(); it != outFromX->end(); it++) {
                Node *yNode = it->second; 
                int yId = yNode->id;
           
                if (yId == uId) { // y¾ÍÊÇuÊ±
                    D[xId].insert(yId);
                    continue;
                } else if (idQ.find(yId) != idQ.end()) { // y ÔÚÂ·¾¶ÉÏÒÑ¾­´æÔÚ£¬Ìø¹ý·ÀÖ¹ÐÎ³É»·
                    continue;
                } else if (D[xId].find(yId) != D[xId].end()) { // y ÒÑ¾­¼ÆËã¹ý£¬Ìø¹ý·ÀÖ¹ÖØ¸´¼ÆËã
                    continue;
				} else if (S.find(yId)!=S.end()){  //y²»ÔÚW£¨V-S£©ÖÐ
					D[xId].insert(yId);
					continue;
				}

                double ppNext = pp * it->first; // pp = pp * b(x,y)
                // ÅÐ¶ÏÊÇ·ñ´ïµ½ãÐÖµ
                if ( ppNext < pruneVal) {
                    inf += ppNext;
                    D[xId].insert(yId); // y is explored
                    continue;
                } 

				//Èô¸ßÓÚãÐÖµÎ´±»¼ôµôÔòËµÃ÷ÕÒµ½Ò»ÌõÐÂÂ·¾¶
                inf += ppNext;  //¸üÐÂuµÄ×ÜÓ°ÏìÖµ
                yNode->pp = ppNext; // ¸üÐÂuÖÁ½áµãyµÄÂ·¾¶¸ÅÂÊÖµ
				yNode->onCurrentPath = true; //´ËÊ±yÔÚµ±Ç°Â·¾¶ÖÐ -->¿É¸Ä

                Q.push_back(yNode); // ½«y½áµã¼ÓÖÁµ±Ç°Â·¾¶½áµãQÖÐ
                idQ.insert(yId); // ½«y½áµãid¼ÓÖÁµ±Ç°Â·¾¶½áµãidQÖÐ
                D[xId].insert(yId); // ¼ÇÂ¼x½áµãµÄÁÚ¾ÓyÒÑ±»·ÃÎÊ¹ý

                break;
            } 

		}

        //½øÐÐBacktrace£¬¼´Ñ°ÕÒÏÂÒ»¸ö¿ÉÌ½Ë÷Â·¾¶µÄµã
        Node *lastNode = Q.back(); 
		int lastId = lastNode->id;
        idQ.erase(lastId);
        D.erase(lastId);
        Q.pop_back();
    }

    return inf;
}
