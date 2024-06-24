#include<iostream>
#include<vector>
#include<set>
#include<queue>
#include<utility>
#include<map>
#include<math.h>
#include<stdio.h>
#include <cstdlib>
#define EPS 0.0001
#define INF 123456789
#define T 1.0
using namespace std;

string inputfile,outputfile;
double ratio;
int cntFace, cntDelFace;//面片数和应该删除的面片数
map < pair<int,int>, bool> isDeleted;

struct Matrix4
{
	double mat[4][4];
	Matrix4(){for(int i=0;i<4;i++)for(int j=0;j<4;j++)mat[i][j]=0.;}
};

struct Vertex
{
	int id;//顶点编号
	double x,y,z;//顶点坐标位置
	set<int> connectV;//邻接点
	Matrix4* Q;
	bool deleted;
};
struct Pair
{
	int v1,v2;
	Vertex v;
	double cost;
};
Vertex vertex[1000000];int tail=1;
//vector<Vertex> vertex;
struct cmp
{//默认是less函数返回true时，a的优先级低于b的优先级（a排在b的后面） 
	bool operator() ( Pair a, Pair b ){ return a.cost> b.cost; }
};
priority_queue <Pair,vector<Pair>,cmp> pairHeap; //这样就可以实现小根堆了 


int  readFile(const char* inputfile)
{cout<<"reading-------------------"<<endl;
	if (freopen(inputfile, "r", stdin) == NULL)
	{
		cout << "Can't open the input file." << endl;
		return 0;
	}
	int cntv = 0, cntf = 0;
	char s[256];
	while (scanf("%s", s) == 1)
	{
		switch (s[0])
		{
			case '#': fgets(s, sizeof s, stdin); cout<<"#####"<<endl;break;//注释部分
			case 'v':
			{
				cntv++;
				double x, y, z;
				scanf("%lf %lf %lf", &x, &y, &z);
				Vertex v;
				v.x=x;v.y=y;v.z=z;
				//v.id=vertex.size();
				v.id=tail;
				v.deleted=false;
				//vertex.push_back(v);
				vertex[tail++]=Vertex(v);
				break;
			}
			case 'f'://建立邻接关系
			{
				cntf++;
				int a, b, c;
				scanf("%d%d%d", &a, &b, &c);
				vertex[a].connectV.insert(b);
				vertex[a].connectV.insert(c);
				vertex[b].connectV.insert(a);
				vertex[b].connectV.insert(c);
				vertex[c].connectV.insert(a);
				vertex[c].connectV.insert(b);
				break;
			}
			default: cout<<"others"<<endl;break;
		}
	}
	cout<<"cntv="<<cntv<<endl;
	fclose(stdin);
	return cntf;
}


void calQ()
{cout<<"----------calQ"<<endl;
	for(int i=1;i<tail;i++)
	{
		Vertex* v=&vertex[i];
		v->Q=new Matrix4();
		for(set<int>::iterator it1 = v->connectV.begin();it1 != v->connectV.end();it1++)	
			for(set<int>::iterator it2 = v->connectV.begin();it2 != v->connectV.end();it2++)

				if((*it1) < (*it2) && (vertex[(*it1)].connectV.count(*it2) >0 )) //
				{
					Vertex* v1 = &(vertex[*it1]);
					Vertex* v2 = &(vertex[*it2]);
					double coefficient[4];//平面方程四个参数a,b,c,d

					coefficient[0]=(v1->y-v->y)*(v2->z-v->z)-(v1->z-v->z)*(v2->y-v->y);
					coefficient[1]=(v1->z-v->z)*(v2->x-v->x)-(v1->x-v->x)*(v2->z-v->z);
					coefficient[2]=(v1->x-v->x)*(v2->y-v->y)-(v1->y-v->y)*(v2->x-v->x);
					coefficient[3]=-coefficient[0]*(v->x)-coefficient[1]*(v->y)-coefficient[2]*(v->z);

					for(int i = 0;i < 4;i++)
						for(int j = 0;j < 4;j++)
							v->Q->mat[i][j] += coefficient[i] * coefficient[j];
				}
	}
	return;
}	
void calculateBest(Pair& pair)
{
	Vertex* u=&vertex[pair.v1];
	Vertex* v=&vertex[pair.v2];
	double Q[4][4];//两个误差矩阵之和
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
		{
			Q[i][j]=u->Q->mat[i][j]+v->Q->mat[i][j];
		}

	Vertex* mid= &pair.v;//收缩节点	
	mid->deleted=false;
	mid->connectV.clear();
	mid->Q=new Matrix4();
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			mid->Q->mat[i][j]=Q[i][j];
//计算收缩点位置
	mid->x=(u->x + v->x)/2;mid->y=(u->y + v->y)/2;mid->z=(u->z + v->z)/2; //中点位置，应对退化情况
	Q[3][0] = 0;Q[3][1] = 0;Q[3][2] = 0;Q[3][3] = 1;
	double Y[4];Y[0]=Y[1]=Y[2]=0;Y[3]=1;

	for (int i = 0; i < 4; i++) 
	{
		int j = 0;
		while (j < 4 && fabs(Q[i][j]) < EPS) j++;
		if (j == 4)continue;

		for (int k = 0; k < 4; k++) 
		{
			if (k != i) 
			{
				double rate = Q[k][j] / Q[i][j];

				for (int l = 0; l < 4; l++)
					Q[k][l] -= Q[i][l] * rate;
				Y[k] -= Y[i] * rate;
			}
		}
	}	
	double X[4];X[3] = 1;
	for (int i = 0; i < 3; i++) 
	{
		int j = 0;
		while (j < 4 && fabs(Q[i][j]) < EPS) j++;
		if (j == 4){X[3]=-1;break;}
		X[i] = Y[i] / Q[i][j];
	}
	
	if(X[3] > EPS)
		mid->x=X[0];mid->y=X[1];mid->z=X[2];


	// if (vGroup->getCommonVertexNum(e.v1, e.v2) != 2) //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// {
	// 	e.deltaV = 0;
	// 	return;
	// }

//计算收缩代价cost 二次型vTQv
	//double X[4];
	X[0]=mid->x;X[1]=mid->y;X[2]=mid->z;X[3]=1;
	double cost = 0;
	for (int i = 0; i < 4; i++) 
	{
		double c = 0;
		for (int j = 0; j < 4; j++)
			c += X[j] * Q[i][j];
		cost += c * X[i];
	}
	pair.cost= cost;
	return;
}
double calDistance(Vertex& a,Vertex& b)
{
	return pow((a.x-b.x),2)+pow((a.y-b.y),2)+pow((a.z-b.z),2);
}
void buildHeap()
{cout<<"---------------------buildHeap"<<endl;
	//下面将边加入到边堆中
	for (int i = 0; i < tail; i++) //遍历所有顶点
	{
		for (set<int>::iterator it = vertex[i].connectV.begin();it != vertex[i].connectV.end(); it++) //邻接顶点
		{
			if(i<*it)
			{
				Pair pair;
				pair.v1= i;pair.v2=*it;
				calculateBest(pair);
				isDeleted.insert(map<std::pair<int,int>, bool >::value_type(make_pair(pair.v1,pair.v2),false));
				pairHeap.push(pair);
			}
		}
	}
}


void simplify()
{cout<<"---------------------simplify"<<endl;
	for (int i = 0; i < cntDelFace; i += 2) //开始删边	
	{cout<<"deleteface "<<i<<endl;
		if(pairHeap.size() <= 0){cout<<"there is no pair in heap"<<endl;return;}
		while(isDeleted[make_pair(pairHeap.top().v1, pairHeap.top().v2)] == true)pairHeap.pop();
		Pair minpair = pairHeap.top();pairHeap.pop();

		Vertex* v1 = &(vertex[minpair.v1]);    
		Vertex* v2 = &(vertex[minpair.v2]);
		Vertex targetv = minpair.v;
		targetv.id=tail;
		vertex[tail++]=Vertex(targetv);//
		isDeleted[make_pair(minpair.v1,minpair.v2)]=true; //打上边已经删除的标记

		int num1,num2;
		for (set<int>::iterator it = v1->connectV.begin(); it != v1->connectV.end(); it++) 
		{
			if ((*it) != v2->id) 
			{
				if((*it)>v1->id){num1=v1->id;num2=(*it);}else{num1=(*it),num2=v1->id;} 
				isDeleted[make_pair(num1,num2)]=true;
				vertex[(*it)].connectV.erase(v1->id);
				vertex[(*it)].connectV.insert(targetv.id);
				vertex[targetv.id].connectV.insert(*it);
			}
		}
		for (set<int>::iterator it = v2->connectV.begin(); it != v2->connectV.end(); it++) 
		{
			if ((*it) != v1->id && vertex[*it].connectV.count(v1->id)==0) /////////////////////
			{
				if((*it)>v2->id){num1=v2->id;num2=(*it);}else{num1=(*it),num2=v2->id;} 
				isDeleted[make_pair(num1,num2)]=true;
				vertex[(*it)].connectV.erase(v2->id);
				vertex[(*it)].connectV.insert(targetv.id);
				vertex[targetv.id].connectV.insert(*it);
			}
		}
		for(set<int>::iterator it = vertex[targetv.id].connectV.begin();it != vertex[targetv.id].connectV.end();it++)
		{cout << "push into heap" << endl;
			Pair newpair;//
			newpair.v1 = *it; newpair.v2 = targetv.id;
			calculateBest(newpair);
			isDeleted.insert(map<pair<int, int>, bool>::value_type(make_pair(newpair.v1, newpair.v2), false));
			pairHeap.push(newpair);
		}
		v1->deleted= true; //标记结点已经被删除
		v2->deleted = true;
		delete v1->Q,v2->Q;
	}
}
void saveFile(const char* outputfile)
{  cout<<"saving-------------------"<<endl;
	freopen(outputfile,"w",stdout);
	int cnt = 1;
	int cntv=0,cntf=0;
	for(int i = 1;i < tail;i++)
	{//输出所有点

		if(vertex[i].deleted==true)continue;           //如果第i个点已经删掉了，就略去
		vertex[i].id = cnt++;
		printf("v %lf %lf %lf\n",vertex[i].x,vertex[i].y,vertex[i].z);

	}	
	for(int i = 1;i <tail;i++)
	{//输出所有面
		if(vertex[i].deleted==true) continue;                     //如果第i个点已经删掉了，就略去
		for(set<int>::iterator it1 = vertex[i].connectV.begin();it1 != vertex[i].connectV.end();it1++)
		{
			if(i >= (*it1))continue;
			for(set<int>::iterator it2 = vertex[i].connectV.begin();it2 != vertex[i].connectV.end();it2++)
			{
				if((*it1) < (*it2) && (vertex[*it1].connectV.count(*it2)>0 ))
				{
					printf("f %d %d %d\n",vertex[i].id,vertex[*it1].id,vertex[*it2].id);
					cntf++;
				}	
			}
		}
	}
	fclose(stdout);
}
int main(int argc,char** argv)
{
	if(argc == 4)
	{
		inputfile=string(argv[1]);
		outputfile=string(argv[2]);
		ratio=(atof(argv[3]));
	}
	else if (argc == 1)
	{
		inputfile = "dinosaur.2k.obj";
		outputfile = "dinosaur.2k2.obj";
		ratio =0.8;
	}
	else
	{
		cout << "Parameters error." << endl;
		cout << "Usage: MeshSimplification [inputPath] [outputPath] [ratio]" << endl;
		return 1;
	}

	// Vertex emptyv;
	// vertex.push_back(emptyv);
	int temp=readFile(inputfile.c_str());
	if(temp!=0)cntFace=temp;else {cout<<"open fail"<<endl;return 0;}
	cntDelFace=(1-ratio)*cntFace;
	cout<<"cntFace="<<cntFace<<"   cntDelFace="<<cntDelFace<<endl;
	calQ();
	buildHeap();//建堆（将每条边作为一个pair）
	simplify();//迭代化简
	saveFile(outputfile.c_str());
	cout<<"write finished-------"<<endl;
	return 0;
}