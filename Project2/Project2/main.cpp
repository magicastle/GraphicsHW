#include<iostream>
#include <opencv2/opencv.hpp>
#include<vector>
#include <math.h> 
#include<random>
#include <fstream> 
#include<time.h>
#define M_PI 3.141592653589793238462643 
using namespace cv;
using namespace std;

//随机数生成函数
double drand()
{
	return (double)(rand() / (double)RAND_MAX);
}


struct Vec3
{
	double x, y, z;
	Vec3() { x = 0; y = 0; z = 0; }
	Vec3(double _x, double _y = 0, double _z = 0.) : x(_x), y(_y), z(_z) {}
	//Vec3(Vec3& vec) { x = vec.x; y = vec.y; z = vec.z; }
	Vec3 operator*(double b) const { return Vec3(x * b, y * b, z * b); }
	Vec3 operator+ (const Vec3& rhs) const { return Vec3(x + rhs.x, y + rhs.y, z + rhs.z); }
	Vec3 operator- (const Vec3& rhs) const { return Vec3(x - rhs.x, y - rhs.y, z - rhs.z); }
	Vec3 mult(const Vec3 &rhs) const { return Vec3(x * rhs.x, y * rhs.y, z * rhs.z); }
	double dot(const Vec3 &rhs) const { return x * rhs.x + y * rhs.y + z * rhs.z; }
	Vec3 cross(const Vec3 &rhs) { return Vec3(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x); }
	Vec3& norm() { return *this = *this * (1 / sqrt(x*x + y * y + z * z)); }//标准化
	Vec3 clip() const { return Vec3(x > 255 ? 255 : x < 0 ? 0 : x, y>255 ? 255 : y < 0 ? 0 : y, z>255 ? 255 : z < 0 ? 0 : z); }
};

struct Ray
{
	Vec3 o, d;
	Ray(Vec3 _o, Vec3 _d) { o = _o; d = _d; }
};


enum Refl_t { DIFF, SPEC, REFR };//反射属性参数

struct Sphere
{
	double radius;
	Vec3 position, emission, color;
	Refl_t refl;
	bool hastexture = 0;
	int rows, cols;

	Mat texture;
	Sphere(double _r, Vec3 _p, Vec3 _e, Vec3 _c, Refl_t _refl, bool hastexture = 0) :radius(_r), position(_p), emission(_e), color(_c), refl(_refl), hastexture(hastexture)
	{
		if (hastexture)
		{
			texture = imread("m.jpg", 1);
			rows = texture.rows;
			cols = texture.cols;
		}

	}
	Vec3 getcolor(double x, double y, double z)const
	{
		if (!hastexture)return color;
		x = abs(x - position.x); y = abs(y - position.y); z = abs(z - position.z);
		Vec3b color;
		if (x == 0 & y == 0)
			color = texture.at<Vec3b>(0, 0);
		else
			color = texture.at<Vec3b>(rows *asin(z / radius) / (2 * M_PI), cols *atan(y / x) / (2 * M_PI));
		return Vec3(color[0], color[1], color[2])*(0.00392);
	}
	double intersect(const Ray &ray) const
	{
		Vec3 op = position - ray.o;
		double t, eps = 1e-4;
		double b = op.dot(ray.d), det = b * b + radius * radius - op.dot(op);
		if (det < 0) return 0;//无交点
		else det = sqrt(det);
		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);//内部 外部
	}
};



//放置球
Sphere spheres[] = {

	Sphere(1e5,Vec3(1e5 + 3, 40.8, 81.6), Vec3(30,30,30),         Vec3(.12, .24, .17), Refl_t::DIFF), //left
	Sphere(1e5,Vec3(-1e5 + 97,40.8,81.6), Vec3(30,30,30),         Vec3(.39,.12,.86),   Refl_t::DIFF), //Rght 

	Sphere(1e5,Vec3(50,40.8, 1e5),        Vec3(30,30,30),         Vec3(.75,.75,.75),   Refl_t::DIFF), //Back 
	Sphere(1e5,Vec3(50,40.8,-1e5 + 170),  Vec3(30,30,30),         Vec3(),              Refl_t::DIFF), //Frnt 
	Sphere(1e5,Vec3(50, 1e5, 81.6),       Vec3(30,30,30),         Vec3(.75,.75,.75),   Refl_t::DIFF), //Botm 
	Sphere(1e5,Vec3(50,-1e5 + 90,81.6), Vec3(30,30,30),         Vec3(.75,.75,.75),   Refl_t::DIFF), //Top 

	Sphere(16.5,Vec3(27,16.5,47),         Vec3(30,30,30),         Vec3(1,1,1)*.999,    Refl_t::SPEC), //Mirr 
	Sphere(16.5,Vec3(75,50,90),         Vec3(30,30,30),         Vec3(1,1,1)*.999,    Refl_t::REFR,1), //Glas 

	Sphere(600,Vec3(50,690 - 0.27,81.6), Vec3(300,300,300), Vec3(),              Refl_t::DIFF)  //Lite 光源

};


inline bool intersect(const Ray &ray, double &t, int &id) //光线和所有球求交点
{
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i >= 0; i--)
	{
		if ((d = spheres[i].intersect(ray)) != 0 && d < t)
		{
			t = d;//相交物体的距离 留下最近的
			id = i;//相交物体的ID
		}
	}
	return t < inf;
}

Vec3 radiance(const Ray &r, int depth)
{
	double t; // distance to intersection 
	int id = 0; // id of intersected object 
	if (!intersect(r, t, id))return Vec3(); // if miss, return black 
	const Sphere &obj = spheres[id]; // the hit object 
	Vec3 x = r.o + r.d*t,//交点 
		n = (x - obj.position).norm(),
		nl = n.dot(r.d) < 0 ? n : n * -1, //球面法向量
		f = obj.getcolor(x.x, x.y, x.z);
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max
	if (++depth > 6)//深度  
		//if (drand() < p)///////////////
		//	f = f * (1 / p);
		//else
		return obj.emission;

	if (obj.refl == DIFF) // 漫反射 随机生成一个方向进行漫反射
	{
		double r1 = 2 * M_PI*drand(), r2 = drand(), r2s = sqrt(r2);//随机数
		Vec3 w = nl, u = ((fabs(w.x) > .1 ? Vec3(0, 1) : Vec3(1)).cross(w)).norm(), v = w.cross(u);  //w，v，u为正交基
		Vec3 d = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm();//方向
		return obj.emission + f.mult(radiance(Ray(x, d), depth));
	}
	else
	{
		Ray reflRay(x, r.d - n * 2 * n.dot(r.d));
		if (obj.refl == SPEC) // 镜面反射
			return obj.emission + f.mult(radiance(reflRay, depth));

		else//折射加反射
		{
			// 由平行四边形的方法求得反射光的direction ？？？
			bool into = n.dot(nl) > 0;                // 判断光线是否进入球体
			double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
			if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0)    // Total internal reflection 
				return obj.emission + f.mult(radiance(reflRay, depth));
			Vec3 tdir = (r.d*nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
			double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
			double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);

			return obj.emission + f.mult(depth > 2 ? (drand() < P ?   // Russian roulette 轮盘赌
				radiance(reflRay, depth)*RP : radiance(Ray(x, tdir), depth)*TP) :
				radiance(reflRay, depth)*Re + radiance(Ray(x, tdir), depth)*Tr);
		}
	}
}


Vec3* bezierCurve(Vec3* points, int n, double precision)     //de Casteljau 
{
	Vec3 *p = new Vec3[n];
	Vec3* curve = new Vec3[int(1. / precision) + 7];
	int k = 0;
	for (double t = 0; t <= 1; t += precision)
	{
		for (int i = 0; i < n; i++)
			p[i] = Vec3(points[i].x, points[i].y);
		for (int i = 1; i < n; ++i)//n-1轮迭代
			for (int j = 0; j < n - i; ++j)//第i次迭代中，计算第i层的控制点(n-i)个
				p[j] = Vec3(p[j].x*(1 - t) + p[j + 1].x*t, p[j].y*(1 - t) + p[j + 1].y*t);
		curve[k++] = Vec3(p[0].x, p[0].y);
	}

	return curve;
}

int main(int argc, char *argv[])//camera at(50,52,295.6) 往Z轴负方向看
{
	int w = 800, h = 600, samps = 1; // # samples 100 200 400   1024 768  
	Ray cam(Vec3(50, 50, 296), Vec3(0, -0.02, -1).norm()); // cam pos, dir Ray cam(Vec3(50, 52, 295.6), Vec3(0, -0.042612, -1).norm()); 
	Vec3 cx = Vec3(w / h * 0.5135),//0.5135
		cy = (cx.cross(cam.d)).norm() * 0.5135,
		r;

	Vec3 *c = new Vec3[w*h]; //the image

//遍历每个像素点，随机采样的方式求得要射出的光线的方向
	for (int y = 0; y < h; y++)  // Loop over image rows
	{
		printf("\rRendering (%d spp) %5.2f%%--------y=%d", samps * 4, 100.*y / (h - 1), y);
		for (int x = 0; x < w; x++)  // Loop cols 
		{
			for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++)  // 2x2 subpixel rows 
			{
				for (int sx = 0; sx < 2; sx++, r = Vec3())  // 2x2 subpixel cols 
				{
					for (int s = 0; s < samps; s++)//default samps=1
					{
						double r1 = 2 * drand(),//(0,2)
							dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1), //(-1,1)
							r2 = 2 * drand(),
							dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

						Vec3 d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;//出射光线方向

						r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0)* (1. / samps);//发出一束光线
					}
					c[i] = c[i] + r.clip()*0.25;
				}
			}
		}
	}

	Mat im(600, 800, CV_8UC3, Scalar(255, 255, 255));
	Vec3b color;
	for (int i = 0; i < w*h; i++)
	{
		//cout << "i/w,i%w : " << i / w << i % w <<"color:"<< c[i].x<< " "<<c[i].y<<" "<< c[i].z<< endl;
		color = Vec3b(c[i].x, c[i].y, c[i].z);
		im.at<Vec3b>(i / w, i%w) = color;
	}
	imshow("im", im);

	imwrite("p15.png", im);/////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!
	waitKey(0);







	/*Mat grayim(600, 800, CV_8UC1, Scalar::all(255));
	Vec3* controlPoints = new Vec3[6];
	int n = 4;
	controlPoints[0] = Vec3(200, 310); controlPoints[1] = Vec3(215, 350); controlPoints[2] = Vec3(230, 340); controlPoints[3] = Vec3(250, 310);
	for (int i = 0; i < 4; i++)
	{
		cout << "i=" << i << "y="<<controlPoints[i].y<<"   x="<< controlPoints[i].x <<endl;
		grayim.at<uchar>(controlPoints[i].y,controlPoints[i].x) = 0;
	}

	Vec3* curve = bezierCurve(controlPoints, n,0.001);
	for (int i = 0; i < 100; i++)
	{
		cout << "i=" << i << "y=" << curve[i].y << "   x=" << curve[i].x << endl;
		grayim.at<uchar>(curve[i].y,curve[i].x) =100;
	}
	imshow("grayim", grayim);
	imwrite("grayim.png", grayim);
	waitKey(0);*/
	return 0;
}