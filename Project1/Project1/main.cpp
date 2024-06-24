#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
#include<stack>
#define xx 240
#define yy 300

void bresenham(Mat im,Point p0,Point p1, int color)//p0->p1
{
	int x, y, dx, dy, e;
	dx = p1.x - p0.x; dy = p1.y - p0.y;
	int ux = dx > 0 ? 1 : -1;//x伸展方向
	int uy = dy > 0 ? 1 : -1;//y伸展方向
	dx = abs(dx); dy = abs(dy);
	x = p0.x; y = p0.y;
	if (dx > dy)
	{//以x为增量方向计算
		e = -dx;
		do{
			im.at<uchar>(y, x) = color;
			x += ux;
			e = e + 2 * dy;
			if (e >= 0)
			{
				y += uy;
				e = e - 2 * dx;
			}
		} while (x != p1.x);
	}
	else
	{//以y为增量方向计算
		e = -dy;
		do {
			im.at<uchar>(y, x) = color;
			y += uy;
			e = e + 2 * dx;
			if (e >= 0)
			{
				x += ux;
				e = e - 2 * dy;
			}
		} while (y != p1.y);
	}
}
void fill(Mat im, Point s, int boundarycolor, int newcolor)
{
	stack <Point> p;
	p.push(s);
	int color;
	Point point;
	while (!p.empty())
	{
		point = p.top(); p.pop();
		color = im.at<uchar>(point);
		if (color != newcolor && color != boundarycolor)
		{
			im.at<uchar>(point) = newcolor;
			p.push(Point(point.x, point.y + 1));
			p.push(Point(point.x, point.y - 1));
			p.push(Point(point.x - 1, point.y));
			p.push(Point(point.x + 1, point.y));
		}
	}
}
void polygon(Mat im,Point *p,int n,int color)
{
	p[n] = p[0];
	for (int i = 0; i < n; i++)
		bresenham(im, p[i], p[i + 1], color);
}
void AApolygon(Mat im, Point *p,int n,int boundarycolor,int innercolor)
{
	int x, y, dx, dy;
	float k, e;
	int q = 0;
	p[n] = p[0];
	for (int i = 0; i < n; i++)
	{
		dx = p[i+1].x - p[i].x; dy = p[i+1].y - p[i].y;
		int ux = dx > 0 ? 1 : -1;//x伸展方向
		int uy = dy > 0 ? 1 : -1;//y伸展方向
		dx = abs(dx); dy = abs(dy);
		e = 0;
		x = p[i].x; y = p[i].y;
		int co = 255 - innercolor;
		if (dx > dy)
		{//以x为增量方向计算
			if (im.at<uchar>((p[i].y+p[i+1].y)/2 + 1, (p[i].x+p[i+1].x)/2) == innercolor) q = 1;
			else if (im.at<uchar>((p[i].y + p[i + 1].y) / 2 - 1, (p[i].x + p[i + 1].x) / 2) == innercolor)q = -1;
			k = (float)dy / (float)dx;//
			do {
				if (q == uy)
					im.at<uchar>(y, x) = innercolor + co * e;
				else im.at<uchar>(y + uy, x) = 255 - co * e;
				e = e + k;
				if (e >= 1)
				{
					y += uy;
					e = e - 1;
				}
				x += ux;
			} while (x != p[i+1].x);
		}
		else
		{//以y为增量方向计算
			if (im.at<uchar>((p[i].y + p[i + 1].y) / 2, (p[i].x + p[i + 1].x) / 2+1) == innercolor) q = 1;
			else if (im.at<uchar>((p[i].y + p[i + 1].y) / 2, (p[i].x + p[i + 1].x) / 2 - 1) == innercolor)q = -1;
			else continue;
			k = ((float)dx / (float)dy);
			do {
				if (q == ux)
					im.at<uchar>(y, x) = innercolor + co * e;
				else im.at<uchar>(y, x+ux) = 255 - co * e;
				e = e + k;
				if (e >= 1)
				{
					x += ux;
					e = e - 1;
				}
				y += uy;
			} while (y != p[i+1].y);
		}
	}
}
void circlePoints(Mat im, int x, int y, int color)
{
	im.at<uchar>(x + yy, y + xx) = im.at<uchar>(y + yy, x + xx) = color;
	im.at<uchar>(-x + yy, y + xx) = im.at<uchar>(y + yy, -x + xx) = color;
	im.at<uchar>(x + yy, -y + xx) = im.at<uchar>(-y + yy, x + xx) = color;
	im.at<uchar>(-x + yy, -y + xx) = im.at<uchar>(-y + yy, -x + xx) = color;
}
void midPointCircle(Mat im, int r, uchar color)
{
	int x, y;
	float d;
	x = 0; y = r; d = 1.25 - r;
	circlePoints(im, x, y, color);
	while (x <= y)
	{
		if (d < 0)
			d += 2 * x + 3;
		else
		{
			d += 2 * (x - y) + 5;
			y--;
		}
		x++;
		circlePoints(im, x, y, color);
	}
}
int main(int argc, char* argv[])
{
	Mat grayim(600, 800, CV_8UC1, Scalar::all(255));
	Point *p=new Point[20];//小鱼
	p[0] = Point(200, 300); p[1] = Point(400, 50); p[2] = Point(400, 295); p[3] = Point(500, 200);
	p[4] = Point(500, 400); p[5] = Point(400, 305); p[6] = Point(400, 550); 

	polygon(grayim,p ,7, 100);
	fill(grayim, Point(300, 300), 100, 100);
	AApolygon(grayim, p, 7, 100, 100);

	midPointCircle(grayim, 5, 20);
	fill(grayim, Point(242, 300), 20, 20);


	
//显示保存结果
	imshow("grayim", grayim);
	//imwrite("grayim.png", grayim);
	waitKey(0);
	return 0;
}