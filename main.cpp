#define _USE_MATH_DEFINES
#include <vector>
#include <random>
#include <array>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cmath>
using namespace std;

double x_min = -100., x_max = 100.;
const int n = 60;
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> dist(x_min, x_max);
uniform_real_distribution<double> vel(-fabs(x_max - x_min), fabs(x_max - x_min));


class Particle
{
private:
	vector<double> x;
	vector<double> p;
	vector<double> v;
public:
	Particle(vector<double> _x);

	void print()
	{
		for (int i = 0; i < n; ++i)
			cout << x[i] << ' ';
	}

	vector<double> getp() { return p; }
	vector<double> getx();
	void setx(vector<double> _x);
	void setv(vector<double> _v) { v = _v; }
	void update_v(double w,
		double fp,
		double fg,
		vector<double> rp,
		vector<double> rg,
		vector<double> g);

	void move();
};



class Swarm
{
private:
	vector<Particle> swarm;
	const int num_swarm;
	const int dim;
	double w, fp, fg;	//	
	vector<double> g; //the best vector

public:
	Swarm(double, double, double, int, int );
	void next_step();
	vector<double> getBest() { return g; }

};


inline double f(const vector<double> &x)
{
	double s = 0.;
	for (int i = 0; i < n-1; ++i)
		s += pow(1 - x[i], 2) + 100 * pow(x[i + 1] - pow(x[i], 2), 2); //Розенброк

	//int A = 10;
	//s = A*n;
	//for (int i = 0; i < n; ++i)
	//	s += pow(x[i], 2) - A*cos(2 * M_PI*x[i]);

	// s = -0.0001 * pow(fabs(sin(x[0])*sin(x[1]) * exp(fabs(100 - sqrt(x[0] * x[0] + x[1] * x[1]) / M_PI)) + 1), 0.1);

	//for (int i = 0; i < n; ++i)
	//	s += x[i] * x[i];
	return s;
}


int main(int argc, char *argv[])
{
	//random_device rd;
	//mt19937 gen(rd());
	//uniform_real_distribution<double> dist(x_min, x_max);
	//uniform_real_distribution<double> vel(-fabs(x_max - x_min), fabs(x_max - x_min));
	//vector<Particle> swarm;
	//vector<double> g(n);
	double w = 0.72984,
		fp = 1.496172,
		fg = fp;
	Swarm sw(w, fp, fg, n, 100);

	for (int i = 0; i < 10000; ++i)
	{
		sw.next_step();
		//cout << endl;
		//auto best = sw.getBest();
		//cout << f(best) << endl;
	}
	
	//for_each(best.begin(), best.end(), [](double x) { cout << x << "  "; });
	

	//int count_par = 1000;
	//cout << "Count particles: " << count_par << endl
	//	<< "w: " << w << endl
	//	<< "fp: " << fp << endl

	//	<< "fg: " << fg << endl;
	//for (int i = 0; i < count_par; ++i)
	//{
	//	vector<double> tmp;
	//	for (int j = 0; j < n; ++j)
	//		tmp.push_back(dist(gen));

	//	if (i == 0) g = tmp;
	//	else if (f(g) > f(tmp))
	//		g = tmp;
	//	
	//	swarm.push_back(Particle(tmp));

	//	vector<double> v;
	//	for (int j = 0; j < n; ++j)
	//		v.push_back(vel(gen));
	//	swarm[i].setv(v);
	//}

	////	for(int i = 0; i < count_par; ++i)
	////	{
	////		swarm[i].print();
	////		cout << endl;
	////	}
	//
	//int k = 100;
	//for (int i = 0; i < k; ++i)
	//{
	//	for (int j = 0; j < count_par; ++j) //для всех частиц
	//	{
	//		vector<double> rp(n), rg(n);
	//		for (int l = 0; l < n; ++l)
	//		{
	//			rp[l] = generate_canonical<double, 10>(gen);
	//			rg[l] = generate_canonical<double, 10>(gen);
	//		}
	//		swarm[j].update_v(w, fp, fg, rp, rg, g);
	//		swarm[j].move();
	//		if (f(g) > f(swarm[j].getp()))
	//			g = swarm[j].getp();
	//	}
	//	//cout << "f(g): " << f(g) << endl;
	//}




	//cout << endl;
	//	for(int i = 0; i < n; ++i)
	//		cout << g[i] << ' ';
	//cout << "f(g): " << f(g);
	cout << endl;
	system("PAUSE");
	return 0;
}

Swarm::Swarm(double _w, double _fp, double _fg, int _dim, int _par): w(_w), fp(_fp), fg(_fg), dim(_dim), num_swarm(_par)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dist(x_min, x_max);
	uniform_real_distribution<double> vel(-fabs(x_max - x_min), fabs(x_max - x_min));

	for (int i = 0; i < num_swarm; ++i)
	{
		vector<double> tmp(dim);
		for (int j = 0; j < n; ++j)
			tmp[j] = dist(gen);

		if (i == 0) // ничего лучше не придумал
			g = tmp;
		else if (f(g) > f(tmp))
			g = tmp;

		swarm.push_back(Particle(tmp));

		vector<double> v(dim);
		for (int j = 0; j < n; ++j)
			v[j] = vel(gen);
		swarm[i].setv(v);
	}
}

void Swarm::next_step()
{
		for (int j = 0; j < num_swarm; ++j) //для всех частиц
		{
			vector<double> rp(dim), rg(dim);
			for (int l = 0; l < n; ++l)
			{
				rp[l] = generate_canonical<double, 10>(gen);
				rg[l] = generate_canonical<double, 10>(gen);
			}
			swarm[j].update_v(w, fp, fg, rp, rg, g);
			swarm[j].move();
			if (f(g) > f(swarm[j].getp()))
				g = swarm[j].getp();
		}
}

void Particle::setx(vector<double> _x)
{

}

void Particle::update_v(double w,
	double fp,
	double fg,
	vector<double> rp,
	vector<double> rg,
	vector<double> g)
{
	for (int i = 0; i < x.size(); ++i)
		v[i] = w*v[i] +
		fp * rp[i] * (p[i] - x[i]) +
		fg * rg[i] * (g[i] - x[i]);

}

void Particle::move()
{
	for (int i = 0; i < x.size(); ++i)
		x[i] += v[i];
	if (f(x) < f(p))
		p = x;
}

Particle::Particle(vector<double> _x)
{
	x = p = _x;
}