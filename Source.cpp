#include <iostream>
#include <cmath>
#include <fstream>
#include <thread>
#include <string>


using namespace std;

double C = 5;
double a = 0.05;

const int iterationsX = 10;
const int iterationsT = 1000;
double dx = 1. / (iterationsX - 1);
double dt = 1. / (iterationsT - 1);


double approximatedSolution(double wMinus1, double w, double wPlus1)
{

	double add1 = pow(w, (double)-2 / 3) * ((wMinus1 - 2 * w + wPlus1) / pow(dx, 2));
	double add2 = (wPlus1 - wMinus1) / (3 * dx * pow(w, (double)5 / 3));

	return w + dt * (add1 - add2);
}

double correctSolution(double x, double t)
{
	double mult1 = pow((C - 4 * a * t), (double)3 / 2);
	double mult2 = mult1 - pow(x, 2);

	return mult1 * pow(mult2, (double)-3 / 2);
}


void writeToFile(string name, double matrix[iterationsT][iterationsX])
{
	ofstream file(name);
	double x = 0.0;
	double t = 0.0;

	file << '{';
	for (int i = 0; i < iterationsT; ++i)
	{
		x = 0.0;
		for (int j = 0; j < iterationsX; ++j)
		{
			file << '{' << x << ", " << t << ", " << matrix[i][j] << '}' << ',';
			x += dx;
		}
		t += dt;
	}
	file << '}';

	file.close();
}



int main()
{
	setlocale(0, "rus");


	double x = 0.0, t = 0.0;
	double correctMatrix[iterationsT][iterationsX];


	for (int i = 0; i < iterationsT; ++i)
	{
		x = 0.0;
		for (int j = 0; j < iterationsX; ++j)
		{
			correctMatrix[i][j] = correctSolution(x, t);
			x += dx;
		}
		t += dt;
	}

	x = 0.0;

	double approximatedMatrix[iterationsT][iterationsX];

	for (int j = 0; j < iterationsX; ++j)
	{
		approximatedMatrix[0][j] = correctSolution(x, 0);
		x += dx;
	}

	t = 0.0;
	for (int i = 0; i < iterationsT; ++i)
	{
		approximatedMatrix[i][0] = correctSolution(0, t);
		approximatedMatrix[i][iterationsX - 1] = correctSolution(1, t);
		t += dt;
	}

	double errorValue = 0.0;
	int errorIPos;
	int errorJPos;
	for (int i = 1; i < iterationsT; ++i)
	{
#pragma omp parallel for
		for (int j = 1; j < iterationsX - 1; ++j)
		{
			approximatedMatrix[i][j] = approximatedSolution(approximatedMatrix[i - 1][j - 1], approximatedMatrix[i - 1][j],
				approximatedMatrix[i - 1][j + 1]);
			if (fabs(approximatedMatrix[i][j] - correctMatrix[i][j]) > errorValue)
			{
				errorValue = fabs(approximatedMatrix[i][j] - correctMatrix[i][j]);
				errorIPos = i;
				errorJPos = j;
			}
		}
	}
	cout << "Значение абсолютной погрешности: " << errorValue << endl;
	cout << "Место максимальной разницы в матрице: " << errorIPos << " по i, " << errorJPos << " по j." << endl;
	cout << "Относительная погрешность: " <<
		fabs((approximatedMatrix[errorIPos][errorJPos] - correctMatrix[errorIPos][errorJPos])
			/ correctMatrix[errorIPos][errorJPos]) * 100 << "%" << endl;

	writeToFile("super_approximated_result.txt", approximatedMatrix);

	writeToFile("super_correct_result.txt", correctMatrix);

	system("pause");
	return 0;
}