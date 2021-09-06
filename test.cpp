// ConsoleApplication3.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"


#include<iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include "p0_starter.h"

using namespace std;
using namespace scudb;

int _tmain(int argc, char * argv[])
{
	RowMatrix<int> a(2, 3), b(2, 3), c(3, 2), d(3, 3);
	RowMatrix<int> *A = &a, *B = &b, *C = &c, *D = &d;
	vector<int> source{ { 1, 2, 3, 4, 5, 6 } };
	a.FillFrom(source);
	cout << "矩阵行数为：" << a.GetRowCount() << endl;
	cout << "矩阵列数为：" << a.GetColumnCount() << endl;
	cout << "一行二列的数为：" << a.GetElement(1, 2) << endl;
	a.SetElement(0, 0, 8);
	RowMatrixOperations<int> opr;
	cout << "a+b" << endl << " ";
	opr.Add(A, B);
	cout << endl;
	cout << "a*c" << endl << " ";
	opr.Multiply(A, C);
	cout << endl;
	cout << "a * c + d" << endl << " ";
	opr.GEMM(A, C, D);
	cout << endl;
	return 0;
}
