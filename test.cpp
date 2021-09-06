// ConsoleApplication3.cpp : �������̨Ӧ�ó������ڵ㡣
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
	cout << "��������Ϊ��" << a.GetRowCount() << endl;
	cout << "��������Ϊ��" << a.GetColumnCount() << endl;
	cout << "һ�ж��е���Ϊ��" << a.GetElement(1, 2) << endl;
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
