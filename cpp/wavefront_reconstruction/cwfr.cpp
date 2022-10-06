// cwfr.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "cwfr.h"


CWFR::~CWFR()
{
}

MatrixXXd CWFR::operator()(const MatrixXXd& Sx, const MatrixXXd& Sy, const MatrixXXd& X, const MatrixXXd& Y, WFR_METHOD method)
{
	switch (method)
	{
	case WFR_METHOD::HFLI:
		return hfli(Sx, Sy, X, Y);
		break;
	case WFR_METHOD::HFLIQ:
		break;
	default:
		return MatrixXXd();
		break;
	}
}

MatrixXXd CWFR::hfli(const MatrixXXd& Sx, const MatrixXXd& Sy, const MatrixXXd& X, const MatrixXXd& Y)
{
	// 0. build the least-squares system
	// 0.0 Pre-allocate or reserve the required space



	// 1. solve the least-squares system

	// 2. build and return the result
	return MatrixXXd();
}
