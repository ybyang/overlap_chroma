#ifndef READWRITE_H
#define READWRITE_H

#include "actions/ferm/fermstates/fermstates.h"
#include <stdio.h>
#include <vector>
#include "util/ferm/transf.h"


namespace Chroma
{

struct io_para
{
	size_t io_num;
	size_t Vvol;
	std::vector<size_t> latsize;
	std::vector<size_t> subsize;
	size_t this_node;
	size_t node_nums;
	size_t io_ratio;
	size_t xinz;
	size_t memsize;
	size_t little_size;
	size_t io_idx;
	size_t vvol;
	size_t vol;
	bool io_flag;

	StopWatch readtime1;
	StopWatch readtime2;
	StopWatch readtime3;
	StopWatch readtime4;
	StopWatch readtime5;
	StopWatch readtime6;
};


struct io_vec
{
	typedef PSpinVector< PColorVector< RComplex<REAL32>, Nc>, Ns> infermF;
	typedef PSpinVector< PColorVector< RComplex<REAL64>, Nc>, Ns> infermD;

	bool single;
	bool endian;
	LatticeFermionF fieldF;
	LatticeFermionD fieldD;
	LatticeRealD realD;
	LatticeRealF realF;
	io_para para;
	SpinMatrixD U;

	io_vec(bool _single, int io_num, bool _endian=true);

	std::vector<size_t> adcrtesn(size_t ipos);

	size_t adnodeNumber(const std::vector<size_t>& coord);

	size_t adlinearSiteIndex(const std::vector<size_t>& coord);

	void read(FILE* filehand);

	void write(FILE* filehand);

	void readD(FILE* filehand);

	void readF(FILE* filehand);

	void writeD(FILE* filehand);

	void writeF(FILE* filehand);

	void readR(FILE* filehand);

	void readRD(FILE* filehand);

	void readRF(FILE* filehand);
};

}
#endif
