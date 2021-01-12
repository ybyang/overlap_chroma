#ifndef READWRITE_H
#define READWRITE_H

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "linearop.h"
#include "meas/inline/io/named_objmap.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermstates/fermstates.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include <stdio.h>
#include "util/ferm/diractodr.h"
#include "util/ferm/transf.h"


namespace Chroma
{

struct io_para
{
	int io_num;
	int Vvol;
	multi1d<int> latsize;
	int this_node;
	int node_nums;
	int io_ratio;
	int xinz;
	int memsize;
	int little_size;
	int io_idx;
	int vvol;
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
	LatticeFermionF fieldF;
	LatticeFermionD fieldD;
	io_para para;

	io_vec(bool _single, int io_num);

	void read(FILE* filehand);

	void write(FILE* filehand);

	void readD(FILE* filehand);

	void readF(FILE* filehand);

	void writeD(FILE* filehand);

	void writeF(FILE* filehand);
};

}
#endif
