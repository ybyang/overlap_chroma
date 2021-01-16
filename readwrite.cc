#include "readwrite.h"
#define max 1024*1024*1024


namespace Chroma
{

void adfseek(FILE *stream, long long size, int fromwhere){
	long long times = size / max;
	long long res = size % max;
	for(long long i=0; i<times; i++)
	fseek(stream, max, SEEK_CUR);
	fseek(stream, res, SEEK_CUR);
}


void adfread(char *buff, int small, long long size, FILE* stream){
	int times = size / max;
	int res = size % max;
	for(int i=0; i<times; i++)
		fread(buff+i*max,small,max,stream);
	fread(buff+times*max,small,res,stream);
}


void adfwrite(char *buff, int small, long long size, FILE* stream){
	int times = size / max;
	int res = size % max;
	for(int i=0; i<times; i++)
		fwrite(buff+i*max,small,max,stream);
	fwrite(buff+times*max,small,res,stream);
}


io_vec::io_vec(bool _single, int io_num){
	single = _single;
	para.readtime1.reset();
	para.readtime2.reset();
	para.readtime3.reset();
	para.readtime4.reset();
	para.readtime5.reset();
	para.readtime6.reset();
	para.io_num = io_num;
	para.Vvol = Layout::vol();
	para.latsize = Layout::lattSize();
	para.this_node = Layout::nodeNumber();
	para.node_nums = Layout::numNodes();
	para.io_ratio = para.node_nums/io_num;
	para.xinz = Layout::subgridLattSize()[0];
	if(single){
		para.memsize = sizeof(infermF);
		para.little_size = sizeof(REAL32);
	}
	else{
		para.memsize = sizeof(infermD);
		para.little_size = sizeof(REAL64);
	}
	para.io_idx = para.this_node/para.io_ratio;
	para.vvol = para.Vvol/io_num;
	para.io_flag = false;
	if(para.this_node%para.io_ratio == 0) para.io_flag = true;
	else para.io_idx = -1;
}


void io_vec::read(FILE* filehand){
	QMP_barrier();
	if(single) readF(filehand);
	else readD(filehand);
	QMP_barrier();
}


void io_vec::write(FILE* filehand){
	QMP_barrier();
	if(single) writeF(filehand);
	else writeD(filehand);
	QMP_barrier();
}


void io_vec::readD(FILE* filehand){

	typedef infermD inferm;
	typedef REAL64 little;
	char* buff = NULL;
	little* tmp = NULL;
	char* buf = NULL;
	inferm* fout=NULL;
	char* fin = NULL;

	fin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag){
		fout = (inferm*) malloc(para.vvol*para.memsize);
		para.readtime1.start();
		buff = (char*) malloc(para.vvol*para.memsize);
		for(int j=0;j<2*Ns*Nc;j++){
			adfseek(filehand, para.io_idx*para.vvol*para.little_size, SEEK_CUR);
			adfread(buff+j*para.vvol*para.little_size, 1, para.vvol*para.little_size, filehand);
			adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.little_size, SEEK_CUR);
		}
		para.readtime1.stop();
		para.readtime2.start();
		QDPUtil::byte_swap(buff, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime3.start();
		buf = buff;
		for(int reim=0;reim<2;reim++)
		for(int spin=0;spin<Ns;spin++)
		for(int color=0;color<Nc;color++){
			tmp = (little*)buf;
			if(reim==0)
				#pragma omp parallel for shared(fout,tmp)
				for(int site=0;site<para.vvol;site++)
					fout[site].elem(spin).elem(color).real() = tmp[site];
			else
				#pragma omp parallel for shared(fout,tmp)
				for(int site=0;site<para.vvol;site++)
					fout[site].elem(spin).elem(color).imag() = tmp[site];
			tmp = NULL;
			buf+=para.little_size*para.vvol;
		}
		para.readtime3.stop();
		buf = NULL;
		free(buff);
		buff = NULL;
	}
	for(int IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(int site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		int node = Layout::nodeNumber(crtesn(IO_idx*para.vvol+site, para.latsize));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(fin,fout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(int j=0;j<para.xinz;j++){
			long long linear = Layout::linearSiteIndex(crtesn(IO_idx*para.vvol+site+j, para.latsize));
			memcpy((char *)&(fieldD.elem(0))+linear*para.memsize, fin+j*para.memsize, para.memsize);
		}
		para.readtime6.stop();
	}
	if(para.io_flag) free(fout);
	fout = NULL;
	free(fin);
	fin = NULL;

}


void io_vec::readF(FILE* filehand){

	typedef infermF inferm;
	typedef REAL32 little;
	inferm* fout=NULL;
	char* fin = NULL;

	fin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag){
		fout = (inferm*) malloc(para.vvol*para.memsize);
		para.readtime1.start();
		adfseek(filehand, para.io_idx*para.vvol*para.memsize, SEEK_CUR);
		adfread((char*)fout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		para.readtime2.start();
		QDPUtil::byte_swap((char*)fout, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
	}
	for(int IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(int site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		int node = Layout::nodeNumber(crtesn(IO_idx*para.vvol+site, para.latsize));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(fin,fout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(int j=0;j<para.xinz;j++){
			long long linear = Layout::linearSiteIndex(crtesn(IO_idx*para.vvol+site+j, para.latsize));
			memcpy((char *)&(fieldF.elem(0))+linear*para.memsize, fin+j*para.memsize, para.memsize);
		}
		para.readtime6.stop();
	}
	if(para.io_flag) free(fout);
	fout = NULL;
	free(fin);
	fin = NULL;

}


void io_vec::writeD(FILE* filehand){

	typedef infermD inferm;
	typedef REAL64 little;
	char* buff = NULL;
	little* tmp = NULL;
	char* buf = NULL;
	inferm* fout=NULL;
	char* fin = NULL;
	fin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag)
		fout = (inferm*) malloc(para.vvol*para.memsize);

	for(int IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(int site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		int node = Layout::nodeNumber(crtesn(IO_idx*para.vvol+site, para.latsize));
		para.readtime4.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(int j=0;j<para.xinz;j++){
			long long linear = Layout::linearSiteIndex(crtesn(IO_idx*para.vvol+site+j, para.latsize));
			memcpy(fin+j*para.memsize, (char *)&(fieldD.elem(0))+linear*para.memsize, para.memsize);
		}
		para.readtime6.stop();
		para.readtime5.start();
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, node, IO_idx*para.io_ratio, para.xinz*para.memsize);
		if(para.io_idx == IO_idx) memcpy(fout+site,fin,para.memsize*para.xinz);
		para.readtime5.stop();
	}
	free(fin);
	fin = NULL;
	if(para.io_flag){
		buff = (char*) malloc(para.vvol*para.memsize);
		buf = buff;
		para.readtime3.start();
		for(int reim=0;reim<2;reim++)
		for(int spin=0;spin<Ns;spin++)
		for(int color=0;color<Nc;color++){
			tmp = (little*)buf;
			if(reim==0)
			#pragma omp parallel for shared(fout,tmp)
				for(int site=0;site<para.vvol;site++)
					tmp[site] = fout[site].elem(spin).elem(color).real();
			else
				#pragma omp parallel for shared(fout,tmp)
				for(int site=0;site<para.vvol;site++)
					tmp[site] = fout[site].elem(spin).elem(color).imag();
			tmp = NULL;
			buf+=para.little_size*para.vvol;
		}
		para.readtime3.stop();
		buf = NULL;
		free(fout);
		fout = NULL;
		para.readtime2.start();
		QDPUtil::byte_swap(buff, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime1.start();
		for(int j=0;j<2*Ns*Nc;j++){
			adfseek(filehand, para.io_idx*para.vvol*para.little_size, SEEK_CUR);
			adfwrite(buff+j*para.vvol*para.little_size, 1, para.vvol*para.little_size, filehand);
			adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.little_size, SEEK_CUR);
		}
		para.readtime1.stop();
		free(buff);
		buff = NULL;
	}

}


void io_vec::writeF(FILE* filehand){

	typedef infermF inferm;
	typedef REAL32 little;
	inferm* fout=NULL;
	char* fin = NULL;
	fin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag)
		fout = (inferm*) malloc(para.vvol*para.memsize);

	for(int IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(int site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		int node = Layout::nodeNumber(crtesn(IO_idx*para.vvol+site, para.latsize));
		para.readtime4.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(int j=0;j<para.xinz;j++){
			long long linear = Layout::linearSiteIndex(crtesn(IO_idx*para.vvol+site+j, para.latsize));
			memcpy(fin+j*para.memsize, (char *)&(fieldF.elem(0))+linear*para.memsize, para.memsize);
		}
		para.readtime6.stop();
		para.readtime5.start();
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, node, IO_idx*para.io_ratio, para.xinz*para.memsize);
		if(para.io_idx == IO_idx) memcpy(fout+site,fin,para.memsize*para.xinz);
		para.readtime5.stop();
	}
	free(fin);
	fin = NULL;
	if(para.io_flag){
		para.readtime2.start();
		QDPUtil::byte_swap((char*)fout, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime1.start();
		adfseek(filehand, para.io_idx*para.vvol*para.memsize, SEEK_CUR);
		adfwrite((char*)fout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		free(fout);
		fout = NULL;
	}

}

}//end namespace
