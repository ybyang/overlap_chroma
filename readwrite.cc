#include "readwrite.h"
#define max (1024*1024*1024)


namespace Chroma
{

void adfseek(FILE *stream, size_t size, int fromwhere){
	size_t times = size / max;
	size_t res = size % max;
	for(size_t i=0; i<times; i++)
		fseek(stream, max, SEEK_CUR);
	fseek(stream, res, SEEK_CUR);
}


std::vector<size_t> io_vec::adcrtesn(size_t ipos){

	std::vector<size_t> coord(Nd);
	for(int i = 0; i < Nd; ++i){
		coord[i] = ipos % para.latsize[i];
		ipos = ipos / para.latsize[i];
	}
	return coord;
}


size_t io_vec::adnodeNumber(const std::vector<size_t>& coord){

	multi1d<int> tmp_coord(Nd);
	for(int i=0; i<Nd; ++i)
		tmp_coord[i] = coord[i] / para.subsize[i];
	return QMP_get_node_number_from(tmp_coord.slice());
}


size_t adlocal_site(const std::vector<size_t>& coord, const std::vector<size_t>& latt_size){

	size_t order = 0;
	for(int mmu=Nd-1; mmu >= 1; --mmu)
		order = latt_size[mmu-1]*(coord[mmu]+order);
	order += coord[0];
	return order;
}


size_t io_vec::adlinearSiteIndex(const std::vector<size_t>& coord){

	size_t subgrid_vol_cb = para.vol >> 1;
	std::vector<size_t> subgrid_cb_nrow = para.subsize;
	subgrid_cb_nrow[0] >>= 1;
	size_t cb = 0;
	for(int m=0; m < Nd; ++m)
		cb += coord[m];
	cb &= 1;

	std::vector<size_t> subgrid_cb_coord(Nd);
	subgrid_cb_coord[0] = (coord[0] >> 1) % subgrid_cb_nrow[0];
	for(int i=1; i < Nd; ++i)
		subgrid_cb_coord[i] = coord[i] % subgrid_cb_nrow[i];

	return adlocal_site(subgrid_cb_coord, subgrid_cb_nrow) + cb*subgrid_vol_cb;
}



io_vec::io_vec(bool _single, int io_num, bool _endian=true){
	single = _single;
	endian = _endian;
	para.readtime1.reset();
	para.readtime2.reset();
	para.readtime3.reset();
	para.readtime4.reset();
	para.readtime5.reset();
	para.readtime6.reset();
	para.io_num = io_num;
	para.Vvol = Layout::vol();
	para.vol = Layout::sitesOnNode();
	para.latsize.resize(Nd);
	para.latsize[0] = Layout::lattSize()[0];
	para.latsize[1] = Layout::lattSize()[1];
	para.latsize[2] = Layout::lattSize()[2];
	para.latsize[3] = Layout::lattSize()[3];
	para.subsize.resize(Nd);
	para.subsize[0] = Layout::subgridLattSize()[0];
	para.subsize[1] = Layout::subgridLattSize()[1];
	para.subsize[2] = Layout::subgridLattSize()[2];
	para.subsize[3] = Layout::subgridLattSize()[3];
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

	U = zero;
	RealD foo = RealD(1) / sqrt(RealD(2));
	ComplexD one = cmplx(foo, RealD(0));
	ComplexD mone = cmplx(-foo, RealD(0));
	pokeSpin(U, one, 0, 1);
	pokeSpin(U, mone, 0, 3);
	pokeSpin(U, mone, 1, 0);
	pokeSpin(U, one, 1, 2);
	pokeSpin(U, one, 2, 1);
	pokeSpin(U, one, 2, 3);
	pokeSpin(U, mone, 3, 0);
	pokeSpin(U, mone, 3, 2);
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
		for(size_t j=0;j<2*Ns*Nc;j++){
			adfseek(filehand, para.io_idx*para.vvol*para.little_size, SEEK_CUR);
			fread(buff+j*para.vvol*para.little_size, 1, para.vvol*para.little_size, filehand);
			adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.little_size, SEEK_CUR);
		}
		para.readtime1.stop();
		para.readtime2.start();
		if(endian) QDPUtil::byte_swap(buff, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime3.start();
		buf = buff;
		for(int reim=0;reim<2;reim++)
		for(int spin=0;spin<Ns;spin++)
		for(int color=0;color<Nc;color++){
			tmp = (little*)buf;
			if(reim==0)
				#pragma omp parallel for shared(fout,tmp)
				for(size_t site=0;site<para.vvol;site++)
					fout[site].elem(spin).elem(color).real() = tmp[site];
			else
				#pragma omp parallel for shared(fout,tmp)
				for(size_t site=0;site<para.vvol;site++)
					fout[site].elem(spin).elem(color).imag() = tmp[site];
			tmp = NULL;
			buf+=para.little_size*para.vvol;
		}
		para.readtime3.stop();
		buf = NULL;
		free(buff);
		buff = NULL;
	}
	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(fin,fout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
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
		fread((char*)fout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		para.readtime2.start();
		if(endian) QDPUtil::byte_swap((char*)fout, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
    }
	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(fin,fout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)fin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
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

	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
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
				for(size_t site=0;site<para.vvol;site++)
					tmp[site] = fout[site].elem(spin).elem(color).real();
			else
				#pragma omp parallel for shared(fout,tmp)
				for(size_t site=0;site<para.vvol;site++)
					tmp[site] = fout[site].elem(spin).elem(color).imag();
			tmp = NULL;
			buf+=para.little_size*para.vvol;
		}
		para.readtime3.stop();
		buf = NULL;
		free(fout);
		fout = NULL;
		para.readtime2.start();
		if(endian) QDPUtil::byte_swap(buff, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime1.start();
		for(size_t j=0;j<2*Ns*Nc;j++){
			adfseek(filehand, para.io_idx*para.vvol*para.little_size, SEEK_CUR);
			fwrite(buff+j*para.vvol*para.little_size, 1, para.vvol*para.little_size, filehand);
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

	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
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
		if(endian) QDPUtil::byte_swap((char*)fout, para.little_size, 2*Ns*Nc*para.vvol);
		para.readtime2.stop();
		para.readtime1.start();
		adfseek(filehand, para.io_idx*para.vvol*para.memsize, SEEK_CUR);
		fwrite((char*)fout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		free(fout);
		fout = NULL;
	}

}


void io_vec::readR(FILE* filehand){
	QMP_barrier();
	if(single) readRF(filehand);
	else readRD(filehand);
	QMP_barrier();
}


void io_vec::readRD(FILE* filehand){
	typedef REAL64 inreal;
	inreal* rout=NULL;
	char* rin=NULL;
	para.memsize = sizeof(inreal);
	para.little_size = sizeof(inreal);

	rin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag){
		rout = (inreal*) malloc(para.vvol*para.memsize);
		para.readtime1.start();
		adfseek(filehand, para.io_idx*para.vvol*para.memsize, SEEK_CUR);
		fread((char*)rout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		para.readtime2.start();
		if(endian) QDPUtil::byte_swap((char*)rout, para.little_size, para.vvol);
		para.readtime2.stop();
	}
	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(rin,rout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)rin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
			memcpy((char *)&(realD.elem(0))+linear*para.memsize, rin+j*para.memsize, para.memsize);
		}
		para.readtime6.stop();
	}
	if(para.io_flag) free(rout);
	rout = NULL;
	free(rin);
	rin = NULL;
}


void io_vec::readRF(FILE* filehand){
	typedef REAL32 inreal;
	inreal* rout=NULL;
	char* rin=NULL;
	para.memsize = sizeof(inreal);
	para.little_size = sizeof(inreal);

	rin = (char*) malloc(para.xinz*para.memsize);
	if(para.io_flag){
		rout = (inreal*) malloc(para.vvol*para.memsize);
		para.readtime1.start();
		adfseek(filehand, para.io_idx*para.vvol*para.memsize, SEEK_CUR);
		fread((char*)rout, 1, para.vvol*para.memsize, filehand);
		adfseek(filehand, (para.io_num-para.io_idx-1)*para.vvol*para.memsize, SEEK_CUR);
		para.readtime1.stop();
		para.readtime2.start();
		if(endian) QDPUtil::byte_swap((char*)rout, para.little_size, para.vvol);
		para.readtime2.stop();
	}
	for(size_t IO_idx=0;IO_idx<para.io_num;IO_idx++)
	for(size_t site=0;site<para.vvol;site+=para.xinz){
		para.readtime4.start();
		size_t node = adnodeNumber(adcrtesn(IO_idx*para.vvol+site));
		para.readtime4.stop();
		para.readtime5.start();
		if(para.io_idx == IO_idx) memcpy(rin,rout+site,para.memsize*para.xinz);
		if(node != IO_idx*para.io_ratio) QDPInternal::route((void *)rin, IO_idx*para.io_ratio, node, para.xinz*para.memsize);
		para.readtime5.stop();
		para.readtime6.start();
		if(para.this_node == node)
		for(size_t j=0;j<para.xinz;j++){
			size_t linear = adlinearSiteIndex(adcrtesn(IO_idx*para.vvol+site+j));
			memcpy((char *)&(realF.elem(0))+linear*para.memsize, rin+j*para.memsize, para.memsize);
		}
		para.readtime6.stop();
	}
	if(para.io_flag) free(rout);
	rout = NULL;
	free(rin);
	rin = NULL;
}


}//end namespace
