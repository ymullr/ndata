#include "narray.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
		size_t  i=0;
		NdShape rank_shp={3, {45,65,89}};
		size_t length_acc=1;

		for(i=0;i<rank_shp.rank;i++){
				length_acc*=rank_shp.shape[i];
		}

		long *multiarray=new long[length_acc];

		for(i=0;i<length_acc;i++){
				multiarray[i]=i;
		}

		long ptr3d=nd_i(rank_shp,20, 15, 45);
		cout<<"Running PointerMathTest....\n"
		<<"multiarray["<<ptr3d<<"]=?"<<multiarray[ptr3d]<<"\n";
		if( multiarray[ptr3d]==ptr3d && ptr3d!=0){
				cout<<"TEST OK\n";

		}else{
				cerr<<"ERREUR\n";
		}
}
