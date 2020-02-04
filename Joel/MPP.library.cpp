#include "MPP.h"

extern "C" {
	Individual *maker(){
		return new MPP();
	}
}
