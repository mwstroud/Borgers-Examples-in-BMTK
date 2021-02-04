#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _k_reg();
extern void _leak_reg();
extern void _na_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," k.mod");
fprintf(stderr," leak.mod");
fprintf(stderr," na.mod");
fprintf(stderr, "\n");
    }
_k_reg();
_leak_reg();
_na_reg();
}
