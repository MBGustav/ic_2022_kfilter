
#include "err_handler.h"

void file_error() {
    fprintf(stderr, "Invalid input file.\n");
    exit(1);
}

void alloc_error(){
	fprintf(stderr, "Invalid Allocation.\n");
    exit(2);
}

void bound_error(){
    fprintf(stderr, "Invalid Memo. Access.\n");
    exit(3);

}

void length_error(){
    fprintf(stderr, "Invalid Sizes.\n");
    exit(4);
}