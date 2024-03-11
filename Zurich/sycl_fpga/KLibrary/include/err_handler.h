#ifndef _ERR_HANDLER_H_
#define _ERR_HANDLER_H_

#include <stdio.h>
#include <stdlib.h>


// Error types - for file reader
void alloc_error();
void file_error();
void bound_error();
void length_error();


#endif //_ERR_HANDLER_H_