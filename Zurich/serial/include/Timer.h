#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>


typedef struct Timer {
    long beg, end;
    char *fname;
}Timer;

//defining to append, in case it doesnt exist
void timer_setconf(char *namefile, Timer *timer ){
    if(access(namefile, 0) != 0){
        FILE *file_time = fopen(namefile, "w");
        fprintf(file_time,"function,time_us\n");
        fclose(file_time);
    }
    timer->fname = namefile;
}



void timer_start(Timer* timer)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    timer->beg = tv.tv_sec*1e6 + tv.tv_usec*1e-0;
}

void timer_stop(Timer *timer, char label[])
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    timer->end = tv.tv_sec*1e6 + tv.tv_usec*1e-0;
    long res = (timer->end - timer->beg);
    
    FILE *file_time = fopen(timer->fname, "a");
    fprintf(file_time,"%s, %ld\n", label, res);
    fclose(file_time);
}

#endif /*_TIMER_H_*/