#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <chrono>


using namespace std::chrono;

typedef struct Timer {
    time_point<high_resolution_clock>  beg, end;
    const char *fname;
}Timer;


//defining to append, in case it doesnt exist
inline void timer_setconf(const char namefile[], Timer *timer,const char header[]){

        FILE *file_time = fopen(namefile, "w");
        if (file_time == nullptr) file_error();

        const char default_header[] = "function, time_us\n";
        const char *custom_header = (header) ? default_header : ""; 
        
        fprintf(file_time, "function,%s,time_us\n", custom_header);
        fclose(file_time);
    
    timer->fname = namefile;
}



inline void timer_start(Timer* timer)
{
    // struct timeval tv;
    // gettimeofday(&tv, NULL);
    // timer->beg = tv.tv_sec*1e6 + tv.tv_usec*1e-0;

    //using lib chrono
    timer->beg = high_resolution_clock::now();
}

inline void timer_stop(Timer *timer,const char label[])
{
    // struct timeval tv;
    // gettimeofday(&tv, NULL);
    // timer->end = tv.tv_sec*1e6 + tv.tv_usec*1e-0;
    timer->end = high_resolution_clock::now();
    
    
    double res = duration<double>(timer->end - timer->beg).count();
    
    FILE *file_time = fopen(timer->fname, "a");
    fprintf(file_time,"%s, %.8lf\n", label, res*1e9);
    fclose(file_time);
}

#endif /*_TIMER_H_*/