#ifndef _OMP_TIMER_H_
#define _OMP_TIMER_H_

#include <cstdio>
#include <omp.h>
#include <chrono>

using namespace std::chrono;

class ParallelTimer {
public:
    ParallelTimer(const char* namefile, const char* header = nullptr) : fname(namefile) {
        setconf(header);
    }

    void setconf(const char* header = nullptr) {
        FILE* file_time = fopen(fname, "a+");

        if (file_time) {
            fclose(file_time);
            file_time = fopen(fname, "a+");
            
            if (!file_time) {
                perror("Error opening file");
                return;
            }
        } else {
            file_time = fopen(fname, "w");
            
            if (!file_time) {
                perror("Error opening file");
                return;
            }

            fprintf(file_time, "function,%s,time_us\n", (header ? header : ""));
        }

        fclose(file_time);
    }

    void start() {
        beg = omp_get_wtime();
    }

    void stop(const char label[]) {
        end = omp_get_wtime();
        double res = (end - beg) * 1e9;

        FILE* file_time = fopen(fname, "a");
        fprintf(file_time, "%s, %.8lf\n", label, res);
        fclose(file_time);
    }

private:
    const char* fname;
    double beg, end;
};

#endif /*_OMP_TIMER_H_*/