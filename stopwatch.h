#include <time.h>

class stopwatch
{
public:
    stopwatch() : start(clock()){} //start counting time
    ~stopwatch();
private:
    clock_t start;
};
