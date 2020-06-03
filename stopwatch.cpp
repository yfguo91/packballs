#include <iostream>
#include "stopwatch.h"

using namespace std;

stopwatch::~stopwatch()
{
    clock_t total = clock()-start; //get elapsed time
    cout << "total of ticks for this activity: " << total << endl;
    cout << "in seconds: " << double(total/CLK_TCK) << endl;
	cout << "Clock ticks per second = " << CLK_TCK << endl;
}
