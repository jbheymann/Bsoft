/**
@file	timer.h 
@brief	Utilities for timing functions 
@author Bernard Heymann 
@date	Created: 20010316
@date	Modified: 20230321
**/
 
//#include <time.h>
//#include <sys/time.h>


// Function prototypes 
tm*			get_localtime();
void		show_localtime(time_t time_sec);
double		getwalltime();
double		getcputime();
double		timer_start();
double		timer_report(double lasttime);

