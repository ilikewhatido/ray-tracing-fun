#ifndef MATHHELPER_H
#define MATHHELPER_H

//PI
const float Pi       =  3.14159265358979323846f;

//PI over 2
const float PiOver2  = Pi / 2.0f;

//2PI
const float TwoPi    = Pi * 2.0f;

// error tolerance for check
const float EPSILON  = 0.005f;


#define	DEGTORAD(x)	( ((x) * Pi) / 180.0 )
#define	RADTODEG(x)	( ((x) * 180.0) / Pi )

#define	SQR(x)		( (x) * (x) )

// limits a value to low and high
#define LIMIT_RANGE(low, value, high)	{	if (value < low)	value = low;	else if(value > high)	value = high;	}

// set float to 0 if within tolerance
#define ZERO_CLAMP(x)	( (EPSILON > fabs(x))?0.0f:(x) )


__inline void clampf(float &v, float min, float max)
{
    if(v < min)
        v = min;
    else if(v > max)
        v = max;
}

__inline void clampd(double &v, double min, double max)
{
    if(v < min)
        v = min;
    else if(v > max)
        v = max;
}

__inline bool FLOAT_EQ(float x, float v)
{
	return ( ((v) - EPSILON) < (x) && (x) < ((v) + EPSILON) );
}

__inline bool FLOAT_EQ(float x, float v, float epsi)
{
	return ( ((v) - epsi) < (x) && (x) < ((v) + epsi) );
}

//Swap 2 floating point numbers
__inline void SWAP(float &x, float &y)
{	
	float temp;	
	temp = x;	
	x = y;	
	y = temp;	
}

__inline void SWAP(int &x, int &y)
{
	int temp; 
	temp = x; 
	x = y;
	y = temp;
}

//Round for Accuracy. The same method we used by AMin
// round a float to a specified degree of accuracy
__inline float ROUND(const float value, const int accuracy)
{
	double integer, fraction;

	fraction = modf(value, &integer);		// get fraction and int components

	return (float(integer + (float(int(fraction*pow(float(10), float(accuracy))))) / pow(float(10), float(accuracy)) ) );
}




#endif
