

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "CVector.h"
#include "CMatrix.h"


#define	DONT_INTERSECT    0
#define	DO_INTERSECT      1
#define COLLINEAR         2

#define SAME_SIGNS( a, b )	\
		(((int) ((unsigned int) a ^ (unsigned int) b)) >= 0 )


int lines_intersect(int x1,int y1,   /* First line segment */
				 int x2,int y2,

				 int x3,int y3,   /* Second line segment */
				 int x4,int y4,

				 int *x,
				 int *y         /* Output value:
										* point of intersection */
							 )
{
		int a1, a2, b1, b2, c1, c2; /* Coefficients of line eqns. */
		int r1, r2, r3, r4;         /* 'Sign' values */
		int denom, offset, num;     /* Intermediate values */

		/* Compute a1, b1, c1, where line joining points 1 and 2
		 * is "a1 x  +  b1 y  +  c1  =  0".
		 */

		a1 = y2 - y1;
		b1 = x1 - x2;
		c1 = x2 * y1 - x1 * y2;

		/* Compute r3 and r4.
		 */


		r3 = a1 * x3 + b1 * y3 + c1;
		r4 = a1 * x4 + b1 * y4 + c1;

		/* Check signs of r3 and r4.  If both point 3 and point 4 lie on
		 * same side of line 1, the line segments do not intersect.
		 */

		if ( r3 != 0 &&
				 r4 != 0 &&
				 SAME_SIGNS( r3, r4 ))
				return ( DONT_INTERSECT );

		/* Compute a2, b2, c2 */

		a2 = y4 - y3;
		b2 = x3 - x4;
		c2 = x4 * y3 - x3 * y4;

		/* Compute r1 and r2 */

		r1 = a2 * x1 + b2 * y1 + c2;
		r2 = a2 * x2 + b2 * y2 + c2;

		/* Check signs of r1 and r2.  If both point 1 and point 2 lie
		 * on same side of second line segment, the line segments do
		 * not intersect.
		 */

		if ( r1 != 0 &&
				 r2 != 0 &&
				 SAME_SIGNS( r1, r2 ))
				return ( DONT_INTERSECT );

		/* Line segments intersect: compute intersection point. 
		 */

		denom = a1 * b2 - a2 * b1;
		if ( denom == 0 )
				return ( COLLINEAR );
		offset = denom < 0 ? - denom / 2 : denom / 2;

		/* The denom/2 is to get rounding instead of truncating.  It
		 * is added or subtracted to the numerator, depending upon the
		 * sign of the numerator.
		 */

		num = b1 * c2 - b2 * c1;
		*x = ( num < 0 ? num - offset : num + offset ) / denom;

		num = a2 * c1 - a1 * c2;
		*y = ( num < 0 ? num - offset : num + offset ) / denom;

		return ( DO_INTERSECT );
		} /* lines_intersect */

/* A main program to test the function.
 */

bool Line::Intersect(const Line& src,Vector3d<float> &sect)
{
	int x,y;

	int ret = lines_intersect( int(mP1.x), int(mP1.y),
														 int(mP2.x), int(mP2.y),
														 int(src.mP1.x), int(src.mP1.y),
														 int(src.mP2.x), int(src.mP2.y),&x,&y );

	if ( ret == DO_INTERSECT )
	{
		sect.x = float(x);
		sect.y = float(y);
		sect.z = 0;
		return true;
	}
	return false;
}



void BoundingBox::TransformBoundAABB(const MyMatrix &transform,const BoundingBox &b)
{
	InitMinMax();
	BoundTest(transform,b.bmin.x,b.bmin.y,b.bmin.z);
	BoundTest(transform,b.bmax.x,b.bmin.y,b.bmin.z);
	BoundTest(transform,b.bmax.x,b.bmax.y,b.bmin.z);
	BoundTest(transform,b.bmin.x,b.bmax.y,b.bmin.z);
	BoundTest(transform,b.bmin.x,b.bmin.y,b.bmax.z);
	BoundTest(transform,b.bmax.x,b.bmin.y,b.bmax.z);
	BoundTest(transform,b.bmax.x,b.bmax.y,b.bmax.z);
	BoundTest(transform,b.bmin.x,b.bmax.y,b.bmax.z);
}


void BoundingBox::BoundTest(const MyMatrix &transform,float x,float y,float z)
{
	Vector3d<float> pos(x,y,z);
	Vector3d<float> t;
	transform.Transform(pos,t);
	MinMax(t);
};

