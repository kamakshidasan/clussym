#include "Utils.hpp"

unsigned int Index(unsigned int x, unsigned int y, unsigned int z, int SIZEX, int SIZEY, int SIZEZ) 
{

	return x*SIZEZ*SIZEY + y*SIZEZ + z;
//	return z*SIZEX*SIZEY + y*SIZEX + x;
}

void DeIndex(unsigned int i, unsigned int & x, unsigned int & y, unsigned int & z, int SIZEX, int SIZEY, int SIZEZ)
{
	
//	z = ((i/SIZEX)/SIZEY)%SIZEZ;
//	y = (i/SIZEX)%SIZEY;
//	x = i%SIZEX;
	x = ((i/SIZEZ)/SIZEY)%SIZEX;
	y = (i/SIZEZ)%SIZEY;
	z = i%SIZEZ;
}


