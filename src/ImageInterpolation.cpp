#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>
#include "NxNDCT.h"



void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	uchar* y_old = (uchar*)malloc(xSize * ySize);
	char* u_old = (char*)malloc(xSize * ySize / 4);
	char* v_old = (char*)malloc(xSize * ySize / 4);

	uchar* y_new = (uchar*)malloc(newXSize * newYSize);
	char* u_new = (char*)malloc(newXSize * newYSize / 4);
	char* v_new = (char*)malloc(newXSize * newYSize / 4);

	double Fh = (double)newXSize / xSize;
	double Fv = (double)newYSize / ySize;

	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	int i_new = 0;
	int j_new = 0;

	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {
			i_new = (i - 1) / Fv;
			j_new = (j - 1) / Fh;

			if (i_new < newYSize - 1) {
				i_new += 1;
			}

			if (j_new < newXSize - 1) {
				j_new += 1;
			}

			y_new[i*newXSize + j] = y_old[i_new * xSize + j_new];
		}
	}

	for (int i = 0; i < newYSize / 2; i++) {
		for (int j = 0; j < newXSize / 2; j++) {
			i_new = (i - 1) / Fv;
			j_new = (j - 1) / Fh;

			if (i_new < newYSize - 1) {
				i_new += 1;
			}

			if (j_new < newXSize - 1) {
				j_new += 1;
			}

			u_new[i*newXSize / 2 + j] = u_old[i_new * xSize / 2 + j_new];
			v_new[i*newXSize / 2 + j] = v_old[i_new * xSize / 2 + j_new];
		}
	}

	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	free(y_old);
	free(u_old);
	free(v_old);
	free(y_new);
	free(u_new);
	free(v_new);

}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	uchar* y_old = (uchar*)malloc(xSize * ySize);
	char* u_old = (char*)malloc(xSize * ySize / 4);
	char* v_old = (char*)malloc(xSize * ySize / 4);

	uchar* y_new = (uchar*)malloc(newXSize * newYSize);
	char* u_new = (char*)malloc(newXSize * newYSize / 4);
	char* v_new = (char*)malloc(newXSize * newYSize / 4);

	double sh = (double)newXSize / xSize;
	double sv = (double)newYSize / ySize;

	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	int i_new = 0;
	int j_new = 0;
	int i_new_2 = 0;
	int j_new_2 = 0;

	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {
			double a = i / sv - floor(i / sv);
			double b = j / sh - floor(j / sh);

			i_new = i / sv;
			j_new = j / sh;

			i_new_2 = i_new;
			j_new_2 = j_new;

			if (i_new_2 < ySize - 1) {
				i_new_2 += 1;
			}

			if (j_new_2 < xSize - 1) {
				j_new_2 += 1;
			}

			y_new[i*newXSize + j] = (1 - a)*(1 - b)*y_old[i_new * xSize + j_new]
				+ (1 - a)*b*y_old[i_new * xSize + j_new_2]
				+ a*(1 - b)*y_old[i_new_2* xSize + j_new]
				+ a*b*y_old[i_new_2*xSize + j_new_2];
		}
	}

	for (int i = 0; i < newYSize / 2; i++) {
		for (int j = 0; j < newXSize / 2; j++) {
			double a = i / sv - floor(i / sv);
			double b = j / sh - floor(j / sh);

			i_new = i / sv;
			j_new = j / sh;

			i_new_2 = i_new;
			j_new_2 = j_new;

			if (i_new_2 < ySize / 2 - 1) {
				i_new_2 += 1;
			}

			if (j_new_2 < xSize / 2 - 1) {
				j_new_2 += 1;
			}

			u_new[i*newXSize / 2 + j] = (1 - a)*(1 - b)*u_old[i_new * xSize / 2 + j_new]
				+ (1 - a)*b*u_old[i_new * xSize / 2 + j_new_2]
				+ a*(1 - b)*u_old[i_new_2* xSize / 2 + j_new]
				+ a*b*u_old[i_new_2*xSize / 2 + j_new_2];

			v_new[i*newXSize / 2 + j] = (1 - a)*(1 - b)*v_old[i_new * xSize / 2 + j_new]
				+ (1 - a)*b*v_old[i_new * xSize / 2 + j_new_2]
				+ a*(1 - b)*v_old[i_new_2* xSize / 2 + j_new]
				+ a*b*v_old[i_new_2*xSize / 2 + j_new_2];
		}
	}


	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	free(y_old);
	free(u_old);
	free(v_old);
	free(y_new);
	free(u_new);
	free(v_new);

}

double wHelperFcn(double d) {
	if (abs(d) < 1) {
		return 1.5 * pow(abs(d), 3) - 2.5 * d * d + 1;
	}
	else if((abs(d) >= 1) && (abs(d) < 2)){
		return -0.5 * pow(abs(d), 3) + 2.5 * d * d - 4 * abs(d) + 2;
	}
	else
	{
		return 0;
	}
}

uchar cubicInterpolate(uchar points[4], double d) {
	double w[4];
	int ret = 0;

	w[0] = wHelperFcn(d + 1);
	w[1] = wHelperFcn(d);
	w[2] = wHelperFcn(1-d);
	w[3] = wHelperFcn(2-d);

	for (int i = 0; i < 4; i++) {
		ret += points[i] * w[i];
	}

	if (ret > 255) {
		ret = 255;
	}
	else  if (ret < 0) {
		ret = 0;
	}

	return (uchar)ret;
}

char cubicInterpolate(char points[4], double d) {
	double w[4];
	int ret = 0;

	w[0] = wHelperFcn(d + 1);
	w[1] = wHelperFcn(d);
	w[2] = wHelperFcn(1 - d);
	w[3] = wHelperFcn(2 - d);

	for (int i = 0; i < 4; i++) {
		ret += points[i] * w[i];
	}

	if (ret > 127) {
		ret = 127;
	}
	else  if (ret < -128) {
		ret = -128;
	}

	return (char)ret;
}


void bicubicInterpolate(uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */

	uchar* y_old = new uchar[xSize * ySize];
	char* u_old = new char[xSize * ySize / 4];
	char* v_old = new char[xSize * ySize / 4];

	int extendedXsizeY;
	int extendedYsizeY;
	int extendedXsizeUV;
	int extendedYsizeUV;

	extendBorders(y_old, xSize, ySize, 2, &y_old, &extendedXsizeY, &extendedYsizeY);
	extendBorders(u_old, xSize / 2, ySize / 2, 2, &u_old, &extendedXsizeUV, &extendedYsizeUV);
	extendBorders(v_old, xSize / 2, ySize / 2, 2, &v_old, &extendedXsizeUV, &extendedYsizeUV);
	
	uchar* y_new = new uchar[newXSize * newYSize];
	char* u_new = new char[newXSize * newYSize / 4];
	char* v_new = new char[newXSize * newYSize / 4];

	double sh = (double)newXSize / xSize;
	double sv = (double)newYSize / ySize;

	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	for (int i = 0; i < newYSize; i++) {
		for (int j = 0; j < newXSize; j++) {
			double dv = i / sv - floor(i / sv);
			double dh = j / sh - floor(j / sh);

			int i_new = i / sv;
			int j_new = j / sh;

			uchar cubic_x[4];
			uchar cubic_y[4];

			int yInd = 0;

			for (int k = i_new - 1; k < i_new + 3; k++, yInd++) {
				cubic_x[0] = y_old[k * xSize + j_new - 1];
				cubic_x[1] = y_old[k * xSize + j_new];
				cubic_x[2] = y_old[k * xSize + j_new + 1];
				cubic_x[3] = y_old[k * xSize + j_new + 2];
				
				cubic_y[yInd] = cubicInterpolate(cubic_x, dh);
			}

			y_new[i * newXSize + j] = cubicInterpolate(cubic_y, dv);
		}
	}

	

	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	delete[] y_old;
	delete[] u_old;
	delete[] v_old;
	delete[] y_new;
	delete[] u_new;
	delete[] v_new;
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
}