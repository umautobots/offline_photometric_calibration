/**
* This file is part of DSO.
* 
* Copyright 2016 Technical University of Munich and Intel.
* Developed by Jakob Engel <engelj at in dot tum dot de>,
* for more information see <http://vision.in.tum.de/dso>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* DSO is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DSO. If not, see <http://www.gnu.org/licenses/>.
*/


#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <iterator>
#include "Undistort.h"

Undistort::~Undistort()
{
	if(remapX != 0) delete[] remapX;
	if(remapY != 0) delete[] remapY;
}

Undistort* Undistort::getUndistorterForFile(std::string configFilename)
{
	printf("Reading Calibration from file %s",configFilename.c_str());

	std::ifstream f(configFilename.c_str());
	if (!f.good())
	{
		f.close();
		printf(" ... not found. Cannot operate without calibration, shutting down.\n");
		f.close();
		return 0;
	}

	printf(" ... found!\n");
	std::string l1;
	std::getline(f,l1);
	f.close();

	float ic[10];

	Undistort* u;

    if(std::sscanf(l1.c_str(), "KannalaBrandt %f %f %f %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4], &ic[5], &ic[6], &ic[7]) == 8)
    {
        u = new UndistortKB(configFilename.c_str());
        if(!u->isValid()) {delete u; return 0; }
    }


    else if(std::sscanf(l1.c_str(), "RadTan %f %f %f %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4], &ic[5], &ic[6], &ic[7]) == 8)
    {
        u = new UndistortRadTan(configFilename.c_str());
        if(!u->isValid()) {delete u; return 0; }
    }


    else if(std::sscanf(l1.c_str(), "EquiDistant %f %f %f %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4], &ic[5], &ic[6], &ic[7]) == 8)
    {
        u = new UndistortEquidistant(configFilename.c_str());
        if(!u->isValid()) {delete u; return 0; }
    }


    else if(std::sscanf(l1.c_str(), "FOV %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4]) == 5)
    {
        u = new UndistortFOV(configFilename.c_str());
        if(!u->isValid()) {delete u; return 0; }
    }


    else if(std::sscanf(l1.c_str(), "Pinhole %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4]) == 5)
    {
        u = new UndistortPinhole(configFilename.c_str());
        if(!u->isValid()) {delete u; return 0; }
    }

    else
    {
        printf("could not read calib file! exit.");
        exit(1);
    }

	return u;
}

template<typename T>
void Undistort::undistort(const T* input, float* output, int nPixIn, int nPixOut) const
{
	if(!valid) return;

	if(nPixIn != wOrg*hOrg)
	{
		printf("ERROR: undistort called with wrong input image dismesions (expected %d pixel, got %d pixel)\n",
				wOrg*hOrg, nPixIn);
		return;
	}
	if(nPixOut != w*h)
	{
		printf("ERROR: undistort called with wrong output image dismesions (expected %d pixel, got %d pixel)\n",
				w*h, nPixOut);
		return;
	}

	if (!passthrough)
	{
		for(int idx = 0; idx < w*h;idx++)
		{
			// get interp. values
			float xx = remapX[idx];
			float yy = remapY[idx];

			if(xx<0)
				output[idx] = 0;
			else
			{
				// get integer and rational parts
				int xxi = xx;
				int yyi = yy;
				xx -= xxi;
				yy -= yyi;
				float xxyy = xx*yy;

				// get array base pointer
				const T* src = input + xxi + yyi * wOrg;

				// interpolate (bilinear)
				output[idx] =  xxyy * src[1+wOrg]
									+ (yy-xxyy) * src[wOrg]
									+ (xx-xxyy) * src[1]
									+ (1-xx-yy+xxyy) * src[0];
			}
		}
	}
	else
	{
		for(int idx = 0; idx < w*h;idx++)
		{
			const T* src = input + idx;
			output[idx] = src[0];
		}
	}
}
template void Undistort::undistort<float>(const float* input, float* output, int nPixIn, int nPixOut) const;
template void Undistort::undistort<unsigned char>(const unsigned char* input, float* output, int nPixIn, int nPixOut) const;
template void Undistort::undistort<unsigned short>(const unsigned short* input, float* output, int nPixIn, int nPixOut) const;

void Undistort::makeOptimalK_crop()
{
	printf("finding CROP optimal new model!\n");
	K.setIdentity();

	// 1. stretch the center lines as far as possible, to get initial coarse quess.
	float* tgX = new float[100000];
	float* tgY = new float[100000];
	float minX = 0;
	float maxX = 0;
	float minY = 0;
	float maxY = 0;

	for(int x=0; x<100000;x++)
	{tgX[x] = (x-50000.0f) / 10000.0f; tgY[x] = 0;}
	distortCoordinates(tgX, tgY,tgX, tgY,100000);
	for(int x=0; x<100000;x++)
	{
		if(tgX[x] > 0 && tgX[x] < wOrg-1)
		{
			if(minX==0) minX = (x-50000.0f) / 10000.0f;
			maxX = (x-50000.0f) / 10000.0f;
		}
	}
	for(int y=0; y<100000;y++)
	{tgY[y] = (y-50000.0f) / 10000.0f; tgX[y] = 0;}
	distortCoordinates(tgX, tgY,tgX, tgY,100000);
	for(int y=0; y<100000;y++)
	{
		if(tgY[y] > 0 && tgY[y] < hOrg-1)
		{
			if(minY==0) minY = (y-50000.0f) / 10000.0f;
			maxY = (y-50000.0f) / 10000.0f;
		}
	}
	delete[] tgX;
	delete[] tgY;

	minX *= 1.01;
	maxX *= 1.01;
	minY *= 1.01;
	maxY *= 1.01;



	printf("initial range: x: %.4f - %.4f; y: %.4f - %.4f!\n", minX, maxX, minY, maxY);



	// 2. while there are invalid pixels at the border: shrink square at the side that has invalid pixels,
	// if several to choose from, shrink the wider dimension.
	bool oobLeft=true, oobRight=true, oobTop=true, oobBottom=true;
	int iteration=0;
	while(oobLeft || oobRight || oobTop || oobBottom)
	{
		oobLeft=oobRight=oobTop=oobBottom=false;
		for(int y=0;y<h;y++)
		{
			remapX[y*2] = minX;
			remapX[y*2+1] = maxX;
			remapY[y*2] = remapY[y*2+1] = minY + (maxY-minY) * (float)y / ((float)h-1.0f);
		}
		distortCoordinates(remapX, remapY,remapX, remapY,2*h);
		for(int y=0;y<h;y++)
		{
			if(!(remapX[2*y] > 0 && remapX[2*y] < wOrg-1))
				oobLeft = true;
			if(!(remapX[2*y+1] > 0 && remapX[2*y+1] < wOrg-1))
				oobRight = true;
		}



		for(int x=0;x<w;x++)
		{
			remapY[x*2] = minY;
			remapY[x*2+1] = maxY;
			remapX[x*2] = remapX[x*2+1] = minX + (maxX-minX) * (float)x / ((float)w-1.0f);
		}
		distortCoordinates(remapX, remapY,remapX, remapY,2*w);


		for(int x=0;x<w;x++)
		{
			if(!(remapY[2*x] > 0 && remapY[2*x] < hOrg-1))
				oobTop = true;
			if(!(remapY[2*x+1] > 0 && remapY[2*x+1] < hOrg-1))
				oobBottom = true;
		}


		if((oobLeft || oobRight) && (oobTop || oobBottom))
		{
			if((maxX-minX) > (maxY-minY))
				oobBottom = oobTop = false;	// only shrink left/right
			else
				oobLeft = oobRight = false; // only shrink top/bottom
		}

		if(oobLeft) minX *= 0.995;
		if(oobRight) maxX *= 0.995;
		if(oobTop) minY *= 0.995;
		if(oobBottom) maxY *= 0.995;

		iteration++;


		printf("iteration %05d: range: x: %.4f - %.4f; y: %.4f - %.4f!\n", iteration,  minX, maxX, minY, maxY);
		if(iteration > 500)
		{
			printf("FAILED TO COMPUTE GOOD CAMERA MATRIX - SOMETHING IS SERIOUSLY WRONG. ABORTING \n");
			exit(1);
		}
	}

	K(0,0) = ((float)w-1.0f)/(maxX-minX);
	K(1,1) = ((float)h-1.0f)/(maxY-minY);
	K(0,2) = -minX*K(0,0);
	K(1,2) = -minY*K(1,1);

}

void Undistort::makeOptimalK_full()
{
	// todo
	assert(false);
}

void Undistort::readFromFile(const char* configFileName, int nPars, std::string prefix)
{
	valid = false;
	passthrough=false;
	remapX = 0;
	remapY = 0;
	
	float outputCalibration[5];

	parsOrg = VecX(nPars);

	// read parameters
	std::ifstream infile(configFileName);
	assert(infile.good());

    std::string l1,l2,l3,l4;

	std::getline(infile,l1);
	std::getline(infile,l2);
    std::getline(infile,l3);
    std::getline(infile,l4);

    // l1 & l2
    if(nPars == 5) // fov model
	{
		char buf[1000];
		snprintf(buf, 1000, "%s%%lf %%lf %%lf %%lf %%lf", prefix.c_str());

		if(std::sscanf(l1.c_str(), buf, &parsOrg[0], &parsOrg[1], &parsOrg[2], &parsOrg[3], &parsOrg[4]) == 5 &&
				std::sscanf(l2.c_str(), "%d %d", &wOrg, &hOrg) == 2)
		{
			printf("Input resolution: %d %d\n",wOrg, hOrg);
			printf("In: %f %f %f %f %f\n",
					parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3], parsOrg[4]);
		}
		else
		{
			printf("Failed to read camera calibration (invalid format?)\nCalibration file: %s\n", configFileName);
			infile.close();
			return;
		}
	}
	else if(nPars == 8) // KB, equi & radtan model
	{
		char buf[1000];
		snprintf(buf, 1000, "%s%%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf", prefix.c_str());

		if(std::sscanf(l1.c_str(), buf,
				&parsOrg[0], &parsOrg[1], &parsOrg[2], &parsOrg[3], &parsOrg[4],
				&parsOrg[5], &parsOrg[6], &parsOrg[7]) == 8 &&
				std::sscanf(l2.c_str(), "%d %d", &wOrg, &hOrg) == 2)
		{
			printf("Input resolution: %d %d\n",wOrg, hOrg);
			printf("In: %s%f %f %f %f %f %f %f %f\n",
					prefix.c_str(),
					parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3], parsOrg[4], parsOrg[5], parsOrg[6], parsOrg[7]);
		}
		else
		{
			printf("Failed to read camera calibration (invalid format?)\nCalibration file: %s\n", configFileName);
			infile.close();
			return;
		}
	}
	else
	{
		printf("called with invalid number of parameters.... forgot to implement me?\n");
		infile.close();
		return;
	}



    if(parsOrg[2] < 1 && parsOrg[3] < 1)
    {
        printf("\n\nFound fx=%f, fy=%f, cx=%f, cy=%f.\n I'm assuming this is the \"relative\" calibration file format,"
               "and will rescale this by image width / height to fx=%f, fy=%f, cx=%f, cy=%f.\n\n",
               parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3],
               parsOrg[0] * wOrg, parsOrg[1] * hOrg, parsOrg[2] * wOrg - 0.5, parsOrg[3] * hOrg - 0.5 );

        // rescale and substract 0.5 offset.
        // the 0.5 is because I'm assuming the calibration is given such that the pixel at (0,0)
        // contains the integral over intensity over [0,0]-[1,1], whereas I assume the pixel (0,0)
        // to contain a sample of the intensity ot [0,0], which is best approximated by the integral over
        // [-0.5,-0.5]-[0.5,0.5]. Thus, the shift by -0.5.
        parsOrg[0] = parsOrg[0] * wOrg;
        parsOrg[1] = parsOrg[1] * hOrg;
        parsOrg[2] = parsOrg[2] * wOrg - 0.5;
        parsOrg[3] = parsOrg[3] * hOrg - 0.5;
    }



	// l3
	if(l3 == "crop")
	{
		outputCalibration[0] = -1;
        printf("Out: Rectify Crop\n");
	}
	else if(l3 == "full")
	{
		outputCalibration[0] = -2;
        printf("Out: Rectify Full\n");
	}
	else if(l3 == "none")
	{
		outputCalibration[0] = -3;
        printf("Out: No Rectification\n");
	}
	else if(std::sscanf(l3.c_str(), "%f %f %f %f %f", &outputCalibration[0], &outputCalibration[1], &outputCalibration[2], &outputCalibration[3], &outputCalibration[4]) == 5)
	{
		printf("Out: %f %f %f %f %f\n",
				outputCalibration[0], outputCalibration[1], outputCalibration[2], outputCalibration[3], outputCalibration[4]);

	}
	else
	{
		printf("Out: Failed to Read Output pars... not rectifying.\n");
		infile.close();
		return;
	}


	// l4
	if(std::sscanf(l4.c_str(), "%d %d", &w, &h) == 2)
	{
		printf("Output resolution: %d %d\n",w, h);
	}
	else
	{
		printf("Out: Failed to Read Output resolution... not rectifying.\n");
		valid = false;
    }

    remapX = new float[w*h];
    remapY = new float[w*h];

	if(outputCalibration[0] == -1)
		makeOptimalK_crop();
	else if(outputCalibration[0] == -2)
		makeOptimalK_full();
	else if(outputCalibration[0] == -3)
	{
		if(w != wOrg || h != hOrg)
		{
			printf("ERROR: rectification mode none requires input and output dimenstions to match!\n\n");
			exit(1);
		}
		K.setIdentity();
        K(0,0) = parsOrg[0];
        K(1,1) = parsOrg[1];
        K(0,2) = parsOrg[2];
        K(1,2) = parsOrg[3];
		passthrough = true;
	}
	else
	{


        if(outputCalibration[2] > 1 || outputCalibration[3] > 1)
        {
            printf("\n\n\nWARNING: given output calibration (%f %f %f %f) seems wrong. It needs to be relative to image width / height!\n\n\n",
                   outputCalibration[0],outputCalibration[1],outputCalibration[2],outputCalibration[3]);
        }


		K.setIdentity();
        K(0,0) = outputCalibration[0] * w;
        K(1,1) = outputCalibration[1] * h;
        K(0,2) = outputCalibration[2] * w - 0.5;
        K(1,2) = outputCalibration[3] * h - 0.5;
	}

	for(int y=0;y<h;y++)
		for(int x=0;x<w;x++)
		{
			remapX[x+y*w] = x;
			remapY[x+y*w] = y;
		}

	distortCoordinates(remapX, remapY, remapX, remapY, h*w);


	for(int y=0;y<h;y++)
		for(int x=0;x<w;x++)
		{
			// make rounding resistant.
			float ix = remapX[x+y*w];
			float iy = remapY[x+y*w];

			if(ix == 0) ix = 0.001;
			if(iy == 0) iy = 0.001;
			if(ix == wOrg-1) ix = wOrg-1.001;
			if(iy == hOrg-1) ix = hOrg-1.001;

			if(ix > 0 && iy > 0 && ix < wOrg-1 &&  iy < wOrg-1)
			{
				remapX[x+y*w] = ix;
				remapY[x+y*w] = iy;
			}
			else
			{
				remapX[x+y*w] = -1;
				remapY[x+y*w] = -1;
			}
		}

	valid = true;




	printf("\nRectified Camera Matrix:\n");
	std::cout << K << "\n\n";

}


UndistortFOV::UndistortFOV(const char* configFileName)
{
    printf("Creating FOV undistorter\n");
    readFromFile(configFileName, 5, "FOV ");
}
UndistortFOV::~UndistortFOV()
{
}

void UndistortFOV::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{
	float dist = parsOrg[4];
	float d2t = 2.0f * tan(dist / 2.0f);



	// current camera parameters
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];



	float ofx = K(0,0);
	float ofy = K(1,1);
	float ocx = K(0,2);
	float ocy = K(1,2);

	for(int i=0;i<n;i++)
	{
		float x = in_x[i];
		float y = in_y[i];
		float ix = (x - ocx) / ofx;
		float iy = (y - ocy) / ofy;

		float r = sqrtf(ix*ix + iy*iy);
		float fac = (r==0 || dist==0) ? 1 : atanf(r * d2t)/(dist*r);

		ix = fx*fac*ix+cx;
		iy = fy*fac*iy+cy;

		out_x[i] = ix;
		out_y[i] = iy;
	}
}






UndistortRadTan::UndistortRadTan(const char* configFileName)
{
    printf("Creating RadTan undistorter\n");
    readFromFile(configFileName, 8,"RadTan ");
}
UndistortRadTan::~UndistortRadTan()
{
}

void UndistortRadTan::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{
    // RADTAN
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];
    float k1 = parsOrg[4];
    float k2 = parsOrg[5];
    float r1 = parsOrg[6];
    float r2 = parsOrg[7];

    float ofx = K(0,0);
    float ofy = K(1,1);
    float ocx = K(0,2);
    float ocy = K(1,2);



    for(int i=0;i<n;i++)
    {
        float x = in_x[i];
        float y = in_y[i];

        // RADTAN
        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;
        float mx2_u = ix * ix;
        float my2_u = iy * iy;
        float mxy_u = ix * iy;
        float rho2_u = mx2_u+my2_u;
        float rad_dist_u = k1 * rho2_u + k2 * rho2_u * rho2_u;
        float x_dist = ix + ix * rad_dist_u + 2.0 * r1 * mxy_u + r2 * (rho2_u + 2.0 * mx2_u);
        float y_dist = iy + iy * rad_dist_u + 2.0 * r2 * mxy_u + r1 * (rho2_u + 2.0 * my2_u);
        float ox = fx*x_dist+cx;
        float oy = fy*y_dist+cy;


        out_x[i] = ox;
        out_y[i] = oy;
    }


}



UndistortEquidistant::UndistortEquidistant(const char* configFileName)
{
    printf("Creating Equidistant undistorter\n");
    readFromFile(configFileName, 8,"EquiDistant ");
}
UndistortEquidistant::~UndistortEquidistant()
{
}

void UndistortEquidistant::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{
    // EQUI
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];
    float k1 = parsOrg[4];
    float k2 = parsOrg[5];
    float k3 = parsOrg[6];
    float k4 = parsOrg[7];


    float ofx = K(0,0);
    float ofy = K(1,1);
    float ocx = K(0,2);
    float ocy = K(1,2);



    for(int i=0;i<n;i++)
    {
        float x = in_x[i];
        float y = in_y[i];

        // EQUI
        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;
        float r = sqrt(ix * ix + iy * iy);
        float theta = atan(r);
        float theta2 = theta * theta;
        float theta4 = theta2 * theta2;
        float theta6 = theta4 * theta2;
        float theta8 = theta4 * theta4;
        float thetad = theta * (1 + k1 * theta2 + k2 * theta4 + k3 * theta6 + k4 * theta8);
        float scaling = (r > 1e-8) ? thetad / r : 1.0;
        float ox = fx*ix*scaling + cx;
        float oy = fy*iy*scaling + cy;

        out_x[i] = ox;
        out_y[i] = oy;
    }
}



UndistortKB::UndistortKB(const char* configFileName)
{
	printf("Creating KannalaBrandt undistorter\n");
    readFromFile(configFileName, 8,"KannalaBrandt ");
}
UndistortKB::~UndistortKB()
{
}

void UndistortKB::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{
    const float fx = parsOrg[0];
	const float fy = parsOrg[1];
	const float cx = parsOrg[2];
	const float cy = parsOrg[3];
    const float k0 = parsOrg[4];
    const float k1 = parsOrg[5];
    const float k2 = parsOrg[6];
    const float k3 = parsOrg[7];

    const float ofx = K(0,0);
    const float ofy = K(1,1);
    const float ocx = K(0,2);
    const float ocy = K(1,2);


	for(int i=0;i<n;i++)
	{
		float x = in_x[i];
		float y = in_y[i];

		// RADTAN
		float ix = (x - ocx) / ofx;
		float iy = (y - ocy) / ofy;

	    const float Xsq_plus_Ysq = ix*ix + iy*iy;
	    const float sqrt_Xsq_Ysq = sqrtf(Xsq_plus_Ysq);
	    const float theta = atan2f( sqrt_Xsq_Ysq, 1 );
	    const float theta2 = theta*theta;
	    const float theta3 = theta2*theta;
	    const float theta5 = theta3*theta2;
	    const float theta7 = theta5*theta2;
	    const float theta9 = theta7*theta2;
	    const float r = theta + k0*theta3 + k1*theta5 + k2*theta7 + k3*theta9;

	    if(sqrt_Xsq_Ysq < 1e-6)
	    {
	    	out_x[i] = fx * ix + cx;
	    	out_y[i] = fy * iy + cy;
	    }
	    else
	    {
	    	out_x[i] = (r / sqrt_Xsq_Ysq) * fx * ix + cx;
	    	out_y[i] = (r / sqrt_Xsq_Ysq) * fy * iy + cy;
	    }
	}
}





UndistortPinhole::UndistortPinhole(const char* configFileName)
{
	printf("Creating FOV undistorter\n");
    readFromFile(configFileName, 5,"Pinhole ");
}
UndistortPinhole::~UndistortPinhole()
{
}

void UndistortPinhole::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{
	// current camera parameters
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];

	float ofx = K(0,0);
	float ofy = K(1,1);
	float ocx = K(0,2);
	float ocy = K(1,2);

	for(int i=0;i<n;i++)
	{
		float x = in_x[i];
		float y = in_y[i];
		float ix = (x - ocx) / ofx;
		float iy = (y - ocy) / ofy;
		ix = fx*ix+cx;
		iy = fy*iy+cy;
		out_x[i] = ix;
		out_y[i] = iy;
	}
}
