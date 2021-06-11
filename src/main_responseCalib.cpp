/* Copyright (c) 2016, Jakob Engel
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimer in the documentation and/or 
 * other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its contributors 
 * may be used to endorse or promote products derived from this software without 
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 */




#include "opencv2/opencv.hpp"
#include "opencv2/video/tracking.hpp"

#include <sstream>
#include <fstream>
#include <dirent.h>
#include <algorithm>
#include <math.h>

#include "BenchmarkDatasetReader.h"


int leakPadding=2;
int nits = 10;
int skipFrames = 1;

// The true (embedded) pixel bit depth
int trueBitDepth = 12;
int saturationVal = pow(2, trueBitDepth) - 1;
int numVals = saturationVal + 1;


template<typename T>
Eigen::Vector2d rmse(double* G, double* E, std::vector<double> &exposureVec, std::vector<T*> &dataVec,  int wh)
{
	long double e=0;		// yeah - these will be sums of a LOT of values, so we need super high precision.
	long double num=0;

	int n = dataVec.size();
	for(int i=0;i<n;i++)
	{
		for(int k=0;k<wh;k++)
		{
			if(dataVec[i][k] == saturationVal) continue;
			double r = G[dataVec[i][k]] - exposureVec[i]*E[k];
			if(!std::isfinite(r)) continue;
			e += r*r*1e-10;
			num++;
		}
	}

	return Eigen::Vector2d(1e5*sqrtl((e/num)), (double)num);
}
template Eigen::Vector2d rmse<unsigned char>(double* G, double* E, std::vector<double> &exposureVec,
	std::vector<unsigned char*> &dataVec,  int wh);
template Eigen::Vector2d rmse<unsigned short>(double* G, double* E, std::vector<double> &exposureVec,
	std::vector<unsigned short*> &dataVec,  int wh);


void plotE(double* E, int w, int h, std::string saveTo="")
{
	double Emin=1e10, Emax=-1e10;

	for(int i=0;i<w*h;i++)
	{
		if(E[i] < Emin) Emin = E[i];
		if(E[i] > Emax) Emax = E[i];
	}

	cv::Mat EImg16 = cv::Mat(h,w,CV_16U);

	for(int i=0;i<w*h;i++)
	{
		EImg16.at<ushort>(i) = 255* 255* (E[i]-Emin) / (Emax-Emin);
	}

	printf("Irradiance %f - %f\n", Emin, Emax);

	if(saveTo != "")
	{
		cv::imwrite(saveTo+".png", EImg16);
	}
}
void plotG(double* G, std::string saveTo="")
{
	cv::Mat GImg = cv::Mat(numVals,numVals,CV_32FC1);
	GImg.setTo(0);

	double min=1e10, max=-1e10;

	for(int i=0;i<numVals;i++)
	{
		if(G[i] < min) min = G[i];
		if(G[i] > max) max = G[i];
	}

	for(int i=0;i<numVals;i++)
	{
		double val = numVals*(G[i]-min) / (max-min);
		for(int k=0;k<numVals;k++)
		{
			if(val < k)
				GImg.at<float>(numVals - k,i) = k-val;
		}
	}

	printf("Inv. Response %f - %f\n", min, max);
	if(saveTo != "") cv::imwrite(saveTo, GImg*255);
}


void parseArgument(char* arg)
{
	int option;

	if(1==sscanf(arg,"leakPadding=%d",&option))
	{
		leakPadding = option;
		printf("leakPadding set to %d!\n", leakPadding);
		return;
	}
	if(1==sscanf(arg,"iterations=%d",&option))
	{
		nits = option;
		printf("nits set to %d!\n", nits);
		return;
	}
	if(1==sscanf(arg,"skip=%d",&option))
	{
		skipFrames = option;
		printf("skipFrames set to %d!\n", skipFrames);
		return;
	}
	if(1==sscanf(arg,"trueBitDepth=%d",&option))
	{
		trueBitDepth = option;
		saturationVal = pow(2, trueBitDepth) - 1;
		numVals = saturationVal + 1;
		printf("trueBitDepth set to %d!\n", trueBitDepth);
		return;
	}

	printf("could not parse argument \"%s\"!!\n", arg);
}


template<typename T>
void loadData(std::vector<T*> &dataVec, std::vector<double> &exposureVec, int &w, int &h, DatasetReader* reader)
{
	w = 0;
	h = 0;

	for(int i=0;i<reader->getNumImages();i+=skipFrames)
	{
		cv::Mat img = reader->getImageRaw_internal(i);
		if(img.rows==0 || img.cols==0) continue;

		if((w!=0 && w != img.cols) || img.cols==0)
		{ printf("width mismatch!\n"); exit(1); };
		if((h!=0 && h != img.rows) || img.rows==0)
		{ printf("height mismatch!\n"); exit(1); };
		w = img.cols;
		h = img.rows;

		T* data = new T[img.rows*img.cols];
		memcpy(data, img.data, sizeof(T)*img.rows*img.cols);
		dataVec.push_back(data);
		exposureVec.push_back((double)(reader->getExposure(i)));
	}

	int maxVal = 0;
	for(int i=0;i<dataVec.size();i++)
	{
		for(int k=0;k<w*h;k++)
		{
			int b = dataVec[i][k];
			if(maxVal < b) maxVal = b;
		}
	}
	if(maxVal < saturationVal)
	{
		printf("The maximum pixel value across all images is %i, this is lower than the expected saturation value %i,"
			" the saturation value will be set to the maximum value.\n", maxVal, saturationVal);
		saturationVal = maxVal;
		numVals = saturationVal + 1;
	}

	for(int i=0;i<reader->getNumImages();i+=skipFrames)
	{
		T* data = dataVec[i];
		T* data2 = new T[w*h];
		for(int it=0;it<leakPadding;it++)
		{
			memcpy(data2, data, sizeof(T)*w*h);
			for(int y=1;y<h-1;y++)
				for(int x=1;x<w-1;x++)
				{
					if(data[x+y*w]==saturationVal)
					{
						data2[x+1 + w*(y+1)] = saturationVal;
						data2[x+1 + w*(y  )] = saturationVal;
						data2[x+1 + w*(y-1)] = saturationVal;

						data2[x   + w*(y+1)] = saturationVal;
						data2[x   + w*(y  )] = saturationVal;
						data2[x   + w*(y-1)] = saturationVal;

						data2[x-1 + w*(y+1)] = saturationVal;
						data2[x-1 + w*(y  )] = saturationVal;
						data2[x-1 + w*(y-1)] = saturationVal;
					}
				}
			memcpy(data, data2, sizeof(T)*w*h);
		}
		delete[] data2;
	}
}
template void loadData<unsigned char>(std::vector<unsigned char*> &dataVec, std::vector<double> &exposureVec, 
	int &w, int &h, DatasetReader* reader);
template void loadData<unsigned short>(std::vector<unsigned short*> &dataVec, std::vector<double> &exposureVec,
	int &w, int &h, DatasetReader* reader);

template<typename T>
void optimize(std::vector<T*> &dataVec, std::vector<double> exposureVec, int w, int h)
{
	int n = dataVec.size();


	printf("loaded %d images\n", n);


	double* E = new double[w*h];		// scene irradiance
	double* En = new double[w*h];		// scene irradiance
	double* G = new double[numVals];		// inverse response function

	// set starting scene irradiance to mean of all images.
	memset(E,0,sizeof(double)*w*h);
	memset(En,0,sizeof(double)*w*h);
	memset(G,0,sizeof(double)*(numVals));
	for(int i=0;i<n;i++)
		for(int k=0;k<w*h;k++)
		{
			E[k] += dataVec[i][k];
			En[k] ++;
		}
	for(int k=0;k<w*h;k++) E[k] = E[k]/En[k];



	if(-1 == system("rm -rf photoCalibResult")) printf("could not delete old photoCalibResult folder!\n");
	if(-1 == system("mkdir photoCalibResult")) printf("could not create photoCalibResult folder!\n");


	std::ofstream logFile;
	logFile.open("photoCalibResult/log.txt", std::ios::trunc | std::ios::out);
	logFile.precision(15);


	printf("init RMSE = %f! \t", rmse(G, E, exposureVec, dataVec, w*h )[0]);
	plotE(E,w,h, "photoCalibResult/E-0");


	bool optE = true;
	bool optG = true;

	Eigen::Vector2d err;

	for(int it=0;it<nits;it++)
	{
		if(optG)
		{
			// optimize log inverse response function.
			double* GSum = new double[numVals];
			double* GNum = new double[numVals];
			memset(GSum,0,numVals*sizeof(double));
			memset(GNum,0,numVals*sizeof(double));
			for(int i=0;i<n;i++)
			{
				for(int k=0;k<w*h;k++)
				{
					int b = dataVec[i][k];
					if(b == saturationVal) continue;
					GNum[b]++;
					GSum[b]+= E[k] * exposureVec[i];
				}
			}
			for(int i=0;i<numVals;i++)
			{
				G[i] = GSum[i] / GNum[i];
				if(!std::isfinite(G[i]) && i > 1) G[i] = G[i-1] + (G[i-1]-G[i-2]);
			}
			delete[] GSum;
			delete[] GNum;
			printf("optG RMSE = %f! \t", rmse(G, E, exposureVec, dataVec, w*h )[0]);

			char buf[1000]; snprintf(buf, 1000, "photoCalibResult/G-%d.png", it+1);
			plotG(G, buf);
		}





		if(optE)
		{
			// optimize scene irradiance function.
			double* ESum = new double[w*h];
			double* ENum = new double[w*h];
			memset(ESum,0,w*h*sizeof(double));
			memset(ENum,0,w*h*sizeof(double));
			for(int i=0;i<n;i++)
			{
				for(int k=0;k<w*h;k++)
				{
					int b = dataVec[i][k];
					if(b == saturationVal) continue;
					ENum[k] += exposureVec[i]*exposureVec[i];
					ESum[k] += (G[b]) * exposureVec[i];
				}
			}
			for(int i=0;i<w*h;i++)
			{
				E[i] = ESum[i] / ENum[i];
				if(E[i] < 0) E[i] = 0;
			}

			delete[] ENum;
			delete[] ESum;
			printf("OptE RMSE = %f!  \t", rmse(G, E, exposureVec, dataVec, w*h )[0]);

			char buf[1000]; snprintf(buf, 1000, "photoCalibResult/E-%d", it+1);
			plotE(E,w,h, buf);
		}


		// rescale such that maximum response is saturationVal (fairly arbitrary choice).
		double rescaleFactor= double(saturationVal) / G[saturationVal];
		for(int i=0;i<w*h;i++)
		{
			E[i] *= rescaleFactor;
			if(i<numVals) G[i] *= rescaleFactor;
		}
		err = rmse(G, E, exposureVec, dataVec, w*h );
		printf("resc RMSE = %f!  \trescale with %f!\n",  err[0], rescaleFactor);

		logFile << it << " " << n << " " << err[1] << " " << err[0] << "\n";

		// Save the current result.
		char buf[1000]; snprintf(buf, 1000, "photoCalibResult/pcalib-%d.txt", it+1);
		std::ofstream lg;
		lg.open(buf, std::ios::trunc | std::ios::out);
		lg.precision(15);
		for(int i=0;i<numVals;i++)
			lg << G[i] << " ";
		lg << "\n";
		lg.flush();
		lg.close();
	}

	// Check if the inverse response is monotonically increasing.
	bool monotonic = true;
	for(int i=0;i<saturationVal;i++)
	{
		if(G[i+1] <= G[i])
		{
			printf("Final G is not monotonically increasing, gaps will be interpolated over.\n");
			monotonic = false;
			break;
		}
	}

	// Interpolate over gaps if G is not mononically increasing.
	if(!monotonic)
	{
		// Interpolate between each point and the next higher point.
		int gap;
		int j = 0;
		for(int i = 1; i <= saturationVal; i++)
		{
			if(G[i] > G[j])
			{
				gap = i - j;
				for(int k = 1; k < gap; k++)
				{
					G[i - k] = G[i] - (G[i] - G[j]) * k / gap;
				}
				j = i;
			}
		}

		// If G[saturationVal] is not the largest output, extrapolate.
		if(j != saturationVal)
		{
			// Extrapolate out from the last increasing segment (assuming max(G) != G[0]).
			double slope = G[j] - G[j - 1];
			for(int i = j + 1; i <= saturationVal; i++)
			{
				G[i] = G[j] + slope * (i - j);
			}

			// Rescale results a final time such that the maximum inverse response is saturationVal.
			double rescaleFactor= double(saturationVal) / G[saturationVal];
			for(int i = 0; i < w * h; i++)
			{
				E[i] *= rescaleFactor;
				if(i<numVals) G[i] *= rescaleFactor;
			}
		}

		// Compute the final RMSE.
		err = rmse(G, E, exposureVec, dataVec, w*h );
		printf("RMSE after monotonic correction = %f!\n",  err[0]);
		logFile << "Monotonic Correction " << n << " " << err[1] << " " << err[0] << "\n";

		// Save the final G plot.
		plotG(G, "photoCalibResult/G-monotonic-correction.png");
	}

	logFile.flush();
	logFile.close();

	// Save the final result.
	std::ofstream lg;
	lg.open("photoCalibResult/pcalib.txt", std::ios::trunc | std::ios::out);
	lg.precision(15);
	for(int i=0;i<numVals;i++)
		lg << G[i] << " ";
	lg << "\n";
	lg.flush();
	lg.close();

	delete[] E;
	delete[] En;
	delete[] G;
	for(int i=0;i<n;i++) delete[] dataVec[i];
}
template void optimize<unsigned char>(std::vector<unsigned char*> &dataVec, std::vector<double> exposureVec, int w, 
	int h);
template void optimize<unsigned short>(std::vector<unsigned short*> &dataVec, std::vector<double> exposureVec, int w, 
	int h);

int main( int argc, char** argv )
{
	// parse arguments
	for(int i=2; i<argc;i++)
		parseArgument(argv[i]);

	std::vector<double> exposureVec;
	int w, h;
	DatasetReader* reader = new DatasetReader(argv[1], trueBitDepth);
	int image_type = reader->getImageType();
	if(image_type == CV_8UC1)
	{
		std::vector<unsigned char*> dataVec;
		loadData(dataVec, exposureVec, w, h, reader);
		optimize(dataVec, exposureVec, w, h);
	}
	if(image_type == CV_16UC1)
	{
		std::vector<unsigned short*> dataVec;
		loadData(dataVec, exposureVec, w, h, reader);
		optimize(dataVec, exposureVec, w, h);
	}
	else
	{
		printf("ERROR: OpenCV image type %i not supported.\n", image_type);
		return -1;
	}

	return 0;
}
