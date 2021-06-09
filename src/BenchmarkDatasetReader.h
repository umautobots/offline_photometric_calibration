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


#pragma once
#include <sstream>
#include <fstream>
#include <dirent.h>
#include <algorithm>
#include <math.h>

#include "opencv2/opencv.hpp"
#include "undistort/Undistort.h"
#include "ExposureImage.h"
#include "PhotometricUndistorter.h"


inline int getdir (std::string dir, std::vector<std::string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	std::string name = std::string(dirp->d_name);

    	if(name != "." && name != "..")
    		files.push_back(name);
    }
    closedir(dp);


    std::sort(files.begin(), files.end());

    if(dir.at( dir.length() - 1 ) != '/') dir = dir+"/";
	for(unsigned int i=0;i<files.size();i++)
	{
		if(files[i].at(0) != '/')
			files[i] = dir + files[i];
	}

    return files.size();
}



/*
 * provides read functionality for one of the dataset sequences.
 * * path is the folder for the sequence, with trailing slash (e.g. /home/Peter/datasets/sequenceX/
 * 		=> if it contains a folder "images", images will be read from there.
 */
class DatasetReader
{
public:
	DatasetReader(std::string folder, int trueBitDepth)
	{
		this->path = folder;
		for(int i=0;i<3;i++)
		{
			undistorter=0;
		}

		getdir (path+"images/", files);
		printf("Load Dataset %s: found %d files in folder /images; assuming that all images are there.\n",
				path.c_str(), (int)files.size());
		loadTimestamps(path+"times.txt");

		// Set the saturation value, image type, number of zero-padded least significant bits
		int saturationVal = pow(2, trueBitDepth) - 1;
		cv::Mat tmpImg = cv::imread(files[0], cv::IMREAD_ANYDEPTH);
		imageType = tmpImg.type();
		if(imageType == CV_8UC1)
		{
			lsbPadNum = 8 - trueBitDepth;
		}
		else if(imageType == CV_16UC1)
		{
			lsbPadNum = 16 - trueBitDepth;
		}
		else
		{
			throw std::runtime_error("ERROR: input image type not supported.\n");
		}
		if(lsbPadNum < 0)
		{
			throw std::runtime_error("ERROR: The trueBitDepth exceeds the bit depth of the input image.\n");
		}

		// create undistorter.
		undistorter = Undistort::getUndistorterForFile((path+"camera.txt").c_str());
		photoUndistorter = new PhotometricUndistorter(path+"pcalib.txt", path+"vignette.png",
			undistorter->getOriginalSize()[0],undistorter->getOriginalSize()[1], saturationVal);


		// get image widths.
		widthOrg = undistorter->getOriginalSize()[0];
		heightOrg = undistorter->getOriginalSize()[1];
		width= undistorter->getSize()[0];
		height= undistorter->getSize()[1];

		internalTempBuffer=new float[widthOrg*heightOrg];

		printf("Dataset %s: Got %d files!\n", path.c_str(), (int)getNumImages());
	}
	~DatasetReader()
	{
		if(undistorter!=0) delete undistorter;
		if(photoUndistorter!=0) delete photoUndistorter;
		delete[] internalTempBuffer;

	}

	Undistort* getUndistorter()
	{
		return undistorter;
	}

	PhotometricUndistorter* getPhotoUndistorter()
	{
		return photoUndistorter;
	}

	int getNumImages()
	{
		return files.size();
	}

	double getTimestamp(int id)
	{
		if(id >= (int)timestamps.size()) return 0;
		if(id < 0) return 0;
		return timestamps[id];
	}

	float getExposure(int id)
	{
		if(id >= (int)exposures.size()) return 0;
		if(id < 0) return 0;
		return exposures[id];
	}

	ExposureImage* getImage(int id, bool rectify, bool removeGamma, bool removeVignette, bool nanOverexposed)
	{
		cv::Mat imageRaw = getImageRaw_internal(id);
		ExposureImage* ret=0;

		if(removeGamma || removeVignette || nanOverexposed)
		{
			if(!rectify)
			{
				// photo undist only.
				ret = new ExposureImage(widthOrg, heightOrg, timestamps[id], exposures[id], id);
				if(imageType == CV_8UC1)
				{
					photoUndistorter->unMapImage(imageRaw.data, ret->image, widthOrg*heightOrg, removeGamma,
						removeVignette, nanOverexposed );
				}
				else if(imageType == CV_16UC1)
				{
					photoUndistorter->unMapImage(imageRaw.ptr<unsigned short>(), ret->image, widthOrg*heightOrg,
						removeGamma, removeVignette, nanOverexposed );
				}
			}
			else
			{
				// photo undist to buffer, then rect
				ret = new ExposureImage(width, height, timestamps[id], exposures[id], id);
				if(imageType == CV_8UC1)
				{
					photoUndistorter->unMapImage(imageRaw.data, internalTempBuffer, widthOrg*heightOrg, removeGamma,
						removeVignette, nanOverexposed );
				}
				else if(imageType == CV_16UC1)
				{
					photoUndistorter->unMapImage(imageRaw.ptr<unsigned short>(), internalTempBuffer, widthOrg*heightOrg,
						removeGamma, removeVignette, nanOverexposed );
				}
				undistorter->undistort(internalTempBuffer, ret->image, widthOrg*heightOrg, width*height);
			}
		}
		else
		{
			if(rectify)
			{
				// rect only.
				ret = new ExposureImage(width, height, timestamps[id], exposures[id], id);
				if(imageType == CV_8UC1)
				{
					undistorter->undistort(imageRaw.data, ret->image, widthOrg*heightOrg, width*height);
				}
				else if(imageType == CV_16UC1)
				{
					undistorter->undistort(imageRaw.ptr<unsigned short>(), ret->image, widthOrg*heightOrg,
						width*height);
				}
			}
			else
			{
				// do nothing.
				ret = new ExposureImage(widthOrg, heightOrg, timestamps[id], exposures[id], id);
				if(imageType == CV_8UC1)
				{
					for(int i=0;i<widthOrg*heightOrg;i++)
						ret->image[i] = imageRaw.at<unsigned char>(i);
				}
				else if(imageType == CV_16UC1)
				{
					for(int i=0;i<widthOrg*heightOrg;i++)
						ret->image[i] = imageRaw.at<unsigned short>(i);
				}
			}
		}
		return ret;
	}



	cv::Mat getImageRaw_internal(int id)
	{
		// Check the index is in range
		if(!(id >= 0 && id < (int)files.size()))
		{
			throw std::runtime_error("ERROR: image index out of range.\n");
		}

		// Import the image
		cv::Mat imageRaw = cv::imread(files[id], cv::IMREAD_ANYDEPTH);

		// Check the image parameters against expected values
		if(!(imageRaw.rows == heightOrg && imageRaw.cols == widthOrg && imageRaw.type() == imageType))
		{
			throw std::runtime_error("ERROR: input image parameters do not match expected values.\n");
		}

		// Remove the zero-padded least significant bits
		imageRaw = imageRaw / pow(2, lsbPadNum);

		return imageRaw;
	}

	int getImageType()
	{
		return imageType;
	}


private:


	inline void loadTimestamps(std::string timesFile)
	{
		std::ifstream tr;
		tr.open(timesFile.c_str());
		timestamps.clear();
		exposures.clear();
		while(!tr.eof() && tr.good())
		{
			std::string line;
			char buf[1000];
			tr.getline(buf, 1000);

			int id;
			double stamp;
			float exposure = 0;

			if(3 == sscanf(buf, "%d %lf %f", &id, &stamp, &exposure))
			{
				timestamps.push_back(stamp);
				exposures.push_back(exposure);
			}

			else if(2 == sscanf(buf, "%d %lf", &id, &stamp))
			{
				timestamps.push_back(stamp);
				exposures.push_back(0);
			}
		}
		tr.close();

		if((int)exposures.size()!=(int)getNumImages())
		{
			printf("DatasetReader: Mismatch between number of images and number of timestamps / exposure times. Set all to zero.");
			timestamps.clear();
			exposures.clear();
			for(int i=0;i<(int)getNumImages();i++)
			{
				timestamps.push_back(0.0);
				exposures.push_back(0);
			}
		}
	}

	// data is here.
	std::vector<std::string> files;
	std::vector<double> timestamps;
	std::vector<float> exposures;

	int width, height;
	int widthOrg, heightOrg;

	std::string path;

	// type of images. assumes every image has the same type
	int imageType;

	// the number of zero-padded least significant bits
	int lsbPadNum;

	// internal structures.
	Undistort* undistorter;
	PhotometricUndistorter* photoUndistorter;

	float* internalTempBuffer;
};

