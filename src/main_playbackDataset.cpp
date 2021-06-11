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


#include <locale.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>


#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>

#include "BenchmarkDatasetReader.h"

int main( int argc, char** argv )
{
	setlocale(LC_ALL, "");

	std::string dataset = argv[1];
	printf("Playback dataset %s!\n", dataset.c_str());

	char* trueBitDepthInput = argv[2];
	int option;
	int trueBitDepth;
	if(1==sscanf(trueBitDepthInput,"trueBitDepth=%d",&option))
	{
		trueBitDepth = option;
		printf("trueBitDepth set to %d!\n", trueBitDepth);
	}
	else
	{
		std::cout << "trueBitDepth not set, exiting.\n";
		return -1;
	}

	DatasetReader* reader = new DatasetReader(dataset, trueBitDepth);
	int saturationVal = reader->getPhotoUndistorter()->getSaturationVal();

	bool autoPlay = false;
	bool rect = false;
	bool removeGamma = false;
	bool removeVignette = false;
	bool killOverexposed = false;

	for(int i=0;i<reader->getNumImages();i++)
	{
		while(true)
		{
			ExposureImage* I = reader->getImage(i, rect, removeGamma, removeVignette, killOverexposed);
			cv::imshow("Image", cv::Mat(I->h, I->w, CV_32F, I->image) * (1 / (float)saturationVal));
			printf("Read image %d, time %.5f, exposure %.5fms. Rect (r): %d, remove gamma (g) %d, remove vignette (v): %d, kill overesposed (o)%d\n",
					I->id, I->timestamp, I->exposure_time,
					(int)rect, (int)removeGamma, (int)removeVignette, (int)killOverexposed);


			char k;
			if(autoPlay) k = cv::waitKey(1);
			else k = cv::waitKey(0);

			delete I;


			if(k==' ') break;
			if(k=='s' || k == 'S') {i+=30; break;};
			if(k=='a' || k == 'A') autoPlay=!autoPlay;
			if(k=='v' || k == 'V') removeVignette=!removeVignette;
			if(k=='g' || k == 'G') removeGamma=!removeGamma;
			if(k=='o' || k == 'O') killOverexposed=!killOverexposed;
			if(k=='r' || k == 'R') rect=!rect;
			if(autoPlay) break;
		}
	}

	delete reader;
}

