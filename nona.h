#pragma once
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include "guidedfilter.h"
#include "interpolationmethods.h"
using namespace std;
using namespace cv;
void demosaic_nona(cv::Mat& Bayer, cv::Mat& Dst, int flag1 = 1, int flag2 = 1, int flag3 = 1, double k = 1.0, float sigma = 1.0);