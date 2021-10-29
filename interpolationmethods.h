#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <vector>
using namespace std;
using namespace cv;
Mat bilinearfilter_method1(Mat residual);
Mat bilinearfilter_method4(Mat residual, double k);
Mat diffilter_H_method1(Mat Src1ch);
Mat diffilter_V_method1(Mat Src1ch);
Mat diffilter_H_method4(Mat Src1ch, double k);
Mat diffilter_V_method4(Mat Src1ch, double k);
Mat method1_h(Mat Src1ch);
Mat method1_v(Mat Src1ch);
Mat method2_h(Mat Src1ch);
Mat method2_v(Mat Src1ch);
Mat method3_h(Mat Src1ch);
Mat method3_v(Mat Src1ch);
Mat method4_h(Mat Src1ch, float k);
Mat method4_v(Mat Src1ch, float k);
