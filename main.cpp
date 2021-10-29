#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <vector>
#include <opencv2/opencv.hpp>
#include "nona.h"
using namespace std;
using namespace cv;

int main() {
	string inFileName = "input.raw";
	int rows = 9000;
	int cols = 12000;
	Mat raw, out;
	raw.create(rows, cols, CV_16UC1);
	ifstream inFile;
	inFile.open(inFileName, ios::binary);
	if (!inFile.is_open()) {
		cout << "Unable to Open File" << endl;
		return 0;
	}
	inFile.unsetf(std::ios::skipws);
	std::streampos fileSize;
	inFile.seekg(0, std::ios::end);
	fileSize = inFile.tellg();
	inFile.seekg(0, std::ios::beg);
	inFile.read((char*)raw.data, rows * cols * sizeof(short));
	inFile.close();
	cout << "CV : " << CV_VERSION << endl;
	cout << "rows : " << raw.rows << endl;
	cout << "cols : " << raw.cols << endl;

	cout << "Preprocessing (black level subtraction + white balance)" << endl;
	float offset = 40.0f;
	float gainR = 1.6484f;
	float gainB = 1.3642f;
	raw.convertTo(raw, CV_32FC1);
	preprocessing(raw, offset, gainR, gainB);
	cout << "Demosaicing" << endl;
	int64 t0 = cv::getTickCount();
	int flag1, flag2, flag3;
	double k;
	cout << "Enter Method flag for Linear Interpolation (1-4), Difference Interpolation(1 OR 4), Bilinear Interpolation( 1 OR 4)" << endl;
	cin >> flag1 >> flag2 >> flag3;
	if (flag1 == 4 || flag2 == 4 || flag3 == 4) {
		cout << "Enter the distance - weightage relation value" << endl;
	}
	demosaic_nona(raw, out, flag1, flag2, flag3, k);

	int64 t1 = cv::getTickCount();
	double secs = (t1 - t0) / cv::getTickFrequency();
	cout << " Demosaicing : " << secs << endl;

	out.convertTo(out, CV_8UC3, 255.0 / 1023.0);
	cv::imwrite("output.bmp", out);
	return 0;
}