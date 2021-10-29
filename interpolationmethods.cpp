#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
using namespace std;
using namespace cv;
Mat bilinearfilter_method1(Mat residual) {
	Mat result = Mat::zeros(residual.rows, residual.cols, CV_32F);
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			result.at<uchar>(row, col) += residual.at<uchar>(row - 3, col - 3) * (0.25);
			result.at<uchar>(row, col) += residual.at<uchar>(row - 3, col) * (0.5);
			result.at<uchar>(row, col) += residual.at<uchar>(row - 3, col + 3) * (0.25);
			result.at<uchar>(row, col) += residual.at<uchar>(row, col - 3) * (0.5);
			result.at<uchar>(row, col) += residual.at<uchar>(row, col);
			result.at<uchar>(row, col) += residual.at<uchar>(row, col + 3) * (0.5);
			result.at<uchar>(row, col) += residual.at<uchar>(row + 3, col - 3) * (0.25);
			result.at<uchar>(row, col) += residual.at<uchar>(row + 3, col) * (0.5);
			result.at<uchar>(row, col) += residual.at<uchar>(row + 3, col + 3) * (0.25);
		}
	}
	return result;
}
void weight_generator_bilinear(vector<vector<double>>& weights, double k) {
	for (int i = 0; i < weights.size(); i++) {
		vector<int> sum(0, 9);
		vector<double> factors(0, 9);
		for (int j = 0; j < 9; j++) {

		}
	}
	return;
}
Mat bilinearfilter_method4(Mat residual, double k) {
	Mat result = Mat::zeros(residual.rows, residual.cols, CV_32F);
	int rowk, colk;
	int top_left_row_start, top_left_col_start;
	int top_middle_row_start, top_middle_col_start;
	int top_right_row_start, top_right_col_start;
	int middle_left_row_start, middle_left_col_start;
	int middle_middle_row_start, middle_middle_col_start;
	int middle_right_row_start, middle_right_col_start;
	int bottom_left_row_start, bottom_left_col_start;
	int bottom_middle_row_start, bottom_middle_col_start;
	int bottom_right_row_start, bottom_right_col_start;
	int w;
	int tl, tm, tr, ml, mm, mr, bl, bm, br;
	vector<vector<double>> distance(9, vector<double>(0, 81));
	vector<vector<double>> weights(9, vector<double>(0, 81));
	weight_generator_bilinear(weights, k);
	//for (int row = 0; row < 3; row++) {//TOP-MOST PIXELS
	//	for (int col = 3; col < residual.cols - 3; col++) {

	//	}
	//}
	//for (int col = 0; col < 3; col++) {//LEFT-MOST PIXELS
	//	for (int row = 3; row < residual.rows - 3; row++) {

	//	}
	//}
	//for (int row = residual.rows - 1; row > residual.rows - 4; row--) {//BOTTOM-MOST PIXELS
	//	for (int col = 3; col < residual.cols - 3; col++) {

	//	}
	//}
	//for (int col = residual.cols - 1; col > residual.cols - 1; col--) {//RIGHT-MOST PIXELS
	//	for (int row = 3; row < residual.rows - 3; row++) {

	//	}
	//}
	for (int row = 3; row < residual.rows - 3; row++) {
		for (int col = 3; col < residual.cols - 3; col++) {
			rowk = row / 3;
			colk = col / 3;
			top_left_row_start = row - rowk - 3;
			top_left_col_start = col - colk - 3;
			top_middle_row_start = row - rowk - 3;
			top_middle_col_start = col - colk;
			top_right_row_start = row - rowk - 3;
			top_right_col_start = col - colk + 3;
			middle_left_row_start = row - rowk;
			middle_left_col_start = col - colk - 3;
			middle_middle_row_start = row - rowk;
			middle_middle_col_start = col - colk;
			middle_right_row_start = row - rowk;
			middle_right_col_start = col - colk + 3;
			bottom_left_row_start = row - rowk + 3;
			bottom_left_col_start = col - colk - 3;
			bottom_middle_row_start = row - rowk + 3;
			bottom_middle_col_start = col - colk;
			bottom_right_row_start = row - rowk + 3;
			bottom_right_col_start = col - colk + 3;
			rowk = row % 3;
			colk = col % 3;
			tl = 0;
			tm = 9;
			tr = 18;
			ml = 27;
			mm = 36;
			mr = 45;
			bl = 54;
			bm = 63;
			br = 72;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((residual.at<uchar>(top_left_row_start + i, top_left_col_start + j)) * (weights[w][tl++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(top_middle_row_start + i, top_middle_col_start + j)) * (weights[w][tm++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(top_right_row_start + i, top_right_col_start + j)) * (weights[w][tr++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(middle_left_row_start + i, middle_left_col_start + j)) * (weights[w][ml++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(middle_middle_row_start + i, middle_middle_col_start + j)) * (weights[w][mm++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(middle_right_row_start + i, middle_right_col_start + j)) * (weights[w][mr++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(bottom_left_row_start + i, top_left_col_start + j)) * (weights[w][bl++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(bottom_middle_row_start + i, bottom_middle_col_start + j)) * (weights[w][bm++]));
					result.at<uchar>(row, col) += ((residual.at<uchar>(bottom_right_row_start + i, bottom_right_col_start + j)) * (weights[w][br++]));
				}
			}

		}
	}
	return result;
}
Mat diffilter_H_method1(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int colend = (Src1ch.cols / 3) - 1;
	for (int row = 0; row < result.rows; row++) {
		for (int col = 0; col < 3; col++) {
			result.at<uchar>(row, col) = 0 - Src1ch.at<uchar>(row, col + 3);
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = result.cols - 1; col > result.cols - 4; col--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row, col - 3);
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col + 3)) * (-1.0)) + ((Src1ch.at<uchar>(row, col - 3)) * (1.0));
		}
	}
	return result;
}
Mat diffilter_V_method1(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowend = (Src1ch.rows / 3) - 1;
	for (int col = 0; col < result.cols; col++) {
		for (int row = 0; row < 3; row++) {
			result.at<uchar>(row, col) = 0 - Src1ch.at<uchar>(row + 3, col);
		}
	}
	for (int col = 0; col < result.cols; col++) {
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row - 3, col);
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + 3, col)) * (-1.0)) + ((Src1ch.at<uchar>(row - 3, col)) * (1.0));
		}
	}
	return result;
}
void weight_generator_diff(vector<vector<double>>& weights, vector<vector<double>>& distance, double k) {
	double sum1, sum2, factor1, factor2;
	for (int i = 0; i < weights.size(); i++) {
		sum1 = 0;
		sum2 = 0;
		factor1 = 1;
		factor2 = -1;
		for (int j = 0; j < 9; j++) {
			sum1 += pow(distance[i][j], k);
		}
		for (int j = 9; j < 18; j++) {
			sum2 += pow(distance[i][j], k);
		}
		factor1 /= sum1;
		factor2 /= sum2;
		for (int j = 0; j < 9; j++) {
			weights[i][j] = factor1 * pow(distance[i][j], k);
		}
		for (int j = 9; j < 18; j++) {
			weights[i][j] = factor2 * pow(distance[i][j], k);
		}
	}
	return;
}
void distance_generator(vector<vector<double>>& distance, bool horizontalflag) {
	int rowk, colk;
	int rowk2, colk2;
	int sum;
	for (int i = 0; i < distance.size(); i++) {
		colk = i % 3;
		rowk = i / 3;
		sum = 0;
		for (int j = 0; j < 9; j++) {
			colk2 = j % 3;
			rowk2 = j / 3;
			if (horizontalflag) {
				if (colk2 <= colk) {
					sum = abs(colk2 - colk) + abs(rowk2 - rowk) + 3;
				}
				else {
					sum = 3 - abs(colk2 - colk) + abs(rowk2 - rowk);
				}
			}
			else {
				sum = abs(colk2 - colk) + abs(rowk2 - rowk);
				if (rowk2 <= rowk) {
					sum = abs(colk2 - colk) + abs(rowk2 - rowk) + 3;
				}
				else {
					sum = 3 - abs(rowk2 - rowk) + abs(colk2 - colk);
				}
			}
			distance[i][j] = sum;
		}
		for (int j = 9; j < 18; j++) {
			colk2 = (j - 9) % 3;
			rowk2 = (j - 9) / 3;
			if (horizontalflag) {
				if (colk2 >= colk) {
					sum = abs(colk2 - colk) + abs(rowk2 - rowk) + 3;
				}
				else {
					sum = 3 - abs(colk2 - colk) + abs(rowk2 - rowk);
				}
			}
			else {
				sum = abs(colk2 - colk) + abs(rowk2 - rowk);
				if (rowk2 >= rowk) {
					sum = abs(colk2 - colk) + abs(rowk2 - rowk) + 3;
				}
				else {
					sum = 3 - abs(rowk2 - rowk) + abs(colk2 - colk);
				}
			}
			distance[i][j] = sum;
		}
	}
}
Mat diffilter_H_method4(Mat Src1ch, double k) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int colend = (Src1ch.cols / 3) - 1;
	int rowk;
	int colk;
	int left_startrow, left_startcol, right_startrow, right_startcol;
	int l, r, w;
	vector<vector<double>> weights_h(9, vector<double>(18, 0));
	vector<vector<double>> dist_h(9, vector<double>(18, 0));
	distance_generator(dist_h, 1);
	weight_generator_diff(weights_h, dist_h, k);
	for (int row = 0; row < result.rows; row++) {//LEFT - MOST PIXELS
		for (int col = 0; col < 3; col++) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col + 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col + 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {//RIGHT - MOST PIXELS
		for (int col = result.cols - 1; col > result.cols - 4; col--) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col - 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col - 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col - 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col + 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	return result;
}
Mat diffilter_V_method4(Mat Src1ch, double k) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowend = (Src1ch.rows / 3) - 1;
	int rowk;
	int colk;
	int top_startrow, top_startcol, bottom_startrow, bottom_startcol;
	int l, r, w;
	vector<vector<double>> weights_v(9, vector<double>(18, 0));
	vector<vector<double>> dist_v(9, vector<double>(18, 0));
	distance_generator(dist_v, 0);
	weight_generator_diff(weights_v, dist_v, k);
	for (int col = 0; col < result.cols; col++) {//TOP - MOST PIXELS
		for (int row = 0; row < 3; row++) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row + 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row + 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	for (int col = 0; col < result.cols; col++) {//BOTTOM-MOST PIXELS
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row - 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row - 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row - 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row + 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	return result;
}
Mat method1_h(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int colend = (Src1ch.cols / 3) - 1;
	for (int row = 0; row < result.rows; row++) {
		for (int col = 0; col < 3; col++) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row, col + 3);
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = result.cols - 1; col > result.cols - 4; col--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row, col - 3);
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col + 3)) * 0.5) + ((Src1ch.at<uchar>(row, col - 3)) * 0.5);
		}
	}
	return result;
}
Mat method1_v(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowend = (Src1ch.rows / 3) - 1;
	for (int col = 0; col < result.cols; col++) {
		for (int row = 0; row < 3; row++) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row + 3, col);
		}
	}
	for (int col = 0; col < result.cols; col++) {
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row - 3, col);
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + 3, col)) * 0.5) + ((Src1ch.at<uchar>(row - 3, col)) * 0.5);
		}
	}
	return result;
}
Mat method2_h(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowk, colk;
	int colend = (Src1ch.cols / 3) - 1;
	for (int row = 0; row < result.rows; row++) {
		for (int col = 0; col < 3; col++) {
			if (row % 3 == 0) {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row, col + 3 - (col % 3)) * 0.5) + (Src1ch.at<uchar>(row + 1, col + 3 - (col % 3)) * 0.3333) + (Src1ch.at<uchar>(row + 2, col + 3 - (col % 3)) * 0.1667);
			}
			else if (row % 3 == 1) {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row - 1, col + 3 - (col % 3)) * 0.2) + (Src1ch.at<uchar>(row, col + 3 - (col % 3)) * 0.6) + (Src1ch.at<uchar>(row - 1, col + 3 - (col % 3)) * 0.2);
			}
			else {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row - 2, col + 3 - (col % 3)) * 0.1667) + (Src1ch.at<uchar>(row - 1, col + 3 - (col % 3)) * 0.3333) + (Src1ch.at<uchar>(row, col + 3 - (col % 3)) * 0.5);
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = result.cols - 1; col > result.cols - 4; col--) {

			if (row % 3 == 0) {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row, col - 1 - (col % 3)) * 0.5) + (Src1ch.at<uchar>(row + 1, col - 1 - (col % 3)) * 0.3333) + (Src1ch.at<uchar>(row + 2, col - 1 - (col % 3)) * 0.1667);
			}
			else if (row % 3 == 1) {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row - 1, col - 1 - (col % 3)) * 0.2) + (Src1ch.at<uchar>(row, col - 1 - (col % 3)) * 0.6) + (Src1ch.at<uchar>(row - 1, col - 1 - (col % 3)) * 0.2);
			}
			else {
				result.at<uchar>(row, col) = (Src1ch.at<uchar>(row - 2, col - 1 - (col % 3)) * 0.1667) + (Src1ch.at<uchar>(row - 1, col - 1 - (col % 3)) * 0.3333) + (Src1ch.at<uchar>(row, col - 1 - (col % 3)) * 0.5);
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			rowk = row % 3;
			colk = col % 3;
			if (rowk == 0 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col - 1)) * 0.375) + ((Src1ch.at<uchar>(row, col + 3)) * 0.125) + ((Src1ch.at<uchar>(row + 1, col - 1)) * 0.25) + ((Src1ch.at<uchar>(row + 1, col + 3)) * 0.0833) + ((Src1ch.at<uchar>(row + 2, col - 1)) * 0.125) + ((Src1ch.at<uchar>(row + 2, col + 3)) * 0.0416);
			}
			else if (rowk == 1 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col - 1)) * 0.15) + ((Src1ch.at<uchar>(row - 1, col + 3)) * 0.05) + ((Src1ch.at<uchar>(row, col - 1)) * 0.45) + ((Src1ch.at<uchar>(row, col + 3)) * 0.15) + ((Src1ch.at<uchar>(row + 1, col - 1)) * 0.15) + ((Src1ch.at<uchar>(row + 1, col + 3)) * 0.05);
			}
			else if (rowk == 2 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col - 1)) * 0.125) + ((Src1ch.at<uchar>(row - 2, col + 3)) * 0.0416) + ((Src1ch.at<uchar>(row - 1, col - 1)) * 0.25) + ((Src1ch.at<uchar>(row - 1, col + 3)) * 0.0833) + ((Src1ch.at<uchar>(row, col - 1)) * 0.375) + ((Src1ch.at<uchar>(row, col + 3)) * 0.125);
			}
			else if (rowk == 0 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col - 2)) * 0.25) + ((Src1ch.at<uchar>(row, col + 2)) * 0.25) + ((Src1ch.at<uchar>(row + 1, col - 2)) * 0.1666) + ((Src1ch.at<uchar>(row + 1, col + 2)) * 0.1666) + ((Src1ch.at<uchar>(row + 2, col - 2)) * 0.0833) + ((Src1ch.at<uchar>(row + 2, col + 2)) * 0.0833);
			}
			else if (rowk == 1 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col - 2)) * 0.1) + ((Src1ch.at<uchar>(row - 1, col + 2)) * 0.1) + ((Src1ch.at<uchar>(row, col - 2)) * 0.3) + ((Src1ch.at<uchar>(row, col + 2)) * 0.3) + ((Src1ch.at<uchar>(row + 1, col - 2)) * 0.1) + ((Src1ch.at<uchar>(row + 1, col + 2)) * 0.1);
			}
			else if (rowk == 2 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col - 2)) * 0.0833) + ((Src1ch.at<uchar>(row - 2, col + 2)) * 0.0833) + ((Src1ch.at<uchar>(row - 1, col - 2)) * 0.1666) + ((Src1ch.at<uchar>(row - 1, col + 2)) * 0.1666) + ((Src1ch.at<uchar>(row, col - 2)) * 0.25) + ((Src1ch.at<uchar>(row, col + 2)) * 0.25);
			}
			else if (rowk == 0 && colk == 2) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col - 3)) * 0.125) + ((Src1ch.at<uchar>(row, col + 1)) * 0.375) + ((Src1ch.at<uchar>(row + 1, col - 3)) * 0.0833) + ((Src1ch.at<uchar>(row + 1, col + 1)) * 0.25) + ((Src1ch.at<uchar>(row + 2, col - 3)) * 0.04166) + ((Src1ch.at<uchar>(row + 2, col + 1)) * 0.125);
			}
			else if (rowk == 1 && colk == 2) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col - 3)) * 0.05) + ((Src1ch.at<uchar>(row - 1, col + 1)) * 0.15) + ((Src1ch.at<uchar>(row, col - 3)) * 0.15) + ((Src1ch.at<uchar>(row, col + 1)) * 0.45) + ((Src1ch.at<uchar>(row + 1, col - 3)) * 0.05) + ((Src1ch.at<uchar>(row + 1, col + 1)) * 0.15);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col - 3)) * 0.04166) + ((Src1ch.at<uchar>(row - 2, col + 1)) * 0.125) + ((Src1ch.at<uchar>(row - 1, col - 3)) * 0.0833) + ((Src1ch.at<uchar>(row - 1, col + 1)) * 0.25) + ((Src1ch.at<uchar>(row, col - 3)) * 0.125) + ((Src1ch.at<uchar>(row, col + 1)) * 0.375);
			}
		}
	}
	return result;
}
Mat method2_v(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowk, colk;
	int rowend = (Src1ch.rows / 3) - 1;
	for (int col = 0; col < result.cols; col++) {
		for (int row = 0; row < 3; row++) {
			if (col % 3 == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + (3 - (row % 3)), col)) * 0.5) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col + 1)) * 0.3333) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col + 2)) * 0.1667);
			}
			else if (col % 3 == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + (3 - (row % 3)), col - 1)) * 0.2) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col)) * 0.6) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col + 1)) * 0.2);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + (3 - (row % 3)), col - 2)) * 0.1667) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col - 1)) * 0.3333) + ((Src1ch.at<uchar>(row + (3 - (row % 3)), col)) * 0.5);
			}
		}
	}
	for (int col = 0; col < result.cols; col++) {
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			if (col % 3 == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1 - (row % 3), col)) * 0.5) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col + 1)) * 0.3333) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col + 2)) * 0.1667);
			}
			else if (col % 3 == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1 - (row % 3), col - 1)) * 0.2) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col)) * 0.6) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col + 1)) * 0.2);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1 - (row % 3), col - 2)) * 0.1667) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col - 1)) * 0.3333) + ((Src1ch.at<uchar>(row - 1 - (row % 3), col)) * 0.5);
			}
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			rowk = row % 3;
			colk = col % 3;
			if (rowk == 0 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col)) * 0.375) + ((Src1ch.at<uchar>(row - 1, col + 1)) * 0.25) + ((Src1ch.at<uchar>(row - 1, col + 2)) * 0.125) + ((Src1ch.at<uchar>(row + 3, col)) * 0.125) + ((Src1ch.at<uchar>(row + 3, col + 1)) * 0.0833) + ((Src1ch.at<uchar>(row + 3, col + 2)) * 0.0416);
			}
			else if (rowk == 1 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col)) * 0.25) + ((Src1ch.at<uchar>(row - 2, col + 1)) * 0.1666) + ((Src1ch.at<uchar>(row - 2, col + 2)) * 0.0833) + ((Src1ch.at<uchar>(row + 2, col)) * 0.25) + ((Src1ch.at<uchar>(row + 2, col + 1)) * 0.1666) + ((Src1ch.at<uchar>(row + 2, col + 2)) * 0.0833);
			}
			else if (rowk == 2 && colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 3, col)) * 0.125) + ((Src1ch.at<uchar>(row - 3, col + 1)) * 0.0833) + ((Src1ch.at<uchar>(row - 3, col + 2)) * 0.0416) + ((Src1ch.at<uchar>(row + 1, col)) * 0.375) + ((Src1ch.at<uchar>(row + 1, col + 1)) * 0.25) + ((Src1ch.at<uchar>(row + 1, col + 2)) * 0.125);
			}
			else if (rowk == 0 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col - 1)) * 0.15) + ((Src1ch.at<uchar>(row - 1, col)) * 0.45) + ((Src1ch.at<uchar>(row - 1, col + 1)) * 0.15) + ((Src1ch.at<uchar>(row + 3, col - 1)) * 0.05) + ((Src1ch.at<uchar>(row + 3, col)) * 0.15) + ((Src1ch.at<uchar>(row + 3, col + 1)) * 0.05);
			}
			else if (rowk == 1 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col - 1)) * 0.1) + ((Src1ch.at<uchar>(row - 2, col)) * 0.3) + ((Src1ch.at<uchar>(row - 2, col + 1)) * 0.1) + ((Src1ch.at<uchar>(row + 2, col - 1)) * 0.1) + ((Src1ch.at<uchar>(row + 2, col)) * 0.3) + ((Src1ch.at<uchar>(row + 2, col + 1)) * 0.1);
			}
			else if (rowk == 2 && colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 3, col - 1)) * 0.05) + ((Src1ch.at<uchar>(row - 3, col)) * 0.15) + ((Src1ch.at<uchar>(row - 3, col + 1)) * 0.05) + ((Src1ch.at<uchar>(row + 1, col - 1)) * 0.15) + ((Src1ch.at<uchar>(row + 1, col)) * 0.45) + ((Src1ch.at<uchar>(row + 1, col + 1)) * 0.15);
			}
			else if (rowk == 0 && colk == 2) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col - 2)) * 0.125) + ((Src1ch.at<uchar>(row - 1, col - 1)) * 0.25) + ((Src1ch.at<uchar>(row - 1, col)) * 0.375) + ((Src1ch.at<uchar>(row + 3, col - 2)) * 0.0416) + ((Src1ch.at<uchar>(row + 3, col - 1)) * 0.0833) + ((Src1ch.at<uchar>(row + 3, col)) * 0.125);
			}
			else if (rowk == 1 && colk == 2) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col - 2)) * 0.0833) + ((Src1ch.at<uchar>(row - 2, col - 1)) * 0.1666) + ((Src1ch.at<uchar>(row - 2, col)) * 0.25) + ((Src1ch.at<uchar>(row + 2, col - 2)) * 0.0833) + ((Src1ch.at<uchar>(row + 2, col - 1)) * 0.1666) + ((Src1ch.at<uchar>(row + 2, col)) * 0.25);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 3, col - 2)) * 0.04166) + ((Src1ch.at<uchar>(row - 3, col - 1)) * 0.0833) + ((Src1ch.at<uchar>(row - 3, col)) * 0.125) + ((Src1ch.at<uchar>(row + 1, col - 2)) * 0.125) + ((Src1ch.at<uchar>(row + 1, col - 1)) * 0.25) + ((Src1ch.at<uchar>(row + 1, col)) * 0.375);
			}
		}
	}
	return result;
}
Mat method3_h(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int colk;
	int colend = (Src1ch.cols / 3) - 1;
	for (int row = 0; row < result.rows; row++) {
		for (int col = 0; col < 3; col++) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row, col + 3 - (col % 3));
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = result.cols - 1; col > result.cols - 4; col--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row, col - 1 - (col % 3));
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			colk = col % 3;
			if (colk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col - 1)) * 0.66667) + ((Src1ch.at<uchar>(row, col + 3)) * 0.33333);
			}
			else if (colk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col - 2)) * 0.5) + ((Src1ch.at<uchar>(row, col + 2)) * 0.5);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row, col + 1)) * 0.66667) + ((Src1ch.at<uchar>(row, col - 3)) * 0.33333);
			}

		}
	}
	return result;
}
Mat method3_v(Mat Src1ch) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowk;
	int rowend = (Src1ch.rows / 3) - 1;
	for (int col = 0; col < result.cols; col++) {
		for (int row = 0; row < 3; row++) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row + 3 - (row % 3), col);
		}
	}
	for (int col = 0; col < result.cols; col++) {
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			result.at<uchar>(row, col) = Src1ch.at<uchar>(row - 1 - (row % 3), col);
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			rowk = row % 3;
			if (rowk == 0) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 1, col)) * 0.66667) + ((Src1ch.at<uchar>(row + 3, col)) * 0.33333);
			}
			else if (rowk == 1) {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row - 2, col)) * 0.5) + ((Src1ch.at<uchar>(row + 2, col)) * 0.5);
			}
			else {
				result.at<uchar>(row, col) = ((Src1ch.at<uchar>(row + 1, col)) * 0.66667) + ((Src1ch.at<uchar>(row - 3, col)) * 0.33333);
			}
		}
	}
	return result;
}
void weight_generator(vector<vector<double>>& weights, vector<vector<double>>& distance, double k) {
	double sum = 0;
	double factor = 1;
	for (int i = 0; i < weights.size(); i++) {
		sum = 0;
		for (int j = 0; j < weights[0].size(); j++) {
			sum += pow(distance[i][j], k);
		}
		factor = 1.00 / sum;
		for (int j = 0; j < weights[0].size(); j++) {
			weights[i][j] = factor * pow(distance[i][j], k);
		}
	}
	return;
}
Mat method4_h(Mat Src1ch, float k) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int colend = (Src1ch.cols / 3) - 1;
	int rowk;
	int colk;
	int left_startrow, left_startcol, right_startrow, right_startcol;
	int l, r, w;
	vector<vector<double>> weights_h(9, vector<double>(18, 0));
	vector<vector<double>> dist_h(9, vector<double>(18, 0));
	distance_generator(dist_h, 1);
	weight_generator(weights_h, dist_h, k);
	for (int row = 0; row < result.rows; row++) {//LEFT-MOST PIXELS
		for (int col = 0; col < 3; col++) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col + 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col + 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {//RIGHT-MOST PIXELS
		for (int col = result.cols - 1; col > result.cols - 4; col--) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col - 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col - 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	for (int row = 0; row < result.rows; row++) {
		for (int col = 3; col < result.cols - 3; col++) {
			rowk = row % 3;
			colk = col % 3;
			left_startrow = row - rowk;
			left_startcol = col - 3 - colk;
			right_startrow = row - rowk;
			right_startcol = col + 3 - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(left_startrow + i, left_startcol + j)) * (weights_h[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(right_startrow + i, right_startcol + j)) * (weights_h[w][r++]));
				}
			}
		}
	}
	return result;
}
Mat method4_v(Mat Src1ch, float k) {
	Mat result = Mat::zeros(Src1ch.rows, Src1ch.cols, CV_32F);
	int rowend = (Src1ch.rows / 3) - 1;
	int rowk;
	int colk;
	int top_startrow, top_startcol, bottom_startrow, bottom_startcol;
	int l, r, w;
	vector<vector<double>> weights_v(9, vector<double>(18, 0));
	vector<vector<double>> dist_v(9, vector<double>(18, 0));
	distance_generator(dist_v, 0);
	weight_generator(weights_v, dist_v, k);
	for (int col = 0; col < result.cols; col++) {//TOP-MOST PIXELS
		for (int row = 0; row < 3; row++) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row + 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row + 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	for (int col = 0; col < result.cols; col++) {//BOTTOM-MOST PIXELS
		for (int row = result.rows - 1; row > result.rows - 4; row--) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row - 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row - 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	for (int row = 3; row < result.rows - 3; row++) {
		for (int col = 0; col < result.cols; col++) {
			rowk = row % 3;
			colk = col % 3;
			top_startrow = row - 3 - rowk;
			top_startcol = col - colk;
			bottom_startrow = row + 3 - rowk;
			bottom_startcol = col - colk;
			l = 0;
			r = 9;
			w = (rowk * 3) + colk;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(top_startrow + i, top_startcol + j)) * (weights_v[w][l++]));
					result.at<uchar>(row, col) += ((Src1ch.at<uchar>(bottom_startrow + i, bottom_startcol + j)) * (weights_v[w][r++]));
				}
			}
		}
	}
	return result;
}