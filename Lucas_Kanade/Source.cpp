// Source.cpp : Defines the entry point for the console application.
//

#include <afxwin.h>  // necessary for MFC to work properly
#include "Header.h"
#include "../../src/blepo.h"
#include <stack>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

using namespace blepo;


float compute(ImgFloat imgFile_temp, float x, float y)
{
	float x0, y0, ax, ay;
	x0 = floor(x); y0 = floor(y); ax = x - x0; ay = y - y0;
	if (x0 + 1<imgFile_temp.Width() && y0 + 1<imgFile_temp.Height())
		return((1 - ax)*(1 - ay)*imgFile_temp(x0, y0) + ax*(1 - ay)*imgFile_temp(x0 + 1, y0) + (1 - ax)*ay*imgFile_temp(x0, y0 + 1) + ax*ay*imgFile_temp(x0 + 1, y0 + 1));
	else if (y0 + 1 >= imgFile_temp.Height() && x0 + 1 >= imgFile_temp.Width())
		return((1 - ax)*(1 - ay)*imgFile_temp(x0, y0));
	else if (x0 + 1 >= imgFile_temp.Width())
		return((1 - ax)*(1 - ay)*imgFile_temp(x0, y0) + (1 - ax)*ay*imgFile_temp(x0, y0 + 1));
	else if (y0 + 1 >= imgFile_temp.Height())
		return((1 - ax)*(1 - ay)*imgFile_temp(x0, y0) + ax*(1 - ay)*imgFile_temp(x0 + 1, y0));

	else
		return 0;

}



float* computeMatrix(ImgFloat Gx, ImgFloat Gy, int i, int j, int size)
{
	float *temp_Matrix;
	int s_temp = floor(size / 2);
	temp_Matrix = new float[4];
	int x, y;
	float c_1, c_2;

	for (x = 0; x < 4; x++)
		temp_Matrix[x] = 0;
	for (x = i - s_temp; x <= i + s_temp; x++)
	{
		for (y = j - s_temp; y <= j + s_temp; y++)
		{
			if (x<0 || y<0 || x >= Gx.Width() || y >= Gy.Height())
			{
				c_1 = 0; c_2 = 0;
			}
			else
			{
				c_1 = Gx(x, y); c_2 = Gy(x, y);
			}
			temp_Matrix[0] = temp_Matrix[0] + c_1*c_1;
			temp_Matrix[1] = temp_Matrix[1] + c_1*c_2;
			temp_Matrix[2] = temp_Matrix[2];
			temp_Matrix[3] = temp_Matrix[3] + c_2*c_2;
		}
	}
	return temp_Matrix;
	delete[] temp_Matrix;
}


void featureTrack(ImgGray imgFile_temp, int t, double sigma, ImgFloat *imgFile_float)
{
	int i, j;
	ImgFloat Gx, Gy, imgFile_temp1, imgFile_output;
	float sum_1, sum_2;
	Convert(imgFile_temp, &imgFile_temp1);
	(*imgFile_float).Reset(imgFile_temp.Width(), imgFile_temp.Height());
	imgFile_output.Reset(imgFile_temp.Width(), imgFile_temp.Height());
	Gradient(imgFile_temp1, sigma, &Gx, &Gy);
	float *temp_Matrix = new float[4];
	for (i = 0; i<imgFile_temp.Width(); i++)
	{
		for (j = 0; j<imgFile_temp.Height(); j++)
		{
			temp_Matrix = computeMatrix(Gx, Gy, i, j, 3);
			sum_1 = (temp_Matrix[0] + temp_Matrix[3] + sqrt((temp_Matrix[0] - temp_Matrix[3])*(temp_Matrix[0] - temp_Matrix[3]) + 4 * temp_Matrix[2] * temp_Matrix[1])) / 2;
			sum_2 = (temp_Matrix[0] + temp_Matrix[3] - sqrt((temp_Matrix[0] - temp_Matrix[3])*(temp_Matrix[0] - temp_Matrix[3]) + 4 * temp_Matrix[2] * temp_Matrix[1])) / 2;
			if (sum_1<sum_2)
				sum_2 = sum_1;
			if (sum_2 >= t)
				imgFile_output(i, j) = sum_2;
			else
				imgFile_output(i, j) = 0;
		}
	}
	NonmaxSuppress4(imgFile_output, &(*imgFile_float));
}


float* calculateError(ImgFloat imgFile_temp1, ImgFloat imgFile_temp2, ImgFloat Gx, ImgFloat Gy, float x, float y, int size, float u, float v)
{
	float i, j;
	int s_temp = floor(size / 2);
	float *temp_Matrix;
	temp_Matrix = new float[2];

	for (int k = 0; k < 2; k++)
		temp_Matrix[k] = 0;
	for (i = x - 1; i <= x + s_temp; i++)
	{
		for (j = y - s_temp; j <= y + s_temp; j++)
		{
			if (i<0 || j<0 || i >= Gx.Width() || j >= Gy.Height() || i + u<0 || i + u >= Gy.Width() || j + v<0 || j + v >= Gy.Height())
			{
				temp_Matrix[1] = temp_Matrix[1] + 0;
				temp_Matrix[0] = temp_Matrix[0] + 0;
			}
			else
			{
				temp_Matrix[0] = temp_Matrix[0] + compute(Gx, i, j)*(compute(imgFile_temp1, i, j) - compute(imgFile_temp2, i + u, j + v));
				temp_Matrix[1] = temp_Matrix[1] + compute(Gy, i, j)*(compute(imgFile_temp1, i, j) - compute(imgFile_temp2, i + u, j + v));
			}
		}
	}
	return temp_Matrix;
	delete[] temp_Matrix;
}

void solve(float* temp_Matrix, float* e, float *udel, float *vdel)
{
	float temp;
	temp = temp_Matrix[0] * temp_Matrix[3] - temp_Matrix[1] * temp_Matrix[2];
	float ud;
	ud = (1 / temp)*(temp_Matrix[3] * e[0] - temp_Matrix[1] * e[1]);
	*udel = (ud);
	*vdel = ((1 / temp)*(temp_Matrix[0] * e[1] - temp_Matrix[1] * e[0]));

}

void LKAlgorithm(ImgGray imgFile_g1, ImgGray imgFile_g2, ImgFloat *temp, int sigma, const int size, ImgBgr *out)
{
	ImgFloat imgFile_f1, imgFile_f2;
	float i, j;
	ImgFloat Gx, Gy, imgFile_temp1, imgFile_temp2;
	imgFile_f1.Reset(imgFile_g1.Width(), imgFile_g1.Height()); imgFile_f2.Reset(imgFile_g1.Width(), imgFile_g1.Height());
	for (i = 0; i<imgFile_g1.Width(); i++)
		for (j = 0; j<imgFile_g1.Height(); j++)
		{
			imgFile_f1(i, j) = 0; imgFile_f2(i, j) = 0;
		}
	Convert(imgFile_g1, &imgFile_temp1);
	Convert(imgFile_g2, &imgFile_temp2);
	Point p; int count;
	Gradient(imgFile_temp1, 0.7, &Gx, &Gy);
	float udel = 0; float vdel = 0;
	float *temp_Matrix = new float[4];
	float *errorVar = new float[2];

	for (i = 0; i<imgFile_g1.Width(); i++)
	{
		for (j = 0; j<imgFile_g1.Height(); j++)
		{
			if ((*temp)(i, j) != 0)
			{
				(*temp)(i, j) = 0; temp_Matrix = computeMatrix(Gx, Gy, i, j, size);
				count = 0;
				while (count<20)
				{
					errorVar = calculateError(imgFile_temp1, imgFile_temp2, Gx, Gy, i, j, size, imgFile_f1(i, j), imgFile_f2(i, j));
					solve(temp_Matrix, errorVar, &udel, &vdel);
					imgFile_f1(i, j) = imgFile_f1(i, j) + udel; imgFile_f2(i, j) = imgFile_f2(i, j) + vdel;
					count++;
				}
				if ((i + imgFile_f1(i, j) >= 0) && (i + imgFile_f1(i, j)<imgFile_g1.Width()) && (j + imgFile_f2(i, j)) >= 0 && (j + imgFile_f2(i, j))<imgFile_g1.Height())
				{
					(*temp)(i + floor(imgFile_f1(i, j)), j + floor(imgFile_f2(i, j))) = 1; (*out)(i + floor(imgFile_f1(i, j)), j + floor(imgFile_f2(i, j))).r = 255;
					(*out)(i + floor(imgFile_f1(i, j)), j + floor(imgFile_f2(i, j))).b = 0; (*out)(i + floor(imgFile_f1(i, j)), j + floor(imgFile_f2(i, j))).g = 0;
				}
			}
		}
	}
}

int main(int argc, const char* argv[], const char* envp[])
{
	// Initialize MFC and return if failure
	HMODULE hModule = ::GetModuleHandle(NULL);
	if (hModule == NULL || !AfxWinInit(hModule, NULL, ::GetCommandLine(), 0))
	{
		printf("Fatal Error: MFC initialization failed (hModule = %x)\n", (int)hModule);
		return 1;
	}

	try
	{
		// Check for command line parameters
		if (argc < 4)
		{
			printf("Insufficient command line parameters!\n");
			exit(0);
		}

		double sigma;
		sigma = atof(argv[4]);

		int frameNumber = atoi(argv[2]);

		ImgGray imgFile_1;
		ImgBgr imgFile_Original;

		CString fString;

		CString temp = argv[1];
		fString.Format(temp, frameNumber);

		CString path("../../images/");
		CString imgPath = path + fString;

		Load(imgPath, &imgFile_1);
		Load(imgPath, &imgFile_Original);

		// Check for proper file size
		if (imgFile_Original.Width() == 0 || imgFile_Original.Height() == 0)
		{
			printf("Unable to load input file. Improper file size!\n");
			exit(0);
		}

		Figure fig1;
		fig1.SetTitle("Input Image");
		fig1.Draw(imgFile_Original);


		ImgFloat corner;
		featureTrack(imgFile_1, 900, sigma, &corner);

		int i, j;
		for (i = 0; i<corner.Width(); i++)
		{
			for (j = 0; j<corner.Height(); j++)
			{
				if (corner(i, j) != 0)
				{
					imgFile_Original(i, j).r = 255; imgFile_Original(i, j).b = 0; imgFile_Original(i, j).g = 0;
				}
			}
		}

		Figure fig2;
		fig2.SetTitle("Features");
		fig2.Draw(imgFile_Original);

		int last = atoi(argv[3]);
		const int size = atoi(argv[5]);

		j = 0;
		vector<ImgBgr> imgVec(30);
		for (i = frameNumber; i<last; i++)
		{
			CString frame1, frame2;
			CString t;

			t.Format(temp, i);
			frame1 = path + t;
			t.Format(temp, i + 1);
			frame2 = path + t;

			ImgGray imgFile_g1, imgFile_g2;
			Load(frame1, &imgFile_g1); Load(frame1, &imgVec[j]);
			Load(frame2, &imgFile_g2);
			LKAlgorithm(imgFile_g1, imgFile_g2, &corner, sigma, size, &imgVec[j]);
			j++;
		}

		i = 0; j = 0;
		Figure fig3;
		while (i == 0)
		{
			Sleep(1000);
			fig3.Draw(imgVec[j]);
			if (j == 29)
				j = 0;
			else
				j++;
		}
		Figure fig4;
		fig4.Draw(imgFile_Original);

		// Loop forever until user presses Ctrl-C in terminal window.
		EventLoop();
	}

	catch (const Exception& e)
	{
		e.Display();    // display exception 
	}
	return 0;
}
