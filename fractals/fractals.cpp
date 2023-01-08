#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <string>
#include <math.h>

#include <chrono>
using namespace std::chrono;

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

//float X_MIN = -2;
//float X_MAX = 0.47;
//float Y_MIN = -1.12;
//float Y_MAX = 1.12;

//float X_MIN = -0.6884971875 -.0001;
//float X_MAX = -0.6884971875 + .0001;
//float Y_MIN = -0.27988465625 - .0001;
//float Y_MAX = -0.27988465625 + .0001;
//#define RADIUS 10

//float X_MIN = -1.65;
//float X_MAX = 1.65;
//float Y_MIN = -0.9539;
//float Y_MAX = 0.9539;

const double pi = std::acos(-1);
const std::complex<double> I(0, 1);

//#define WIDTH 2200
//#define HEIGHT 2000

//#define WIDTH 2000
//#define HEIGHT 1500

#define WIDTH 2000
#define HEIGHT 2000

// square
/*
float X_MIN = -1.25 - 0.001;
float X_MAX = -1.25+0.001;
float Y_MIN = -0.035-0.001;
float Y_MAX = -0.035+0.001;
*/


// seahorse tail
//double X_MIN = -0.7435669 - 0.0022878;
//double X_MAX = -0.7435669 + 0.0022878;
//double Y_MIN = -0.1314023 - 0.0022878;
//double Y_MAX = -0.1314023 + 0.0022878;

/*
float X_MIN = -0.7435669 - 0.00022878;
float X_MAX = -0.7435669 + 0.00022878;
float Y_MIN = -0.1314023 - 0.00022878;
float Y_MAX = -0.1314023 + 0.00022878;
*/

float X_MIN = -2;
float X_MAX = 2;
float Y_MIN = -2;
float Y_MAX = 2;

//float X_MIN = -2;
//float X_MAX = 2;
//float Y_MIN = -1.5;
//float Y_MAX = 1.5;

//float X_MIN = -1;
//float X_MAX = 1;
//float Y_MIN = -0.75;
//float Y_MAX = 0.75;

//#define WIDTH 320
//#define HEIGHT 185

//int radius = pow(2, 8);
//double radius = DBL_MAX;
double radius = 2;

int maxIterations = 1000;
int kIterations = 1;
int k = 2;

const double epsilon = pow(10, -6);

using namespace std;

#define BLACK 0
#define WHITE 1
#define RED 2
#define BLUE 3

void spectral_color(double& r, double& g, double& b, double l) // RGB <0,1> <- lambda l <400,700> [nm]
{
    double t;  r = 0.0; g = 0.0; b = 0.0;
    if ((l >= 400.0) && (l < 410.0)) { t = (l - 400.0) / (410.0 - 400.0); r = +(0.33 * t) - (0.20 * t * t); }
    else if ((l >= 410.0) && (l < 475.0)) { t = (l - 410.0) / (475.0 - 410.0); r = 0.14 - (0.13 * t * t); }
    else if ((l >= 545.0) && (l < 595.0)) { t = (l - 545.0) / (595.0 - 545.0); r = +(1.98 * t) - (t * t); }
    else if ((l >= 595.0) && (l < 650.0)) { t = (l - 595.0) / (650.0 - 595.0); r = 0.98 + (0.06 * t) - (0.40 * t * t); }
    else if ((l >= 650.0) && (l < 700.0)) { t = (l - 650.0) / (700.0 - 650.0); r = 0.65 - (0.84 * t) + (0.20 * t * t); }
    if ((l >= 415.0) && (l < 475.0)) { t = (l - 415.0) / (475.0 - 415.0); g = +(0.80 * t * t); }
    else if ((l >= 475.0) && (l < 590.0)) { t = (l - 475.0) / (590.0 - 475.0); g = 0.8 + (0.76 * t) - (0.80 * t * t); }
    else if ((l >= 585.0) && (l < 639.0)) { t = (l - 585.0) / (639.0 - 585.0); g = 0.84 - (0.84 * t); }
    if ((l >= 400.0) && (l < 475.0)) { t = (l - 400.0) / (475.0 - 400.0); b = +(2.20 * t) - (1.50 * t * t); }
    else if ((l >= 475.0) && (l < 560.0)) { t = (l - 475.0) / (560.0 - 475.0); b = 0.7 - (t)+(0.30 * t * t); }
}

double lerp(double t, double a, double b)
{
    return a + t * (b - a);
}

void ColorMap(int iteration, double& r, double& g, double& b)
{
    double t = (double)iteration / maxIterations;
    if (t >= 0 && t < 0.16)
    {
        double r1 = 0;
        double g1 = 7.0 / 255;
        double b1 = 100.0 / 255;
        double r2 = 32.0 / 255;
        double g2 = 107.0 / 255;
        double b2 = 203.0 / 255;

        //t = 0.16;
        r = r1 + (t/0.16) * (r2 - r1);
        g = g1 + (t/0.16) * (g2 - g1);
        b = b1 + (t/0.16) * (b2 - b1);
    }
    else if (t >= 0.16 && t < 0.42)
    {
        double r1 = 32.0 / 255;
        double g1 = 107.0 / 255;
        double b1 = 203.0 / 255;
        double r2 = 237.0 / 255;
        double g2 = 1;
        double b2 = 1;
        //t = 0.42;
        r = r1 + ((t - 0.16) / (0.42 - 0.16)) * (r2 - r1);
        g = g1 + ((t - 0.16) / (0.42 - 0.16)) * (g2 - g1);
        b = b1 + ((t - 0.16) / (0.42 - 0.16)) * (b2 - b1);
    }
    else if (t >= 0.42 && t < 0.6425)
    {
        double r1 = 237.0 / 255;
        double g1 = 1;
        double b1 = 1;
        double r2 = 1;
        double g2 = 170.0 / 255;
        double b2 = 0;

        //t = 0.6425;
        r = r1 + ((t - 0.42) / (0.6425 - 0.42)) * (r2 - r1);
        g = g1 + ((t - 0.42) / (0.6425 - 0.42)) * (g2 - g1);
        b = b1 + ((t - 0.42) / (0.6425 - 0.42)) * (b2 - b1);
    }
    else if (t >= 0.6425 && t < 0.8575)
    {
        double r1 = 1;
        double g1 = 170.0 / 255;
        double b1 = 0;
        double r2 = 0;
        double g2 = 2.0 / 255;
        double b2 = 0;
       // t = 0.8575;
        r = r1 + ((t - 0.6425) / (0.8575 - 0.6425)) * (r2 - r1);
        g = g1 + ((t - 0.6425) / (0.8575 - 0.6425)) * (g2 - g1);
        b = b1 + ((t - 0.6425) / (0.8575 - 0.6425)) * (b2 - b1);
    }
    else
    {
        double r1 = 0;
        double g1 = 2.0 / 255;
        double b1 = 0;
        double r2 = 0;
        double g2 = 7.0 / 255;
        double b2 = 100.0 / 255;
        //t = 1;
        r = r1 + ((t - 0.8575) / (1 - 0.8575)) * (r2 - r1);
        g = g1 + ((t - 0.8575) / (1 - 0.8575)) * (g2 - g1);
        b = b1 + ((t - 0.8575) / (1 - 0.8575)) * (b2 - b1);
    }

    return;
}

void PreComputeColors(vector<vector<double>>& colors)
{
    vector<double> color = { 0,0, 0 };
    for (int i = 0; i < maxIterations; i++)
    {
        ColorMap(i, color[0], color[1], color[2]);
        colors[i] = color;
    }
}

void HSVtoRGB(double& fR, double& fG, double& fB, double& fH, double& fS, double& fV) {
    double fC = fV * fS; // Chroma
    //double fHPrime = fmod(fH / 60.0, 6);
    double fX = fC * (1 - fabs(fmod(fH/60, 2) - 1));
    double fM = fV - fC;

    if (0 <= fH && fH < 60){
    //if (0 <= fHPrime && fHPrime < 1) {
        fR = fC;
        fG = fX;
        fB = 0;
    }
    else if (60 <= fH && fH < 120) {
    //else if (1 <= fHPrime && fHPrime < 2) {
        fR = fX;
        fG = fC;
        fB = 0;
    }
    else if (120 <= fH && fH < 180) {
    //else if (2 <= fHPrime && fHPrime < 3) {
        fR = 0;
        fG = fC;
        fB = fX;
    }
    else if (180 <= fH && fH < 240) {
    //else if (3 <= fHPrime && fHPrime < 4) {
        fR = 0;
        fG = fX;
        fB = fC;
    }
    else if (240 <= fH && fH < 300) {
    //else if (4 <= fHPrime && fHPrime < 5) {
        fR = fX;
        fG = 0;
        fB = fC;
    }
    else if (300 <= fH && fH < 360) {
    //else if (5 <= fHPrime && fHPrime < 6) {
        fR = fC;
        fG = 0;
        fB = fX;
    }
    else {
        fR = fC;
        fG = 0;
        fB = fX;
    }

    fR += fM;
    fG += fM;
    fB += fM;
}


double Length(complex<double> number)
{
    return number.imag() * number.imag() + number.real() * number.real();
}

vector<double> GetVecColor(int iteration)
{
    vector<double> color;
    vector<double> a{ 0.5, 0.5, 0.5 };
    vector<double> b{ 0.5, 0.5, 0.5 };
    vector<double> c{ 1, 1, 0.5 };
    vector<double> d{ 0.8, 0.9, 0.3 };
    //vector<double> t = { (double)iteration / maxIterations, (double)iteration / maxIterations, (double)iteration / maxIterations };
    double t = (double)iteration / maxIterations;
    c[0] *= t;
    c[1] *= t;
    c[2] *= t;
    c[0] += d[0];
    c[1] += d[1];
    c[2] += d[2];
    c[0] *= 2 * pi;
    c[1] *= 2 * pi;
    c[2] *= 2 * pi;
    b[0] *= cos(c[0]);
    b[1] *= cos(c[1]);
    b[2] *= cos(c[2]);
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
    //color = a + b * cos(6.28318 * (c + d));
    //if (t < 0.01)
    //    std::cout << t << std::endl;
    return a;
}

vector<double> GetColor(int iteration, complex<double> z0, complex<double> start, complex<double> p, complex<double> zn)
{
    vector<double> color;
    complex<double> z0n = z0;
    if (iteration > maxIterations)
    {
        int i;
        for (i = 0; i < 100000; i++)
        {
            //std::cout << "z0n: " << z0n.real() << " " << z0n.imag() << std::endl;
            complex<double> z_k = z0n;
            for (int j = 0; j < k; j++)
            {
                z_k *= z_k;
                z_k += p;
            }
            //std::cout << "z_k: " << z_k.real() << " " << z_k.imag() << std::endl;
            double len = abs(z_k - z0n);
            //std::cout << "len: " << len << std::endl;
            if (len < epsilon)
            {
                //std::cout << "episilon" << std::endl;
                iteration = i;
                break;
            }
            z0n *= z0n;
            z0n += p;
        }
        if (iteration % 2 == 0)
        {
            color = { 1, 0, 0 };
        }
        else
        {
            color = { 0, 0, 1 };
        }
        //color = GetVecColor(iteration);
        return color;
    }

    double smooth = iteration + 1 - log(log(abs(zn))) / log(2);
    smooth /= maxIterations;
    
    
    if (iteration == maxIterations)
        color = { 0, 0, 0 };
    else if (iteration % 2 == 0)
        color = { 1 ,0, 0 };
    else
        color = { 1, 1, 1 };
    //color = GetVecColor(iteration);
    //color = { smooth, smooth, smooth };
    //double c = 1 - iteration / maxIterations;
    //color = { c, c, c };
    return color;
}

int GetIterations(complex<double> num, complex<double> start)
{
    int iteration;
    for (iteration = 0; iteration < maxIterations; iteration++)
    {
        if (num.imag() * num.imag() + num.real() * num.real() > 4)
        {
            break;
        }
        num *= num;
        num += start;
    }
    return iteration;
}

int GetIterations(complex<double> num, double power, double radius)
{
    int iteration;
    complex<double> current(0, 0);
    for (iteration = 0; iteration < maxIterations; iteration++)
    {
        if (current.imag() * current.imag() + current.real() * current.real() > radius * radius)
        {
            break;
        }
        //current = pow(current, power);
        current *= current;
        current += num;
    }
    return iteration;
}

int Mandelbrot(complex<double> pixelC, complex<double>* zn)
{
    complex<double> z0 = { 0, 0 };
    complex<double> z_k = { 0, 0 };
    int iteration;
    for (iteration = 1; iteration < maxIterations; iteration++)
    {
        z0 *= z0;
        z0 += pixelC;

        /*z_k = z0;
        for (int j = 0; j < k; j++)
        {
            z_k *= z_k;
            z_k += pixelC;
        }

        if (abs(z_k - z0) < epsilon)
        {
            std::cout << "epsilon" << std::endl;
        }*/

        if (Length(z0) > radius * radius)
        //if (abs(z_k-z0) < epsilon)
        {
            break;
        }
    }
    *zn = z0;
    return iteration;
}

int JuliaExponent(complex<double> pixelC, double a)
{
    complex<double> z0 = pixelC;
    int iteration;
    for (iteration = 1; iteration < maxIterations; iteration++)
    {
        z0 *= z0;
        z0 += 0.7885 * exp(I*a);

        if (Length(z0) > radius * radius)
        {
            break;
        }
    }
    return iteration;
}

int Julia(complex<double> pixelC, complex<double> p)
{
    complex<double> z0 = pixelC;
    int iteration;
    for (iteration = 1; iteration < maxIterations; iteration++)
    {
        z0 *= z0;
        z0 += p;

        if (Length(z0) > radius * radius)
        {
            break;
        }
    }
    return iteration;
}

int main()
{
    vector<char> data;
    vector<vector<int>> nums;
    vector<vector<complex<double>>> zns;
    nums.resize(HEIGHT);
    zns.resize(HEIGHT);
    vector<int> NumIterationsPerPixel;
    NumIterationsPerPixel.resize(maxIterations+1);
    data.resize(WIDTH * HEIGHT * 3);
    int dataPos = 0;
    float x_pos = -1 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (1 - -1)));
    float y_pos = -1 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (1 - -1)));
    //complex<double> start = complex<double>(x_pos, y_pos);
    complex<double> start = { -0.77146, 0.10119 };
    int totalIterations = 0;

    vector<vector<double>> colors;
    colors.resize(maxIterations);
    PreComputeColors(colors);

    auto startTime = high_resolution_clock::now();
    for (int k = 0; k < kIterations; k++)
    {
        data.clear();
        data.resize(WIDTH * HEIGHT * 3);
        dataPos = 0;
        for (int i = 0; i < HEIGHT; i++)
        {
            std::cout << "Starting row: " << i << std::endl;
            nums[i].resize(WIDTH);
            zns[i].resize(WIDTH);
            for (int j = 0; j < WIDTH; j++)
            {
                complex<double> num = complex<double>(X_MIN + (double)j / WIDTH * (X_MAX - X_MIN), Y_MIN + (double)i / HEIGHT * (Y_MAX - Y_MIN));
                
                //int iterations = GetIterations(num, start);
                complex<double> zn = (0, 0);
                int iterations = Mandelbrot(num, &zn);
                //int iterations = Julia(num, complex<double>(0.281215625, -0.0113825));
                complex<double> p = { -0.77146, 0.10119 };
                //int iterations = Julia(num, p);
                //int iterations = JuliaExponent(num, (double)k/kIterations*2*pi);
                //vector<double> color = GetColor(iterations, num, start, p, zn);
                nums[i][j] = iterations;
                //nums[i][j] = j;
                zns[i][j] = zn;
                //NumIterationsPerPixel[iterations]++;
                //vector<double> colorVec = GetVecColor(iterations);
                //vector<double> color = GetColor(nums[i][j], num, start, p, zn);
                //data[dataPos++] = color[0]*255;
                //data[dataPos++] = color[1]*255;
                //data[dataPos++] = color[2]*255;
            }
        }
        for (int i = 0; i <= maxIterations; i++)
        {
            totalIterations += NumIterationsPerPixel[i];
        }
    

        for (int i = 0; i < HEIGHT; i++)
        {
            for (int j = 0; j < WIDTH; j++)
            {
                
                int iteration = nums[i][j];
                if (iteration == maxIterations)
                {
                    data[dataPos++] = 0;
                    data[dataPos++] = 0;
                    data[dataPos++] = 0;
                    continue;
                }
                /* historgram
                
                auto zn = zns[i][j];
                float log_zn = log(zn.real() * zn.real() + zn.imag() * zn.imag()) / 2;
                float nu = log(log_zn / log(2)) / log(2);
                //iteration = floor(iteration + 1 - nu);
                float tIter = 0;
                for (int k = 0; k <= iteration; k++)
                    tIter += (float)NumIterationsPerPixel[k] / totalIterations;
                */
                
                auto zn = zns[i][j];
                double log_zn = log(zn.real() * zn.real() + zn.imag() * zn.imag()) / 2;
                double nu = log(log_zn / log(2)) / log(2);
                double fIt = iteration + 1 - nu;
                double s = (double)iteration / maxIterations;
                //double s1 = (double)floor(fIt) / maxIterations;
                //double s2 = (double)(floor(fIt)+1) / maxIterations;

                //vector<double> rgbColor1 = { 0, 0, 0 };
                //vector<double> rgbColor2 = { 0, 0, 0 };
                //double h1 = fmod(powf(s1 * 360, 1.5), 360);
                //double h2 = fmod(powf(s2 * 360, 1.5), 360);
                //vector<double> hsvColor1 = { h1, 1, s1*100 };
                //vector<double> hsvColor2 = { h2, 1, s2* 100 };
                //HSVtoRGB(rgbColor1[0], rgbColor1[1], rgbColor1[2], hsvColor1[0], hsvColor1[1], hsvColor1[2]);
                //HSVtoRGB(rgbColor2[0], rgbColor2[1], rgbColor2[2], hsvColor2[0], hsvColor2[1], hsvColor2[2]);
                //double redLerp = rgbColor1[0] + fmod(fIt, 1) * (rgbColor2[0] - rgbColor1[0]);
                //double greenLerp = rgbColor1[1] + fmod(fIt, 1) * (rgbColor2[1] - rgbColor1[1]);
                //double blueLerp = rgbColor1[2] + fmod(fIt, 1) * (rgbColor2[2] - rgbColor1[2]);
                //vector<double> finalColor = { redLerp, greenLerp, blueLerp };

                vector<double> finalColor = { 0, 0, 0 };
                double l = 400 + s * (700 - 400);
                //spectral_color(finalColor[0], finalColor[1], finalColor[2], l);

                iteration = floor(fIt);
                vector<double> color1 = { 0, 0, 0 };
                vector<double> color2 = { 0, 0, 0 };
                //ColorMap(iteration, color1[0], color1[1], color1[2]);
                //ColorMap(iteration+1, color2[0], color2[1], color2[2]);
                color1 = colors[iteration];
                color2 = colors[iteration + 1];
                //fIt = floor(fIt);
                finalColor[0] = lerp(fmod(fIt, 1), color1[0], color2[0]);
                finalColor[1] = lerp(fmod(fIt, 1), color1[1], color2[1]);
                finalColor[2] = lerp(fmod(fIt, 1), color1[2], color2[2]);

                //tIter = fmod(pow((pow(((float)iteration / maxIterations), 2) * 255), 1.5), 255);
               // tIter = 0 + fmod(fiteration, 1) * (1 - 0);
                //data[dataPos++] = (char)(tIter * tIter * tIter * 255);
                //data[dataPos++] = (char)(tIter * tIter * tIter * 255);
                //data[dataPos++] = (char)((0.8 + 0.2*tIter * tIter * tIter) * 255);
                data[dataPos++] = (char)(finalColor[0]*255);
                data[dataPos++] = (char)(finalColor[1]*255);
                data[dataPos++] = (char)(finalColor[2]*255);
            }
        }
        string extension = ".png";
        string fileName = "pics/image";
        fileName.append(std::to_string(k).c_str());
        fileName.append(extension.c_str());
        stbi_write_png(fileName.c_str(), WIDTH, HEIGHT, 3, data.data(), WIDTH * 3);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - startTime);

    cout << duration.count() << endl;

}
