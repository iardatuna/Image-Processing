//I. Arda Tuna
//240201031
// ------------------------------
// Written by Mustafa Ozuysal
// Contact <mustafaozuysal@iyte.edu.tr> for comments and bug reports
// ------------------------------
// Copyright (c) 2018, Mustafa Ozuysal
// All rights reserved.
// ------------------------------
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the copyright holders nor the
//       names of his/its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// ------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ------------------------------
#include "image.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ostream>
#include <memory>

using std::cerr;
using std::clog;
using std::cos;
using std::min;
using std::endl;
using std::exp;
using std::ifstream;
using std::ios;
using std::memset;
using std::ofstream;
using std::sin;
using std::string;
using std::unique_ptr;
using std::cout;

namespace ceng391 {

Image::Image(int width, int height, int n_channels, int step)
{
        m_width = width;
        m_height = height;
        m_n_channels = n_channels;

        m_step = m_width*m_n_channels;
        if (m_step < step)
                m_step = step;
        m_data = new uchar[m_step*height];
}

Image::~Image()
{
        delete [] m_data;
}

Image* Image::new_gray(int width, int height)
{
        return new Image(width, height, 1);
}

Image* Image::new_rgb(int width, int height)
{
        return new Image(width, height, 3);
}

Image *Image::new_copy(Image *img)
{
        Image *cpy = new Image(img->w(), img->h(), img->n_ch());
        for (int y = 0; y < img->h(); ++y)
                memcpy(cpy->data(y), img->data(y), img->w() * img->n_ch());
        return cpy;
}

void Image::set_rect(int x, int y, int width, int height, uchar red, uchar green, uchar blue)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        if (m_n_channels == 1) {
                int value = 0.3*red + 0.59*green + 0.11*blue;
                if (value > 255)
                        value = 255;
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i] = value;
                        }
                }
        } else if (m_n_channels == 3) {
                for (int j = y; j < y+height; ++j) {
                        if (j >= m_height)
                                break;
                        uchar* row_data = data(j);
                        for (int i = x; i < x+width; ++i) {
                                if (i >= m_width)
                                        break;
                                row_data[i*3]     = red;
                                row_data[i*3 + 1] = green;
                                row_data[i*3 + 2] = blue;
                        }
                }
        }
}

void Image::set_rect(int x, int y, int width, int height, uchar value)
{
        if (x < 0) {
                width += x;
                x = 0;
        }

        if (y < 0) {
                height += y;
                y = 0;
        }

        for (int j = y; j < y+height; ++j) {
                if (j >= m_height)
                        break;
                uchar* row_data = data(j);
                for (int i = x; i < x+width; ++i) {
                        if (i >= m_width)
                                break;
                        for (int c = 0; c < m_n_channels; ++c)
                                row_data[i*m_n_channels + c] = value;
                }
        }
}

void Image::to_grayscale()
{
        if (m_n_channels == 1) {
                return;
        } else if (m_n_channels == 3) {
                int new_step = m_width;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar red = row_old[3*x];
                                uchar green = row_old[3*x + 1];
                                uchar blue = row_old[3*x + 2];
                                int value = 0.3*red + 0.59*green + 0.11*blue;
                                if (value > 255)
                                        value = 255;
                                row_new[x] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 1;
        }
}

void Image::to_rgb()
{
        if (m_n_channels == 3) {
                return;
        } else if (m_n_channels == 1) {
                int new_step = m_width * 3;
                uchar *new_data = new uchar[new_step * m_height];
                for (int y = 0; y < m_height; ++y) {
                        uchar *row_old = m_data + m_step * y;
                        uchar *row_new = new_data + new_step * y;
                        for (int x = 0; x < m_width; ++x) {
                                uchar value = row_old[x];
                                row_new[3*x]     = value;
                                row_new[3*x + 1] = value;
                                row_new[3*x + 2] = value;
                        }
                }

                delete [] m_data;
                m_data = new_data;
                m_step = new_step;
                m_n_channels = 3;
        }
}

bool Image::write_pnm(const std::string& filename) const
{
        ofstream fout;

        string magic_head = "P5";
        string extended_name = filename + ".pgm";
        if (m_n_channels == 3) {
                magic_head = "P6";
                extended_name = filename + ".ppm";
        }

        fout.open(extended_name.c_str(), ios::out | ios::binary);
        if (!fout.good()) {
                cerr << "Error opening file " << extended_name << " for output!" << endl;
                return false;
        }

        fout << magic_head << "\n";
        fout << m_width << " " << m_height << " 255\n";
        for (int y = 0; y < m_height; ++y) {
                const uchar *row_data = data(y);
                fout.write(reinterpret_cast<const char*>(row_data), m_width*m_n_channels*sizeof(uchar));
        }
        fout.close();

        return true;
}

bool Image::read_pnm(const std::string& filename)
{
        ifstream fin(filename.c_str(), ios::in | ios::binary);
        if (!fin.good()) {
                cerr << "Error opening PNM file " << filename << endl;
                return false;
        }

        int width;
        int height;
        int max_val;
        int n_channels = 1;
        string head = "00";
        head[0] = fin.get();
        head[1] = fin.get();
        if (head == "P5") {
                clog << "Loading PGM Binary" << endl;
                n_channels = 1;
        } else if (head == "P6") {
                clog << "Loading PPM Binary" << endl;
                n_channels = 3;
        } else {
                cerr << "File " << filename << " is not a Binary PGM or PPM!" << endl;
                return false;
        }

        fin >> width;
        fin >> height;
        fin >> max_val;
        if (fin.peek() == '\n')
                fin.get();

        int step = width * n_channels;
        uchar *new_data = new uchar[step*height];
        for (int y = 0; y < height; ++y) {
                fin.read(reinterpret_cast<char*>(new_data + y*step), step*sizeof(uchar));
        }
        fin.close();

        delete [] m_data;
        m_data = new_data;
        m_width = width;
        m_height = height;
        m_step = step;
        m_n_channels = n_channels;

        return true;
}

short *Image::deriv_x() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dx = new short[m_width * m_height];
        for (int y = 0; y < m_height; ++y) {
                const uchar *row = this->data(y);
                short *drow = dx + y * m_width;
                drow[0] = 0;
                for (int x = 1; x < m_width - 1; ++x) {
                        drow[x] = row[x + 1] - row[x - 1];
                }
                drow[m_width - 1] = 0;
        }

        return dx;
}

short *Image::deriv_y() const
{
        if (m_n_channels == 3) {
                cerr << "Image derivatives only implemented for grayscale images!" << endl;
                return nullptr;
        }

        short *dy = new short[m_width * m_height];

        memset(dy, 0, m_width * sizeof(*dy));
        for (int y = 1; y < m_height - 1; ++y) {
                const uchar *rowm = this->data(y - 1);
                const uchar *rowp = this->data(y + 1);
                short *drow = dy + y * m_width;
                for (int x = 0; x < m_width; ++x) {
                        drow[x] = rowp[x] - rowm[x];
                }
        }
        memset(dy + (m_height - 1) * m_width, 0, m_width * sizeof(*dy));

        return dy;
}

void Image::rotate(Image *rotated, double theta, double tx, double ty) const
{
        if (m_n_channels != 1) {
                cerr << "Rotate only works on grayscale images!" << endl;
                return;
        }
        rotated->to_grayscale();

        double ct = cos(theta);
        double st = sin(theta);
        double tx_inv = -ct * tx + st * ty;
        double ty_inv = -st * tx - ct * ty;

        int wp = rotated->w();
        int hp = rotated->h();

        for (int yp = 0; yp < hp; ++yp) {
                uchar *row_p = rotated->data(yp);
                for (int xp = 0; xp < wp; ++xp) {
                        double x = ct * xp - st * yp + tx_inv;
                        double y = st * xp + ct * yp + ty_inv;

                        int x0 = static_cast<int>(x);
                        int y0 = static_cast<int>(y);

                        int value = 0;
                        if (x0 < 0 || y0 < 0 || x0 >= m_width || y0 >= m_height) {
                                value = 0;
                        } else {
                                const uchar *row = this->data(y0);
                                value = row[x0];
                        }

                        row_p[xp] = value;
                }
        }
}

void Image::rotate_centered(Image *rotated, double theta) const
{
        double ct = cos(theta);
        double st = sin(theta);
        double hw = m_width / 2.0;
        double hh = m_height / 2.0;
        double hwp = rotated->w() / 2.0;
        double hhp = rotated->h() / 2.0;

        double tx_cap = -ct * hw - st * hh + hwp;
        double ty_cap =  st * hw - ct * hh + hhp;
        this->rotate(rotated, theta, tx_cap, ty_cap);
}

void Image::smooth_x(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_width + 2 * l]);

        for (int y = 0; y < m_height - 1; ++y) {
                copy_to_buffer(buffer.get(), this->data(y), m_width, l, 1);
                convolve_buffer(buffer.get(), m_width, kernel.get(), k);
                copy_from_buffer(this->data(y), buffer.get(), m_width, 1);
        }
}

void Image::smooth_y(float sigma)
{
        if (m_n_channels != 1) {
                cerr << "Smooth-x only works on grayscale images!" << endl;
                return;
        }

        int k = 0;
        unique_ptr<float []> kernel(gaussian_kernel(sigma, &k));

        int l = k / 2;
        unique_ptr<float []>  buffer(new float[m_height + 2 * l]);

        for (int x = 0; x < m_width - 1; ++x) {
                copy_to_buffer(buffer.get(), m_data + x, m_height, l, m_step);
                convolve_buffer(buffer.get(), m_height, kernel.get(), k);
                copy_from_buffer(m_data + x, buffer.get(), m_height, m_step);
        }
}

void Image::smooth(float sigma_x, float sigma_y)
{
        smooth_x(sigma_x);
        smooth_y(sigma_y);
}

std::vector<Keypoint> Image::harris_corners(float threshold, float k,
                                            float sigma)
{
        short *Ix = this->deriv_x();
        short *Iy = this->deriv_y();

        int n = m_width * m_height;
        short *Ix2  = vec_mul(n, Ix, Ix);
        short *Iy2  = vec_mul(n, Iy, Iy);
        short *IxIy = vec_mul(n, Ix, Iy);

        smooth_short_buffer(m_width, m_height, Ix2, sigma);
        smooth_short_buffer(m_width, m_height, Iy2, sigma);
        smooth_short_buffer(m_width, m_height, IxIy, sigma);

        const float *score = harris_corner_score(m_width, m_height,
                                           Ix2, Iy2, IxIy, k);

        std::vector<Keypoint> keys;
        int border = min(1, static_cast<int>(2.0f * sigma));
        for (int y = border; y < m_height - border; ++y) {
                const float *Rm = score + (y - 1) * m_width;
                const float *R  = score + y * m_width;
                const float *Rp = score + (y + 1) * m_width;
                for (int x = border; x < m_width - border; ++x) {
                        if (R[x] > threshold && R[x] > R[x - 1]
                            && R[x] > R[x + 1] && R[x] > Rm[x]
                            && R[x] > Rp[x]) {
                                Keypoint key;
                                key.x = x;
                                key.y = y;
                                key.score = R[x];
                                keys.push_back(key);
                        }
                }
        }

        delete [] score;
        delete [] IxIy;
        delete [] Iy2;
        delete [] Ix2;
        delete [] Iy;
        delete [] Ix;

        return keys;
}
std::vector<Descriptor> Image::compute_brief(std::vector<Keypoint> keypointVector){
    //a
    //make copy
    Image *keypointCopy = Image::new_copy(this);
    //gaussian filter around 2.5 
    keypointCopy->smooth(2.5f,2.5f);
    std::vector<Descriptor> descriptorVector;
    static int randomOffsets[256][4] = {{2, -5, -4, 1}, {-6, 7, 6, -7}, {0, 0, 3, 6}, {-6, 7, 5, -3}, {-4, -3, -5, 7}, {-6, -5, -3, 1}, 
                        {3, 4, 4, -4}, {-3, -1, 0, 5}, {-5, 1, -3, 8}, {1, -6, -3, 5}, {2, -1, -5, -6}, {-1, 0, -5, -3}, {0, 5, 0, -5}, 
                        {2, 1, -8, 8}, {-2, -5, -8, 8}, {-7, -3, -3, 3}, {4, -3, 1, -8}, {-6, 5, 5, -6}, {1, -7, -3, -2}, {6, -5, -8, -2}, 
                        {-6, -4, 3, 3}, {-3, -3, -6, 0}, {-1, 8, 4, -3}, {-3, -5, 7, -4}, {-5, 1, 6, 7}, {-1, 7, -1, 8}, {-5, -2, 1, -1}, 
                        {-7, 0, 0, 8}, {3, 2, 6, -6}, {2, 0, 8, 1}, {6, 4, -4, 8}, {-4, 4, 8, -2}, {-6, 8, -1, 0}, {-3, -8, 3, -1}, {3, 3, -3, -7}, 
                        {1, 5, 6, 6}, {-3, 0, -2, 8}, {-1, 4, 0, 8}, {-2, 3, 0, 6}, {7, 6, 2, -3}, {5, -1, 4, 7}, {-6, 7, -7, 6}, {-3, 2, 3, 1}, 
                        {8, -4, 7, 4}, {-3, -2, -6, -2}, {-3, 4, -5, -8}, {-3, -8, 7, 6}, {5, -8, 8, 0}, {-6, 6, -3, -1}, {8, 3, 5, 5}, {-7, 7, 3, 8}, 
                        {-4, -2, 2, 5}, {-3, -4, -6, -8}, {-8, -5, -8, 1}, {-4, -5, 0, -2}, {8, -4, 8, -3}, {-8, -6, 6, -1}, {8, -4, 6, -1}, {3, 5, -2, -6}, 
                        {-3, 4, 7, 0}, {-2, -4, -4, -7}, {-4, -4, 3, -8}, {-3, 4, 4, -5}, {-7, -7, 7, 4}, {3, 4, -4, -1}, {-5, 2, -3, 7}, {-7, 2, 3, -2}, 
                        {3, -8, 3, 3}, {-8, -2, -1, -7}, {-6, -1, -3, -5}, {7, 7, -6, -2}, {0, 8, 7, 0}, {1, -5, 4, -8}, {-2, -4, 2, -6}, {-8, -3, 8, 3}, 
                        {-5, 8, -2, -3}, {3, 1, -6, 0}, {5, 5, -2, -5}, {3, 0, -3, -4}, {6, -1, 5, 7}, {-7, -6, 0, 2}, {4, 1, -1, 6}, {1, -6, -8, -1}, 
                        {-3, 8, -5, -3}, {2, -3, -4, 7}, {8, 1, -8, 4}, {7, -7, 5, -2}, {-1, 4, -1, -6}, {8, -7, -3, 1}, {0, -7, -4, 7}, {7, -4, 0, 3}, 
                        {-2, 1, -8, 1}, {1, 0, -3, -1}, {-4, 4, -3, -3}, {-6, 6, -2, 2}, {-7, 1, -6, -2}, {7, -7, -8, -5}, {8, -6, 6, 8}, {4, 1, -3, 2}, 
                        {1, 2, -3, -2}, {8, 7, 5, -1}, {-6, 7, 0, 8}, {-5, -3, 2, 8}, {0, -4, 0, -7}, {7, 1, -1, -6}, {1, 5, -1, 4}, {0, -3, -5, 2}, 
                        {-7, 2, -8, -7}, {-3, -6, 5, 6}, {-2, 5, -3, -6}, {-8, -6, 4, 1}, {-6, 4, 3, -5}, {4, -6, -6, -4}, {1, 6, -5, -8}, {2, 0, 7, -5}, 
                        {7, 0, -4, 7}, {-2, -8, 0, 3}, {-5, 5, -7, -5}, {8, -5, 5, -8}, {8, 3, 4, 7}, {-3, -3, -7, 3}, {-8, -1, 0, 3}, {6, 3, -7, -8}, 
                        {3, -5, -5, 4}, {5, 6, 4, 8}, {-3, 7, -1, 1}, {8, -1, -2, 7}, {-8, -8, -7, -8}, {8, 2, -6, 8}, {3, -7, -1, -4}, {5, -2, 4, -5}, 
                        {2, 3, -3, 4}, {0, -2, 3, -2}, {4, 8, -5, 3}, {-8, -8, 6, -6}, {-3, -3, -8, 3}, {-2, -8, -6, 5}, {4, -4, 6, -7}, {8, 4, 5, 0}, 
                        {-3, 8, -4, 1}, {-7, -3, 1, -1}, {2, 1, 7, 1}, {4, -4, -1, 2}, {-7, 8, 0, -2}, {4, 1, 8, -2}, {-1, 4, 3, 6}, {-8, 5, -6, 0}, 
                        {-4, 4, 4, -4}, {4, 6, -2, 8}, {2, -8, -8, -5}, {0, 1, -2, 6}, {5, 5, -1, 3}, {-1, -7, -2, 0}, {3, -7, -4, 2}, {8, 6, 0, 4}, 
                        {8, 1, 5, 8}, {7, 7, -6, -7}, {-6, 0, 6, -7}, {6, -3, 3, -8}, {4, 7, 3, -7}, {-3, 3, -1, 5}, {4, -1, 6, 5}, {5, -1, 1, 3}, 
                        {-5, -8, -3, 5}, {1, -6, -1, -3}, {6, -3, -4, 6}, {1, 1, -3, -4}, {-2, 4, -6, -6}, {-5, 5, -5, -6}, {1, 8, 5, -1}, {0, -5, -6, 7}, 
                        {2, 1, 3, 7}, {-8, -3, -3, -3}, {8, -8, -4, -3}, {0, -8, 7, -5}, {1, 1, -6, -8}, {-3, -6, -1, -6}, {7, -5, -2, 5}, {-1, 0, 4, -3}, 
                        {8, -1, 7, 8}, {7, 2, 1, -4}, {-2, -6, 8, -3}, {5, 6, 4, -5}, {5, -2, -2, 2}, {-8, -1, -1, 8}, {1, 5, -5, 1}, {-5, -3, 8, 8}, {0, 1, 0, 2}, 
                        {-8, -2, 4, 5}, {-5, 3, -7, -2}, {-2, -4, -4, 0}, {-5, -7, -6, -4}, {4, -5, -8, 3}, {-5, 8, 4, 2}, {-6, -5, -3, -7}, {-2, -8, 7, -3}, 
                        {0, 5, -3, 5}, {-5, 7, -4, -3}, {7, -4, -3, 1}, {8, 6, 2, -8}, {-5, 6, -1, 7}, {1, 0, 7, -5}, {4, 3, -1, -8}, {-4, -1, -1, -8}, {4, 3, -6, 5}, 
                        {-5, -8, -5, 8}, {-6, 5, 0, -2}, {2, 3, -1, 0}, {-1, -6, 0, 1}, {2, -6, 0, 2}, {-2, 6, 2, -8}, {-4, -3, -5, -6}, {-5, 7, 8, 8}, {7, 5, -7, 4}, 
                        {-2, 4, -8, 3}, {8, -1, 8, 1}, {-1, 7, 0, 4}, {-4, -8, -2, -3}, {-7, 5, -8, -3}, {5, 7, 0, 0}, {5, 1, -8, 6}, {-8, 3, 1, 3}, {-1, 0, 2, 7}, {2, 2, -1, 4}, 
                        {-5, 6, 3, 2}, {-8, -6, -2, 2}, {3, -2, -1, -5}, {3, 7, 8, -7}, {4, 3, 5, -3}, {3, -8, 6, 7}, {-4, -6, 4, -6}, {-7, -3, -8, -2}, {8, 7, 3, -7}, {1, 8, -4, -2},
                        {1, -2, -4, -4}, {-3, -8, 2, 5}, {3, -6, -4, -4}, {7, -4, 8, 5}, {6, 7, -5, 3}, {-2, 4, 6, 4}, {4, 8, 7, 8}, {8, 6, -5, -7}, {0, -8, -1, 0}, {-4, -4, -8, 2}, 
                        {7, -3, 6, 8}, {-7, -1, 0, -3}, {-1, -1, 8, -2}, {-1, 0, 8, 2}, {4, 1, 7, -6}, {3, -6, 4, -6}, {8, -2, -7, 2}, {4, -4, -4, 0}, {-7, 1, -4, 4}, {0, 6, 5, -3}, {7, 8, 8, 7}, {-7, 5, 3, 5}};
    
    for (int i=0; i<keypointVector.size();i++){
        Keypoint currentKey = keypointVector.at(i);
        int currentX= currentKey.x;
        int currentY= currentKey.y;
        if(currentX < 8 || currentY < 8 || (keypointCopy->m_width-8) <= currentX || (keypointCopy->m_height-8) <= currentY){ 
            continue;
        }
        Descriptor descs;
        for(int j=0; j<32; j++){
            uchar bits=0;
            for(int k=0; k<8; k++){
                int currentOffsetX0=randomOffsets[32*k+j][0];
                int currentOffsetY0=randomOffsets[32*k+j][1];
                int currentOffsetX1=randomOffsets[32*k+j][2];
                int currentOffsetY1=randomOffsets[32*k+j][3];                
                uchar I0=keypointCopy->m_data[currentX+currentOffsetX0+(currentY+currentOffsetY0)*keypointCopy->m_width];
                uchar I1=keypointCopy->m_data[currentX+currentOffsetX1+(currentY+currentOffsetY1)*keypointCopy->m_width];                
                if(I0>I1){
                    bits+=1;
                }
                bits<<=1; //sll
            }
            descs.desc[j]=bits;
            descs.key_id=i;
        }
        descriptorVector.push_back(descs);
        

    }
    delete keypointCopy;                    
    return descriptorVector;
}
int hammingDist(uchar* desc1, uchar* desc2)
{
    int result;
    int setbits;
    for (int i=0; i<32;i++){
        result = desc1[i]^desc2[i];
        setbits=0;
    }
    while (result>0){
        setbits += result & 1;
        result>>=1;
    }
    return setbits;
}
std::vector<Match> Image::match_brief(std::vector<Descriptor> descriptorVector1,std::vector<Descriptor> descriptorVector2) {
    
    std::vector<Match> matchVector;
    Match match;
    int compareHammingDist;
    int initialHammingDist;
    int minHammingDist;
    for (int i=0; i<descriptorVector1.size();i++){
        Descriptor currentDi = descriptorVector1.at(i);
        int currentDiId = currentDi.key_id;
        uchar* currentDiDesc = currentDi.desc;
        initialHammingDist = hammingDist(currentDiDesc,descriptorVector2.at(0).desc);
        for(int j=0; j<descriptorVector2.size();j++){
            Descriptor currentDj= descriptorVector2.at(j);
            int currentDjId = currentDj.key_id;
            uchar* currentDjDesc = currentDj.desc;
            compareHammingDist = hammingDist(currentDiDesc,currentDjDesc); 
            if(compareHammingDist<initialHammingDist){
                minHammingDist = compareHammingDist;
            }   
            match.key_id1= currentDjId;        
        }
        match.key_id0= currentDiId;
        match.distance=minHammingDist;
        matchVector.push_back(match);
    }
    return matchVector;
}
}
