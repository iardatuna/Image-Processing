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
#include <cstdlib>
#include <iostream>

#include "image.hpp"

using ceng391::Image;
using ceng391::short_to_image;
using ceng391::Keypoint;
using ceng391::Descriptor;
using ceng391::Match;
using std::vector;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
        Image img(4, 4, 1);
        img.read_pnm("/tmp/small_city.pgm");
        Image rotated(img.w(), img.h(), 1);
        double theta = 1.0 * 3.1415926 / 180;
        img.rotate_centered(&rotated, theta);
        rotated.write_pnm("/tmp/small_city_crotated_1");

        float threshold = 1000.0f;
        float k = 0.06f;
        float sigma = 2.5f;
        
        //Harris corner
        vector<Keypoint> keys = img.harris_corners(threshold, k, sigma);
        cout << "Detected " << keys.size()
             << " keypoints on small_city.pgm" << endl;

        //Rotated Harris corner     
        vector<Keypoint> rotatedKeys = rotated.harris_corners(threshold, k, sigma);
        cout << "Detected " << rotatedKeys.size()
             << " keypoints on rotated small_city.pgm" << endl;
        rotated.write_pnm("/tmp/small_city_keypointed_crotated_1");

        //Descriptor 
        img.write_pnm("/tmp/keys");
        vector<Descriptor> descs = img.compute_brief(keys);
        cout << "Detected " << descs.size()
             << " descriptors on small_city.pgm" << endl;

        //Rotated Descriptor
        vector<Descriptor> rotateddescs = rotated.compute_brief(rotatedKeys);
        cout << "Detected " << rotatedKeys.size()
             << " descriptors on rotated small_city.pgm" << endl;   

        //Match Descriptor
        vector<Match> matchVec = Image::match_brief(descs,rotateddescs);
        cout << "Detected " << matchVec.size()
             << " matched descriptors on small_city.pgm and rotated small_city.pgm" << endl;        

        return EXIT_SUCCESS;
}

