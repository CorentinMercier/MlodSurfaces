// --------------------------------------------------------------------------
// This file is part of the reference implementation for the paper
//    Moving Level-of-Detail Surfaces.
//    C. Mercier, T. Lescoat, P. Roussillon, T. Boubekeur, and J-M. Thiery
//    ACM Transaction On Graphics 2022
//    DOI: 10.1145/3528223.3530151
//
// All rights reserved. Use of this source code is governed by a
// MIT license that can be found in the LICENSE file.
// --------------------------------------------------------------------------

#version 450

layout (location = 0) in vec4 fPosition;

uniform vec3 pointsColor;
uniform vec3 ambient_color;

layout (location = 0) out vec4 color;

void main(void)
{
	color = vec4(pointsColor, 1.0);
	//color = vec4(normalize(fNormal.xyz), 1.0);//vec4(0.5, 0.1, 0.7, 1.0);
}
