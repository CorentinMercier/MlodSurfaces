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

uniform mat4 camera_modelview;
uniform mat4 camera_projection;
uniform mat4 normal_matrix;

//uniform vec3 camPos;
uniform int pointSize;

layout (location = 0) in vec4 position;
layout (location = 1) in vec4 normal;
layout (location = 2) in vec3 perPointColor;

layout (location = 0) out vec4 fPosition;
layout (location = 1) out vec4 fNormal;
layout (location = 2) out vec3 fPerPointColor;

void main(void)
{
	fPerPointColor = perPointColor;
    fNormal = normalize(normal_matrix * normal);
    fPosition = camera_modelview * position;
    gl_PointSize = pointSize;// max(1.0,pow(20.0/dot(abs(position.xyz-camPos.xyz),vec3(1.0,1.0,1.0)),2.0));
    gl_Position = camera_projection * fPosition;
}
