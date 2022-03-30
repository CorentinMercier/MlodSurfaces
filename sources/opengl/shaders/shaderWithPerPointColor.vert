// --------------------------------------------------------------------------
// Source code provided FOR REVIEW ONLY, as part of the submission entitled
// "Moving Level-of-Detail Surfaces".
//
// A proper version of this code will be released if the paper is accepted
// with the proper licence, documentation and bug fix.
// Currently, this material has to be considered confidential and shall not
// be used outside the review process.
//
// All right reserved. The Authors
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
