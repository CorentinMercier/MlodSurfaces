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

layout (location = 0) in vec4 fPosition;

uniform vec3 pointsColor;
uniform vec3 ambient_color;

layout (location = 0) out vec4 color;

void main(void)
{
	color = vec4(pointsColor, 1.0);
	//color = vec4(normalize(fNormal.xyz), 1.0);//vec4(0.5, 0.1, 0.7, 1.0);
}
