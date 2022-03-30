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
layout (location = 1) in vec4 fNormal;
layout (location = 2) in vec3 fPerPointColor;

struct LightSource {
	vec3 position;
	vec3 color;
	float intensity;
};

uniform LightSource lightSource;
uniform vec3 ambient_color;

layout (location = 0) out vec4 color;

void main(void)
{
	vec3 lp = lightSource.position;
	vec3 n = normalize(fNormal).xyz;
	vec3 wi = normalize (lp - fPosition.xyz);
	vec3 wo = normalize (-fPosition).xyz;
	vec3 wh = normalize (wi + wo);
	vec4 fs = 0.3 * pow(max(0.0, dot(n, wh)), 5.0) * vec4(lightSource.color, 1.0);
	vec4 fd = 0.8 * vec4(ambient_color, 0.5);
	color = vec4(fPerPointColor, 1.0) * lightSource.intensity * (fd + fs) * max (0.0, dot (n, wi));
	//color = vec4(normalize(fNormal.xyz), 1.0);//vec4(0.5, 0.1, 0.7, 1.0);
}
