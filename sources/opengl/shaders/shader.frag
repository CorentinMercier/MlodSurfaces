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
layout (location = 1) in vec4 fNormal;

struct LightSource {
	vec3 position;
	vec3 color;
	float intensity;
};

uniform LightSource lightSource1;
uniform LightSource lightSource2;
uniform vec3 ambient_color;
uniform vec4 pointsColor;
uniform mat4 camera_modelview;
uniform bool normal_mode;

layout (location = 0) out vec4 color;

void main(void)
{
	vec3 n = normalize(fNormal).xyz;
	vec3 wo = normalize (-fPosition).xyz;
	//Lightsource1
	vec3 lp = lightSource1.position;//(camera_modelview * vec4(lightSource1.position, 1.0)).xyz;
	vec3 wi = normalize (lp - fPosition.xyz);
	vec3 wh = normalize (wi + wo);
	vec4 fs = 0.4 * pow(max(0.0, dot(n, wh)), 5.0) * vec4(lightSource1.color, 1.0);
	//Lightsource2
	vec3 lp2 = lightSource2.position;//(camera_modelview * vec4(lightSource2.position, 1.0)).xyz;
	vec3 wi2 = normalize (lp2 - fPosition.xyz);
	vec3 wh2 = normalize (wi2 + wo);
	vec4 fs2 = 0.4 * pow(max(0.0, dot(n, wh2)), 5.0) * vec4(lightSource2.color, 1.0);

	vec4 fd = 0.3 * vec4(ambient_color, 0.5);
	if (normal_mode)
		color = normalize(fNormal);
	else
		color = pointsColor * (lightSource1.intensity * (fd + fs) * max (0.0, dot (n, wi)) + lightSource2.intensity * (fd + fs2) * max (0.0, dot (n, wi2)));
}
