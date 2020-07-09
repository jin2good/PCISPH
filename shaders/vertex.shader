#version 330 core

layout(location = 0) in vec4 position;

out float fragDepth;

uniform mat4 worldMat, viewMat, projMat; // Model View Projection Matrix


void main()
{
	gl_Position = projMat * viewMat * worldMat * position;
	fragDepth = gl_Position.z / gl_Position.w;
	gl_PointSize = (1.0 - fragDepth) * 50.0;
}
