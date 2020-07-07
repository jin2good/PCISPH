#version 330 core

in float fragDepth;

layout(location = 0) out vec4 vFragColor;

uniform vec4 u_Color;

vec3 lightDir = vec3(1.0);
vec4 mat_specular = vec4(1);
vec4 light_specular = vec4(1);

void main()
{

    // calculate normal from texture coordinates
    vec3 N;
    N.xy = gl_PointCoord - vec2(0.5);
    float mag = dot(N.xy, N.xy);
    if (mag > 0.25)
        discard;   // kill pixels outside circle
    N.z = sqrt(1.0 - mag); //why minus?

    // calcuate depth


    // calculate lighting
    float diffuse = max(0.0, dot(lightDir, N));

    vFragColor = vec4(vec3(u_Color), 1.0f) * diffuse;
}
