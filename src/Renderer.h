#pragma once

#include "Global.h"

#include "VertexArray.h"
#include "IndexBuffer.h"
#include "shader.h"


class Renderer
{
public:
    void Clear() const;
    void Draw(const VertexArray& va, const IndexBuffer& ib, const Shader& shader ) const; // Needs 1. vertex array, 2. index buffer, 3. shader
    void DrawPoints(const VertexArray& va, const Shader& shader, unsigned int count) const;

};