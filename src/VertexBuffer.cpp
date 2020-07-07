#include "VertexBuffer.h"


VertexBuffer::VertexBuffer(const void* data, unsigned int size)
{
    GLCall(glGenBuffers(1, &m_RendererID));               // id assigned for one buffer, later i can use this buffer using id 
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, m_RendererID));
    GLCall(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STREAM_DRAW));
}

VertexBuffer::VertexBuffer(const void* data, unsigned int size, unsigned int RendererID)
{
    m_RendererID = RendererID;            // id assigned for one buffer, later i can use this buffer using id 
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, m_RendererID));
    GLCall(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STREAM_DRAW));
}

VertexBuffer::~VertexBuffer()
{
    GLCall(glDeleteBuffers(1, &m_RendererID));
}

void VertexBuffer::Bind() const
{
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, m_RendererID));
}

void VertexBuffer::Unbind() const
{
    GLCall(glBindBuffer(GL_ARRAY_BUFFER, 0));
}

unsigned int VertexBuffer::getRendererID()
{
    return this->m_RendererID;
}

void VertexBuffer::Update(const void* data, unsigned int size) {
    GLCall(glBufferData(GL_ARRAY_BUFFER, size, data, GL_STREAM_DRAW));
}