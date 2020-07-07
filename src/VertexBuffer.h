#pragma once
#include "Global.h"

class VertexBuffer
{
private:
	unsigned int m_RendererID;
public:
	VertexBuffer(const void* data, unsigned int size);
	VertexBuffer(const void* data, unsigned int size,unsigned int RendererID);
	~VertexBuffer();

	void Bind() const;
	void Unbind() const;
	unsigned int getRendererID();
	void Update(const void* data, unsigned int size);
};