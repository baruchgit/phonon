#pragma once
template <typename T, size_t N>
class List
{
	T* buffer = nullptr;
public:
	List()
	{
		buffer = (T*) malloc(sizeof(T) * N);
	}
	List(List&& rhs) noexcept
	{
		rhs.buffer = this.buffer;
		this.buffer = nulltpr;
	}
	~List()
	{
		if (buffer != nullptr)
		{
			free(buffer);
		}
	}

	T& operator[](std::size_t idx) 
	{ 
		return buffer[idx]; 
	}
	const T& operator[](std::size_t idx) const { return buffer[idx]; }

	T* data() { return buffer; }
};

