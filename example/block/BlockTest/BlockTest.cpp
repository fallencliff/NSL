// BlockTest.cpp : Defines the entry point for the console application.
//
//The following program shows how to allocate a block
#include <stdio.h>
#include <block.h>

using namespace gslcpp;

int main(int argc, char* argv[])
{
	CBlock b(100);
	
	printf ("length of block = %u\n", b.Size());
	printf ("block data address = %#x\n",b.GetDataPtr());

	return 0;
}

