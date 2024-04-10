/*
*
* Yasol: nomp.h -- Copyright (c) 2012-2017 Ulf Lorenz
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _NOMP
#define _NOMP

#define MPI_COMM_WORLD 0
#define MPI_BYTE 0
#define MPI_CHAR 0
#define MPI_ANY_TAG 0
#define MPI_ANY_SOURCE 0

typedef struct MPI_Status_t 
{
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
} MPI_Status;

typedef struct MPI_Request_t 
{
  int MPI_SOURCE;
} MPI_Request;

void MPI_Init(int*, char***);
void MPI_Buffer_attach(void*,int);
void MPI_Comm_rank(int, int*);
void MPI_Comm_size(int, int*);
void MPI_Wait(MPI_Request *, MPI_Status*);
void MPI_Barrier(int);
void MPI_Finalize(void);
void MPI_Ibsend(void *, int, int , int, int, int, MPI_Request *);
void MPI_Isend(void *, int, int , int, int, int, MPI_Request *);
void MPI_Irecv(void *, int , int , int , int , int, MPI_Request*);
void MPI_Iprobe(int , int, int, int*, MPI_Status*);
void MPI_Test(MPI_Request*, int *, MPI_Status*);
void MPI_Recv(void *buf, int count,int datatype, int source, int tag, int* comm, MPI_Status *status);
void MPI_Send(const void *buf, int count, int datatype, int dest, int tag, int* comm);
void MPI_Get_count(MPI_Status*, int, int*);

#endif
