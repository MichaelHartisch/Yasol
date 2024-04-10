/*
*
* Yasol: commprint.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#include "commprint.h"


void CommPrint::reverse(char s[])
{
  char c;
  size_t i,j;
  for (i = 0, j = strlen(s)-1;i < j;i++,j--)
    {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
    }
}

void CommPrint::itoa(int n, char s[])
{
  int i,sign;
  if ((sign = n) < 0) n = -n;
  i = 0;
  do {
   s[i++] = n % 10 + '0';
  } while ((n /= 10) > 0);
  if (sign < 0) s[i++] = '-';
  s[i] = 0;
  reverse(s);
}


int CommPrint::mefprint(int ID, char *fmt ...)
{
   char fname[200];
   char ending[10];
   char IDStr[10];
   //return 0; 
   va_list argpoint;
   FILE *file=NULL;

   strcpy(fname,"./Solutions_");
   strcpy(ending,".log");

   sprintf(IDStr,"%d",ID);
   //fprintf(stderr,"%s\n",fname);
   //fprintf(stderr,"%s\n",IDStr);
   strcat(fname,IDStr);
   strcat(fname,ending);
   if(!(file=fopen( fname, "a" ))){
     puts("!!!!!!!!!!! Konnte mefprint-File nicht oeffnen!\n");
     return 0;
   }

   va_start(argpoint, fmt );
   vfprintf( file, fmt, argpoint );
   fflush(file);
   fclose(file);
   va_end(argpoint);
   return 0;
}

void CommPrint::solprint(std::string toWrite, const char *FileName){

  FILE *file = NULL;

  const char *test = toWrite.c_str();

  //strcopy(test, toWrite);

  file = fopen(FileName, "a");
  fprintf(file,"%s", test);

  fclose(file);
}
