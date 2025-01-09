# ifndef FILE_IO_H_
# define FILE_IO_H_

using namespace std;

//----------------------------------------//
//----------File IO Manipulation----------//
//----------------------------------------//
int file_length(string);                    // Determine length (number of rows) in input file
int file_width(string);                     // Determine width (number of columns) in input file
void file_2d_dims(string, int &, int &);    // Determine dimensions of a 2D field stored as x : y : f(x,y)
bool string2bool (std::string v);           // Determine if string is true or false

#endif /* FILE_IO_H_ */