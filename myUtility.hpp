/* Author: Ming Du
 * Feb.16 2016
 * 
 * Contains function that can be used to coorporate with alignment result and 3D coordinates 
 * that stored in files.
 */




#include "DataObject/Alignment.h"
#include "DataObject/Point.h"
#include "json/json.h"
#include "Alignment/BlastParser.h"

#include <vector>

void test();


/*
 * Read blast alignment results file into a vector of alignment class
 * my version using Json library
 */
bool ReadAligns(const char* filename, vector<Alignment>& AlVct);




/* Read .coords file into a array of Point
 * if p == null then function will allocate a new array for it
 * based on the size of the file
 */
bool ReadCoord(const char* filename, Point*& p);
