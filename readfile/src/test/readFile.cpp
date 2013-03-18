/*
 * readFile.cpp
 *
 *  Created on: Mar 14, 2013
 *      Author: didi
 */

#include "readFile.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>

using namespace std;

bool ReadFile::compareFiles()	{
	bool result = true;
		ifstream file1;
		ifstream file2;
		string line_file1;
		string line_file2;
		string name_file1(("/home/didi/mytesting/readfile/src/test/test1.txt"));
		string name_file2(("/home/didi/mytesting/readfile/src/test/test2.txt"));
		char *endptr_line1;
		char *endptr_line2;
		double out_double1;
		double out_double2;

		int cpt=1;

		// open 2 files

		file1.open(name_file1.c_str());
		file2.open(name_file2.c_str());



		while(result && (file1.good() && file2.good())){
			// read line from files
			getline(file1,line_file1);
			getline(file2,line_file2);

			out_double1 = strtod(line_file1.c_str(),&endptr_line1);
			out_double2 = strtod(line_file2.c_str(),&endptr_line2);

			while(result) {

			//	cout << cpt++ << setw(5) << " : "    <<  " mutation : " << setw(15) << out_double1 <<  "\t mutation++ : " << out_double2 << endl;
			//	break;
				while (result && (endptr_line1 != '\0' && (out_double1 != 0)) && (endptr_line2 != '\0'&& (out_double2	 != 0))){
					cout << cpt++ << setw(5) << " : "    <<  " mutation : " << setw(15) << out_double1 <<  "\t mutation++ : " << out_double2 << endl;
//break;					
						out_double1 = strtod(endptr_line1,&endptr_line1);
						out_double2 = strtod(endptr_line2,&endptr_line2);
				//		cout << cpt++ << setw(5) << " : "  << " mutation : " << setw(15) << out_double1 <<  "\t mutation++ : " << out_double2 << endl;
						result = (out_double1 == out_double2);

				}
break;
			}
			

		}
			cout <<  "AVANT result  : " << result << endl;
			if(!file1.eof() || !file2.eof())	{
				result = false;
			}
			cout <<  "APRES result  : " << result << endl;

		file1.close();
		file2.close();



		return result;


}


