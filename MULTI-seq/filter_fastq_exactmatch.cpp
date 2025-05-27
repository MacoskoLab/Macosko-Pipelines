#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <time.h>

using namespace std;

// Check if file exists
bool is_file_exist(const string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

// Split a string into an array of substrings
vector<string> splitstr2string(const string str, const char delim)
{
	stringstream ss(str);
	string s;
	vector<string> res;
	while (getline(ss, s, delim))
	{
		res.push_back(s);
	}
	return res;
}

int main(int argc, char const *argv[]) 
{	
	time_t my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
	if (argc != 9)
	{
		cout << "Please provide eight parameters" << endl;
		return 0;
	}
	
    string read1file = argv[1];
	string read2file = argv[2];
	string barcodelist = argv[3];
	string outputFile1 = argv[4];
	string outputFile2 = argv[5];
	int dist_threshold = stoi(argv[6]);
	string correctbarcode = argv[7];
	string uniquematch = argv[8];
	
	if (!is_file_exist(read1file))
	{
		cout << read1file + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(read2file))
	{
		cout << read2file + " not found" << endl;
		return 0;
	}
	
	if (!is_file_exist(barcodelist))
	{
		cout << barcodelist + " not found" << endl;
		return 0;
	}
	
	std::ifstream infile(barcodelist);
	std::string line;
	unordered_set<string> barcodes;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			vector<string> items = splitstr2string(line, '-');
			barcodes.insert(items[0]);
		}
		infile.close();
	}
	
	unordered_set<string> selectedReads;
	
	ofstream outfile1;
	outfile1.open (outputFile1);

	//@NB501583:616:HLKTHAFXY:1:11101:14435:1045 1:N:0:ACTATGGA+NGAGC
	//TACTANACGAGTGTTCAACAAGGCGGCAATTGGTTGATCGAGGCCAGTCCTCATT
	//+
	//AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    std::ifstream infile2(read1file);
	string line1, line2, line3;
    int j = 0, k = 0;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line))
		{
			j++;
			
			if (j % 4 == 0)
			{
				string bc = line2.substr(0, 16);
				if (barcodes.find(bc) != barcodes.end())
				{
					outfile1 << line1 << endl;
					outfile1 << line2 << endl;
					outfile1 << line3 << endl;
					outfile1 << line << endl;
					selectedReads.insert(splitstr2string(line1, ' ')[0]);
					k++;
				}
			}
			else if (j % 4 == 1) {
				line1 = line;
			}
			else if (j % 4 == 2) {
				line2 = line;
			}
			else if (j % 4 == 3) {
				line3 = line;
			}
		}
		infile2.close();
	}
	outfile1.close();
	
	j = j / 4;
	cout << "Total reads: " << to_string(j) << endl;
	cout << "Unique matched reads: " << to_string(k) << endl;
	cout << "Percentage of unique matched reads: " << to_string(((float)k / j) * 100) << endl;
	
	ofstream outfile2;
	outfile2.open (outputFile2);

    std::ifstream infile3(read2file);
    j = 0;
	if (infile3.is_open())
	{
		while (std::getline(infile3, line))
		{
			j++;
			
			if (j % 4 == 0) {
				string read = splitstr2string(line1, ' ')[0];
				if (selectedReads.find(read) != selectedReads.end()) {
					outfile2 << line1 << endl;
					outfile2 << line2 << endl;
					outfile2 << line3 << endl;
					outfile2 << line << endl;
				}
			}
			else if (j % 4 == 1) {
				line1 = line;
			}
			else if (j % 4 == 2) {
				line2 = line;
			}
			else if (j % 4 == 3) {
				line3 = line;
			}
		}
		infile3.close();
	}
	outfile2.close();
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0; 
}

