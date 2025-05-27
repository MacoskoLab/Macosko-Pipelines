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

// Calculate hamming distance between two strings
int hammingDist(string seq1, string seq2)
{
    int cou = 0;
    for (int i = 0; i < min(seq1.size(), seq2.size()); i++)
	{
        if (seq1[i] != seq2[i])
            cou++;
	}
    return cou;
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
	vector<string> barcodes;
	if (infile.is_open())
	{
		while (std::getline(infile, line)) 
		{
			vector<string> items = splitstr2string(line, '-');
			barcodes.push_back(items[0]);
		}
		infile.close();
	}
	
	unordered_map<string, string> mp;
	unordered_map<string, int> mpcou;
	unordered_map<int, int> distcou;
	unordered_set<string> selectedReads;
	
	ofstream outfile1;
	outfile1.open (outputFile1);

	//@NB501583:616:HLKTHAFXY:1:11101:14435:1045 1:N:0:ACTATGGA+NGAGC
	//TACTANACGAGTGTTCAACAAGGCGGCAATTGGTTGATCGAGGCCAGTCCTCATT
	//+
	//AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

    std::ifstream infile2(read1file);
	string line1, line2, line3;
    int j = 0, k = 0, k1 = 0;
	if (infile2.is_open())
	{
		while (std::getline(infile2, line))
		{
			j++;
			
			if (j % 4 == 0)
			{
				string bc = line2.substr(0, 32);

				if (mp.find(bc) == mp.end())
				{
					// calculate hamming distance		
					int minDist = 100;
					int minCou = 0;
					string matchedb = "";
					for (int i = 0; i < barcodes.size(); i++)
					{
						int dist = hammingDist(bc, barcodes[i]);
						if (dist == minDist)
							minCou++;
						if (dist < minDist)
						{
							minDist = dist;
							minCou = 1;
							matchedb = barcodes[i];
							
							if (minDist == 0)
								break;
						}
					}
					
					distcou[minDist]++;
					if (minDist <= dist_threshold && matchedb != "")
					{
						if ((uniquematch == "true" && minCou == 1) || uniquematch != "true")
						{
							mp[bc] = matchedb;
							mpcou[bc] = minCou;
						}
						else
						{
							mp[bc] = "NA";
							mpcou[bc] = 0;
						}
					}
					else
					{
						mp[bc] = "NA";
						mpcou[bc] = 0;
					}
				}
				if (mp[bc] != "NA")
				{
					outfile1 << line1 << endl;
					if (correctbarcode == "true")
						outfile1 << mp[bc] << line2.substr(32) << endl;
					else
						outfile1 << line2 << endl;
					outfile1 << line3 << endl;
					outfile1 << line << endl;
					selectedReads.insert(splitstr2string(line1, ' ')[0]);
					k++;
					if (mpcou[bc] > 1)
						k1++;
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
	cout << "Unique matched reads: " << to_string(k - k1) << endl;
	cout << "Percentage of unique matched reads: " << to_string(((float)(k - k1) / j) * 100) << endl;
	cout << "Multiple matched reads: " << to_string(k1) << endl;
	cout << "Percentage of multiple matched reads: " << to_string(((float)k1 / j) * 100) << endl;
	
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
	
	for (auto it = distcou.begin(); it != distcou.end(); it++)
		cout << to_string(it->first) << ": " << to_string(it->second) << endl;
	
	my_time = time(NULL); 
    printf("%s", ctime(&my_time));
	
    return 0; 
}

