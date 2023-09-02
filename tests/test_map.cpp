#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

int main() {
	cout << "teste" << endl;
	ifstream infile("C:\\Users\\andriy\\Desktop\\monomers-git\\monomers\\params.in");
	
	string strKey, strValue;
	
	map<string, string> config;
	
	while (infile >> strKey >> strValue) {
		config[strKey] = strValue;
	}
	for (map<string, string>::iterator it=config.begin(); it!=config.end(); ++it)
		cout << it->first << " => " << it->second << '\n';
	
	return 0;
}



	
