#include <random>
#include <cstdint>
#include <iostream>
#include <fstream>
 
int main()
{
    std::seed_seq seq{1,2,3,4,5};
    std::vector<std::uint32_t> seeds(10);
    seq.generate(seeds.begin(), seeds.end());
    
    std::ofstream myfile;
    myfile.open ("seeds.txt");
    
    for (std::uint32_t n : seeds) {
        myfile << n << '\n';
    }
    myfile.close();

}
