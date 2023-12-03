#ifndef IOInterface_hpp__
#define IOInterface_hpp__

#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <cstring>
#include <string>
#include <sstream>

namespace evogen
{
    class IOInterface
    {
    public:

        IOInterface();

        void set_fname(std::string file);                                 /* setting IO file name */        
        void set_fname(std::string file, bool is_genotypes);              /* setting IO file name, and explicitely define whether we use genotypes or haplotypes */
        
        void fgetdata( size_t samples,
                       size_t variants,
                       std::vector<std::vector<int>> &out );              /* Reads a binary SNP data stored in the PLINK bed format */
        
        void fgetdata( std::vector< std::vector<unsigned long> > &id_map,
                       std::vector<std::vector<bool>> &out );            /* Reads ASCII SNP data with extra info (the very first three columns -> outputs in the additional data structure) */

        void fgetdata( std::vector<std::vector<short>> &out );           /* Reads ASCII SNP data; this is for markers */
        void fgetdata( std::vector<std::vector<bool>> &out );            /* Reads ASCII SNP data; this is for markers */
        
        template <typename T>
        void fgetdata( std::vector<std::vector<T>> &out );               /* this is for data tables of type: float, double, long int, size_t; we are not parsing 1,0 as markers */                                   /* Reads general ASCII formated data */
        
        bool is_genotype_data();                                         /* returns true if data extracted is genotypes; otherwise it is haplotypes */

    private:

        template <typename T>
        void str_parse( std::string& snpStr,
                        std::vector<T>& markers );                      /* Parsing a SNP string */
        
        std::string io_file;                                                                  /* IO file name */
        bool use_genotypes;                                             /* determines if we read genotypes instead of haplotypes */
    };

} // end of namespace evogen

#endif // IOInterface_hpp__
