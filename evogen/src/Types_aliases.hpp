#ifndef typesaliases_hpp__
#define typesaliases_hpp__

typedef unsigned short popid_t; // population ID
typedef unsigned long genlen_t; // distance in a genome
typedef int poplen_t; // distance in a genome
typedef std::tuple<genlen_t, popid_t> ancestry_segment; // pair that defines an origin of genomic segment
typedef std::tuple<float, size_t> selection_candidate; // selection candidate: <sel.value, id_pos>

#endif // typesaliases_hpp__