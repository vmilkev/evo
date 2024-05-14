#include <iostream>
#include <plinkio/plinkio.h>

int
main(int argc, char *argv[])
{   
    if ( argc < 2 )
        std::cerr<<"There is not enough input arguments!"<<"\n";
        
    struct pio_file_t plink_file;
    snp_t *snp_buffer;
    size_t sample_id;
    int locus_id;

    if( pio_open( &plink_file, argv[ 1 ] ) != PIO_OK )
    {
        printf( "Error: Could not open %s\n", argv[ 1 ] );
        return EXIT_FAILURE;
    }

    if( !pio_one_locus_per_row( &plink_file ) )
    {
        printf( "This script requires that snps are rows and samples columns.\n" );
        return EXIT_FAILURE;
    }

    // get number of samples and SNPs:
    size_t n_samples = pio_num_samples(&plink_file);
    size_t n_variants = pio_num_loci(&plink_file);

    std::cout<<"samples: "<<n_samples<<", SNPs: "<<n_variants<<"\n";

    locus_id = 0;
    snp_buffer = (snp_t *) malloc( pio_row_size( &plink_file ) );
    while( pio_next_row( &plink_file, snp_buffer ) == PIO_OK )
    {
        for( sample_id = 0; sample_id < n_samples; sample_id++)
        {
            struct pio_sample_t *sample = pio_get_sample( &plink_file, sample_id );
            struct pio_locus_t *locus = pio_get_locus( &plink_file, locus_id );
            
            if ( (sample_id < 10) && (locus_id < 10) )
            {
                std::cout<<"Individual: "<<sample->iid<<" genotype: "<< (int)snp_buffer[ sample_id ] <<" at locus: "<<locus->name<<" at distance: "<<locus->bp_position<<"\n";
                //printf( "Individual %s has genotype %d for snp %s.\n", sample->iid, snp_buffer[ sample_id ], locus->name );
            }
        }

        locus_id++;
    }

    std::cout<<"samples: "<<pio_num_samples( &plink_file )<<", SNPs: "<<locus_id<<"\n";

    free( snp_buffer );
    pio_close( &plink_file );
    
    return EXIT_SUCCESS;
}