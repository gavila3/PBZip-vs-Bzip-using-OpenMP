#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <omp.h>

#include <sys/stat.h>
// check if file already exists (w/o using c++17)
bool file_exists( const std::string& filename )
{
    struct stat f;
    if ( stat( filename.c_str(), &f ) != -1 )
        return true;
    else
        return false;
}

#include <bzlib.h>

static int bz_compression_level = 9;
static int bz_workFactor = 0; // default
static int bz_verbose    = 0;

int max_chunks = 1; 
static int verbose = 0;

std::vector<char>
compress_chunk( const std::vector<char>& idata )
{
    using uint = unsigned int;
    const uint isize = idata.size();
          uint zsize = uint( isize * 1.02 ) + 600;

    std::vector<char> zdata( zsize );

    char* iptr = const_cast<char*>( idata.data() );
    char* zptr = const_cast<char*>( zdata.data() );

    auto ierr = BZ2_bzBuffToBuffCompress( zptr, &zsize,
                                          iptr, isize,
                                          bz_compression_level,
                                          bz_verbose,
                                          bz_workFactor );

    if ( ierr != BZ_OK ) {
        fprintf(stderr,"bzBuffToBuffCompress failed with %d\n", ierr);
        exit(1);
    }

    // Reclaim any extra memory.
    zdata.resize( zsize );

    return std::move( zdata );
}


std::array<size_t,3>
compress_stream( FILE* istream, FILE* ostream )
{
    // How big are the input chunks?
    // The compression level switch is [1-9] which means 100-900k.
    size_t chunk_size = bz_compression_level * 100 * 1024;

    int n_chunks = 0;
    size_t ibytes = 0;
    size_t obytes = 0;

    while ( not(feof(istream)) )
    {
        // Space for the input stream data.
        std::vector<std::vector<char>> chunks( max_chunks );
        int n_chunks_read = 0;
        
        // Read a chunk of data from the input stream.
        for (int i = 0; i < max_chunks; i++)  {
            chunks [i].resize(chunk_size); 
            auto bytes_read = fread( chunks[i].data(), sizeof(char), chunk_size, istream );

            // for debugging

            ibytes += bytes_read;
            n_chunks_read++;
            if (feof(istream))
              break;
          }
        std::vector<std::vector<char>> ochunks( n_chunks_read );
        #pragma omp parallel for
        for (int i = 0; i < n_chunks_read; i++) {
          ochunks[i] = compress_chunk( chunks[i] );
          }
        for (int i = 0; i < n_chunks_read; i++) {
          obytes += ochunks[i].size();
          
          fwrite( ochunks[i].data(), sizeof(char), ochunks[i].size(), ostream );
          }
        if (verbose > 1)
            fprintf(stderr,"chunk: %d, idata: %lu, odata: %lu\n", n_chunks, ibytes, obytes);

        n_chunks += n_chunks_read;
    }

    return { size_t(n_chunks), ibytes, obytes };
}


void usage(FILE *os)
{
    fprintf(os, "Usage: pbzip2 <options> <files>\n");
    fprintf(os, " -h | --help           print this message\n");
    fprintf(os, " -l | --level=<#>      compression level [1-9] (default: 9)\n");
    fprintf(os, " -c | --stdout         write the compressed data to stdout instead of <input>.bz2\n");
    fprintf(os, " -i | --stdin          read the input data from stdin instead of a file.\n");
    fprintf(os, " -v | --verbose        increment the verbosity level.\n");
    fprintf(os, " -f | --force          force overwrite of exiting .bz2 file if present.\n");
}

std::string to_string( const bool value ) { return (value) ? "true" : "false"; }

template <typename T>
std::string to_string (const std::vector<T>& v)
{
    std::stringstream ss;

    ss << "[";
    for (int i = 0; i < v.size(); ++i)
    {
        ss << v[i];
        if (i == v.size()-1)
            ss << "]";
        else
            ss << ", ";
    }

    return ss.str();
}


int main (int argc, char* argv[])
{
    bool from_stdin = false;
    bool to_stdout = false;
    bool force_overwrite = false;
    std::vector< std::string > in_files;

    for (int i = 1; i < argc; )
    {
        std::string key = argv[i++];

        if ( key == "-h" or key == "--help")
        {
            usage(stdout);
            return 0;
        }
        else if (key == "-l" or key == "--level")
            bz_compression_level = atoi( argv[i++] );
        else if (key == "-v" or key == "--verbose")
            verbose++;
        else if (key == "-c" or key == "--stdout")
            to_stdout = true;
        else if (key == "-i" or key == "--stdin")
            from_stdin = true;
        else if (key == "-f" or key == "--force")
            force_overwrite = true;
        else {
            in_files.push_back( key );
        }
    }

    if (from_stdin)
    {
        if (not(to_stdout) and verbose)
            fprintf(stderr,"Warning: writing to stdout by default when reading from stdin.\n");

        to_stdout = true;

        if ( in_files.size() != 0 ) {
            fprintf(stderr, "Error: cannot read from stdin and a file list\n");
            return 1;
        }
    }
    else {
        if ( in_files.size() == 0 ) {
            fprintf(stderr, "Error: no input files given\n");
            return 1;
        }
    }

    if (verbose > 1) {
        std::cerr << "Options:"
                  << " compression_level= " << bz_compression_level
                  << " verbose=" << to_string(verbose )
                  << " force="   << to_string(force_overwrite )
                  << " stdout="  << to_string(to_stdout )
                  << " stdin="   << to_string(from_stdin )
                  << " files="   << to_string( in_files ) << "\n";
    }

    if (from_stdin)
        in_files.push_back( "(stdin)" );

    size_t max_filename_length = 0;
    for (int i = 0; i < in_files.size(); ++i)
        max_filename_length = std::max( max_filename_length, in_files[i].length() );
#ifdef _OPENMP
    max_chunks = omp_get_max_threads();
#endif    
    for (auto ifile: in_files)
    {
        FILE *istream = stdin;
        if (not(from_stdin)) {
            istream = fopen(ifile.c_str(), "rb");
            if ( istream == nullptr ) {
                fprintf(stderr,"Error opening input file %s\n", ifile.c_str());
                return 1;
            }
        }

        FILE *ostream = stdout;
        if ( not(to_stdout) ) {
            std::string ofile = ifile + ".bz2";

            if ( not(force_overwrite) and file_exists(ofile) ) {
                fprintf(stderr, "Error: will not overwrite exiting .bz2 file %s. To overwrite, specify --force", ofile.c_str());
                return 2;
            }

            ostream = fopen( ofile.c_str(), "wb" );
            if ( ostream == nullptr ) {
                fprintf(stderr, "Error: failed to open %d for writing\n", ofile.c_str());
                return 3;
            }
        }

        if (verbose) {
            std::stringstream ss;
            ss << "  %" << max_filename_length <<  "s: ";
            fprintf(stderr, ss.str().c_str(), ifile.c_str() );
        }

        auto res = compress_stream( istream, ostream );

        if ( not(from_stdin) )
            fclose( istream );

        if ( not(to_stdout) )
            fclose( ostream );

        if (verbose)
        {
            size_t ibytes = res[1];
            size_t obytes = res[2];

            double ratio = double(ibytes) / double(obytes);
            double bits = 8. * double(obytes) / double(ibytes);
            double percent = 100. * (double(ibytes) - double(obytes)) / double(ibytes);

            fprintf(stderr, "%6.3f:1, %6.3f bits/byte, %5.2f%% saved, %lu in, %lu out.\n", ratio, bits, percent, ibytes, obytes);
            //"""   openmpi/openmpi-4.1.4.tar: 11.454:1,  0.698 bits/byte, 91.27% saved, 115025920 in, 10042839 out."""
        }
    }

    return 0;
}
