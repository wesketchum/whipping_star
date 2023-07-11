
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>


#include <h5cpp/all>

#include <chrono>
#include <vector>
#include <algorithm>

#include <cstddef>

using namespace std;

#include "opts.h"

typedef diy::DiscreteBounds Bounds;


struct Block
{
    static void*    create()            { return new Block; }
    static void     destroy(void* b)    { delete static_cast<Block*>(b); }

    void show_link(const diy::Master::ProxyWithLink& cp)
    {
      diy::RegularLink<Bounds>* link = static_cast<diy::RegularLink<Bounds>*>(cp.link());
      std::cout << "Block (" << cp.gid() << "): "
                << link->core().min[0]   << ' ' << link->core().min[1]   << ' ' << link->core().min[2] << " - "
                << link->core().max[0]   << ' ' << link->core().max[1]   << ' ' << link->core().max[2] << " : "
                << link->bounds().min[0] << ' ' << link->bounds().min[1] << ' ' << link->bounds().min[2] << " - "
                << link->bounds().max[0] << ' ' << link->bounds().max[1] << ' ' << link->bounds().max[2] << " : "
                << link->size()   << ' ' //<< std::endl
                << std::dec
                << std::endl;
    }
    //int nUni;
    //int nPts;
    ////std::vector<int> myUnis;
    ////std::vector<int> myPts;
    //int dsIndex(int uni, int pt) {
};


void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, size_t nbins, size_t nrows, bool verbose, const h5::ds_t& ds_1, const h5::ds_t& ds_2, const vector<double> data, const vector<double> specdata)
{
   //auto ds = h5::open(f_out, "/dataset");
   h5::write( ds_1, data, h5::current_dims{nrows,size}, h5::offset{cp.gid(),0}, h5::count{1,1}, h5::collective );

   //auto ds2d = h5::open(f_out, "/dataset2d");
   h5::write( ds_2, specdata, h5::current_dims{nrows,size}, h5::offset{cp.gid(),0}, h5::count{1, nbins}, h5::collective );
   
}


// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBlocks = 0;
    size_t nBins=10;
    size_t chunksize= 0;
    int nPoints=1000;
    int nUniverses=1000;
    std::string out_file="test.hdf5";
    // get command line arguments
    using namespace opts;
    Options ops(argc, argv);
    ops >> Option('p', "npoints",   nPoints,   "Number of parameter points");
    ops >> Option('u', "nuniverse", nUniverses, "Number of universes");
    ops >> Option('b', "nblocks",   nBlocks,   "Number of blocks");
    ops >> Option('n', "nbins",   nBins,   "Number of bins in 2d dataset");
    ops >> Option('c', "chunksize",   chunksize,   "chunksize");
    ops >> Option('o', "output",    out_file,  "Output filename.");
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    MPI_Info info  = MPI_INFO_NULL;


    size_t nRows = nPoints*nUniverses;
    if (chunksize==0) chunksize=nRows;

    h5::fd_t f_out = h5::create(out_file, H5F_ACC_TRUNC, h5::fcpl,  h5::mpiio({MPI_COMM_WORLD, info}));
    // NOTE: the chunk size must be >= nRows!!!
    // 1D data set
    h5::ds_t ds  = h5::create<double>(f_out,"dataset",   h5::max_dims{nRows, 1},     h5::chunk{chunksize,1}      | h5::alloc_time_early );
    // 2D data set
    h5::ds_t ds2 = h5::create<double>(f_out,"dataset2d", h5::max_dims{nRows, nBins}, h5::chunk{chunksize, nBins} | h5::alloc_time_early );

    std::vector<double> data;
    data.push_back(world.rank());

    std::vector<double> specdata;
    for (size_t i=0; i<nBins;++i) { specdata.push_back(1.1*world.rank()); }
    // diy initialization
    int dim(1);
    
    size_t blocks;
    if (nBlocks==0) blocks= nPoints*nUniverses;//world.size() * threads;
    else blocks=nBlocks;

    diy::FileStorage storage("./DIY.XXXXXX"); // used for blocks moved out of core
    Bounds domain;
    for (int i = 0; i < dim; ++i) {
      domain.min[i] = 0;
      domain.max[i] = blocks-1;
    }
    ////// choice of contiguous or round robin assigner
    diy::ContiguousAssigner   assigner(world.size(), blocks);
    //// decompose the domain into blocks
    //// This is a DIY regular way to assign neighbors. You can do this manually.
    diy::RegularDecomposer<Bounds> decomposer(dim, domain, blocks);

    int k = 2;       // the radix of the k-ary reduction tree
    diy::RegularBroadcastPartners comm(decomposer, k, true);

    diy::RegularMergePartners  partners(decomposer,  // domain decomposition
                                        k,           // radix of k-ary reduction
                                        true); // contiguous = true: distance doubling

    diy::Master master(world, 1, -1, &Block::create, &Block::destroy);
    diy::decompose(dim, world.rank(), domain, assigner, master);//, share_face, wrap, ghosts);

    //AddBlock create(master);
    //decomposer.decompose(world.rank(), assigner, create); 

    if( world.rank()==0 ) {
      fmt::print(stderr, "\n*** This is diy running h5cppwrite ***\n");
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    DIY blocks:              {}\n", blocks);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
      fmt::print(stderr, "***********************************\n");
    }


    master.foreach([world, verbose, nBins, nRows, ds, ds2, data, specdata](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block(b, cp, world.size(), world.rank(), nBins, nRows, verbose, ds, ds2, data, specdata ); });

    return 0;
}

