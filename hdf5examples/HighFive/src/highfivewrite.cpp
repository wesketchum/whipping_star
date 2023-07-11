
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

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>


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


void process_block(Block* b, diy::Master::ProxyWithLink const& cp, int size, int rank, size_t nbins, bool verbose, HighFive::File* f_out, const vector<double> data, const vector<double> specdata)
{
   HighFive::DataSet test = f_out->getDataSet("dataset");
   test.select({std::size_t(cp.gid())}, {1}).write(data);
   

   HighFive::DataSet spec = f_out->getDataSet("dataset2d");
   spec.select(   {std::size_t(    cp.gid()), 0}, {1, nbins}).write(specdata);
}


// --- main program ---//
int main(int argc, char* argv[])
{
    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    size_t nBlocks = 0;
    size_t nBins=10;
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
    ops >> Option('o', "output",    out_file,  "Output filename.");
    bool verbose     = ops >> Present('v', "verbose", "verbose output");
    if (ops >> Present('h', "help", "Show help"))
    {
        std::cout << "Usage:  [OPTIONS]\n";
        std::cout << ops;
        return 1;
    }

    
    // Create hdf5 file structure here 
    HighFive::File* f_out  = new HighFive::File(out_file,
			HighFive::File::ReadWrite|HighFive::File::Create|HighFive::File::Truncate,
			HighFive::MPIOFileDriver(MPI_COMM_WORLD,MPI_INFO_NULL));


    std::vector<size_t> ds_dims(1);
    ds_dims[0] = nPoints*nUniverses;
    ds_dims[1] = 1;
    f_out->createDataSet<double>("dataset", HighFive::DataSpace(ds_dims));

    std::vector<size_t> spec_dims(2);
    spec_dims[0] = nPoints*nUniverses;
    spec_dims[1] = nBins;
    f_out->createDataSet<double>("dataset2d", HighFive::DataSpace(spec_dims));

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
      fmt::print(stderr, "\n*** This is diy running highfivewrite ***\n");
      fmt::print(stderr, "\n    Output will be written to {}\n", out_file);
      fmt::print(stderr, "\n    DIY blocks:              {}\n", blocks);
      fmt::print(stderr, "\n    Points:  {}\n", nPoints);
      fmt::print(stderr, "\n    Universes: {}\n", nUniverses);
      fmt::print(stderr, "\n    Total size of dataset:  {}\n", nPoints*nUniverses);
      fmt::print(stderr, "***********************************\n");
    }


    master.foreach([world, verbose, nBins, f_out, data, specdata](Block* b, const diy::Master::ProxyWithLink& cp)
                           {process_block(b, cp, world.size(), world.rank(), nBins, verbose, f_out, data, specdata); });

    delete f_out;
    return 0;
}

