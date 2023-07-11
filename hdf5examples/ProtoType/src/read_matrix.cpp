//#ifdef H5_USE_EIGEN
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Easy.hpp>
//#endif
#include <iostream>
#include <iomanip>
#include <TFile.h>
#include <TMatrixT.h>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <file_name> <matrix_name> <grid_index>" << std::endl;
        return 1;
    }

    const char* fileName = argv[1];
    const char* matrixName = argv[2];
    Int_t igrid = std::stoi(argv[3]);

    // Open HDF5
    std::string fname = "/Users/wospakrk/testfile_80dm2_s2-2Tue_s2-T24.hdf5";
    H5Easy::File h5file(fname, H5Easy::File::ReadOnly);
    Eigen::MatrixXf spectrums = H5Easy::load<Eigen::MatrixXf>(h5file, "tree_spectrum/vec_energy_spectrum");
    Eigen::VectorXf spectrum = spectrums.row(igrid);

    // Open the ROOT file
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open the file '" << fileName << "'" << std::endl;
        return 1;
    }

    // Get the TMatrixT object from the file
    TMatrixT<Double_t>* matrix = dynamic_cast<TMatrixT<Double_t>*>(file->Get(matrixName));
    if (!matrix) {
        std::cerr << "Failed to retrieve the matrix '" << matrixName << "' from the file" << std::endl;
        file->Close();
        return 1;
    }

    // Get the number of bins in the z axis
    Int_t numBinsZ = matrix->GetNrows();
    Int_t numBinsX = matrix->GetNcols();
    Eigen::VectorXf rootspectrum(numBinsX);
    // Extract the specified column as a TMatrixT object
    //TMatrixT<Double_t> columnMatrix = matrix->GetSub(0, numBinsZ - 1, column, column);
    TMatrixT<Double_t> columnMatrix = matrix->GetSub(0, numBinsZ - 1, 0, numBinsX - 1);
    // Iterate over the elements in the column
    for (Int_t i = 0; i < numBinsZ; i++) {
        for (Int_t j = 0; j < numBinsX; j++) {
            Double_t element = columnMatrix(i, j);  // Access the element in the column
            rootspectrum(j) = element;
            // Double_t element = matrix(i, j);  // Access the element in the column
            // Do something with the element
            // std::cout << "Column " << column << ", Element " << i << ": " << element << std::endl;
            // std::cout << std::setprecision(16) << element << ", ";
        }
    }
    //std::cout << ";" << std::endl;
    std::cout << "diff: " << spectrum - rootspectrum << std::endl;
    // Close the file
    file->Close();

    return 0;
}
