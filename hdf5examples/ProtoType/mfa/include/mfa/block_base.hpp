//--------------------------------------------------------------
// base class for one diy block
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
//--------------------------------------------------------------

#include    <mfa/mfa.hpp>

#include    <diy/master.hpp>
#include    <diy/reduce-operations.hpp>
#include    <diy/decomposition.hpp>
#include    <diy/assigner.hpp>
#include    <diy/io/block.hpp>
#include    <diy/pick.hpp>
#include    <diy/thirdparty/fmt/format.h>
#include    <diy/reduce.hpp>
#include    <diy/partners/merge.hpp>

#include    <stdio.h>

#include    <Eigen/Dense>

#include    <random>

using namespace std;

using Index = MatrixXf::Index;

template <typename T>
using Bounds = diy::Bounds<T>;
template <typename T>
using RCLink = diy::RegularLink<diy::Bounds<T>>;

template <typename T>
using Decomposer = diy::RegularDecomposer<Bounds<T>>;



// struct ModelInfo
// {
//     ModelInfo(int dom_dim_, int pt_dim_) :
//         dom_dim(dom_dim_),
//         pt_dim(pt_dim_),
//         weighted(true),
//         local(false),
//         verbose(1)
//     {
//         geom_p.resize(dom_dim);
//         vars_p.resize(pt_dim - dom_dim);
//         for (auto i = 0; i < vars_p.size(); i++)
//             vars_p[i].resize(dom_dim);
//         geom_nctrl_pts.resize(dom_dim);
//         vars_nctrl_pts.resize(pt_dim - dom_dim);
//         for (auto i = 0; i < vars_nctrl_pts.size(); i++)
//             vars_nctrl_pts[i].resize(dom_dim);

//         regularization = 0;
//     }
//     virtual ~ModelInfo()                        {}

//     int                 dom_dim;                // domain dimensionality
//     int                 pt_dim;                 // point dimensionality (> dom_dim)
//     VectorXi            model_dims;             // (currently unused) pt dim for each mfa_data
//     vector<int>         geom_p;                 // degree in each dimension of geometry
//     vector<vector<int>> vars_p;                 // degree in each dimension of each science variable vars_p[var][dim]
//     vector<int>         geom_nctrl_pts;         // number of input points in each dimension of geometry
//     vector<vector<int>> vars_nctrl_pts;         // number of input pts in each dim of each science variable vars_nctrl_pts[var][dim]
//     bool                weighted;               // solve for and use weights (default = true)
//     bool                local;                  // solve locally (with constraints) each round (default = false)
//     float               regularization;         // smoothing parameter for unstructured data with nonuniform point density (value of 0 does nothing)
//     bool                reg1and2;               // regularize with 1st and 2nd derivatives (false --> 2nd derivs only)
//     int                 verbose;                // debug level (default = 1)
// };


// block
template <typename T>
struct BlockBase
{
    // info for DIY
    int                 dom_dim;                // dimensionality of domain (geometry)
    int                 pt_dim;                 // dimensionality of full point (geometry + science vars)
    VectorX<T>          bounds_mins;            // local domain minimum corner
    VectorX<T>          bounds_maxs;            // local domain maximum corner
    VectorX<T>          core_mins;              // local domain minimum corner w/o ghost
    VectorX<T>          core_maxs;              // local domain maximum corner w/o ghost

    bool is_ray_model{false};
    VectorX<T>  box_mins;   // TODO remove this. This is solely for line integrals but should not be a permanent member
    VectorX<T>  box_maxs;   

    // data sets
    mfa::PointSet<T>    *input;                 // input data
    mfa::PointSet<T>    *approx;                // output data
    mfa::PointSet<T>    *errs;                  // error field

    // MFA object
    mfa::MFA<T>         *mfa;

    // errors for each science variable
    vector<T>           max_errs;               // maximum (abs value) distance from input points to curve
    vector<T>           sum_sq_errs;            // sum of squared errors

    VectorX<T>          overlaps;               // local domain overlaps in each direction

    // blending output
    // will be computed only in core, for each block, using neighboring approximations if necessary
    //  needed only at decoding stage
    VectorXi            ndom_outpts;            // number of output points in each dimension
    MatrixX<T>          blend;                  // output result of blending (lexicographic order)
    vector<Bounds<T>>   neighOverlaps;          // only the overlapping part of the neighbors to current core is important
                                                // will store the actual overlap of neighbors
    vector<Bounds<T>>   overNeighCore;          // part of current domain over neighbor core;
                                                // will send decoded values for this patch, in the resolution
                                                //  needed according to that core; symmetric with neighOverlaps
    vector<int>         map_dir;                // will map current directions with global directions
    vector<T>           max_errs_reduce;        // max_errs used in the reduce operations, plus location (2 T vals per entry)

    // zero-initialize pointers during default construction
    BlockBase() : 
        mfa(nullptr), 
        input(nullptr), 
        approx(nullptr), 
        errs(nullptr) { }

    ~BlockBase()
    {
        delete mfa;
        delete input;
        delete approx;
        delete errs;
    }

    // initialize an empty block that was previously added
    void init_block(
            const Bounds<T>&    core,               // block bounds without any ghost added
            const Bounds<T>&    domain,             // global data bounds
            int                 dom_dim_,           // domain dimensionality
            int                 pt_dim_,            // point dimensionality
            T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
    {
        dom_dim = dom_dim_;
        pt_dim  = pt_dim_;

        // NB: using bounds to hold full point dimensionality, but using core to hold only domain dimensionality
        bounds_mins.resize(pt_dim);
        bounds_maxs.resize(pt_dim);
        core_mins.resize(dom_dim);
        core_maxs.resize(dom_dim);
        // blending
        overlaps.resize(dom_dim);

        // manually set ghosted block bounds as a factor increase of original core bounds
        for (int i = 0; i < dom_dim; i++)
        {
            T ghost_amount = ghost_factor * (core.max[i] - core.min[i]);
            if (core.min[i] > domain.min[i])
                bounds_mins(i) = core.min[i] - ghost_amount;
            else
                bounds_mins(i)= core.min[i];

            if (core.max[i] < domain.max[i])
                bounds_maxs(i) = core.max[i] + ghost_amount;
            else
                bounds_maxs(i) = core.max[i];
            core_mins(i) = core.min[i];
            core_maxs(i) = core.max[i];
        }

        // debug
//         cerr << "core_mins: " << core_mins.transpose() << endl;
//         cerr << "core_maxs: " << core_maxs.transpose() << endl;
//         cerr << "bounds_mins: " << bounds_mins.transpose() << endl;
//         cerr << "bounds_maxs: " << bounds_maxs.transpose() << endl;
    }

    void setup_MFA(
            const       diy::Master::ProxyWithLink& cp,
            MFAInfo&    info)
    {
        if (mfa != nullptr)
        {
            cerr << "\nWarning: Overwriting existing MFA in setup_MFA!\n" << endl;
            delete mfa;
        }

        // Silence verbose output if not on rank 0
        info.verbose = info.verbose && cp.master()->communicator().rank() == 0; 

        // Construct MFA from MFAInfo
        this->mfa = new mfa::MFA<T>(info);
    }


    // fixed number of control points encode block
    void fixed_encode_block(
            const       diy::Master::ProxyWithLink& cp,
            MFAInfo&  info)
    {
        mfa->FixedEncode(*input, info.regularization, info.reg1and2, info.weighted, false);
    }

    // adaptively encode block to desired error limit
    void adaptive_encode_block(
            const diy::Master::ProxyWithLink& cp,
            T                                 err_limit,
            int                               max_rounds,
            MFAInfo&                          info)
    {
        if (!input->structured)
        {
            // Can't use adaptive encoding until support for unstructured data is
            // added to NewKnots and tmesh.all_knot_param_idxs
            cerr << "ERROR: Adaptive encoding not currently supported for unstructured data" << endl;
            exit(1);
        }

        VectorX<T> extents = bounds_maxs - bounds_mins;

        mfa->AdaptiveEncode(*input, err_limit, info.weighted, extents, max_rounds);
    }

    // decode entire block at the same parameter locations as 'input'
    void decode_block(
            const   diy::Master::ProxyWithLink& cp,
            int                                 verbose,        // debug level
            bool                                saved_basis)    // whether basis functions were saved and can be reused
    {
        if (approx)
        {
            cerr << "WARNING: Overwriting \"approx\" pointset in BlockBase::decode_block" << endl;
            delete approx;
        }
        approx = new mfa::PointSet<T>(input->params, input->model_dims());  // Set decode params from input params

        mfa->Decode(*approx, saved_basis);
    }

    // decode entire block over a regular grid
    void decode_block_grid(
        const   diy::Master::ProxyWithLink& cp,
        int                                 verbose,
        vector<int>&                        grid_size)
    {
        VectorXi grid_npts(dom_dim);
        for (int k = 0; k < dom_dim; k++)
        {
            grid_npts(k) = grid_size[k];
        }

        shared_ptr<mfa::Param<T>> grid_params = make_shared<mfa::Param<T>>(grid_npts, VectorX<T>::Zero(dom_dim), VectorX<T>::Ones(dom_dim));

        if (approx)
        {
            cerr << "WARNING: Overwriting \"approx\" pointset in BlockBase::decode_block_grid" << endl;
            delete approx;
        }
        approx = new mfa::PointSet<T>(grid_params, input->model_dims());

        mfa->Decode(*approx, false);
    }

    // decode one point
    void decode_point(
            const   diy::Master::ProxyWithLink& cp,
            const VectorX<T>&                   param,          // parameters of point to decode
            Eigen::Ref<VectorX<T>>              cpt)            // (output) decoded point
                                                                // using Eigen::Ref instead of C++ reference so that pybind11 can pass by reference
    {
        mfa->Decode(param, cpt);
    }

    // differentiate one point
    void differentiate_point(
            const diy::Master::ProxyWithLink&   cp,
            const VectorX<T>&                   param,      // parameters of point to decode
            int                                 deriv,      // which derivative to take (1 = 1st, 2 = 2nd, ...) in each domain dim.
            int                                 partial,    // limit to partial derivative in just this dimension (-1 = no limit)
            int                                 var,        // differentiate only this one science variable (0 to nvars -1, -1 = all vars)
            Eigen::Ref<VectorX<T>>              cpt)        // (output) decoded point
                                                            // using Eigen::Ref instead of C++ reference so that pybind11 can pass by reference
    {
        VectorXi derivs(dom_dim);                           // degree of derivative in each domain dimension

        // degree of derivative is same for all dimensions
        for (auto i = 0; i < derivs.size(); i++)
            derivs(i) = deriv;

        // optional limit to one partial derivative
        if (deriv && dom_dim > 1 && partial >= 0)
        {
            for (auto i = 0; i < dom_dim; i++)
            {
                if (i != partial)
                    derivs(i) = 0;
            }
        }

        mfa->Decode(param, cpt, derivs);
    }

    void definite_integral(
            const diy::Master::ProxyWithLink&   cp,
            int                                 verbose,
                  VectorX<T>&                   output,
            const VectorX<T>&                   lim_a,
            const VectorX<T>&                   lim_b)
    {
        for (auto k = 0; k < mfa->nvars(); k++)
        {
            mfa->DefiniteIntegral(mfa->var(k), output, verbose, lim_a, lim_b);
        }

        T scale = (core_maxs - core_mins).prod();
        output *= scale;
    }

    void integrate_axis_ray(
        const diy::Master::ProxyWithLink&   cp,
        T                                   alpha,
        T                                   rho,
        T                                   u0,
        T                                   u1,
        T                                   scale,
        VectorX<T>&                         output)
    {
        // TODO: Create a subclass RayMFA from MFA which has methods like this in it?
        if (!is_ray_model)
        {
            cerr << "ERROR: Attempting to call BlockBase::integrate_axis_ray but is_ray_model=false" << endl;
            exit(1);
        }

        T alpha_param = (alpha - bounds_mins(2)) / (bounds_maxs(2) - bounds_mins(2));
        T rho_param = (rho - bounds_mins(1)) / (bounds_maxs(1) - bounds_mins(1));
        
        // TODO: this is first science variable only
        mfa->IntegrateAxisRay(mfa->var(0), alpha_param, rho_param, u0, u1, output);

        output *= scale;
    }

    void integrate_block(
            const diy::Master::ProxyWithLink&   cp,
            int                                 int_dim,
            int                                 verbose)
    {
        if (approx)
        {
            cerr << "WARNING: Overwriting \"approx\" pointset in BlockBase::integrate_block" << endl;
            delete approx;
        }
        approx = new mfa::PointSet<T>(input->params, input->model_dim());

        approx->domain.leftCols(dom_dim) = input->domain.leftCols(dom_dim);

        for (auto k = 0; k < mfa->nvars(); k++)
        {
            mfa->IntegratePointSet(mfa->var(k), *approx, int_dim, verbose, dom_dim + k, dom_dim + k);
        }
        
        // scale the integral
        // MFA computes integral with respect to parameter space
        // This scaling factor returns the integral wrt physical space (assuming domain parameterization)
        T scale = core_maxs(int_dim) - core_mins(int_dim);
        approx->domain.rightCols(mfa->pt_dim - mfa->geom_dim) *= scale;
    }

    // differentiate entire block at the same parameter locations as 'input'
    void differentiate_block(
            const diy::Master::ProxyWithLink& cp,
            int                               verbose,  // debug level
            int                               deriv,    // which derivative to take (1 = 1st, 2 = 2nd, ...) in each domain dim.
            int                               partial,  // limit to partial derivative in just this dimension (-1 = no limit)
            int                               var)      // differentiate only this one science variable (0 to nvars -1, -1 = all vars)
    {
        if (approx)
        {
            cerr << "WARNING: Overwriting \"approx\" pointset in BlockBase::differentiate_block" << endl;
            delete approx;
        }
        approx = new mfa::PointSet<T>(input->params, input->model_dim());  // Set decode params from input params
        
        VectorXi derivs = deriv * VectorXi::Ones(dom_dim);
        if (deriv && dom_dim > 1 && partial >= 0)   // optional limit to one partial derivative
        {
            for (auto i = 0; i < dom_dim; i++)
            {
                if (i != partial)
                    derivs(i) = 0;
            }
        }

        mfa->Decode(*approx, false, derivs);

        // the derivative is a vector of same dimensionality as domain
        // derivative needs to be scaled by domain extent because u,v,... are in [0.0, 1.0]
        if (deriv)
        {
            if (dom_dim == 1 || partial >= 0) // TODO: not for mixed partials
            {
                if (dom_dim == 1)
                    partial = 0;
                for (auto j = 0; j < approx->pt_dim; j++)
                    // scale once for each derivative
                    for (auto i = 0; i < deriv; i++)
                        approx->domain.col(j) /= (bounds_maxs(partial) - bounds_mins(partial));
            }
        }

        // for plotting, set the geometry coordinates to be the same as the input
        if (deriv)
            for (auto i = 0; i < dom_dim; i++)
                approx->domain.col(i) = input->domain.col(i);
    }

    // compute error field and maximum error in the block
    // uses coordinate-wise difference between values
    void range_error(
            const   diy::Master::ProxyWithLink& cp,
            int     verbose,                                // output level
            bool    decode_block_,                          // decode entire block first
            bool    saved_basis)                            // whether basis functions were saved and can be reused
    {
        if (input == nullptr)
        {
            cerr << "ERROR: Cannot compute range_error; no valid input data" << endl;
            exit(1);
        }
        if (errs != nullptr)
        {
            cerr << "Warning: Overwriting existing error field in BlockBase::range_error" << endl;
            delete errs;
        }
        errs = new mfa::PointSet<T>(input->params, input->model_dims());

        // saved_basis only applies when not using tmesh
#ifdef MFA_TMESH
        saved_basis = false;
#endif
        // Decode entire block and then compare to input
        if (decode_block_)
        {
            decode_block(cp, verbose, saved_basis);

            input->abs_diff(*approx, *errs, verbose);
        }
        else // Compute error at each input point on the fly
        {
            mfa->AbsPointSetError(*input, *errs, verbose);
        }

        // Compute error metrics
        for (int k = 0; k < mfa->nvars(); k++)
        {
            T l2sq = 0, linf = 0;
            sum_sq_errs[k] = 0;
            max_errs[k] = 0;
            VectorX<T> var_err(errs->var_dim(k));

            for (int i = 0; i < errs->npts; i++)
            {
                errs->var_coords(i, k, var_err);    // Fill var_err

                l2sq = var_err.squaredNorm();
                // linf = var_err.cwiseAbs().maxCoeff();

                sum_sq_errs[k] += l2sq;
                if (l2sq > max_errs[k])
                    max_errs[k] = l2sq;
            }

            // Important! Take square root of max squared norm
            max_errs[k] = sqrt(max_errs[k]);
        }

        // copy the computed error in a new array for reduce operations
        max_errs_reduce.resize(max_errs.size() * 2);
        for (auto i = 0; i < max_errs.size(); i++)
        {
            max_errs_reduce[2 * i] = max_errs[i];
            max_errs_reduce[2 * i + 1] = cp.gid(); // use converter from type T to integer
        }
     }

    void print_block(const diy::Master::ProxyWithLink& cp,
            bool                              error)       // error was computed
    {
        //fprintf(stderr, "gid = %d\n", cp.gid());
//         cerr << "domain\n" << domain << endl;

        // max errors over all science variables
        // Initializing to 0 so we don't have to initialize with first var's error below
        T all_max_err = 0, all_max_norm_err = 0, all_max_sum_sq_err = 0, all_max_rms_err = 0, all_max_norm_rms_err = 0;
        // variables where the max error occurs
        int all_max_var, all_max_norm_var, all_max_sum_sq_var, all_max_rms_var, all_max_norm_rms_var;   

        VectorXi tot_nctrl_pts_dim = VectorXi::Zero(dom_dim);        // total num. ctrl. pts. per dim.
        size_t tot_nctrl_pts = 0;                                                       // total number of control points

        // geometry
        //cerr << "\n------- geometry model -------" << endl;
        for (auto j = 0; j < mfa->geom().tmesh.tensor_prods.size(); j++)
        {
            tot_nctrl_pts_dim += mfa->geom().tmesh.tensor_prods[j].nctrl_pts;
            tot_nctrl_pts += mfa->geom().tmesh.tensor_prods[j].nctrl_pts.prod();
        }
        // print number of control points per dimension only if there is one tensor
        //if (mfa->geom().tmesh.tensor_prods.size() == 1)
        //    cerr << "# output ctrl pts     = [ " << tot_nctrl_pts_dim.transpose() << " ]" << endl;
        //cerr << "tot # output ctrl pts = " << tot_nctrl_pts << endl;

        //  debug: print control points and weights
//         print_ctrl_weights(geometry.mfa_data->tmesh);
        // debug: print knots
//         print_knots(geometry.mfa_data->tmesh);

        //fprintf(stderr, "# output knots        = [ ");
        //for (auto j = 0 ; j < mfa->geom().tmesh.all_knots.size(); j++)
        //{
        //    fprintf(stderr, "%ld ", mfa->geom().tmesh.all_knots[j].size());
        //}
        //fprintf(stderr, "]\n");

        //cerr << "-----------------------------" << endl;

        // science variables
        //cerr << "\n----- science variable models -----" << endl;
        for (auto i = 0; i < mfa->nvars(); i++)
        {
            int min_dim = mfa->var(i).min_dim;
            int max_dim = mfa->var(i).max_dim;
            int vardim  = mfa->var_dim(i);
            MatrixX<T> varcoords = input->domain.middleCols(min_dim, vardim);
            // T range_extent = input->domain.col(dom_dim + i).maxCoeff() - input->domain.col(dom_dim + i).minCoeff();

            // range_extents_max is a vector containing the range extent in each component of the science variable
            // So, the size of 'range_extents_max' is the dimension of the science variable
            // VectorX<T> range_extents_max = varcoords.colwise().maxCoeff() - varcoords.colwise().minCoeff();

            // 'range_extents' is the difference between the largest and smallest vector norms of the input data
            // So, for each point in the input, we compute the norm of that point; then we take the difference of the 
            // largest and smallest norms. 
            //
            // NOTE! For scalar science variables, this is different from our older definition. Essentially, we 
            // are taking the absolute value of every point before taking the max/min. Thus this could vary from
            // the older definition of range_extent by up to a factor of 2.
            T range_extent  = varcoords.rowwise().norm().maxCoeff() - varcoords.rowwise().norm().minCoeff();

            //cerr << "\n---------- var " << i << " ----------" << endl;
            tot_nctrl_pts_dim   = VectorXi::Zero(dom_dim);
            tot_nctrl_pts       = 0;
            for (auto j = 0; j < mfa->var(i).tmesh.tensor_prods.size(); j++)
            {
                tot_nctrl_pts_dim   += mfa->var(i).tmesh.tensor_prods[j].nctrl_pts;
                tot_nctrl_pts       += mfa->var(i).tmesh.tensor_prods[j].nctrl_pts.prod();
            }
            // print number of control points per dimension only if there is one tensor
            //if (mfa->var(i).tmesh.tensor_prods.size() == 1)
              //  cerr << "# ouput ctrl pts      = [ " << tot_nctrl_pts_dim.transpose() << " ]" << endl;
            //cerr << "tot # output ctrl pts = " << tot_nctrl_pts << endl;

            //  debug: print control points and weights
//             print_ctrl_weights(mfa->var(i)->tmesh);
            // debug: print knots
//             print_knots(mfa->var(i)->tmesh);


            //fprintf(stderr, "# output knots        = [ ");
            //for (auto j = 0 ; j < mfa->var(i).tmesh.all_knots.size(); j++)
            //{
            //    fprintf(stderr, "%ld ", mfa->var(i).tmesh.all_knots[j].size());
            //}
            //fprintf(stderr, "]\n");

            //cerr << "-----------------------------" << endl;

            if (error)
            {
                T rms_err = sqrt(sum_sq_errs[i] / (input->npts));
                fprintf(stderr, "range extent          = %e\n",  range_extent);
                fprintf(stderr, "max_err               = %e\n",  max_errs[i]);
                fprintf(stderr, "normalized max_err    = %e\n",  max_errs[i] / range_extent);
                fprintf(stderr, "sum of squared errors = %e\n",  sum_sq_errs[i]);
                fprintf(stderr, "RMS error             = %e\n",  rms_err);
                fprintf(stderr, "normalized RMS error  = %e\n",  rms_err / range_extent);

                // find max over all science variables
                if (max_errs[i] > all_max_err)
                {
                    all_max_err = max_errs[i];
                    all_max_var = i;
                }
                if (max_errs[i] / range_extent > all_max_norm_err)
                {
                    all_max_norm_err = max_errs[i] / range_extent;
                    all_max_norm_var = i;
                }
                if (sum_sq_errs[i] > all_max_sum_sq_err)
                {
                    all_max_sum_sq_err = sum_sq_errs[i];
                    all_max_sum_sq_var = i;
                }
                if (rms_err > all_max_rms_err)
                {
                    all_max_rms_err = rms_err;
                    all_max_rms_var = i;
                }
                if (rms_err / range_extent > all_max_norm_rms_err)
                {
                    all_max_norm_rms_err = rms_err / range_extent;
                    all_max_norm_rms_var = i;
                }
            }
            //cerr << "-----------------------------" << endl;
        }

        if (error)
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "Maximum errors over all science variables:\n");
            fprintf(stderr, "max_err                (var %d)    = %e\n",  all_max_var,          all_max_err);
            fprintf(stderr, "normalized max_err     (var %d)    = %e\n",  all_max_norm_var,     all_max_norm_err);
            fprintf(stderr, "sum of squared errors  (var %d)    = %e\n",  all_max_sum_sq_var,   all_max_sum_sq_err);
            fprintf(stderr, "RMS error              (var %d)    = %e\n",  all_max_rms_var,      all_max_rms_err);
            fprintf(stderr, "normalized RMS error   (var %d)    = %e\n",  all_max_norm_rms_var, all_max_norm_rms_err);
        }

       //cerr << "\n-----------------------------------" << endl;

        //  debug: print approximated points
//         cerr << approx->npts << " approximated points\n" << approx->domain << endl;

        //fprintf(stderr, "# input points        = %ld\n", input->npts);
        //fprintf(stderr, "compression ratio     = %.2f\n", compute_compression());
    }

    // compute compression ratio
    float compute_compression()
    {
        float in_coords = (input->npts) * (input->pt_dim);
        float out_coords = 0.0;
        for (auto j = 0; j < mfa->geom().tmesh.tensor_prods.size(); j++)
            out_coords += mfa->geom().tmesh.tensor_prods[j].ctrl_pts.rows() *
                mfa->geom().tmesh.tensor_prods[j].ctrl_pts.cols();
        for (auto j = 0; j < mfa->geom().tmesh.all_knots.size(); j++)
            out_coords += mfa->geom().tmesh.all_knots[j].size();
        for (auto i = 0; i < mfa->nvars(); i++)
        {
            for (auto j = 0; j < mfa->var(i).tmesh.tensor_prods.size(); j++)
                out_coords += mfa->var(i).tmesh.tensor_prods[j].ctrl_pts.rows() *
                    mfa->var(i).tmesh.tensor_prods[j].ctrl_pts.cols();
            for (auto j = 0; j < mfa->var(i).tmesh.all_knots.size(); j++)
                out_coords += mfa->var(i).tmesh.all_knots[j].size();
        }
        return(in_coords / out_coords);
    }

    //     DEPRECATE
//     // compute compression ratio
//     float compute_compression()
//     {
//         // TODO: hard-coded for one tensor product
//         float in_coords = (input->npts) * (input->pt_dim);
//         float out_coords = geometry.mfa_data->tmesh.tensor_prods[0].ctrl_pts.rows() *
//             geometry.mfa_data->tmesh.tensor_prods[0].ctrl_pts.cols();
//         for (auto j = 0; j < geometry.mfa_data->tmesh.all_knots.size(); j++)
//             out_coords += geometry.mfa_data->tmesh.all_knots[j].size();
//         for (auto i = 0; i < vars.size(); i++)
//         {
//             out_coords += (vars[i].mfa_data->tmesh.tensor_prods[0].ctrl_pts.rows() *
//                     vars[i].mfa_data->tmesh.tensor_prods[0].ctrl_pts.cols());
//             for (auto j = 0; j < vars[i].mfa_data->tmesh.all_knots.size(); j++)
//                 out_coords += vars[i].mfa_data->tmesh.all_knots[j].size();
//         }
//         return(in_coords / out_coords);
//     }

    //  debug: print control points and weights in all tensor products of a tmesh
    void print_ctrl_weights(mfa::Tmesh<T>& tmesh)
    {
        for (auto i = 0; i < tmesh.tensor_prods.size(); i++)
        {
            cerr << "tensor_prods[" << i << "]:\n" << endl;
            cerr << tmesh.tensor_prods[i].ctrl_pts.rows() <<
                " final control points\n" << tmesh.tensor_prods[i].ctrl_pts << endl;
            cerr << tmesh.tensor_prods[i].weights.size()  <<
                " final weights\n" << tmesh.tensor_prods[i].weights << endl;
        }
    }

    // debug: print knots in a tmesh
    void print_knots(mfa::Tmesh<T>& tmesh)
    {
        for (auto j = 0 ; j < tmesh.all_knots.size(); j++)
        {
            fprintf(stderr, "%ld knots[%d]: [ ", tmesh.all_knots[j].size(), j);
            for (auto k = 0; k < tmesh.all_knots[j].size(); k++)
                fprintf(stderr, "%.3lf ", tmesh.all_knots[j][k]);
            fprintf(stderr, " ]\n");
        }
        fprintf(stderr, "\n");
    }

    void print_deriv(const diy::Master::ProxyWithLink& cp)
    {
        fprintf(stderr, "gid = %d\n", cp.gid());
        cerr << "domain\n" << input->domain << endl;
        cerr << approx->npts << " derivatives\n" << approx->domain << endl;
        fprintf(stderr, "\n");
    }

    // write original and approximated data in raw format
    // only for one block (one file name used, ie, last block will overwrite earlier ones)
    void write_raw(const diy::Master::ProxyWithLink& cp)
    {
        int last = input->domain.cols() - 1;           // last column in domain points

        // write original points
        ofstream domain_outfile;
        domain_outfile.open("orig.raw", ios::binary);
        vector<T> out_domain(input->domain.rows());
        for (auto i = 0; i < input->domain.rows(); i++)
            out_domain[i] = input->domain(i, last);
        domain_outfile.write((char*)(&out_domain[0]), input->domain.rows() * sizeof(T));
        domain_outfile.close();

#if 0
        // debug: read back original points
        ifstream domain_infile;
        vector<T> in_domain(input->domain.rows());
        domain_infile.open("orig.raw", ios::binary);
        domain_infile.read((char*)(&in_domain[0]), input->domain.rows() * sizeof(T));
        domain_infile.close();
        for (auto i = 0; i < input->domain.rows(); i++)
            if (in_domain[i] != input->domain(i, last))
                fprintf(stderr, "Error writing raw data: original data does match writen/read back data\n");
#endif

        // write approximated points
        ofstream approx_outfile;
        approx_outfile.open("approx.raw", ios::binary);
        vector<T> out_approx(approx->npts);
        for (auto i = 0; i < approx->npts; i++)
            out_approx[i] = approx->domain(i, last);
        approx_outfile.write((char*)(&out_approx[0]), approx->npts * sizeof(T));
        approx_outfile.close();

#if 0
        // debug: read back original points
        ifstream approx_infile;
        vector<T> in_approx(approx->npts);
        approx_infile.open("approx.raw", ios::binary);
        approx_infile.read((char*)(&in_approx[0]), approx->npts * sizeof(T));
        approx_infile.close();
        for (auto i = 0; i < approx->npts; i++)
            if (in_approx[i] != approx->domain(i, last))
                fprintf(stderr, "Error writing raw data: approximated data does match writen/read back data\n");
#endif
    }

//#define BLEND_VERBOSE

    void DecodeRequestGrid(int verbose,
                           MatrixX<T> &localBlock) // local block is in a grid, in lexicographic order
    {
        VectorXi &ndpts = input->ndom_pts;

        VectorX<T> min_bd(dom_dim);
        VectorX<T> max_bd(dom_dim);
        size_t cs = 1;
        for (int i = 0; i < dom_dim; i++) {
            min_bd(i) = input->domain(0, i);
            max_bd(i) = input->domain(cs * (ndpts(i) - 1), i);
            cs *= ndpts(i);
        }
        if (verbose) {
            std::cerr << " actual bounds, for full local decoded domain:\n";
            for (int i = 0; i < dom_dim; i++) {
                std::cerr << "   " << min_bd(i) << " " << max_bd(i) << "\n";
            }

            std:: cout << localBlock << "\n";
        }
        // determine the max and min extension of localBlock, and dimensions of localBlock
        // it is an inverse problem, from the block extensions, find out the maximum and index
        //
        VectorXi grid_counts(dom_dim);
        VectorX<T> min_block(dom_dim);
        VectorX<T> max_block(dom_dim);
        VectorXi index_max(dom_dim);
        for (int i = 0; i < dom_dim; i++)
        {
            min_block[i] = localBlock(0, i);
            max_block[i] = localBlock(0, i); // this will be modified in the next loop
        }
        int num_rows = (int) localBlock.rows();
        for (int k = 0; k < num_rows; k++)
        {
            for (int i = 0; i < dom_dim; i++)
            {
                if (max_block[i] < localBlock(k,i))
                {
                    max_block[i] = localBlock(k,i);
                    index_max[i] = k; // first index to this extent !
                }
            }
        }
        cs = 1;
        for (int i = 0; i < dom_dim; i++)
        {
            grid_counts[i] = index_max[i] / cs + 1;
            cs = cs * grid_counts[i];
        }
        // assert product (grid_counts[i] ) == num_rows ?
        //use domain parametrization logic to compute actual par for block within bounds,
        // to be able to use DecodeAtGrid
        VectorX<T> p_min(dom_dim); // initially, these are for the whole decoded domain
        VectorX<T> p_max(dom_dim);
#ifdef CURVE_PARAMS
                 // it is not possible
#else
        // use same logic as DomainParams(domain_, params);          // params spaced according to domain spacing
        // compute the par min and par max from bounds and core (core is embedded in bounds)

        for (int i = 0; i < dom_dim; i++) {
            T diff = bounds_maxs[i] - bounds_mins[i];
            p_min[i] = (min_block[i] - bounds_mins[i]) / diff;
            p_max[i] = (max_block[i] - bounds_mins[i]) / diff;
        }
#endif
        // geometry first
        /*mfa->DecodeAtGrid(*geometry.mfa_data, 0, dom_dim - 1, p_min, p_max,
                grid_counts, griddom); */
        for (size_t j = 0; j < mfa->nvars(); j++) {
            mfa->DecodeAtGrid(mfa->var(j), p_min, p_max, grid_counts, localBlock);
        }
    }
//#define BLEND_VERBOSE
    void decode_core_ures(const diy::Master::ProxyWithLink &cp,
            std::vector<int> &resolutions) // we already know the core, decode at desired resolution, without any requests sent to neighbors
            {
        // dom_dim  actual geo dimension for output points
        int tot_core_pts = 1;

        ndom_outpts.resize(dom_dim);
        for (int i = 0; i < dom_dim; i++) {
            ndom_outpts[i] = resolutions[map_dir[i]];
            tot_core_pts *= resolutions[map_dir[i]];
        }
        // get local core bounds are decided already:  core_mins, core_maxs
        blend.resize(tot_core_pts, dom_dim + 1); // these are just the evaluation points (core domain pts)

        // deltas in each direction
        VectorX<T> d(dom_dim);
        // starting point in each dimension will be core_mins
        int nghost_pts;           // number of ghost points in current dimension
        for (int i = 0; i < dom_dim; i++)
            d(i) = (core_maxs(i) - core_mins(i)) / (ndom_outpts(i) - 1);

        // assign values to the domain (geometry)
        int cs = 1;                        // stride of a coordinate in this dim
        T eps = 1.0e-10;                  // floating point roundoff error
        for (int i = 0; i < dom_dim; i++)  // all dimensions in the domain
                {
            int k = 0;
            int co = 0;            // j index of start of a new coordinate value
            for (int j = 0; j < tot_core_pts; j++) {
                if (core_mins(i) + k * d(i) > core_maxs(i) + eps)
                    k = 0;
                blend(j, i) = core_mins(i) + k * d(i);
                if (j + 1 - co >= cs) {
                    k++;
                    co = j + 1;
                }
            }
            cs *= ndom_outpts(i);
        }
        // now decode at resolution
        this->DecodeRequestGrid(0, blend); // so just compute the value, without any further blending
        //mfa->DecodeDomain(*(vars[0].mfa_data), 0, blend, dom_dim, dom_dim, false);
#ifdef BLEND_VERBOSE
       cerr << " block: " << cp.gid() << " decode no blend:" << blend << "\n";
#endif
    }

#undef BLEND_VERBOSE
    // each block will store the relevant neighbors' encroachment in its core domain
    // will be used to see if a point should be decoded by the neighbor block
    // we can also compute the bounds of the patch that will be sent to each neighbor
    void compute_neighbor_overlaps(const diy::Master::ProxyWithLink &cp) {
        RCLink<T> *l = static_cast<RCLink<T>*>(cp.link());
        // compute the overlaps of the neighbors outside main loop, and keep them in a vector
        Bounds<T> neigh_overlaps(dom_dim); // neighbor overlaps have the same dimension as domain , core
        Bounds<T> overNeigh(dom_dim); // over neigh core have the same dimension as domain , core
                                      // these bounds are symmetric to neigh overlaps
                                      // these cover the core of neigh [k] , and will be decoded and sent for blending
                                      // local grid is decided by that core [k] and resolutions
        diy::Direction dirc(dom_dim, 0); // neighbor direction
        diy::BlockID bid;  // this is stored for debugging purposes
        // if negative direction, correct min, if positive, correct max
        for (auto k = 0; k < l->size(); k++) {
            // start from current core, and restrict it
            //neigh_overlaps = l->bounds(k);
            dirc = l->direction(k);
            bid = l->target(k);
            //Bounds<T> boundsk = l->bounds(k);
            //Bounds<T> corek = l->core(k);
            for (int di = 0; di < dom_dim; di++) {
                int gdir = map_dir[di];
                // start with full core, then restrict it based on direction
                neigh_overlaps.max[di] = core_maxs[di];
                neigh_overlaps.min[di] = core_mins[di];
                // start with full bounds, then restrict to neigh core, based on direction
                overNeigh.max[di] = bounds_maxs[di];
                overNeigh.min[di] = bounds_mins[di];
                if (dirc[gdir] < 0) {
                    // it is reciprocity; if my bounds extended, neighbor bounds extended too in opposite direction
                    neigh_overlaps.max[di] = core_mins[di]
                            + (core_mins[di] - bounds_mins[di]);
                    overNeigh.max[di] = core_mins[di];
                }
                else if (dirc[gdir] > 0) {
                    neigh_overlaps.min[di] = core_maxs[di]
                            - (bounds_maxs[di] - core_maxs[di]);
                    overNeigh.min[di] = core_maxs[di];
                }
                else // if (dirc[gdir] == 0 )
                {
                    overNeigh.max[di] = core_maxs[di];
                    overNeigh.min[di] = core_mins[di];

                }
#ifdef BLEND_VERBOSE
          cerr << "overlaps block:" << cp.gid() << " nb:" << bid.gid << " di:" << di << " min:" <<
              neigh_overlaps.min[di] << " max:"<< neigh_overlaps.max[di] <<"\n";

#endif
            }
            neighOverlaps.push_back(neigh_overlaps);
            overNeighCore.push_back(overNeigh);
        }

    }

    void affected_neighbors(vector<T> &dom_pt, T &eps, vector<int> &dests) {
        Bounds<T> overlap(dom_dim); // neighbor overlap
        // for all neighbors of this block
        for (size_t i = 0; i < neighOverlaps.size(); i++) {
            overlap = neighOverlaps[i];
            bool insideDomain = true;
            for (int j = 0; j < dom_dim; j++) {
                T val = dom_pt[j];
                if (val + eps <= overlap.min[j]
                        || val - eps >= overlap.max[j]) {
                    insideDomain = false;
                    break;
                }
            }
            if (insideDomain)
                dests.push_back(i);
        }

    }

    // receive requests , decode, and send back to the requester only the decoded value
    // entire block was already encoded, we just need to decode at a grid of points
    void decode_patches(const diy::Master::ProxyWithLink &cp,
                        std::vector<int> &resolutions)
    {
        T eps = 1.0e-10;
        int verbose = 0;
        diy::Direction dirc(dom_dim, 0); // neighbor direction

        RCLink<T> *l = static_cast<RCLink<T>*>(cp.link());
        for (auto k = 0; k < l->size(); k++) {
            diy::BlockID bid = l->target(k);

            // alternatively, build local block for neighbor k using its core, global resolutions and overlap over core
            MatrixX<T> localBlockOverCoreK;
            Bounds<T>  domainK = overNeighCore[k]; //
            auto  coreK = l->core(k);
            dirc = l->direction(k);
            VectorXi counts(dom_dim);      // how many core points in overNeigh, per direction
            int  sizeBlock = 1;
            Bounds<T> resol(dom_dim);
            VectorX<T> d(dom_dim);                 // step in domain points in each dimension
            for (int j = 0; j < dom_dim; j++)
            {
                int dirGlob = map_dir[j];
                T deltaCore = coreK.max[dirGlob] - coreK.min[dirGlob];
                d(j) = deltaCore / ( resolutions[j] - 1 );
                T deltaOver = domainK.max[j] - domainK.min[j];
                counts[j]  = (int)( (deltaOver + 2 * eps) / d(j) ) + 1; // to account for some roundoff errors at the boundaries

                if (dirc[dirGlob] < 0 )
                {
                    resol.min[j] = coreK.max[dirGlob] - (counts[j] - 1) * d(j);
                    resol.max[j] = coreK.max[dirGlob];
                }
                else if (dirc[dirGlob] > 0 )
                {
                    resol.min[j] = coreK.min[dirGlob];
                    resol.max[j] = coreK.min[dirGlob] + (counts[j] - 1) * d(j);
                }
                else // dirc[dirGlob] == 0
                {
                    resol.min[j] = coreK.min[dirGlob]; // full extents
                    resol.max[j] = coreK.max[dirGlob];
                    counts[j] = resolutions[j]; // there should be no question about it
                }
                sizeBlock *= counts[j];
            }


            localBlockOverCoreK.resize(sizeBlock, dom_dim + 1);
            //assert(sizeBlock == num_points_received);
            // we know the bounds (resol) and counts, we can form the grid
            // assign values to the domain (geometry)
            mfa::VolIterator vol_it(counts);
            // current index of domain point in each dim, initialized to 0s
            // flattened loop over all the points in a domain
            while (!vol_it.done())
            {
                int j = (int)vol_it.cur_iter();
                // compute geometry coordinates of domain point
                for (auto i = 0; i < dom_dim; i++)
                    localBlockOverCoreK(j, i) = resol.min[i] + vol_it.idx_dim(i) * d(i);

                vol_it.incr_iter();
            }

            //this->DecodeRequestGrid(0, localBlock); //
            this->DecodeRequestGrid(0, localBlockOverCoreK);
            //mfa->DecodeDomain(*(vars[0].mfa_data), 0, localBlock, dom_dim, dom_dim, false);
            // local block localBlock(:, dom_dim) will have the decoded  values for science variable
            vector<T> computed_values(sizeBlock);
            for (int k = 0; k < sizeBlock; k++) {
                computed_values[k] = localBlockOverCoreK(k, dom_dim); // just last dimension
            }
#ifdef BLEND_VERBOSE
        cerr << "   block gid " << cp.gid() << " for block " << bid.gid << "\n" << localBlockOverCoreK << endl;
#endif
            cp.enqueue(bid, computed_values);
        }
    }
    // receive requests , decode, and send back to the requester only the decoded value
    // entire block was already encoded, we just need to decode at a grid of points
    void decode_patches_discrete(const diy::Master::ProxyWithLink &cp,
                            std::vector<int> &resolutions)
    {
        T eps = 1.0e-10;
        int verbose = 0;
        diy::Direction dirc(dom_dim, 0); // neighbor direction

        RCLink<int> *l = static_cast<RCLink<int>*>(cp.link());
        for (auto k = 0; k < l->size(); k++) {
            diy::BlockID bid = l->target(k);

            // alternatively, build local block for neighbor k using its core, global resolutions and overlap over core
            MatrixX<T> localBlockOverCoreK;
            Bounds<T>  domainK = overNeighCore[k]; //
            auto  coreK = l->core(k);
            dirc = l->direction(k);
            VectorXi counts(dom_dim);      // how many core points in overNeigh, per direction
            int  sizeBlock = 1;
            Bounds<T> resol(dom_dim);
            VectorX<T> d(dom_dim);                 // step in domain points in each dimension
            for (int j = 0; j < dom_dim; j++)
            {
                int dirGlob = map_dir[j];
                T deltaCore = coreK.max[dirGlob] - coreK.min[dirGlob];
                d(j) = deltaCore / ( resolutions[j] - 1 );
                T deltaOver = domainK.max[j] - domainK.min[j];
                counts[j]  = (int)( (deltaOver + 2 * eps) / d(j) ) + 1; // to account for some roundoff errors at the boundaries

                if (dirc[dirGlob] < 0 )
                {
                    resol.min[j] = coreK.max[dirGlob] - (counts[j] - 1) * d(j);
                    resol.max[j] = coreK.max[dirGlob];
                }
                else if (dirc[dirGlob] > 0 )
                {
                    resol.min[j] = coreK.min[dirGlob];
                    resol.max[j] = coreK.min[dirGlob] + (counts[j] - 1) * d(j);
                }
                else // dirc[dirGlob] == 0
                {
                    resol.min[j] = coreK.min[dirGlob]; // full extents
                    resol.max[j] = coreK.max[dirGlob];
                    counts[j] = resolutions[j]; // there should be no question about it
                }
                sizeBlock *= counts[j];
            }


            localBlockOverCoreK.resize(sizeBlock, dom_dim + 1);
            //assert(sizeBlock == num_points_received);
            // we know the bounds (resol) and counts, we can form the grid
            // assign values to the domain (geometry)
            mfa::VolIterator vol_it(counts);
            // current index of domain point in each dim, initialized to 0s
            // flattened loop over all the points in a domain
            while (!vol_it.done())
            {
                int j = (int)vol_it.cur_iter();
                // compute geometry coordinates of domain point
                for (auto i = 0; i < dom_dim; i++)
                    localBlockOverCoreK(j, i) = resol.min[i] + vol_it.idx_dim(i) * d(i);

                vol_it.incr_iter();
            }

            //this->DecodeRequestGrid(0, localBlock); //
            this->DecodeRequestGrid(0, localBlockOverCoreK);
            //mfa->DecodeDomain(*(vars[0].mfa_data), 0, localBlock, dom_dim, dom_dim, false);
            // local block localBlock(:, dom_dim) will have the decoded  values for science variable
            vector<T> computed_values(sizeBlock);
            for (int k = 0; k < sizeBlock; k++) {
                computed_values[k] = localBlockOverCoreK(k, dom_dim); // just last dimension
            }
#ifdef BLEND_VERBOSE
        cerr << "   block gid " << cp.gid() << " for block " << bid.gid << "\n" << localBlockOverCoreK << endl;
#endif
            cp.enqueue(bid, computed_values);
        }
    }
    // receive requests , decode, and send back


    // blending in 1d, on interval -inf, [a , b] +inf, -inf < a < b < +inf
    T blend_1d(T t, T a, T b) {
        if (a == b)
            return 1; // assume t is also a and b
        assert(a < b);
        T x = (t - a) / (b - a);
        if (x < 0)
            return (T) 1.;
        if (x >= 1.)
            return (T) 0.;
        T tmp = x * x * x;
        tmp = 1 - tmp * (6 * x * x - 15 * x + 10);
        return tmp;
    }

    T alfa(T x, T border, T overlap) {
        return blend_1d(x, border - overlap, border + overlap);
    }

    T one_alfa(T x, T border, T overlap) {
        return 1 - alfa(x, border, overlap);
    }

    // now we know we need to receive requested values decoded by neighbors
    void recv_and_blend(const diy::Master::ProxyWithLink &cp) {

        T eps = 1.0e-10;
        RCLink<T> *l = static_cast<RCLink<T>*>(cp.link());
        vector<T> dom_pt(dom_dim);   // only domain coords of point
        std::vector<diy::Direction> directions;
        diy::BlockID bid;
        std::vector<int> in;
        cp.incoming(in);
        // receive all values ,
        std::map<int, std::vector<T>> incoming_values;
        std::map<int, int> pointer_in_incoming; // values will be processed one by one in the incoming vector,
        /// in the same order they were first processed, then sent, then received
        for (auto k = 0; k < l->size(); k++) {
            bid = l->target(k);
            cp.dequeue(bid.gid, incoming_values[bid.gid]);
            pointer_in_incoming[bid.gid] = 0;
        }
#ifdef BLEND_VERBOSE
      for (typename std::map<int, std::vector<T>>::iterator mit= incoming_values.begin(); mit!=incoming_values.end(); mit++)
      {
        int from_gid = mit->first;
        std::vector<T> & vectVals = mit->second;
        cerr << "on block" << cp.gid() << " receive from block " << from_gid << " values: \n";
        for (size_t j=0; j<vectVals.size(); j++)
        {
          cerr << " " << vectVals[j] ;
        }
        cerr << "\n";
      }
#endif

        Bounds<T> neigh_bounds(dom_dim); // neighbor block bounds
        diy::Direction dirc(dom_dim, 0); // neighbor direction

        for (auto i = 0; i < (size_t) blend.rows(); i++) {
            vector<int> dests;               // link neighbor targets (not gids)
            for (int j = 0; j < dom_dim; j++)
                dom_pt[j] = blend(i, j);
            // decide if a point is inside a neighbor domain;
            // return a list of destinations , similar to diy::near(*l, dom_pt, eps, insert_it, decomposer.domain);
            // this logic will be used at decoding too
            affected_neighbors(dom_pt, eps, dests);
            // skip blending if there are no affected neighbors! blend value will stay unmodified
            if (dests.empty())
                continue;
            diy::Direction genDir(dom_dim, 0);
            for (size_t k = 0; k < dests.size(); k++) {
                dirc = l->direction(dests[k]);
                for (int kk = 0; kk < dom_dim; kk++) {
                    int gdir = map_dir[kk];
                    if (dirc[gdir] != 0)
                        genDir[kk] = dirc[gdir];
                }
            }
            VectorX<T> border(dom_dim); // one variable only right now
            VectorX<T> myBlend(dom_dim); // should be 1 in general

            for (int kk = 0; kk < dom_dim; kk++) {
                if (genDir[kk] < 0) {
                    border(kk) = core_mins(kk);
                    myBlend(kk) = one_alfa(dom_pt[kk], border(kk),
                            overlaps(kk));
                } else if (genDir[kk] > 0) {
                    border(kk) = core_maxs(kk);
                    myBlend(kk) = alfa(dom_pt[kk], border(kk), overlaps(kk));
                } else {
                    border(kk) = -1.e20; // better to be unset, actually, or nan, to be sure we will not use it
                    myBlend(kk) = 1.;
                }
            }
            // now, finally, add contribution from each neighbor
            // my contribution is already set, by product of blends
            // other points will be the same product of blends, but with (1-blend) if the direction is
            // active, otherwise it is just 1

            // current point value, assume 1 variable
            // this needs to be blended
            T val = blend(i, dom_dim);
#ifdef BLEND_VERBOSE
        cerr << " blend point " << dom_pt[0] << " " << dom_pt[1] << " val:" << val << " factors: " << myBlend(0) <<
            " " << myBlend(1) << "\n";
        cerr << " destinations: " << dests.size() << " " ;
        for (size_t i1=0; i1<dests.size(); i1++)
          cerr << " " << dests[i1];
        cerr << "\n";
#endif
            for (int kk = 0; kk < dom_dim; kk++)
                val = val * myBlend(kk); // some might be 1, but all should be greater than 1/2., because local contribution is
                                         //  the biggest
                // now contributions from the neighbors
            for (size_t k = 0; k < dests.size(); k++) {
                T factor = 1;
                dirc = l->direction(dests[k]);
                for (int kk = 0; kk < dom_dim; kk++) {
                    int gdir = map_dir[kk];
                    if (dirc[gdir] != 0)
                        factor = factor * (1 - myBlend(kk)); // we know that the other blend should be less than .5
                    else
                        factor = factor * myBlend(kk); // similar with my own contribution,
                                                       // because it is too far from the boundary in this direction
                }
                // now process the value!!
                bid = l->target(dests[k]);
                int index = pointer_in_incoming[bid.gid];
                T incomingValue = incoming_values[bid.gid][index];
                pointer_in_incoming[bid.gid]++;
#ifdef BLEND_VERBOSE
          cerr << "  incoming value " << incomingValue << " from block: " << bid.gid << " index: " <<index <<
             " factor: " << factor << endl;
#endif
                val = val + incomingValue * factor;
            }
            // finally, new corrected value for ptc:
            blend(i, dom_dim) = val;
#ifdef BLEND_VERBOSE
        cerr << "     general direction: " << genDir << " border: " << border.transpose() <<
                         "  blending: "<< myBlend.transpose() << " new val:" << val << endl;
#endif
        }

    }
    void print_brief_block(const diy::Master::ProxyWithLink &cp, bool error) // error was computed
            {

        for (auto i = 0; i < mfa->nvars(); i++) {
            cerr << "block gid = " << cp.gid() << " var: " << i
                    << " max error: " << max_errs[i] << "\n";
        }

        if (0 == cp.gid()) {
            for (auto i = 0; i < mfa->nvars(); i++) {
                cerr << " \n Maximum error:" << max_errs_reduce[2 * i]
                        << " at block id::" << max_errs_reduce[2 * i + 1]
                        << "\n";
            }
        }
    }
};

namespace mfa
{
    template<typename B>                        // B = block object
        void* create()          { return new B; }

    template<typename B>                        // B = block object
    void destroy(void* b)       { delete static_cast<B*>(b); }

    template<typename B, typename T>                // B = block object,  T = float or double
        void add(                                       // add the block to the decomposition
                int                 gid,                // block global id
                const Bounds<T>&    core,               // block bounds without any ghost added
                const Bounds<T>&    bounds,             // block bounds including any ghost region added
                const Bounds<T>&    domain,             // global data bounds
                const RCLink<T>&    link,               // neighborhood
                diy::Master&        master,             // diy master
                int                 dom_dim,            // domain dimensionality
                int                 pt_dim,             // point dimensionality
                T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
        {
            B*              b   = new B;
            RCLink<T>*      l   = new RCLink<T>(link);
            diy::Master&    m   = const_cast<diy::Master&>(master);
            m.add(gid, b, l);

            b->init_block(core, domain, dom_dim, pt_dim, ghost_factor);
        }

    template<typename B, typename T>                // B = block object,  T = float or double
        void save(
                const void*        b_,
                diy::BinaryBuffer& bb)
        {
            B* b = (B*)b_;

            // top-level mfa data
            diy::save(bb, b->dom_dim);
            diy::save(bb, b->pt_dim);

            // block bounds
            diy::save(bb, b->bounds_mins);
            diy::save(bb, b->bounds_maxs);
            diy::save(bb, b->core_mins);
            diy::save(bb, b->core_maxs);

            // TODO: don't save data sets in practice
            diy::save(bb, b->input);
            diy::save(bb, b->approx);
            diy::save(bb, b->errs);

            // save mfa
            diy::save(bb, b->mfa);

            // output for blending
            diy::save(bb, b->ndom_outpts);
            diy::save(bb, b->blend);
        }

    template<typename B, typename T>                // B = block object, T = float or double
        void load(
                void*              b_,
                diy::BinaryBuffer& bb)
        {
            B* b = (B*)b_;

            // top-level mfa data
            diy::load(bb, b->dom_dim);
            diy::load(bb, b->pt_dim);

            // block bounds
            diy::load(bb, b->bounds_mins);
            diy::load(bb, b->bounds_maxs);
            diy::load(bb, b->core_mins);
            diy::load(bb, b->core_maxs);

            // TODO: Don't load data sets in practice
            diy::load(bb, b->input);
            diy::load(bb, b->approx);
            diy::load(bb, b->errs);

            // load mfa
            diy::load(bb, b->mfa);

            if (b->pt_dim != b->mfa->pt_dim)
                cerr << "WARNING: Block::pt_dim and MFA::pt_dim do not match!" << endl;

            // output for blending
            diy::load(bb, b->ndom_outpts);
            diy::load(bb, b->blend);
        }
}                       // namespace

namespace diy
{
    template <typename T>
        struct Serialization<MatrixX<T>>
        {
            static
                void save(diy::BinaryBuffer& bb, const MatrixX<T>& m)
                {
                    diy::save(bb, m.rows());
                    diy::save(bb, m.cols());
                    for (size_t i = 0; i < m.rows(); ++i)
                        for (size_t j = 0; j < m.cols(); ++j)
                            diy::save(bb, m(i, j));
                }
            static
                void load(diy::BinaryBuffer& bb, MatrixX<T>& m)
                {
                    Index rows, cols;
                    diy::load(bb, rows);
                    diy::load(bb, cols);
                    m.resize(rows, cols);
                    for (size_t i = 0; i < m.rows(); ++i)
                        for (size_t j = 0; j < m.cols(); ++j)
                            diy::load(bb, m(i, j));
                }
        };

    template <typename T>
        struct Serialization<VectorX<T>>
        {
            static
                void save(diy::BinaryBuffer& bb, const VectorX<T>& v)
                {
                    diy::save(bb, v.size());
                    for (size_t i = 0; i < v.size(); ++i)
                        diy::save(bb, v(i));
                }
            static
                void load(diy::BinaryBuffer& bb, VectorX<T>& v)
                {
                    Index size;
                    diy::load(bb, size);
                    v.resize(size);
                    for (size_t i = 0; i < size; ++i)
                        diy::load(bb, v(i));
                }
        };

    template<>
        struct Serialization<VectorXi>
        {
            static
                void save(diy::BinaryBuffer& bb, const VectorXi& v)
                {
                    diy::save(bb, v.size());
                    for (size_t i = 0; i < v.size(); ++i)
                        diy::save(bb, v(i));
                }
            static
                void load(diy::BinaryBuffer& bb, VectorXi& v)
                {
                    Index size;
                    diy::load(bb, size);
                    v.resize(size);
                    for (size_t i = 0; i < size; ++i)
                        diy::load(bb, v.data()[i]);
                }
        };

    template <typename T>
        struct Serialization<mfa::Param<T>>
        {
            static
                void save(diy::BinaryBuffer& bb, const mfa::Param<T>& p)
                {
                    diy::save(bb, p.structured);
                    diy::save(bb, p.dom_dim);
                    diy::save(bb, p.ndom_pts);

                    if (p.structured)  // save param_grid
                    {
                        for (int k = 0; k < p.dom_dim; k++)
                        {
                            diy::save(bb, p.param_grid[k]);
                        }
                    }
                    else    // save param_list
                    {
                        diy::save(bb, p.param_list);
                    }
                }
            static
                void load(diy::BinaryBuffer& bb, mfa::Param<T>& p)
                {
                    diy::load(bb, p.structured);
                    diy::load(bb, p.dom_dim);
                    diy::load(bb, p.ndom_pts);

                    if (p.structured)
                    {
                        p.param_grid.resize(p.dom_dim);
                        for (int k = 0; k < p.dom_dim; k++)
                        {
                            diy::load(bb, p.param_grid[k]);
                        }
                    }
                    else
                    {
                        diy::load(bb, p.param_list);
                    }
                }
        };

    template <typename T>
        struct Serialization<mfa::PointSet<T>*>
        {
            static
                void save(diy::BinaryBuffer& bb, mfa::PointSet<T>* const & ps)
                {
                    if (ps == nullptr)
                        diy::save(bb, (int) 0);
                    else
                    {
                        diy::save(bb, (int) 1); // indicate that there is a nonempty PointSet to load

                        diy::save(bb, ps->dom_dim);
                        diy::save(bb, ps->model_dims());
                        diy::save(bb, ps->npts);
                        diy::save(bb, ps->ndom_pts);

                        if (ps->params.get() == nullptr)
                        {
                            diy::save(bb, false);
                        }
                        else
                        {
                            diy::save(bb, true);
                            diy::save(bb, *(ps->params));
                        }
                        
                        diy::save(bb, ps->dom_mins);
                        diy::save(bb, ps->dom_maxs);

                        diy::save(bb, ps->domain);
                    }
                }
            static
                // NOTE: we use reference to pointer so the pointer is updated after a call to new
                void load(diy::BinaryBuffer& bb, mfa::PointSet<T>*& ps) 
                {
                    int load_flag;
                    diy::load(bb, load_flag);

                    if (load_flag == 0)
                    {
                        ps = nullptr;
                    }
                    else
                    {
                        int dom_dim, npts = 0;
                        VectorXi model_dims, ndom_pts;
                        VectorX<T> dom_mins, dom_maxs;

                        // Basic domain structure
                        diy::load(bb, dom_dim);
                        diy::load(bb, model_dims);
                        diy::load(bb, npts);
                        diy::load(bb, ndom_pts);
                        ps = new mfa::PointSet<T>(dom_dim, model_dims, npts, ndom_pts);

                        // Param info. We need to handle cases where Params may or may not be initialized
                        bool load_param;
                        diy::load(bb, load_param);
                        if (load_param)
                        {
                            ps->params = make_shared<mfa::Param<T>>();
                            diy::load(bb, *(ps->params));
                        }
                        diy::load(bb, dom_mins);
                        diy::load(bb, dom_maxs);

                        // Point info
                        diy::load(bb, ps->domain);
                    }   
                }
        };

        // Template for pointers to MFA_Data so that we can read/write uninitialized models
        template <typename T>
        struct Serialization<mfa::MFA_Data<T>*>
        {
            static
                void save(diy::BinaryBuffer& bb, mfa::MFA_Data<T>* const & md)
                {
                    if (md == nullptr)
                        diy::save(bb, (int) 0);
                    else
                    {
                        diy::save(bb, (int) 1); // indicate that there is a MFA_Data to load

                        diy::save(bb, md->p);
                        diy::save(bb, md->tmesh.tensor_prods.size());
                        diy::save(bb, md->min_dim);
                        diy::save(bb, md->max_dim);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.nctrl_pts);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.ctrl_pts);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.weights);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.knot_mins);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.knot_maxs);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::save(bb, t.level);
                        diy::save(bb, md->tmesh.all_knots);
                        diy::save(bb, md->tmesh.all_knot_levels);
                    }
                }
            static
                // NOTE: we use reference to pointer so the pointer is updated after a call to new
                void load(diy::BinaryBuffer& bb, mfa::MFA_Data<T>*& md) 
                {
                    int load_flag;
                    diy::load(bb, load_flag);

                    if (load_flag == 0)
                    {
                        md = nullptr;
                    }
                    else
                    {
                        VectorXi    p;                  // degree of the mfa
                        size_t      ntensor_prods;      // number of tensor products in the tmesh
                        int         min_dim;
                        int         max_dim;

                        // geometry
                        diy::load(bb, p);
                        diy::load(bb, ntensor_prods);
                        diy::load(bb, min_dim);
                        diy::load(bb, max_dim);
                        md = new mfa::MFA_Data<T>(p, ntensor_prods, min_dim, max_dim);

                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.nctrl_pts);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.ctrl_pts);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.weights);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.knot_mins);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.knot_maxs);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            diy::load(bb, t.level);
                        diy::load(bb, md->tmesh.all_knots);
                        diy::load(bb, md->tmesh.all_knot_levels);
                        for (TensorProduct<T>& t: md->tmesh.tensor_prods)
                            md->tmesh.tensor_knot_idxs(t);
                    }   
                }
        };

        template <typename T>
        struct Serialization<mfa::MFA<T>*>
        {
            static
                void save(diy::BinaryBuffer& bb, mfa::MFA<T>* const & mfa)
                {
                    if (mfa == nullptr)
                    {
                        diy::save(bb, (int) 0);
                    }
                    else
                    {
                        diy::save(bb, (int) 1);

                        diy::save(bb, mfa->dom_dim);
                        diy::save(bb, mfa->verbose);
                        diy::save(bb, mfa->nvars());
                        diy::save(bb, mfa->geometry.get());
                        for (int i = 0; i < mfa->nvars(); i++)
                        {
                            diy::save(bb, mfa->vars[i].get());
                        }
                    }
                }
            static
                void load(diy::BinaryBuffer& bb, mfa::MFA<T>*& mfa)
                {
                    int load_flag = 0;
                    diy::load(bb, load_flag);
                    
                    if (load_flag == 0)
                    {
                        mfa = nullptr;
                    }
                    else
                    {
                        int dom_dim = 0;
                        int verbose = 0;
                        int nvars = 0;
                        mfa::MFA_Data<T>* geometry;
                        vector<mfa::MFA_Data<T>*> vars;

                        diy::load(bb, dom_dim);
                        diy::load(bb, verbose);
                        diy::load(bb, nvars);
                        vars.resize(nvars);
                        diy::load(bb, geometry);
                        for (int i = 0; i < vars.size(); i++)
                        {
                            diy::load(bb, vars[i]);
                        }

                        // n.b. this constructor takes ownership of MFA_Data pointers
                        mfa = new mfa::MFA<T>(dom_dim, verbose, geometry, vars);

                        // zero-out raw pointers; the resources are now managed by MFA
                        // not strictly necessary, but just good practice
                        geometry = nullptr;
                        for (int i = 0; i < vars.size(); i++)
                        {
                            vars[i] = nullptr;
                        }
                    }
                }
        };
}                       // namespace
