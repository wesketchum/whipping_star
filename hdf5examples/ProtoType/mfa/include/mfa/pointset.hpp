//--------------------------------------------------------------
// mfa point set data structure
//
// Tom Peterka
// Argonne National Laboratory
// tpeterka@mcs.anl.gov
// 
// David Lenz
// Argonne National Laboratory
// dlenz@anl.gov
//--------------------------------------------------------------
#ifndef _POINTSET_HPP
#define _POINTSET_HPP

#include <memory>

namespace mfa
{
    template <typename T>
    struct PointSet
    {
        // Defined during construction
        int         dom_dim;
        int         pt_dim;
        int         npts;

        // Members that track how the columns of domain correspond to different models
        VectorXi    mdims;      // Dimensionality of each model (including geometry!)
        VectorXi    dim_mins;   // Start column index for each science variable model
        VectorXi    dim_maxs;   // End column index for each science variable model

        // List of points
        MatrixX<T>  domain;

        // Parameter values corresponding to domain
        shared_ptr<Param<T>>    params;
        VectorX<T>              dom_mins;
        VectorX<T>              dom_maxs;

        // Optional grid data
        bool        structured;
        VectorXi    ndom_pts;   // TODO: remove since this is contained in GridInfo?  
        GridInfo    g;

        // Generic constructor with only dimensionality and total points
        PointSet(
                size_t          dom_dim_,
                const VectorXi& mdims_,
                size_t          npts_) :
            dom_dim(dom_dim_),
            pt_dim(mdims_.sum()),
            npts(npts_),
            mdims(mdims_),
            structured(false)
        {
            domain.resize(npts, pt_dim);

            // Fill dim_mins/maxs
            dim_mins.resize(nvars());
            dim_maxs.resize(nvars());
            dim_mins[0] = geom_dim();
            dim_maxs[0] = dim_mins[0] + var_dim(0) - 1;
            for (int k = 1; k < nvars(); k++)
            {
                dim_mins[k] = dim_maxs[k-1] + 1;
                dim_maxs[k] = dim_mins[k] + var_dim(k) - 1;
            }

            validate();
        }   

        // Constructor for (possibly) grid-aligned data
        // If ndom_pts_ is empty, PointSet is unstructured
        PointSet(
                size_t          dom_dim_,
                const VectorXi& mdims_,
                size_t          npts_,
                const VectorXi& ndom_pts_) :
            PointSet(dom_dim_, mdims_, npts_)
        {
            add_grid(ndom_pts_);
        }

        // Constructor for a PointSet mapped from existing parameters
        //   N.B. this is useful when decoding at the same params as encoding,
        //        or when Params are constructed ahead of time
        PointSet(
                shared_ptr<Param<T>>    params_,
                const VectorXi&         mdims_) :
            dom_dim(params_->dom_dim),
            pt_dim(mdims_.sum()),
            npts(params_->npts()),
            mdims(mdims_),
            params(params_), 
            structured(params->structured)
        {
            domain.resize(npts, pt_dim);

            // Fill dim_mins/maxs
            dim_mins.resize(nvars());
            dim_maxs.resize(nvars());
            dim_mins[0] = geom_dim();
            dim_maxs[0] = dim_mins[0] + var_dim(0) - 1;
            for (int k = 1; k < nvars(); k++)
            {
                dim_mins[k] = dim_maxs[k-1] + 1;
                dim_maxs[k] = dim_mins[k] + var_dim(k) - 1;
            }

            if (params_->structured)
                add_grid(params_->ndom_pts);
        }

        void set_bounds(const VectorX<T>& mins_, const VectorX<T>& maxs_)
        {
            if ( (mins_.size() > 0 && mins_.size() != dom_dim) ||
                 (mins_.size() != maxs_.size()) )
            {
                cerr << "ERROR: Invalid bounds passed to PointSet" << endl;
                cerr << "  dom_mins: " << mins_ << endl;
                cerr << "  dom_maxs: " << maxs_ << endl;
                exit(1);
            }

            if (params != nullptr)
            {
                cerr << "Warning: Setting PointSet bounds after parameter initialization, this likely will have no effect" << endl;
            }

            dom_mins = mins_;
            dom_maxs = maxs_;
        }

        int nvars() const {return mdims.size() - 1;};

        int geom_dim() const {return mdims[0];}

        int var_dim(int k) const {return mdims[k+1];}

        int var_min(int k) const {return dim_mins[k];}

        int var_max(int k) const {return dim_maxs[k];}

        VectorXi model_dims() const {return mdims;}

        // Initialize parameter values from domain points and bounding box (dom_mins/maxs)
        //   N.B. If dom_mins/maxs are empty, they are set to be the min/max values of the
        //        domain points during Param construction.
        void init_params()
        {
            // If mins/maxs are NULL, set them to be the min/max coordinates of the pointset
            if (dom_mins.size() == 0 || dom_maxs.size() == 0)
            {
                dom_mins = domain.leftCols(dom_dim).colwise().minCoeff();
                dom_maxs = domain.leftCols(dom_dim).colwise().maxCoeff();
            }

            // n.b. A structured grid which has been rotated will still have its parameters computed correctly.
            //      dom mins/maxs are not used in the computation of structured parameters, so the parameters
            //      are computed to be the correct "rotated" grid
            params = make_shared<Param<T>>(dom_dim, dom_mins, dom_maxs, ndom_pts, domain, structured);
        }

        // Specify bounding box and initialize parameters simultaneously
        void init_params(const VectorX<T>& dom_mins_, const VectorX<T>& dom_maxs_)
        {
            set_bounds(dom_mins_, dom_maxs_);
            init_params();
        }

        void set_params(shared_ptr<Param<T>> params_)
        {
            if (params != nullptr)
            {
                cerr << "Warning: Overwriting existing parameters in PointSet" << endl;
            }

            if (!check_param_domain_agreement(*params_))
            {
                cerr << "ERROR: Attempted to add mismatched Params to PointSet" << endl;
                exit(1);
            }
 
            params = params_;

            // If Param is structured and *this does not yet have a grid structure, inherit the 
            // grid structure from Param
            if (params_->structured && !this->structured)
                add_grid(params_->ndom_pts);
        }


        // Checks that the Param object does not contradict existing members of PointSet
        //   N.B. If *this is structured, then Param object must be too.
        //        However, if *this is not structured, Param is allowed to be structured;
        //        this allows a new grid structure to be imposed from the Param grid structure.
        bool check_param_domain_agreement(const Param<T>& params_)
        {
            bool agreement =    (dom_dim == params_.dom_dim)
                            &&  (npts == params_.npts())
                            &&  (structured ? params_.structured == true : true)
                            &&  (structured ? ndom_pts == params_.ndom_pts : true)
                            ;

            return agreement;
        }

        // Add a grid structure to point set
        // Does nothing if ndom_pts_ is empty
        void add_grid(VectorXi ndom_pts_)
        {
            if (ndom_pts_.size() == 0)
            {
                return;
            }
            else
            {
                if (npts != ndom_pts_.prod())
                {
                    cerr << "ERROR: Invalid grid added to PointSet. Total points do not match." << endl;
                    cerr << "  npts = " << npts << endl;
                    cerr << "  ndom_pts = " << ndom_pts_ << endl;
                    exit(1);
                }

                structured = true;
                ndom_pts = ndom_pts_;
                g.init(dom_dim, ndom_pts);

                validate(); // Prints warning and aborts if invalid
            }
        }

        // Test that user-provided data meets basic sanity checks
        bool validate() const
        { 
            bool is_valid =     (dom_dim > 0)
                            &&  (geom_dim() > 0)
                            &&  (pt_dim > dom_dim)
                            &&  (pt_dim == domain.cols())
                            &&  (npts == domain.rows())
                            &&  (structured ? ndom_pts.size() == dom_dim : true)
                            &&  (structured ? ndom_pts.prod() == domain.rows() : true)
                            ;

            for (int k = 0; k < nvars(); k++)
            {
                is_valid = is_valid && (mdims[k+1] > 0);
                is_valid = is_valid && (dim_maxs[k] - dim_mins[k] + 1 == mdims[k+1]);
            }
            is_valid = is_valid && (dim_mins[0] == geom_dim()) && (dim_maxs[nvars()-1] == pt_dim-1);

            if (is_valid) return is_valid;
            else 
            {
                cerr << "ERROR: PointSet initialized with incompatible data" << endl;
                cerr << "  structured: " << boolalpha << structured << endl;
                cerr << "  dom_dim: " << dom_dim << ", geom_dim: " << geom_dim() << ",  pt_dim: " << pt_dim << endl;
                cerr << "  npts: " << npts << endl;
                cerr << "  ndom_pts: ";
                for (size_t k=0; k < ndom_pts.size(); k++) 
                    cerr << ndom_pts(k) << " ";
                cerr << endl;
                cerr << "  domain matrix dims: " << domain.rows() << " x " << domain.cols() << endl;

                cerr << "  nvars: " << nvars() << endl;
                cerr << "  model_dims:";
                for (int k = 0; k < nvars() + 1; k++)
                    cerr << mdims[k] << " ";
                cerr << endl;
                cerr << "  dim_mins:";
                for (int k = 0; k < nvars(); k++)
                    cerr << dim_mins[k] << " ";
                cerr << endl;
                cerr << "  dim_maxs:";
                for (int k = 0; k < nvars(); k++)
                    cerr << dim_maxs[k] << " ";
                cerr << endl;
                exit(1);
                
                return is_valid;
            }
        }

        bool is_same_layout(const PointSet<T>& ps, int verbose = 1) const
        {
            bool is_same =      (dom_dim    == ps.dom_dim)
                            &&  (pt_dim     == ps.pt_dim)
                            &&  (npts       == ps.npts)
                            &&  (structured == ps.structured)
                            &&  (mdims == ps.mdims);
            
            if (structured)
            {
                is_same = is_same && (ndom_pts == ps.ndom_pts);
            }

            if (is_same) return is_same;
            else
            {
                if (verbose)
                {
                    cerr << "Pair of PointSets do not have matching layout" << endl;
                    cerr << "  dom_dim    = " << dom_dim << ",\t" << ps.dom_dim << endl;
                    cerr << "  pt_dim     = " << pt_dim << ",\t" << ps.pt_dim << endl;
                    cerr << "  npts       = " << npts << ",\t" << ps.npts << endl;
                    cerr << "  structured = " << boolalpha << structured << ",\t" << ps.structured << endl;
                    if (structured || ps.structured)
                    {
                        cerr << "  ndom_pts: " << ndom_pts << ps.ndom_pts << endl;
                    }
                    cerr << "  model_dims = " << mdims << "\n" << ps.mdims << endl;
                }

                return is_same;
            }
        }

        void abs_diff(
            const   mfa::PointSet<T>& other,
                    mfa::PointSet<T>& diff,
                    int               verbose) const
        {
            if (!this->is_same_layout(other) || !this->is_same_layout(diff))
            {
                cerr << "ERROR: Incompatible PointSets in PointSet::abs_diff" << endl;
                exit(1);
            }

#ifdef MFA_SERIAL
            diff.domain.leftCols(geom_dim()) = this->domain.leftCols(geom_dim());
            diff.domain.rightCols(pt_dim-geom_dim()) = (this->domain.rightCols(pt_dim-geom_dim()) - other.domain.rightCols(pt_dim-geom_dim())).cwiseAbs();
#endif // MFA_SERIAL
#ifdef MFA_TBB
            parallel_for (size_t(0), (size_t)diff.npts, [&] (size_t i)
                {
                    for (auto j = 0; j < geom_dim(); j++)
                    {
                        diff.domain(i,j) = this->domain(i,j); // copy the geometric location of each point
                    }
                });

            parallel_for (size_t(0), (size_t)npts, [&] (size_t i)
                {
                    for (auto j = geom_dim(); j < pt_dim; j++)
                    {
                        diff.domain(i,j) = fabs(this->domain(i,j) - other.domain(i,j)); // compute distance between each science value
                    }
                });
#endif // MFA_TBB
        }

        PointSet(const PointSet&) = delete;
        PointSet(PointSet&&) = delete;
        PointSet& operator=(const PointSet&) = delete;
        PointSet& operator=(PointSet&&) = delete;

        // PointSet& operator=(PointSet&& other)
        // {
        //     swap(*this, other);
        //     return *this;
        // }

        // friend void swap(PointSet& first, PointSet& second)
        // {
        //     swap(first.dom_dim, second.dom_dim);
        //     swap(first.pt_dim, second.pt_dim);
        //     swap(first.structured, second.structured);
        //     first.ndom_pts.swap(second.ndom_pts);
        //     first.dom_mins.swap(second.dom_mins);
        //     first.dom_maxs.swap(second.dom_maxs);
        //     swap(first.g, second.g);
        //     first.domain.swap(second.domain);
        //     swap(first.params, second.params);
        // }

        class PtIterator
        {
            const bool  structured;
            size_t      lin_idx;
            mfa::VolIterator vol_it;
            const PointSet&  pset;

        public:
            PtIterator(const PointSet& pset_, size_t idx_) :
                structured(pset_.structured),
                lin_idx(idx_),
                vol_it(structured ? VolIterator(pset_.ndom_pts, idx_) : VolIterator()),
                pset(pset_)
            { }

            // prefix increment
            PtIterator operator++()
            {
                if(structured)
                    vol_it.incr_iter();
                
                lin_idx++;

                return *this;
            }


            bool operator!=(const PtIterator& other)
            {
                return (structured != other.structured) || (lin_idx != other.lin_idx);
            }

            bool operator==(const PtIterator& other)
            {
                return !(*this!=other);
            }

            // Full set of coordinates at this point
            void coords(VectorX<T>& coord_vec)
            {
                coord_vec = pset.domain.row(lin_idx);
            }

            // Geometric coordinates at this point
            void geom_coords(VectorX<T>& coord_vec)
            {
                coord_vec = pset.domain.row(lin_idx).head(geom_dim());
            }

            // Coordinates for variable k at this point
            void var_coords(int k, VectorX<T>& coord_vec)
            {
                coord_vec = pset.domain.row(lin_idx).segment(dim_mins[k], var_dim(k));
            }

            void coords(VectorX<T>& coord_vec, size_t min_dim, size_t max_dim)
            {
                coord_vec = pset.domain.block(idx(), min_dim, 1, max_dim - min_dim + 1).transpose();
            }

            void params(VectorX<T>& param_vec)
            {
                if(structured)
                    param_vec = pset.params->pt_params(vol_it);
                else
                    param_vec = pset.params->pt_params(lin_idx);
            }

            void ijk(VectorXi& ijk_vec)
            {
                if (!structured)
                {
                    cerr << "ERROR: No ijk values in PtIterator for unstructured input" << endl;
                    exit(1);
                }

                ijk_vec = vol_it.idx_dim();
            }

            int ijk(int k)
            {
                if (!structured)
                {
                    cerr << "ERROR: No ijk values in PtIterator for unstructured input" << endl;
                    exit(1);
                }

                return vol_it.idx_dim(k);
            }

            int idx()
            {
                return lin_idx;
            }
        };  // PtIterator

        PtIterator iterator(size_t idx) const
        {
            return PtIterator(*this, idx);
        }

        PtIterator begin() const
        {
            return PtIterator(*this, 0);
        }

        PtIterator end() const
        {
            return PtIterator(*this, npts);
        }

        void pt_coords(size_t idx, VectorX<T>& coord_vec) const
        {
            coord_vec = domain.row(idx);
        }

        void geom_coords(size_t idx, VectorX<T>& coord_vec) const
        {
            coord_vec = domain.row(idx).head(geom_dim());
        }

        void var_coords(size_t idx, size_t k, VectorX<T>& coord_vec) const
        {
            coord_vec = domain.row(idx).segment(var_min(k), var_dim(k));
        }

        void pt_params(size_t idx, VectorX<T>& param_vec) const
        {
            if(structured)
            {
                VectorXi ijk(dom_dim);
                g.idx2ijk(idx, ijk);
                param_vec = params->pt_params(ijk);
            }
            else
            {
                param_vec = params->pt_params(idx);   
            }
        }
    };
}   // namespace mfa

#endif // _POINTSET_HPP