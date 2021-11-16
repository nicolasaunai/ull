
#include <array>
#include <cstddef>
#include <random>
#include <iostream>
#include <stdexcept>

template<size_t dim>
struct Particle
{
    static_assert(dim > 0 and dim < 4, "Only dimensions 1,2,3 are supported.");
    static const size_t dimension = dim;

    double weight;
    double charge;

    std::array<int, dim> iCell    = std::array<int, dim>{};
    std::array<double, dim> delta = std::array<double, dim>{};
    std::array<double, 3> v       = std::array<double, 3>{};

    double Ex = 0, Ey = 0, Ez = 0;
    double Bx = 0, By = 0, Bz = 0;
};




template<std::size_t bucket_size, typename T>
class ull
{
    class iterator : public std::iterator<std::forward_iterator_tag, T>
    {
    public:
        iterator(ull* u, std::size_t curr_bucket = 0, std::size_t curr_pos = 0)
            : curr_bucket_{curr_bucket}
            , curr_pos_{curr_pos}
            , ull_{u}
        {
        }

    public:
        const T* operator*() { return ull_->buckets_[curr_bucket_][curr_pos_]; }

        iterator operator++()
        {
            curr_pos_++;
            if (curr_pos_ == bucket_size)
            {
                curr_bucket_++;
                curr_pos_ = 0;
            }
            // std::cout << "++ : " << curr_bucket_ << " " << curr_pos_ << "\n";
            return *this;
        }


        bool operator!=(iterator const& other) const
        {
            return (other.curr_bucket_ != curr_bucket_ or other.curr_pos_ != curr_pos_)
                   or ull_ != other.ull_;
        }


    private:
        std::size_t curr_bucket_ = 0, curr_pos_ = 0;
        ull<bucket_size, T>* ull_;
    };

public:
    ull()
        : buckets_(1)
    {
    }

    void add(T const& t)
    {
        if (curr != bucket_size)
        {
            buckets_[bucket_idx][curr++] = &t;
        }
        else
        {
            buckets_.push_back(bucket{});
            bucket_idx++;
            curr = 0;
            add(t);
        }
    }

    std::size_t total() const { return bucket_size * (bucket_idx) + curr; }

    auto begin() { return iterator{this}; }
    auto end()
    {
        auto it = iterator{this, bucket_idx, curr};
        return it;
    }



private:
    using bucket           = std::array<const T*, bucket_size>;
    std::size_t bucket_idx = 0;
    std::size_t curr       = 0;
    std::vector<bucket> buckets_;
};




template<std::size_t dim>
struct Box
{
    class iterator
    {
    public:
        iterator(Box* b, std::array<std::size_t, dim> index = std::array<std::size_t, dim>{})
            : box_{b}
            , index_{index}
        {
        }

    public:
        std::array<std::size_t, dim> operator*() { return index_; }

        iterator operator++()
        {
            index_[dim - 1]++;
            for (auto idim = dim - 1; idim > 0; idim--)
            {
                if (index_[idim] == box_->upper[idim] and index_[idim - 1] != box_->upper[idim - 1])
                {
                    index_[idim] = 0;
                    index_[idim - 1]++;
                }
            }
            return *this;
        }


        bool operator!=(iterator const& other) const
        {
            return box_ != other.box_ or index_ != other.index_;
        }


    private:
        Box* box_;
        std::array<std::size_t, dim> index_;
    };


    auto begin() { return iterator{this, lower}; }
    auto end() { return iterator{this, upper}; }

    std::array<std::size_t, dim> lower;
    std::array<std::size_t, dim> upper;
};




template<std::size_t dim, typename T>
class grid
{
public:
    grid(std::size_t nx, std::size_t ny)
        : nx_{ny}
        , ny_{ny}
    {
        ulls_.resize(nx_ * ny_);
    }

    void addToCell(std::array<int, dim> cell, T const& obj)
    {
        ulls_[cell[1] + ny_ * cell[0]].add(obj);
    }

    std::size_t total(std::size_t ix, std::size_t iy)
    {
        auto icell = iy + ix * ny_;
        return ulls_[icell].total();
    }


    auto select(std::vector<Particle<dim>> const& particles, Box<dim> const& box)
    {
        // count all particles
        std::size_t nbrTot = 0;
        for (auto ix = box.lower[0]; ix <= box.upper[0]; ++ix)
        {
            for (auto iy = box.lower[1]; iy <= box.upper[1]; ++iy)
            {
                nbrTot += total(ix, iy);
            }
        }

        std::vector<Particle<dim>> selection(nbrTot);

        for (auto ix = box.lower[0]; ix <= box.upper[0]; ++ix)
        {
            for (auto iy = box.lower[1]; iy <= box.upper[1]; ++iy)
            {
                for (auto const* p : ulls_[iy + ix * ny_])
                {
                    // selection[ip++] = *p;
                }
            }
        }
    }


    ull<200, T>& cell(std::size_t ix, std::size_t iy)
    {
        auto& c = ulls_[iy + ix * ny_];
        return c;
    }

    ull<200, T>& cell(std::array<std::size_t, dim> cell)
    {
        auto& c = ulls_[cell[1] + cell[0] * ny_];
        return c;
    }

private:
    std::size_t nx_, ny_;
    std::vector<ull<200, T>> ulls_;
};




int main()
{ //
    std::size_t nx = 100, ny = 150;
    std::size_t nppc   = 4;
    constexpr auto dim = 2u;
    std::vector<Particle<dim>> particles;
    particles.reserve(nx * ny * nppc);

    for (auto ix = 0u; ix < nx; ++ix)
    {
        for (auto iy = 0u; iy < ny; ++iy)
        {
            for (auto ip = 0; ip < nppc; ++ip)
            {
                auto p     = Particle<dim>{};
                p.iCell[0] = ix;
                p.iCell[1] = iy;
                particles.push_back(std::move(p));
            }
        }
    }


    std::cout << "making the grid...\n";
    grid<dim, Particle<dim>> myGrid(nx, ny);

    std::cout << "pushing...\n";
    // pretend pushing
    for (auto ip = 0; ip < nppc * nx * ny; ++ip)
    {
        myGrid.addToCell(particles[ip].iCell, particles[ip]);
    }


    std::cout << "testing nbr of particles per cell registered\n";
    for (auto ix = 0u; ix < nx; ++ix)
    {
        for (auto iy = 0u; iy < ny; ++iy)
        {
            if (myGrid.total(ix, iy) != nppc)
            {
                std::runtime_error("invalid number of particles in cell ");
            }
        }
    }
    std::cout << "nbr of particles ok.\n";


    auto& cell  = myGrid.cell(10, 12);
    auto it_beg = cell.begin();
    auto it_end = cell.end();

    Box<dim> selection_box;
    selection_box.lower[0] = 20;
    selection_box.lower[1] = 15;
    selection_box.upper[0] = 22;
    selection_box.upper[1] = 17;
    std::cout << "cells in box : [(" << selection_box.lower[0] << "," << selection_box.lower[1]
              << "),(" << selection_box.upper[0] << "," << selection_box.upper[1] << ")]"
              << "\n";
    for (auto c : selection_box)
    {
        std::cout << " " << c[0] << " " << c[1] << "\n";
    }
    std::cout << "now listing particles for each cell\n";
    for (auto const& c : selection_box)
    {
        auto plist = myGrid.cell(c);
        std::cout << "cell : " << c[0] << " " << c[1] << "\n";
        for (Particle<dim> const* p : plist)
        {
            std::cout << "particle : " << p->iCell[0] << ", " << p->iCell[1] << "\n";
        }
    }
    Box<dim> selection_box2;
    selection_box.lower[0] = 20;
    selection_box.lower[1] = 15;
    selection_box.upper[0] = 38;
    selection_box.upper[1] = 44;
    std::cout << "particle in cell 10,12 : " << (*it_beg)->iCell[0] << "\n";




    return 0;
}

