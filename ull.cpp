#include <chrono>
#include <array>
#include <cstddef>
#include <random>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <fstream>

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



class Timer
{
public:
    Timer(std::size_t& save, bool print = false)
        : tstart_{std::chrono::high_resolution_clock::now()}
        , save_{save}
        , print_{print}
    {
    }

    ~Timer()
    {
        tstop_ = std::chrono::high_resolution_clock::now();
        std::size_t duration
            = std::chrono::duration_cast<std::chrono::microseconds>(tstop_ - tstart_).count();
        save_ = duration;
        if (print_)
            std::cout << duration << " us\n";
    }

private:
    std::chrono::high_resolution_clock::time_point tstart_, tstop_;
    std::size_t& save_;
    bool print_;
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
        const T* operator*()
        {
            // std::cout << curr_bucket_ << " " << curr_pos_ << "\n";
            return ull_->buckets_[curr_bucket_][curr_pos_];
        }

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
            // std::cout << "adding..." << bucket_idx << " " << curr << "\n";
            buckets_[bucket_idx][curr++] = &t;
        }
        else
        {
            buckets_.push_back(bucket{});
            bucket_idx++;
            curr = 0;
            add(t);
            std::cout << "buckets" << buckets_.size() << " " << buckets_.capacity() << "\n";
        }
    }

    std::size_t total() const { return bucket_size * (bucket_idx) + curr; }

    auto begin() { return iterator{this}; }
    const auto begin() const { return iterator{this}; }
    auto end()
    {
        if (curr != bucket_size)
        {
            auto it = iterator{this, bucket_idx, curr};
            return it;
        }
        else
        {
            auto it = iterator{this, bucket_idx + 1, 0};
            return it;
        }
    }
    const auto end() const
    {
        if (curr != bucket_size)
        {
            auto it = iterator{this, bucket_idx, curr};
            return it;
        }
        else
        {
            auto it = iterator{this, bucket_idx + 1, 0};
            return it;
        }
    }


    void empty()
    {
        bucket_idx = 0;
        curr       = 0;
    }

    bool is_empty() const { return bucket_idx == 0 and curr == 0; }

    std::size_t capacity() const { return buckets_.size() * bucket_size; }


private:
    using bucket           = std::array<const T*, bucket_size>;
    std::size_t bucket_idx = 0;
    std::size_t curr       = 0;
    std::vector<bucket> buckets_;
};




template<std::size_t dim>
struct Box
{
    Box(std::array<std::size_t, dim> lower_cell, std::array<std::size_t, dim> upper_cell)
        : lower{lower_cell}
        , upper{upper_cell}
    {
    }
    Box operator*(Box const& other) const
    {
        Box intersection{other};
        for (auto idim = 0u; idim < dim; ++idim)
        {
            intersection.lower[idim] = std::max(lower[idim], other.lower[idim]);
            intersection.upper[idim] = std::min(upper[idim], other.upper[idim]);
            if (intersection.lower[idim] > intersection.upper[idim])
                throw std::runtime_error("invalid intersection");
        }
        return intersection;
    }



    class iterator
    {
    public:
        iterator(Box const* b, std::array<std::size_t, dim> index = std::array<std::size_t, dim>{})
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
                if (index_[idim] == box_->upper[idim] + 1
                    and index_[idim - 1] < box_->upper[idim - 1])
                {
                    index_[idim] = box_->lower[idim];
                    index_[idim - 1]++;
                }
            }
            return *this;
        }
        iterator operator++() const
        {
            index_[dim - 1]++;
            for (auto idim = dim - 1; idim > 0; idim--)
            {
                if (index_[idim] == box_->upper[idim] + 1
                    and index_[idim - 1] < box_->upper[idim - 1])
                {
                    index_[idim] = box_->lower[idim];
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
        Box const* box_;
        std::array<std::size_t, dim> index_;
    };


    auto size() const { return (upper[1] - lower[1] + 1) * (upper[0] - lower[0] + 1); }
    std::array<std::size_t, dim> shape() const
    {
        return {upper[0] - lower[0] + 1, upper[1] - lower[1] + 1};
    }
    auto begin() { return iterator{this, lower}; }
    auto end() { return iterator{this, {upper[0], upper[1] + 1}}; }

    const auto begin() const { return iterator{this, lower}; }
    const auto end() const { return iterator{this, {upper[0], upper[1] + 1}}; }
    std::array<std::size_t, dim> lower;
    std::array<std::size_t, dim> upper;
};




template<std::size_t dim, typename T, std::size_t bucket_size>
class grid
{
public:
    grid(std::size_t nx, std::size_t ny)
        : nx_{ny}
        , ny_{ny}
    {
        ulls_.resize(nx_ * ny_);
    }

    grid(std::array<std::size_t, dim> shape)
        : nx_{shape[0]}
        , ny_{shape[1]}
    {
        ulls_.resize(nx_ * ny_);
    }

    void addToCell(std::array<int, dim> cell, T const& obj)
    {
        ulls_[cell[1] + ny_ * cell[0]].add(obj);
    }


    void add(std::vector<Particle<dim>> const& particles)
    {
        for (auto const& particle : particles)
        {
            addToCell(particle.iCell, particle);
        }
    }


    std::size_t total(std::size_t ix, std::size_t iy)
    {
        auto icell = iy + ix * ny_;
        return ulls_[icell].total();
    }


    auto select(std::vector<Particle<dim>>& particles, Box<dim> const& box) const
    {
        // count all particles
        std::size_t nbrTot = 0;
        for (auto const& c : box)
        {
            auto cc = cell(c);
            nbrTot += cc.total();
        }
        std::vector<Particle<dim>> selection(nbrTot);

        std::size_t ipart = 0;
        for (auto const& c : box)
        {
            auto plist = cell(c);
            for (Particle<dim> const* p : plist)
            {
                selection[ipart++] = *p;
            }
        }
        return selection;
    }


    ull<bucket_size, T>& cell(std::size_t ix, std::size_t iy)
    {
        auto& c = ulls_[iy + ix * ny_];
        return c;
    }
    ull<bucket_size, T> const& cell(std::size_t ix, std::size_t iy) const
    {
        auto& c = ulls_[iy + ix * ny_];
        return c;
    }

    ull<bucket_size, T>& cell(std::array<std::size_t, dim> cell)
    {
        auto& c = ulls_[cell[1] + cell[0] * ny_];
        return c;
    }

    ull<bucket_size, T> const& cell(std::array<std::size_t, dim> cell) const
    {
        auto& c = ulls_[cell[1] + cell[0] * ny_];
        return c;
    }


    auto capacity() const
    {
        auto tot = 0u;
        for (auto const& ull : ulls_)
        {
            tot += ull.capacity();
        }
        return tot;
    }


    void reset()
    {
        for (auto& ull : ulls_)
        {
            ull.empty();
        }
    }




private:
    std::size_t nx_, ny_;
    std::vector<ull<bucket_size, T>> ulls_;
};



template<std::size_t dim>
auto make_particles_in(Box<dim> box, std::size_t nppc)
{
    std::vector<Particle<dim>> particles;
    particles.reserve(box.size() * nppc);

    for (auto ix = 0u; ix <= box.upper[0]; ++ix)
    {
        for (auto iy = 0u; iy <= box.upper[1]; ++iy)
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
    return particles;
}

void test()
{
    //
    std::cout << "---------------START TESTS ------------------------------\n";
    constexpr auto dim = 2u;
    Box<dim> domain{{0, 0}, {9, 19}};
    std::size_t nppc = 4;
    auto particles   = make_particles_in(domain, nppc);

    std::cout << "making the grid...\n";
    grid<dim, Particle<dim>, 200> myGrid(domain.shape());

    std::cout << "pushing...\n";
    // pretend pushing
    for (auto ip = 0; ip < nppc * domain.size(); ++ip)
    {
        myGrid.addToCell(particles[ip].iCell, particles[ip]);
    }

    std::cout << "total grid capacity : " << myGrid.capacity() << "\n";


    std::cout << "testing nbr of particles per cell registered\n";
    for (auto ix = 0u; ix <= domain.upper[0]; ++ix)
    {
        for (auto iy = 0u; iy <= domain.upper[1]; ++iy)
        {
            if (myGrid.total(ix, iy) != nppc)
            {
                throw std::runtime_error("invalid number of particles in cell "
                                         + std::to_string(myGrid.total(ix, iy)));
            }
        }
    }
    std::cout << "nbr of particles ok.\n";


    auto& cell  = myGrid.cell(2, 1);
    auto it_beg = cell.begin();
    auto it_end = cell.end();
    std::cout << "particle in cell 10,12 : " << (*it_beg)->iCell[0] << "\n";

    Box<dim> selection_box({4, 3}, {6, 4});
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

    std::cout << "selecting particle...\n";
    auto selected = myGrid.select(particles, selection_box);
    std::cout << "nbr of particles selected : " << selected.size() << "\n";
    std::cout << "nbr of particles expected : " << selection_box.size() * nppc << "\n";
    if (selected.size() != selection_box.size() * nppc)
        throw std::runtime_error("invalid number of found particles");

    std::cout << "---------------END OF TESTS ------------------------------\n";
}

template<std::size_t dim>
auto box_generator(Box<dim> const& domain, std::size_t lower_size, std::size_t upper_size,
                   std::size_t nbr_boxes)
{
    std::vector<Box<dim>> boxes;
    boxes.reserve(nbr_boxes);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> size_dist(lower_size, upper_size);
    std::array<std::uniform_int_distribution<std::size_t>, dim> pos_dist;
    for (auto idim = 0u; idim < dim; ++idim)
        pos_dist[idim]
            = std::uniform_int_distribution<std::size_t>(domain.lower[idim], domain.upper[idim]);

    for (auto ibox = 0; ibox < nbr_boxes; ++ibox)
    {
        auto size = size_dist(gen);
        auto ipos = pos_dist[0](gen);
        auto jpos = pos_dist[1](gen);
        boxes.push_back(Box<dim>{{ipos, jpos}, {ipos + size - 1, jpos + size - 1}});
    }
    return boxes;
}



template<std::size_t dim>
inline bool isInBox(Box<dim> const& box, Particle<dim> const& particle)
{
    auto const& iCell = particle.iCell;

    auto const& lower = box.lower;
    auto const& upper = box.upper;


    if (iCell[0] >= lower[0] && iCell[0] <= upper[0])
    {
        if constexpr (dim > 1)
        {
            if (iCell[1] >= lower[1] && iCell[1] <= upper[1])
            {
                if constexpr (dim > 2)
                {
                    if (iCell[2] >= lower[2] && iCell[2] <= upper[2])
                    {
                        return true;
                    }
                }
                else
                {
                    return true;
                }
            }
        }
        else
        {
            return true;
        }
    }
    return false;
}

int main()
{
    test();
    constexpr auto dim = 2u;
    Box<dim> domain{{0, 0}, {30, 40}};
    std::size_t nppc = 100;

    auto constexpr bucket_size = 64u;

    auto boxes = box_generator(domain, 5, 10, 10);


    std::array<std::size_t, 100> nppcs;
    std::vector<std::size_t> m1times(nppcs.size());
    std::vector<std::size_t> m2times(nppcs.size());
    std::generate(std::begin(nppcs), std::end(nppcs), []() {
        static std::size_t nppc = 1;
        return nppc++;
    });


    std::cout << "first method\n";
    auto m1time = std::begin(m1times);
    for (auto const& nppc : nppcs)
    {
        std::cout << nppc << " particles per cell\n";
        auto particles = make_particles_in(domain, nppc);
        grid<dim, Particle<dim>, bucket_size> myGrid(domain.shape());
        {
            auto _ = Timer(*(m1time));
            myGrid.add(particles);
            for (auto const& box : boxes)
            {
                auto intersection = domain * box;
                auto found        = myGrid.select(particles, intersection);
            }
        }
        m1time++;
    }

    std::cout << "second method\n";
    auto m2time = std::begin(m2times);
    for (auto const& nppc : nppcs)
    {
        std::cout << nppc << " particles per cell\n";
        auto particles = make_particles_in(domain, nppc);
        {
            auto _ = Timer(*(m2time));
            std::vector<Particle<dim>> selected;
            for (auto const& box : boxes)
            {
                for (auto const& p : particles)
                {
                    auto intersection = domain * box;
                    if (isInBox(intersection, p))
                    {
                        selected.push_back(p);
                    }
                }
            }
        }
        m2time++;
    }


    std::cout << "saving results....\n";
    std::ofstream output;
    output.open("test_" + std::to_string(bucket_size) + ".txt");
    for (auto i = 0u; i < nppcs.size(); ++i)
    {
        output << nppcs[i] << " " << m1times[i] << " " << m2times[i] << "\n";
    }
    output.close();


    return 0;
}

