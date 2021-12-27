/*
 * author: max klein
 */

#ifndef EXP_HIST_H
#define EXP_HIST_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>
#if __cplusplus <= 199711L
#include <map>
#else
#include <unordered_map>
#endif

#include "lm/Iterator.h"
#include "lm/Print.h"
//#include "lm/protowrap/NDArray.h"
#include "lm/Types.h"
#include "robertslab/pbuf/Sparse.pb.h"

template <typename Container>
struct CHash {
    std::size_t operator()(Container const& c) const {
        std::size_t seed = c.size();
        for (typename Container::const_iterator it=c.begin(); it!=c.end(); it++) {
            seed ^= *it + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

namespace experimental
{

template <class H>
void printHist(H& h)
{
    for (typename H::const_iterator it=h.begin(); it!=h.end(); it++)
    {
        std::cout << it->first << " " << it->second << std::endl;
    }
}

template <class H>
void printHistObs(H& h)
{
    for (KeyIter<H> keys(h.begin()); keys!=h.end(); keys++)
    {
        std::cout << *keys << std::endl;
    }
}

template <class H>
void printHistVals(H& h)
{
    for (ValIter<H> values(h.begin()); values!=h.end(); values++)
    {
        std::cout << *values << std::endl;
    }
}


#if __cplusplus <= 199711L
template<class Obs, class Val=uint64_t, class=void>
struct _DefaultHistMap
{
    typedef std::map<Obs, Val> type;
};
#else
template<class Obs, class Val=uint64_t, class=void>
struct _DefaultHistMap
{
    typedef std::unordered_map<Obs, Val> type;
};

template<class Obs, class Val>
struct _DefaultHistMap<Obs, Val, typename EnableIf<IsVector<Obs>::value>::type>
{
    typedef std::unordered_map<Obs, Val, CHash<Obs>> type;
};
#endif

template <class Obs, class Val=uint64_t, class Map=typename _DefaultHistMap<Obs, Val>::type>
class Sparse
{
public:
    Sparse(): map() {}

    typedef robertslab::pbuf::Sparse Msg;
    typedef KeyIter<Map> KeyIter;
    typedef ValIter<Map> ValIter;
    typedef flattening_iterator<KeyIter> KeyFlatIter;

    inline void addObs(const Obs& obs)
    {
        map[obs]++;
    }

    template <class C, class=typename EnableIfNot<IsSame<C, Obs>::value>::type>
    inline void addObs(const C& container)
    {
        for (typename C::const_iterator it = container.begin(); it!=container.end(); it++)
        {
            addObs(*it);
        }
    }

    KeyIter beginObs() {return KeyIter(map).begin();}
    KeyIter endObs() {return KeyIter(map).end();}

    ValIter beginVal() {return ValIter(map).begin();}
    ValIter endVal() {return ValIter(map).end();}

    KeyFlatIter beginObsFlat() {return flatten(beginObs(), endObs());}
    KeyFlatIter endObsFlat() {return flatten(endObs());}

    void printObs()
    {
        printIterable(std::cout, beginObs());
    }

    void printVals()
    {
        printIterable(std::cout, beginVal());
    }

    void printObsFlat()
    {
        printIterable(std::cout, beginObsFlat());
    }

public:
    Map map;
};

}

#endif //EXP_HIST_H
