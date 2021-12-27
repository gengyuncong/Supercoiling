/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Max Klein
 */
#ifndef LM_ITERATOR_H_
#define LM_ITERATOR_H_

#include <iterator>
#include <map>
#include <memory>
#include <string>

#include "lm/Traits.h"
#include "lm/Types.h"

// Forward declarations
template <typename Iter> Iter peek(Iter iter);

template <typename Iter, typename Cont> bool isLast(Iter iter, const Cont& cont);

template<class T> struct OutputIteratorTraits;

template <class Derived, class WrappedIter, class value_type, class IterTag=std::forward_iterator_tag> class WrapIter;
template <class Derived, class OuterIter, class InnerIter=typename OuterIter::value_type::iterator, class IterTag=std::forward_iterator_tag> class FlatIter;

template<class Map> class KeyIter;
template<class Map> class KeyConstIter;
template<class Map> class ValIter;
template<class Map> class ValConstIter;

template<typename Map> class KeyFlatIter;
template<typename Map> class KeyFlatConstIter;

template <class Iter, class IterTag=std::forward_iterator_tag, bool=IsVector<typename Iter::value_type>::value> struct FlatIterPolicy;
template <class Map, bool=IsVector<typename Map::key_type>::value> struct KeyFlatIterPolicy;
template <class Map, bool=IsVector<typename Map::key_type>::value> struct KeyFlatConstIterPolicy;


template <typename Iter>
Iter peek(Iter iter)
{
    return ++iter;
}

/*
 * - function to check if an iterator points to the last element in a container
 *     - modified from http://stackoverflow.com/a/3516224/425458
 */
template <typename Iter, typename Cont>
bool isLast(Iter iter, const Cont& cont)
{
    // if the iterator points to the end, then return true only if the container is zero sized
    if (cont.end()==iter)
    {
        return (cont.size()==0);
    }
    // otherwise, check if the next iterator value brings us to the container end
    else
    {
        return (cont.end()==peek(iter));
    }
}

/*
 * - iterator traits specialized by tag (output, forward, etc). more generic than the STL version
 *     - OutputIteratorTraits modified from http://stackoverflow.com/a/29084919/425458
 */
template<class T>
struct OutputIteratorTraits
:std::iterator_traits<T> {};

template< class OutputIt, class T>
struct OutputIteratorTraits<std::raw_storage_iterator<OutputIt, T> >
:std::iterator<std::output_iterator_tag, T> {};

template<class Container>
struct OutputIteratorTraits<std::back_insert_iterator<Container> >
:std::iterator<std::output_iterator_tag, typename Container::value_type> {};

template<class Container>
struct OutputIteratorTraits<std::front_insert_iterator<Container> >
:std::iterator<std::output_iterator_tag, typename Container::value_type> {};

template<class Container>
struct OutputIteratorTraits<std::insert_iterator<Container> >
:std::iterator<std::output_iterator_tag, typename Container::value_type> {};

#if __cplusplus > 199711L

//template <class T, class charT = char, class traits = std::char_traits<charT> >
//struct OutputIteratorTraits<std::ostream_iterator<T, charT, traits> >
//:std::iterator<std::output_iterator_tag, T> {};
//
//template <class charT, class traits = std::char_traits<charT> >
//struct OutputIteratorTraits<std::ostreambuf_iterator<charT, traits> >
//:std::iterator<std::output_iterator_tag, charT> {};

#endif

template <class Derived, class WrappedIter, class ValueType, class IterTag>
class WrapIter : public std::iterator<IterTag, ValueType>
{
public:
    typedef Derived     this_type;
    typedef WrappedIter wrapped_iterator;

    typedef ValueType                                  value_type;
    typedef typename wrapped_iterator::difference_type difference_type;
    typedef typename wrapped_iterator::pointer         pointer;
    typedef typename wrapped_iterator::reference       reference;

    WrapIter(): wrapit() {};
    WrapIter(this_type begin, this_type end): wrapit(begin.wrapit), wrapit_begin(begin.wrapit), wrapit_end(end.wrapit) {};
    WrapIter(wrapped_iterator begin, wrapped_iterator end): wrapit(begin), wrapit_begin(begin), wrapit_end(end) {};

    // explicit operator wrapped_iterator() const {return wrapit;}

    this_type& operator++() {++wrapit; return *static_cast<this_type*>(this);}
    this_type operator++(int) {this_type tmp(*static_cast<this_type*>(this)); operator++(); return tmp;}

    bool operator==(const this_type& rhs) const {return wrapit==rhs.wrapit;}
    bool operator!=(const this_type& rhs) const {return wrapit!=rhs.wrapit;}

    bool operator==(const wrapped_iterator& rhs) const {return wrapit==rhs;}
    bool operator!=(const wrapped_iterator& rhs) const {return wrapit!=rhs;}

    this_type begin() const {return this_type(wrapit_begin, wrapit_end);}
    this_type end() const {return this_type(wrapit_end, wrapit_end);}
    
public:
    wrapped_iterator wrapit;

protected:
    wrapped_iterator wrapit_begin;
    wrapped_iterator wrapit_end;
};

/** Base class for forward iterators that "flattens" a container of containers. For example,
 * a vector<vector<int> > containing {{1, 2, 3}, {4, 5, 6}} is iterated as
 * a single range, {1, 2, 3, 4, 5, 6}.
 *
 * Heavily modified from https://stackoverflow.com/a/3623597/425458 */
template <class Derived, class OuterIter, class InnerIter, class IterTag>
class FlatIter : public std::iterator<IterTag, typename InnerIter::value_type>
{
public:
    typedef Derived this_type;
    
    typedef OuterIter outer_iterator;
    typedef InnerIter inner_iterator;

    typedef typename inner_iterator::value_type      value_type;
    typedef typename inner_iterator::difference_type difference_type;
    typedef typename inner_iterator::pointer         pointer;
    typedef typename inner_iterator::reference       reference;


    FlatIter(): outerit() {};
    FlatIter(this_type begin, this_type end): outerit(begin.outerit), outerit_begin(begin.outerit), outerit_end(end.outerit) {init();}
    FlatIter(outer_iterator begin, outer_iterator end): outerit(begin), outerit_begin(begin), outerit_end(end) {init();}

    void init()
    {
        if (outerit == outerit_end) {return;}

        innerit = outerit->begin();
        advance_past_empty_inner_containers();
    }

    // explicit operator inner_iterator() const {return innerit;}

    this_type& operator++() {++innerit; if (innerit == outerit->end()) advance_past_empty_inner_containers(); return *static_cast<this_type*>(this);}
    this_type operator++(int) {this_type tmp(*static_cast<this_type*>(this)); operator++(); return tmp;}

    friend bool operator==(const FlatIter& a, const FlatIter& b)
    {
        if (a.outerit != b.outerit) return false;
        if (a.outerit != a.outerit_end && b.outerit != b.outerit_end && a.innerit != b.innerit) return false;
        return true;
    }
    friend bool operator!=(const FlatIter& a, const FlatIter& b) {return !(a == b);}

    this_type begin() const {return this_type(outerit_begin, outerit_end);}
    this_type end() const {return this_type(outerit_end, outerit_end);}

    outer_iterator outerit;
    
protected:
    void advance_past_empty_inner_containers()
    {
        while (outerit != outerit_end && innerit == outerit->end())
        {
            ++outerit;
            if (outerit != outerit_end) innerit = outerit->begin();
        }
    }

    outer_iterator outerit_begin;
    outer_iterator outerit_end;
    inner_iterator innerit;
};

/*
 * Key/Value iterators for map types, a la the Python map.keys() and map.values() methods.
 * Initialize them from the normal map iterators (eg KeyIter(map.begin()))
 */

template <class Map>
class KeyIter : public WrapIter<KeyIter<Map>, typename Map::iterator, typename Map::key_type>
{
public:
    typedef typename KeyIter::WrapIter Base;
    typedef typename Base::wrapped_iterator wrapped_iterator;
    typedef typename Base::value_type value_type;

    typedef KeyIter<Map> iterator;
    typedef KeyConstIter<Map> const_iterator;
    
    KeyIter(): Base() {}
    KeyIter(wrapped_iterator begin, wrapped_iterator end): Base(begin, end) {}
    KeyIter(Map& map): Base(map.begin(), map.end()) {}

    value_type* operator->() {return (value_type*)&(Base::wrapit->first);}
    value_type& operator*() {return Base::wrapit->first;}
};

template<class Map>
class KeyConstIter : public WrapIter<KeyConstIter<Map>, typename Map::const_iterator, typename Map::key_type>
{
public:
    typedef typename KeyConstIter::WrapIter Base;
    typedef typename Base::wrapped_iterator wrapped_iterator;
    typedef typename Base::value_type value_type;

    typedef KeyIter<Map> iterator;
    typedef KeyConstIter<Map> const_iterator;

    KeyConstIter(): Base() {}
    KeyConstIter(wrapped_iterator begin, wrapped_iterator end): Base(begin, end) {}
    KeyConstIter(const Map& map): Base(map.begin(), map.end()) {}

    const value_type* operator->() const {return (const value_type*)&(Base::wrapit->first);}
    const value_type& operator*() const {return Base::wrapit->first;}
};

template<class Map>
class ValIter: public WrapIter<ValIter<Map>, typename Map::iterator, typename Map::mapped_type>
{
public:
    typedef typename ValIter::WrapIter Base;
    typedef typename Base::wrapped_iterator wrapped_iterator;
    typedef typename Base::value_type value_type;

    typedef ValIter<Map> iterator;
    typedef ValConstIter<Map> const_iterator;

    ValIter(): Base() {}
    ValIter(wrapped_iterator begin, wrapped_iterator end): Base(begin, end) {}
    ValIter(Map& map): Base(map.begin(), map.end()) {}

    value_type* operator->() {return (value_type*)&(Base::wrapit->second);}
    value_type& operator*() {return Base::wrapit->second;}
};

template<class Map>
class ValConstIter: public WrapIter<ValConstIter<Map>, typename Map::const_iterator, typename Map::mapped_type>
{
public:
    typedef typename ValConstIter::WrapIter Base;
    typedef typename Base::wrapped_iterator wrapped_iterator;
    typedef typename Base::value_type value_type;

    typedef ValIter<Map> iterator;
    typedef ValConstIter<Map> const_iterator;

    ValConstIter(): Base() {}
    ValConstIter(wrapped_iterator begin, wrapped_iterator end): Base(begin, end) {}
    ValConstIter(const Map& map): Base(map.begin(), map.end()) {}

    const value_type* operator->() const {return (const value_type*)&(Base::wrapit->second);}
    const value_type& operator*() const {return Base::wrapit->second;}
};

/*
 * Flattening Key/Value iterators for map types. For example, given a map of type Map -> std::map<std::vector<int>, int>,
 * the iterator KeyFlatIter<Map>(map) will iterate over each int in the key vectors.
 */

template <class Map>
class KeyFlatIter : public FlatIter<KeyFlatIter<Map>, KeyIter<Map> >
{
public:
    typedef typename KeyFlatIter::FlatIter Base;
    typedef typename Base::outer_iterator outer_iterator;
    typedef typename Base::value_type value_type;

    typedef KeyFlatIter<Map> iterator;
    typedef KeyFlatConstIter<Map> const_iterator;

    KeyFlatIter(): Base() {}
    KeyFlatIter(outer_iterator begin, outer_iterator end): Base(begin, end) {}
    KeyFlatIter(Map& map): Base(outer_iterator(map).begin(), outer_iterator(map).end()) {}

    value_type* operator->() {return (value_type*)&(*Base::innerit);}
    value_type& operator*() {return *Base::innerit;}
};

template<class Map>
class KeyFlatConstIter : public FlatIter<KeyFlatConstIter<Map>, KeyConstIter<Map>, typename KeyConstIter<Map>::value_type::const_iterator >
{
public:
    typedef typename KeyFlatConstIter::FlatIter Base;
    typedef typename Base::outer_iterator outer_iterator;
    typedef typename Base::value_type value_type;

    typedef KeyFlatIter<Map> iterator;
    typedef KeyFlatConstIter<Map> const_iterator;

    KeyFlatConstIter(): Base() {}
    KeyFlatConstIter(outer_iterator begin, outer_iterator end): Base(begin, end) {}
    KeyFlatConstIter(const Map& map): Base(outer_iterator(map).begin(), outer_iterator(map).end()) {}

    const value_type* operator->() const {return (const value_type*)&(*Base::innerit);}
    const value_type& operator*() const {return *Base::innerit;}
};

/* - policies that can be used to infer the appropriate type of flat iterator  */

template <class Iter, class IterTag, bool> struct FlatIterPolicy {typedef Iter type;};
template <class Iter, class IterTag> struct FlatIterPolicy<Iter, IterTag, true> {typedef FlatIter<Iter, IterTag> type;};

template <class Map, bool> struct KeyFlatIterPolicy {typedef KeyIter<Map> type;};
template <class Map> struct KeyFlatIterPolicy<Map, true> {typedef KeyFlatIter<Map> type;};

template <class Map, bool> struct KeyFlatConstIterPolicy {typedef KeyConstIter<Map> type;};
template <class Map> struct KeyFlatConstIterPolicy<Map, true> {typedef KeyFlatConstIter<Map> type;};

/*
 * - simplified version of transform_iterator from boost
 */
/*
#include <boost/iterator/iterator_adaptor.hpp>

template <class UnaryFunction, class Iterator, class Reference = boost::use_default, class Value = boost::use_default>
class transform_iterator;

namespace internal
{
// Compute the iterator_adaptor instantiation to be used for transform_iterator
template <class UnaryFunc, class Iterator, class Reference, class Value>
struct transform_iterator_base
{
private:

    typedef typename boost::detail::ia_dflt_help<Reference, result_of<const UnaryFunc(typename std::iterator_traits<Iterator>::reference)> >::type reference;
    typedef typename boost::detail::ia_dflt_help<Value, remove_reference<reference> >::type cv_value_type;

public:
    typedef boost::iterator_adaptor<
    transform_iterator<UnaryFunc, Iterator, Reference, Value>, Iterator, cv_value_type,
    boost::use_default,    // Leave the traversal category alone
    reference> type;
};
}

template <class UnaryFunc, class Iterator, class Reference, class Value>
class transform_iterator : public internal::transform_iterator_base<UnaryFunc, Iterator, Reference, Value>::type
{
    typedef typename internal::transform_iterator_base<UnaryFunc, Iterator, Reference, Value>::type super_t;

public:
    transform_iterator() { }

    transform_iterator(Iterator const& x, UnaryFunc f)
    :super_t(x), m_f(f)
    {
    }

    explicit transform_iterator(Iterator const& x)
    :super_t(x)
    {
    }

    template <class OtherUnaryFunction, class OtherIterator, class OtherReference, class OtherValue>
    transform_iterator(transform_iterator<OtherUnaryFunction, OtherIterator, OtherReference, OtherValue> const& t,
                       typename EnableIfConvertible<OtherIterator, Iterator>::type* = 0,
                       typename EnableIfConvertible<OtherUnaryFunction, UnaryFunc>::type* = 0)
    :super_t(t.base()), m_f(t.functor())
    {
    }

    UnaryFunc functor() const { return m_f; }

private:
    typename super_t::reference dereference() const { return m_f(*this->base()); }

    UnaryFunc m_f;
};

template <class UnaryFunc, class Iterator>
inline transform_iterator<UnaryFunc, Iterator> make_transform_iterator(Iterator it, UnaryFunc fun)
{
    return transform_iterator<UnaryFunc, Iterator>(it, fun);
}

template <class Return, class Argument, class Iterator>
inline transform_iterator< Return (*)(Argument), Iterator, Return>
make_transform_iterator(Iterator it, Return (*fun)(Argument))
{
    return transform_iterator<Return (*)(Argument), Iterator, Return>(it, fun);
}
*/

#endif /* LM_ITERATOR_H_ */