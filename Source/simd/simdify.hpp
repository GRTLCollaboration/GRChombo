/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMDIFY_HPP_
#define SIMDIFY_HPP_

#if !defined(SIMD_HPP_)
#error "This file should only be included through simd.hpp"
#endif

#include <type_traits>

/// This class is used to load simd<t> objects from a number (given by the simd
/// width) of t's
/** simd_proxy containts a pointer to the memory location from where the simd<t>
 * is to be loaded. It defines a conversion to simd<t>, which is a load (gather)
 * operation, and an operator =simd<t> which does the store (scatter).
 * Furthermore, several basic arithmetic operators are defined directly on the
 * simd_proxy so that we don't have to load, do the operation, and store back
 * again.
 */
template <typename data_t> struct simd_proxy
{
    typedef simd<std::remove_cv_t<data_t>> simd_t;
    data_t *const m_ptr;

    ALWAYS_INLINE
    simd_proxy(data_t *m_ptr) : m_ptr(m_ptr) {}

    ALWAYS_INLINE
    data_t *operator&() const { return m_ptr; }

    ALWAYS_INLINE
    operator simd_t() const { return simd_t::load(m_ptr); }

    ALWAYS_INLINE
    simd_proxy &operator=(const simd_t &rhs)
    {
        simd_t::store(m_ptr, rhs);
        return *this;
    }

    /*ALWAYS_INLINE
    const simd_proxy& operator=(const simd_t& rhs) const
    {
        const_cast<simd_proxy&>(*this) = rhs;
    }*/

    // The operator definitions below are there just in case we want to do
    // arithmetic directly  with a simd_proxy object. e.g. if we write
    // 3*SIMDIFY<simd<t>>(ptr_t). The behaviour mimics what we  would get if we
    // loaded a simd<t> object, did the operation, and then stored it back to
    // the location that the  simd_proxy points to.
    ALWAYS_INLINE
    simd_t operator+(const simd_t &other)
    {
        return static_cast<simd_t>(*this) + other;
    }

    ALWAYS_INLINE
    simd_t operator-(const simd_t &other)
    {
        return static_cast<simd_t>(*this) - other;
    }

    ALWAYS_INLINE
    simd_t operator*(const simd_t &other)
    {
        return static_cast<simd_t>(*this) * other;
    }

    ALWAYS_INLINE
    simd_t operator/(const simd_t &other)
    {
        return static_cast<simd_t>(*this) / other;
    }

    friend ALWAYS_INLINE simd_t operator+(const simd_t &a, const simd_proxy &b)
    {
        return a + static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE simd_t operator-(const simd_t &a, const simd_proxy &b)
    {
        return a - static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE simd_t operator*(const simd_t &a, const simd_proxy &b)
    {
        return a * static_cast<simd_t>(b);
    }

    friend ALWAYS_INLINE simd_t operator/(const simd_t &a, const simd_proxy &b)
    {
        return a / static_cast<simd_t>(b);
    }

    ALWAYS_INLINE
    simd_proxy &operator+=(const simd_t &rhs)
    {
        (*this) = (*this) + rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy &operator-=(const simd_t &rhs)
    {
        (*this) = (*this) - rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy &operator*=(const simd_t &rhs)
    {
        (*this) = (*this) * rhs;
        return *this;
    }

    ALWAYS_INLINE
    simd_proxy &operator/=(const simd_t &rhs)
    {
        (*this) = (*this) / rhs;
        return *this;
    }
};

// This class can be viewed as an array but with special property that the
// dereferencing operator  and operator [] do not return an element (say type
// data_t) but a simd_proxy<data_t>. In practice,  this means that [i] doesn't
// give the ith element but a simd_proxy which can be used to load a  vector of
// simd width starting at the ith element.
template <typename data_t> struct simd_array_wrapper
{
    data_t *m_ptr;

    ALWAYS_INLINE
    simd_array_wrapper(data_t *m_ptr) : m_ptr(m_ptr) {}

    ALWAYS_INLINE
    simd_proxy<data_t> operator[](int i) const { return &m_ptr[i]; }

    ALWAYS_INLINE
    simd_proxy<data_t> operator*() const { return *this[0]; }
};

// Begin: structs that help create the return type of SIMDIFY and make sure it
// is only called on valid input types. -->

// If the first and second template parameter are the same then this struct will
// have a member of type  given by the third template parameter and called
// "type"
template <typename q1, typename q2, typename t> struct _simd_enable_if_same
{
};
template <typename q, typename t> struct _simd_enable_if_same<q, q, t>
{
    using type = t;
};

// This struct is used to create the right return types for the SIMDIFY function
template <typename t, typename ptr_t> struct _simdify
{
    typedef
        typename _simd_enable_if_same<t, std::remove_cv_t<ptr_t>, ptr_t *>::type
            type;
};

template <typename t, typename ptr_t> struct _simdify<simd<t>, ptr_t>
{
    typedef typename _simd_enable_if_same<t, std::remove_cv_t<ptr_t>,
                                          simd_array_wrapper<ptr_t>>::type type;
};
//<--End: structs that help create the return type of SIMDIFY and make sure it
// is only called on valid input types.

/// SIMDIFY<simd<base_t>>(ptr) is used to load simd data from an array of
/// base_t's pointed to by ptr.  SIMDIFY<t>(ptr) for t not a simd object just
/// returns the ptr.
/** The usage of SIMDIFY<simd<base_t>>(ptr) is as follows:
 * SIMDIY<simd<base_t>>(ptr) returns a simd_array_wrapper which, when supplied
 * with operator [i], yields a simd_proxy at the ith element in the array. This
 * simd_proxy can then be used to load a simd<base_t> object with data starting
 * at the ith position and with width given by the simd length.
 */
// Given SIMDIFY<t>(ptr_t) the return types are:
//- If t=ptr_t (up to const and volatile) and t not simd<base_t>: return type
// just ptr_t* - If t=simd<base_t> (up to const and volatile) where base_t =
// ptr_t: return type simd_array_wrapper<ptr_t> - Else: substitution failure
template <typename t, typename ptr_t>
ALWAYS_INLINE typename _simdify<t, ptr_t>::type // See comments at _simdify for
                                                // explanation of return types
SIMDIFY(ptr_t *ptr)
{
    return ptr;
}

#endif /* SIMDIFY_HPP_ */
