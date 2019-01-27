/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRPARMPARSE_HPP_
#define GRPARMPARSE_HPP_

#include "ParmParse.H"
#include "parstream.H" //Gives us pout()
#include <algorithm>
#include <memory>

/// Helper structs to translate a dataype into a Chombo ParmParse data type
template <class T> struct ParmParseTranslator;

template <> struct ParmParseTranslator<double>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppDouble;
};

template <> struct ParmParseTranslator<int>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppInt;
};

template <> struct ParmParseTranslator<bool>
{
    static constexpr ParmParse::PPType pp_type = ParmParse::ppBool;
};
// End: Helper structs to translate a dataype into a Chombo ParmParse data type

class GRParmParse : public ParmParse
{
  public:
    using ParmParse::ParmParse; // Just use ParmParse's constructor

    // (MK): I called the functions below "load" rather than "get" to avoid
    // clashes with the many  different overloads of "get" in ParmParse. Also, I
    // think load is a more intuitive name.

    /// Loads an array from the parameter file
    template <class data_t, long unsigned int n_comp>
    void load(const char *name, std::array<data_t, n_comp> &array) const
    {
        getarr(name, ParmParseTranslator<data_t>::pp_type, array.data(), 0,
               n_comp, -1);
    }

    /// Loads a vector with num_comp components from the parameter file
    template <class data_t>
    void load(const char *name, std::vector<data_t> &vector,
              const int num_comp) const
    {
        getarr(name, vector, 0, num_comp);
    }

    /// Loads a value from the parameter file
    template <class data_t>
    typename std::enable_if<
        !std::is_enum<data_t>::value>::type // Can't use for enum types
    load(const char *name, data_t &parameter) const
    {
        get(name, parameter);
    }

    /// Loads an enum value from the parameter file
    template <typename enum_type>
    typename std::enable_if<
        std::is_enum<enum_type>::value>::type // Only enabled for enum types
    load(const char *name, enum_type &parameter) const
    {
        int iparam;
        get(name, iparam);
        parameter = static_cast<enum_type>(iparam);
    }

    /// Loads a value from the parameter file, if the value isn't defined it
    /// sets to the supplied default
    template <class data_t>
    void load(const char *name, data_t &parameter,
              const data_t default_value) const
    {
        if (contains(name))
        {
            load(name, parameter);
        }
        else
        {
            parameter = default_value;
            pout() << "Parameter: " << name << " not found in parameter file. "
                   << "It has been set to its default value." << std::endl;
        }
    }

    /// Loads a vector with num_comp components from the parameter file, if the
    /// vector isn't defined it sets all components to the supplied default
    template <class data_t>
    void load(const char *name, std::vector<data_t> &vector, const int num_comp,
              const data_t default_value) const
    {
        if(contains(name))
        {
            load(name, vector, num_comp);
        }
        else
        {
            vector.resize(num_comp);
            vector.assign(num_comp, default_value);
            pout() << "Parameter: " << name << " not found in parameter file. "
                   << "It has been set to its default value." << std::endl;
        }
    }
};

#endif /* GRPARMPARSE_HPP_ */
