/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef FILESYSTEMTOOLS_HPP_
#define FILESYSTEMTOOLS_HPP_

// Chombo includes
#include "SPMD.H" // gives procID()

// Other includes
#include <sys/stat.h> // gives 'stat' and 'S_ISDIR'
#include <unistd.h>   // gives 'mkdir'

// Chombo namespace
#include "UsingNamespace.H"

// Some filesystem useful functions

// This should be changed in the future by the actual C++17 filesystem library
// when we decide to only support from C++17

namespace FilesystemTools
{

static bool directory_exists(const std::string &path)
{
    struct stat stat_struct;
    return (path == "") || (stat(path.c_str(), &stat_struct) == 0 &&
                            S_ISDIR(stat_struct.st_mode));
}

static bool mkdir_recursive(const std::string &path)
{
#ifdef CH_MPI
    // all processes should get here at the same time
    // e.g. if doing:
    //     if (!directory_exists(path))
    //         mkdir_recursive(path);
    // We want to make sure the directory is not created until all
    // ranks have passed the 'directory_exists'
    MPI_Barrier(Chombo_MPI::comm);
#endif

    bool success = true;
    if (procID() == 0)
    {
        // For Windows backslash, use delimeters = "/\\"
        static const std::string delimeter = "/";
        std::size_t found = path.find_first_of(delimeter);

        // Even though this might look like RWX permissions for everyone,
        // this will be masked by the system user mask (see 'man umask' or
        // https://man7.org/linux/man-pages/man2/umask.2.html)
        mode_t permissions = S_IRWXU | S_IRWXG | S_IRWXO;

        // NB: this would be very beautiful recursively, but let's not do it
        // because the function involves MPI_Barrier's
        while (success && found != std::string::npos)
        {
            std::string subpath;
            subpath = path.substr(0, found);

            if (subpath != "." && subpath != "")
            {
                // success if created or if not created because it already
                // exists
                success &= (mkdir(subpath.c_str(), permissions) == 0 ||
                            errno == EEXIST);
            }
            found = path.find_first_of(delimeter, found + 1);
        }
        // if path doesn't finish with "/", one more to do
        if (found != path.size() - 1)
            success &=
                (mkdir(path.c_str(), permissions) == 0 || errno == EEXIST);
    }

#ifdef CH_MPI
    // all processes should wait to make sure directory structure is well set
    // for everyone
    MPI_Barrier(Chombo_MPI::comm);
#endif

    return success;
}
} // namespace FilesystemTools

#endif /* FILESYSTEMTOOLS_HPP_ */
