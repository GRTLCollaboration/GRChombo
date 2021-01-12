/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MULTI_LEVEL_TASK
#define MULTI_LEVEL_TASK

#include "AMR.H"
#include "GRAMRLevel.hpp"
#include "Scheduler.H"

#include <limits> // std::numeric_limits

//! This is just an interface for the AMR scheduler to call some GRAMRLevel (or
//! any other Example specific level) function on every AMRLevel
//! Satisfies syntax of Chmbo's Scheduler such that it can be passed to GRAMR
//! and be scheduled
template <class Level = GRAMRLevel>
class MultiLevelTask : public Scheduler::PeriodicFunction
{
    bool m_reverse_levels;
    std::function<void(Level *)> m_func;

    // Use default condstructor of PeriodicFunction
    using PeriodicFunction::PeriodicFunction;

    AMR *m_amr_ptr; //! pointer to AMR object

  public:
    MultiLevelTask(std::function<void(Level *)> a_func,
                   bool a_reverse_levels = true)
        : m_func(a_func), m_reverse_levels(a_reverse_levels)
    {
    }

    // required from Scheduler::PeriodicFunction
    virtual void setUp(AMR &a_AMR, int a_interval = -1) override
    {
        m_amr_ptr = &a_AMR;
    }

    // required from Scheduler::PeriodicFunction
    virtual void operator()(int a_step = 0, Real a_time = 0.) override
    {
        auto amr_level_ptrs = m_amr_ptr->getAMRLevels().stdVector();

        // need to reverse this vector so that m_func is called in order of
        // finest level to coarsest. This is important for example for
        // 'specificPostTimeStep', which is always run in reverse order of
        // levels
        if (m_reverse_levels)
            std::reverse(std::begin(amr_level_ptrs), std::end(amr_level_ptrs));

        for (AMRLevel *amr_level_ptr : amr_level_ptrs)
            m_func(Level::gr_cast(amr_level_ptr));
    }
};

//! This is just an interface for the AMR scheduler to call some Level
//! function on every AMRLevel
//! This can either be called directly by calling execute, or passed to an AMR
//! (as GRAMR) by doing gr_amr.schedule(me) (this version will make it be called
//! only after plot files are written, if that is ever an interest)
template <class Level = GRAMRLevel>
class MultiLevelTaskPtr : public RefCountedPtr<Scheduler>
{
    RefCountedPtr<MultiLevelTask<Level>> m_ptr_to_func;

  public:
    //! interval defines the frequency with which the scheduler will be called
    //! if added to an AMR
    MultiLevelTaskPtr(std::function<void(Level *)> a_func,
                      bool a_reverse_levels = true,
                      int a_interval = std::numeric_limits<int>::max())
        : RefCountedPtr<Scheduler>(new Scheduler),
          m_ptr_to_func(new MultiLevelTask<Level>(a_func, a_reverse_levels))
    // the two 'new' pointers are deleted by RefCountedPtr, no memory leak
    {
        if (a_interval <= 0) // the user probably means "never again"
            a_interval = std::numeric_limits<int>::max();
        (*this)->schedule(m_ptr_to_func, a_interval);
    }

    // run immediately!
    void execute(AMR &amr)
    {
        m_ptr_to_func->setUp(amr);
        (*m_ptr_to_func)();
    }
};

#endif /* MULTI_LEVEL_TASK */
