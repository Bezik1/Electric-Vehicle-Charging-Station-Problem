#pragma once

#include <utility>

/**
 * @brief This interface purpose is to be used, as a template for
 * creating new optimzation problems using an elegant, objective approach.
 * It is designed to be used for gurobi using class instances.
 * 
 */
class GurobiModel
{
public:

    /**
     * @brief Builds deciding variables.
     * 
     */
    virtual void build_variables() = 0;

    /**
     * @brief Creates criterion / objective function.
     * 
     */
    virtual void build_criterion() = 0;

    /**
     * @brief Build constraints, used inside optimziation process.
     * 
     */
    virtual void build_constraints() = 0;

    virtual ~GurobiModel() = default;
};