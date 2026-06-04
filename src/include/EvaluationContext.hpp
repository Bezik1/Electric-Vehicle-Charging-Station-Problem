#pragma once
#include <vector>

class EvaluationContext {
public:
    virtual ~EvaluationContext() = default;

    virtual int getWidth() const = 0;
    virtual int getHeight() const = 0;

    virtual int getPoiAt(int i, int j) const = 0;
    virtual double getDemandAt(int i, int j) const = 0;
    virtual double getRentalCostAt(int i, int j) const = 0;
    virtual double getDistance(int m, int n, int i, int j) const = 0;

    virtual int getMaxStationsPerCell() const = 0;
    virtual double getBudget() const = 0;
    
    virtual std::pair<double, double> getStationPowers() const = 0;
    virtual std::pair<double, double> getInitialCosts() const = 0;
    virtual std::pair<double, double> getMaintenanceCosts() const = 0;
};