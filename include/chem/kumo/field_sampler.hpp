#pragma once

#include "chem/kumo/density.hpp"

#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace chem::kumo {

struct FieldSample {
    double value = 0.0;
    bool inside = false;
};

/// Trilinear interpolation on a general affine grid, returning zero outside.
/// The referenced Field must outlive this sampler and remain unchanged.
class FieldSampler {
public:
    explicit FieldSampler(const Field& field) : field_(field) {
        const Grid& grid = field_.grid();
        if (!grid.origin.allFinite() || !grid.steps.allFinite())
            throw std::invalid_argument("nonfinite field grid");
        for (int axis = 0; axis < 3; ++axis) {
            if (grid.counts[axis] < 2)
                throw std::invalid_argument("trilinear sampling requires two nodes per axis");
            const double length = grid.steps.col(axis).norm();
            if (!(length > 0.0) || !std::isfinite(length))
                throw std::invalid_argument("invalid grid step");
            lengths_[axis] = length;
            normalized_.col(axis) = grid.steps.col(axis) / length;
        }
        const double determinant = normalized_.determinant();
        if (!std::isfinite(determinant) || std::abs(determinant) < 1e-12)
            throw std::invalid_argument("singular or ill-conditioned field grid");
        inverse_ = normalized_.inverse();
        for (int axis = 0; axis < 3; ++axis)
            inverse_.row(axis) /= lengths_[axis];
        if (!inverse_.allFinite())
            throw std::invalid_argument("invalid inverse grid transform");
    }

    [[nodiscard]] FieldSample sample(const Eigen::Vector3d& point_bohr) const {
        if (!point_bohr.allFinite())
            throw std::invalid_argument("nonfinite sampling point");
        const Eigen::Vector3d delta = point_bohr - field_.grid().origin;
        if (!delta.allFinite())
            throw std::invalid_argument("nonfinite relative position");
        const Eigen::Vector3d fractional = inverse_ * delta;
        if (!fractional.allFinite())
            throw std::invalid_argument("nonfinite grid coordinate");
        Eigen::Vector3i lower;
        Eigen::Vector3d fraction;
        for (int axis = 0; axis < 3; ++axis) {
            const double upper = static_cast<double>(field_.grid().counts[axis] - 1);
            const double epsilon = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, upper);
            if (fractional[axis] < -epsilon || fractional[axis] > upper + epsilon)
                return {0.0, false};
            const double q = std::clamp(fractional[axis], 0.0, upper);
            lower[axis] = std::min(static_cast<int>(std::floor(q)), field_.grid().counts[axis] - 2);
            fraction[axis] = q - static_cast<double>(lower[axis]);
        }
        double value = 0.0;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const double weight = (i ? fraction[0] : 1.0 - fraction[0]) *
                                          (j ? fraction[1] : 1.0 - fraction[1]) *
                                          (k ? fraction[2] : 1.0 - fraction[2]);
                    const double sample_value = field_(lower[0] + i, lower[1] + j, lower[2] + k);
                    if (!std::isfinite(sample_value))
                        throw std::invalid_argument("nonfinite field value");
                    value += weight * sample_value;
                }
        if (!std::isfinite(value))
            throw std::overflow_error("field interpolation overflow");
        return {value, true};
    }

private:
    const Field& field_;
    Eigen::Vector3d lengths_ = Eigen::Vector3d::Zero();
    Eigen::Matrix3d normalized_ = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d inverse_ = Eigen::Matrix3d::Zero();
};

} // namespace chem::kumo
