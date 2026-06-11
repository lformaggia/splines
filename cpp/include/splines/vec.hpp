/**
 * @file vec.hpp
 * @brief Fixed-size spatial vector `Vec<Dim, Real>`.
 *
 * Wraps `std::array<Real, Dim>` with arithmetic operators, dot/cross
 * products, norm utilities, and convenience type aliases.  The type is
 * `constexpr`-friendly: all arithmetic operators and `dot` are `constexpr`.
 *
 * @par Supported operations
 * - Element-wise: `+`, `-`, unary `-`, `*` (scalar), `/` (scalar)
 * - Compound assignment: `+=`, `-=`, `*=`, `/=`
 * - Products: `dot(a,b)`, `cross(a,b)` (3-D only)
 * - Norms: `norm_sq(v)`, `norm(v)`, `normalized(v)`
 * - Comparison: `operator==` (defaulted, exact equality)
 * - Range iteration: `begin()` / `end()` via underlying `std::array`
 *
 * @par Convenience aliases
 * | Alias    | Expansion              |
 * |----------|------------------------|
 * | `Vec2d`  | `Vec<2, double>`       |
 * | `Vec3d`  | `Vec<3, double>`       |
 * | `Vec2f`  | `Vec<2, float>`        |
 * | `Vec3f`  | `Vec<3, float>`        |
 */
#pragma once
#include <array>
#include <cmath>
#include <concepts>
#include "concepts.hpp"

namespace splines {

/**
 * @brief Fixed-size spatial vector wrapping `std::array<Real, Dim>`.
 *
 * @tparam Dim   Number of components.  Must be positive.
 * @tparam Real  Floating-point scalar type (default `double`).
 */
template<int Dim, RealType Real = double>
struct Vec {
    static_assert(Dim > 0, "Dimension must be positive");

    std::array<Real, Dim> data{}; ///< Underlying storage (zero-initialised by default).

    /// @brief Default-construct, zero-initialising all components.
    constexpr Vec() = default;

    /**
     * @brief Construct from individual component values.
     *
     * The number of arguments must exactly equal `Dim`.
     *
     * @par Example
     * @code
     * Vec<3, double> v{1.0, 2.0, 3.0};
     * @endcode
     *
     * @tparam Args  Argument types (each convertible to `Real`).
     */
    template<typename... Args>
        requires (sizeof...(Args) == Dim)
    constexpr Vec(Args&&... args)
        : data{static_cast<Real>(std::forward<Args>(args))...} {}

    /// @brief Mutable element access by index.
    [[nodiscard]] constexpr Real&       operator[](int i)       { return data[i]; }
    /// @brief Const element access by index.
    [[nodiscard]] constexpr const Real& operator[](int i) const { return data[i]; }

    /// @brief Component-wise addition.
    [[nodiscard]] constexpr Vec operator+(const Vec& rhs) const noexcept {
        Vec r;
        for (int i = 0; i < Dim; ++i) r.data[i] = data[i] + rhs.data[i];
        return r;
    }
    /// @brief Component-wise subtraction.
    [[nodiscard]] constexpr Vec operator-(const Vec& rhs) const noexcept {
        Vec r;
        for (int i = 0; i < Dim; ++i) r.data[i] = data[i] - rhs.data[i];
        return r;
    }
    /// @brief Scalar multiplication (vector × scalar).
    [[nodiscard]] constexpr Vec operator*(Real s) const noexcept {
        Vec r;
        for (int i = 0; i < Dim; ++i) r.data[i] = data[i] * s;
        return r;
    }
    /// @brief Scalar division.
    [[nodiscard]] constexpr Vec operator/(Real s) const noexcept {
        return *this * (Real(1) / s);
    }
    /// @brief Unary negation.
    [[nodiscard]] constexpr Vec operator-() const noexcept {
        Vec r;
        for (int i = 0; i < Dim; ++i) r.data[i] = -data[i];
        return r;
    }

    constexpr Vec& operator+=(const Vec& rhs) noexcept { *this = *this + rhs; return *this; }
    constexpr Vec& operator-=(const Vec& rhs) noexcept { *this = *this - rhs; return *this; }
    constexpr Vec& operator*=(Real s)         noexcept { *this = *this * s;   return *this; }
    constexpr Vec& operator/=(Real s)         noexcept { *this = *this / s;   return *this; }

    /// @brief Exact equality comparison (component-wise).
    [[nodiscard]] bool operator==(const Vec&) const = default;

    /// @name Range support (allows range-for)
    ///@{
    auto begin()        { return data.begin(); }
    auto end()          { return data.end();   }
    auto begin()  const { return data.begin(); }
    auto end()    const { return data.end();   }
    ///@}
};

/**
 * @brief Scalar × Vec (commutative scalar multiplication).
 *
 * @tparam Dim   Vector dimension.
 * @tparam Real  Scalar type.
 */
template<int Dim, RealType Real>
[[nodiscard]] constexpr Vec<Dim,Real> operator*(Real s, const Vec<Dim,Real>& v) noexcept {
    return v * s;
}

/**
 * @brief Inner (dot) product of two vectors.
 *
 * @tparam Dim   Vector dimension.
 * @tparam Real  Scalar type.
 * @return       \f$ \sum_{i=0}^{Dim-1} a_i \, b_i \f$
 */
template<int Dim, RealType Real>
[[nodiscard]] constexpr Real dot(const Vec<Dim,Real>& a, const Vec<Dim,Real>& b) noexcept {
    Real r{};
    for (int i = 0; i < Dim; ++i) r += a[i] * b[i];
    return r;
}

/**
 * @brief Squared Euclidean norm.
 * @return `dot(v, v)`
 */
template<int Dim, RealType Real>
[[nodiscard]] constexpr Real norm_sq(const Vec<Dim,Real>& v) noexcept { return dot(v, v); }

/**
 * @brief Euclidean norm (L² norm).
 * @return `std::sqrt(norm_sq(v))`
 */
template<int Dim, RealType Real>
[[nodiscard]] Real norm(const Vec<Dim,Real>& v) { return std::sqrt(norm_sq(v)); }

/**
 * @brief Unit vector in the direction of @p v.
 * @return `v / norm(v)`
 * @warning Behaviour is undefined if `norm(v) == 0`.
 */
template<int Dim, RealType Real>
[[nodiscard]] Vec<Dim,Real> normalized(const Vec<Dim,Real>& v) { return v / norm(v); }

/**
 * @brief Cross product (defined only for 3-D vectors).
 *
 * @tparam Real  Scalar type.
 * @return       \f$ \mathbf{a} \times \mathbf{b} \f$
 */
template<RealType Real>
[[nodiscard]] constexpr Vec<3,Real> cross(const Vec<3,Real>& a, const Vec<3,Real>& b) noexcept {
    return Vec<3,Real>{
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

/// @name Convenience type aliases
///@{
using Vec2d = Vec<2, double>; ///< 2-D double-precision vector.
using Vec3d = Vec<3, double>; ///< 3-D double-precision vector.
using Vec2f = Vec<2, float>;  ///< 2-D single-precision vector.
using Vec3f = Vec<3, float>;  ///< 3-D single-precision vector.
///@}

} // namespace splines
